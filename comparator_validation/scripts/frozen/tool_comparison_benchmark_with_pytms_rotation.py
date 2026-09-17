#!/usr/bin/env python3
"""PhosphoFill tool comparison benchmark with optional real PyTMs run.

Compares:
  1. PhosphoFill Stage 0 internal-coordinate seed from the benchmark TSV
  2. PhosphoFill top-1 from benchmark TSV (S3)
  3. PhosphoFill best-of-3 from benchmark TSV (top3_best_rmsd)
  4. PyTMs actual output, if --run-pytms is enabled

For PyTMs, the script:
  - loads the original phosphorylated crystal/reference structure
  - saves original phosphate coordinates
  - strips SEP/TPO/PTR back to SER/THR/TYR
  - runs PyTMs phosphorylation in headless PyMOL
  - compares the resulting phosphate to the original crystal phosphate
  - uses symmetric phosphate oxygen RMSD

Usage, summary only:
    python3 tool_comparison_benchmark_with_pytms.py \
        --phosphofill-tsv TPO_final_res.tsv SEP_final_res.tsv PTR_final_res.tsv \
        --output tool_comparison_summary.tsv

Usage, real PyTMs on a subset:
    python3 tool_comparison_benchmark_with_pytms.py \
        --phosphofill-tsv TPO_final_res.tsv SEP_final_res.tsv PTR_final_res.tsv \
        --run-pytms \
        --pdb-cache-dir pdb_cache \
        --pymol-path pymol \
        --pytms-workdir pytms_work \
        --max-sites-per-restype 50 \
        --output tool_comparison_summary.tsv \
        --pytms-detail-output pytms_detail.tsv

PyTMs must be available inside the PyMOL Python environment.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import subprocess
import sys
import time
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select
    from Bio.PDB.Atom import Atom, DisorderedAtom
except ImportError:
    print("Required: pip install biopython pandas numpy")
    sys.exit(1)

CIF_PARSER = MMCIFParser(QUIET=True)
PDB_PARSER = PDBParser(QUIET=True)

RESTYPE_MAP = {"T": "TPO", "S": "SEP", "Y": "PTR", "TPO": "TPO", "SEP": "SEP", "PTR": "PTR"}
PARENT_MAP = {"SEP": "SER", "TPO": "THR", "PTR": "TYR"}
PHOSPHO_ATOMS = {"P", "O1P", "O2P", "O3P"}
PHOSPHO_ONLY_ATOMS = {"P", "O1P", "O2P", "O3P", "OP1", "OP2", "OP3", "H1P", "H2P", "H3P"}

COL_ALIASES = {
    "pdb_id": ["ref_pdb_id", "pdb_id", "PDBID", "pdb", "PDB"],
    "chain": ["ref_chain", "chain", "chain_id", "CHAINID"],
    "position": ["position", "target_position", "resseq", "Modified_Sites_Auth", "site_position"],
    "restype": ["restype", "new_resname", "modres"],
    "acc_id": ["acc_id", "ACC_ID", "uniprot", "uniprot_id"],
    "ref_path": ["ref_path", "reference_path", "pdb_path", "structure_path", "input_structure"],
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compare PhosphoFill vs naive S0 and optional real PyTMs.")
    p.add_argument("--phosphofill-tsv", nargs="+", required=True, help="PhosphoFill final result TSVs")
    p.add_argument("--output", default="tool_comparison_summary.tsv", help="Summary output TSV")
    p.add_argument("--detail-output", default="tool_comparison_detail.tsv", help="Merged per-site detail output TSV")
    p.add_argument("--run-pytms", action="store_true", help="Actually run PyTMs through headless PyMOL")
    p.add_argument("--pymol-path", default="pymol", help="PyMOL executable")
    p.add_argument("--pdb-cache-dir", default="pdb_cache", help="PDB/mmCIF cache dir")
    p.add_argument("--pytms-workdir", default="pytms_work", help="Working dir for stripped/PyTMs outputs")
    p.add_argument("--max-sites-per-restype", type=int, default=None, help="Optional balanced subset per residue type")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--timeout", type=int, default=90, help="PyMOL timeout per site in seconds")
    p.add_argument("--keep-workdir", action="store_true", help="Keep PyTMs intermediate files")
    p.add_argument("--pytms-detail-output", default="pytms_detail.tsv", help="PyTMs per-site detail TSV")
    return p.parse_args()


def find_col(df: pd.DataFrame, key: str) -> Optional[str]:
    for c in COL_ALIASES[key]:
        if c in df.columns:
            return c
    return None


def normalize_pos(x) -> Optional[int]:
    if pd.isna(x):
        return None
    text = str(x).strip()
    if not text:
        return None
    try:
        return int(float(text))
    except Exception:
        return None


def normalize_restype(x) -> str:
    val = str(x).strip().upper()
    if val in RESTYPE_MAP:
        return RESTYPE_MAP[val]
    if val in {"SER", "S"}:
        return "SEP"
    if val in {"THR", "T"}:
        return "TPO"
    if val in {"TYR", "Y"}:
        return "PTR"
    return val


def load_structure(path: str, sid: str = "s"):
    suffix = Path(path).suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return CIF_PARSER.get_structure(sid, path)
    return PDB_PARSER.get_structure(sid, path)


def download_pdb(pdb_id: str, cache_dir: str) -> Optional[str]:
    if not pdb_id or str(pdb_id).lower() == "nan":
        return None
    pdb_id = str(pdb_id).strip().lower()
    os.makedirs(cache_dir, exist_ok=True)
    cif_path = os.path.join(cache_dir, f"{pdb_id}.cif")
    if os.path.exists(cif_path):
        return cif_path
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    try:
        urllib.request.urlretrieve(url, cif_path)
        return cif_path
    except Exception as exc:
        print(f"  WARNING: failed to download {pdb_id}: {exc}")
        return None


def find_residue(model, chain_id: str, resseq: int):
    if chain_id not in model:
        return None
    chain = model[chain_id]
    for hetfield in [" ", "H_TPO", "H_SEP", "H_PTR"]:
        try:
            return chain[(hetfield, int(resseq), " ")]
        except KeyError:
            pass
    for res in chain:
        if res.id[1] == int(resseq):
            return res
    return None


def norm_atom_name(name: str) -> str:
    n = str(name).strip().upper()
    if n == "OP1":
        return "O1P"
    if n == "OP2":
        return "O2P"
    if n == "OP3":
        return "O3P"
    return n


def get_phosphate_coords(residue) -> Dict[str, np.ndarray]:
    coords = {}
    if residue is None:
        return coords
    for atom in residue.get_atoms():
        n = norm_atom_name(atom.get_name())
        if n in PHOSPHO_ATOMS:
            coords[n] = np.array(atom.coord, dtype=float).copy()
    return coords


# Dihedral atoms per residue type (atom1-atom2-atom3-atom4)
DIHEDRAL_ATOMS = {
    "TPO": ("CG2", "CB", "OG1", "P"),
    "SEP": ("CA", "CB", "OG", "P"),
    "PTR": ("CE1", "CZ", "OH", "P"),
}


def calc_dihedral(p1, p2, p3, p4):
    """Calculate dihedral angle in degrees from four 3D points."""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    if n1_norm < 1e-8 or n2_norm < 1e-8:
        return float("nan")
    n1 = n1 / n1_norm
    n2 = n2 / n2_norm
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return float(np.degrees(np.arctan2(-y, x)))


def measure_phosphate_dihedral(residue, restype_3: str) -> Optional[float]:
    """Measure the key phosphate dihedral angle for a residue.

    Returns the dihedral in degrees, or None if atoms are missing.
    """
    if residue is None or restype_3 not in DIHEDRAL_ATOMS:
        return None
    a1_name, a2_name, a3_name, a4_name = DIHEDRAL_ATOMS[restype_3]

    # Build atom coord lookup (handle alt naming)
    atom_coords = {}
    for atom in residue.get_atoms():
        n = atom.get_name().strip().upper()
        # Normalize phosphate oxygen names
        if n == "OP1": n = "O1P"
        elif n == "OP2": n = "O2P"
        elif n == "OP3": n = "O3P"
        atom_coords[n] = np.array(atom.coord, dtype=float)

    for name in [a1_name, a2_name, a3_name, a4_name]:
        if name not in atom_coords:
            return None

    return calc_dihedral(
        atom_coords[a1_name], atom_coords[a2_name],
        atom_coords[a3_name], atom_coords[a4_name]
    )


N_FRAMES = 12
STEP_DEG = 30.0
FRAME_ANGLES = [i * STEP_DEG for i in range(N_FRAMES)]


def nearest_frame_index(angle_deg: float) -> Optional[int]:
    """Map an angle to the nearest 30° frame index (0-11)."""
    if angle_deg is None or np.isnan(angle_deg):
        return None
    angle_deg = angle_deg % 360
    dists = [min(abs(angle_deg - fa), 360 - abs(angle_deg - fa)) for fa in FRAME_ANGLES]
    return int(np.argmin(dists))


def phosphate_sym_rmsd(ref: Dict[str, np.ndarray], test: Dict[str, np.ndarray]) -> Optional[float]:
    if "P" not in ref or "P" not in test:
        return None
    ref_ox = [k for k in ["O1P", "O2P", "O3P"] if k in ref]
    test_ox = [k for k in ["O1P", "O2P", "O3P"] if k in test]
    n_ox = min(len(ref_ox), len(test_ox))
    if n_ox == 0:
        return float(np.linalg.norm(ref["P"] - test["P"]))
    best = None
    for perm in itertools.permutations(test_ox, n_ox):
        sq = [np.sum((ref["P"] - test["P"]) ** 2)]
        for rn, tn in zip(ref_ox[:n_ox], perm):
            sq.append(np.sum((ref[rn] - test[tn]) ** 2))
        rmsd = float(np.sqrt(np.mean(sq)))
        if best is None or rmsd < best:
            best = rmsd
    return best


def p_only_dist(ref: Dict[str, np.ndarray], test: Dict[str, np.ndarray]) -> Optional[float]:
    if "P" not in ref or "P" not in test:
        return None
    return float(np.linalg.norm(ref["P"] - test["P"]))


class StripPhosphoSelect(Select):
    def __init__(self, chain_id: str, resseq: int):
        super().__init__()
        self.chain_id = chain_id
        self.resseq = int(resseq)

    def accept_atom(self, atom):
        res = atom.get_parent()
        chain = res.get_parent()
        if chain.id == self.chain_id and res.id[1] == self.resseq:
            an = norm_atom_name(atom.get_name())
            if an in PHOSPHO_ONLY_ATOMS:
                return False
        return True


def strip_phospho_to_parent(input_path: str, chain_id: str, resseq: int, phospho_resname: str, output_pdb: str) -> bool:
    structure = load_structure(input_path, "strip")
    model = list(structure.get_models())[0]
    res = find_residue(model, chain_id, resseq)
    if res is None:
        return False
    parent = PARENT_MAP.get(phospho_resname)
    if parent is None:
        return False

    # Directly remove phosphate atoms from the residue
    atoms_to_remove = []
    for atom in list(res.get_atoms()):
        aname = norm_atom_name(atom.get_name())
        if aname in PHOSPHO_ONLY_ATOMS:
            atoms_to_remove.append(atom.get_id())
    for aid in atoms_to_remove:
        res.detach_child(aid)

    # Change residue name to parent
    res.resname = parent

    # Fix hetfield to standard residue
    old_id = res.id
    if old_id[0] != " ":
        new_id = (" ", old_id[1], old_id[2])
        chain = res.get_parent()
        if old_id in chain.child_dict:
            del chain.child_dict[old_id]
        chain.child_dict[new_id] = res
        res._id = new_id

    # Handle DisorderedAtom containers (altlocs)
    for child in list(res.get_list()):
        if isinstance(child, DisorderedAtom):
            sel_child = child.selected_child
            new_atom = Atom(
                name=sel_child.get_name(),
                coord=sel_child.coord.copy(),
                bfactor=sel_child.get_bfactor(),
                occupancy=1.0,
                altloc=" ",
                fullname=sel_child.get_fullname(),
                serial_number=sel_child.get_serial_number(),
                element=sel_child.element,
            )
            atom_id = child.get_id()
            res.detach_child(atom_id)
            res.add(new_atom)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    return os.path.exists(output_pdb)


PYMOL_PYTMS_SCRIPT = r'''
from pymol import cmd
import json
import traceback

input_pdb = r"{input_pdb}"
chain = "{chain}"
resseq = "{resseq}"
output_pdb = r"{output_pdb}"
status_json = output_pdb + ".status.json"

cmd.reinitialize()
cmd.load(input_pdb, "protein")
sel = f"protein and chain {chain} and resi {resseq}"
n_before = cmd.count_atoms(sel)
status = "OK"
message = ""
try:
    pytms = None
    try:
        import pmg_tk.startup.pytms as pytms
    except Exception:
        try:
            import pytms
        except Exception:
            pytms = None
    if pytms is None:
        raise RuntimeError("Could not import PyTMs in PyMOL Python environment")
    if hasattr(pytms, "phosphorylation"):
        pytms.phosphorylation(selection=sel, optimize=1, interval=30, remove_radius=0)
    elif hasattr(pytms, "phosphorylate"):
        pytms.phosphorylate(selection=sel, optimize=1, interval=30, remove_radius=0)
    elif hasattr(pytms, "add_phosphorylation"):
        pytms.add_phosphorylation(selection=sel, optimize=1, interval=30, remove_radius=0)
    else:
        raise RuntimeError("PyTMs imported, but no phosphorylation function found")
    n_after = cmd.count_atoms(sel)
    cmd.save(output_pdb, "protein")
except Exception as exc:
    status = "FAILED"
    message = str(exc) + "\n" + traceback.format_exc()
    n_after = cmd.count_atoms(sel)
    try:
        cmd.save(output_pdb, "protein")
    except Exception:
        pass
with open(status_json, "w") as fh:
    json.dump({{"status": status, "message": message, "n_atoms_before": n_before, "n_atoms_after": n_after, "selection": sel}}, fh)
cmd.quit()
'''


def run_pytms(input_pdb: str, chain: str, resseq: int, output_pdb: str, pymol_path: str, timeout: int) -> Tuple[bool, str]:
    script_path = output_pdb + ".pml.py"
    status_path = output_pdb + ".status.json"
    deadline = time.monotonic() + timeout

    for stale_path in (output_pdb, status_path):
        try:
            os.remove(stale_path)
        except FileNotFoundError:
            pass
        except OSError as exc:
            return False, f"Could not remove stale PyTMs output {stale_path}: {exc}"

    Path(script_path).write_text(PYMOL_PYTMS_SCRIPT.format(input_pdb=input_pdb, chain=chain, resseq=resseq, output_pdb=output_pdb), encoding="utf-8")
    try:
        proc = subprocess.run([pymol_path, "-cq", "-r", script_path], capture_output=True, text=True, timeout=timeout)
    except Exception as exc:
        return False, f"PyMOL failed: {exc}"

    while not os.path.exists(status_path) and time.monotonic() < deadline:
        time.sleep(0.2)

    if not os.path.exists(status_path):
        return False, f"No status JSON after waiting {timeout}s. PyMOL may not have run the script or may still be detached. stdout={proc.stdout[-500:]} stderr={proc.stderr[-500:]}"
    try:
        status = json.loads(Path(status_path).read_text())
    except Exception as exc:
        return False, f"Could not read status JSON: {exc}"
    if status.get("status") != "OK":
        return False, status.get("message", "PyTMs failed")
    if not os.path.exists(output_pdb):
        return False, "PyTMs did not write output PDB"
    return True, ""


def load_phosphofill_results(tsv_paths: List[str]) -> pd.DataFrame:
    frames = []
    for path in tsv_paths:
        df = pd.read_csv(path, sep="\t")
        if "status" in df.columns:
            df = df[df["status"] == "OK"].copy()
        else:
            df = df.copy()
        basename = Path(path).stem.upper()
        primary_one = None
        if "TPO" in basename or basename.startswith("T_"):
            primary_one = "T"
        elif "SEP" in basename or basename.startswith("S_"):
            primary_one = "S"
        elif "PTR" in basename or basename.startswith("Y_"):
            primary_one = "Y"
        if "restype" in df.columns:
            df["restype_3"] = df["restype"].apply(normalize_restype)
            if primary_one:
                df = df[df["restype_3"] == RESTYPE_MAP[primary_one]].copy()
        else:
            df["restype"] = primary_one if primary_one else ""
            df["restype_3"] = RESTYPE_MAP.get(primary_one, "UNKNOWN")
        df["_source_tsv"] = path
        frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def balanced_subset(df: pd.DataFrame, max_per_restype: Optional[int], seed: int) -> pd.DataFrame:
    if max_per_restype is None:
        return df
    out = []
    rng = np.random.RandomState(seed)
    for _, sub in df.groupby("restype_3"):
        if len(sub) > max_per_restype:
            out.append(sub.sample(n=max_per_restype, random_state=rng))
        else:
            out.append(sub)
    return pd.concat(out, ignore_index=True) if out else df.iloc[0:0].copy()


def resolve_reference_path(row: pd.Series, df: pd.DataFrame, pdb_cache_dir: str) -> Optional[str]:
    ref_col = find_col(df, "ref_path")
    if ref_col is not None:
        p = str(row.get(ref_col, "")).strip()
        if p and p.lower() != "nan" and os.path.exists(p):
            return p
    pdb_col = find_col(df, "pdb_id")
    if pdb_col is None:
        return None
    return download_pdb(str(row.get(pdb_col, "")).strip(), pdb_cache_dir)


def run_pytms_benchmark(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    os.makedirs(args.pytms_workdir, exist_ok=True)
    rows = []
    work = df.copy()
    if "context_n_sites" in work.columns:
        work = work[work["context_n_sites"] == 1].copy()
    work = balanced_subset(work, args.max_sites_per_restype, args.seed)
    print(f"\nRunning PyTMs on {len(work)} single-site rows...")
    chain_col = find_col(work, "chain")
    pos_col = find_col(work, "position")
    if chain_col is None or pos_col is None:
        raise ValueError(f"Cannot run PyTMs. Need chain/position columns. Columns: {list(work.columns)}")
    for _, row in work.iterrows():
        restype = normalize_restype(row.get("restype_3", row.get("restype", "")))
        if restype not in PARENT_MAP:
            continue
        chain = str(row.get(chain_col, "")).strip()
        pos = normalize_pos(row.get(pos_col))
        if not chain or pos is None:
            continue
        ref_path = resolve_reference_path(row, work, args.pdb_cache_dir)
        site_id = f"{row.get('acc_id', row.get('ACC_ID', 'NA'))}_{restype}_{Path(str(ref_path)).stem if ref_path else 'NA'}_{chain}_{pos}"
        try:
            if ref_path is None or not os.path.exists(ref_path):
                raise RuntimeError("missing_reference")
            ref_structure = load_structure(ref_path, "ref")
            ref_model = list(ref_structure.get_models())[0]
            ref_res = find_residue(ref_model, chain, pos)
            ref_phos = get_phosphate_coords(ref_res)
            if "P" not in ref_phos:
                raise RuntimeError("reference phosphate not found")
            stripped_pdb = os.path.join(args.pytms_workdir, f"{site_id}_stripped.pdb")
            pytms_pdb = os.path.join(args.pytms_workdir, f"{site_id}_pytms.pdb")
            if not strip_phospho_to_parent(ref_path, chain, pos, restype, stripped_pdb):
                raise RuntimeError("failed to strip phospho residue")
            ok_run, msg = run_pytms(stripped_pdb, chain, pos, pytms_pdb, args.pymol_path, args.timeout)
            if not ok_run:
                raise RuntimeError(msg)
            out_structure = load_structure(pytms_pdb, "pytms")
            out_model = list(out_structure.get_models())[0]
            out_res = find_residue(out_model, chain, pos)
            out_phos = get_phosphate_coords(out_res)
            rmsd = phosphate_sym_rmsd(ref_phos, out_phos)
            pdist = p_only_dist(ref_phos, out_phos)

            # Measure dihedral angles
            ref_dihedral = measure_phosphate_dihedral(ref_res, restype)
            pytms_dihedral = measure_phosphate_dihedral(out_res, restype)

            # Map dihedrals to nearest 30° frame index
            ref_frame = nearest_frame_index(ref_dihedral) if ref_dihedral is not None else None
            pytms_frame = nearest_frame_index(pytms_dihedral) if pytms_dihedral is not None else None

            # PhosphoFill's selected rotation from benchmark result
            # stage2_selected_angle = the rotation PhosphoFill's scoring picked (top-1)
            # stage1_best_of_12_angle = the oracle best rotation
            pf_dihedral = None
            pf_selected_angle = row.get("stage2_selected_angle")
            if pd.notna(pf_selected_angle):
                try:
                    pf_dihedral = float(pf_selected_angle)
                except (ValueError, TypeError):
                    pass

            oracle_angle = None
            oracle_best = row.get("stage1_best_of_12_angle")
            if pd.notna(oracle_best):
                try:
                    oracle_angle = float(oracle_best)
                except (ValueError, TypeError):
                    pass

            rows.append({
                "status": "OK", "message": "", "restype_3": restype, "chain": chain, "position": pos,
                "reference_path": ref_path, "stripped_input": stripped_pdb, "pytms_output": pytms_pdb,
                "pytms_phosphate_sym_rmsd": round(rmsd, 4) if rmsd is not None else np.nan,
                "pytms_p_only_dist": round(pdist, 4) if pdist is not None else np.nan,
                "stage0_kabsch_rmsd": row.get("stage0_kabsch_rmsd", np.nan),
                "stage3_post_minimization_rmsd": row.get("stage3_post_minimization_rmsd", np.nan),
                "top3_best_rmsd": row.get("top3_best_rmsd", np.nan),
                # Dihedral angles
                "crystal_dihedral_deg": round(ref_dihedral, 1) if ref_dihedral is not None else np.nan,
                "crystal_frame": ref_frame,
                "pytms_dihedral_deg": round(pytms_dihedral, 1) if pytms_dihedral is not None else np.nan,
                "pytms_frame": pytms_frame,
                "pf_selected_angle": round(pf_dihedral, 1) if pf_dihedral is not None else np.nan,
                "pf_frame": nearest_frame_index(pf_dihedral) if pf_dihedral is not None else None,
                "oracle_angle": round(oracle_angle, 1) if oracle_angle is not None else np.nan,
                "oracle_frame": nearest_frame_index(oracle_angle) if oracle_angle is not None else None,
                # Frame match comparisons
                # PyTMs vs crystal: both absolute dihedrals, directly comparable
                "pytms_matches_crystal_frame": (pytms_frame == ref_frame) if (pytms_frame is not None and ref_frame is not None) else np.nan,
                # PF vs oracle: both relative rotations from Kabsch start, directly comparable
                "pf_matches_oracle_frame": (nearest_frame_index(pf_dihedral) == nearest_frame_index(oracle_angle)) if (pf_dihedral is not None and oracle_angle is not None) else np.nan,
                "oracle_matches_crystal_frame": (nearest_frame_index(oracle_angle) == ref_frame) if (oracle_angle is not None and ref_frame is not None) else np.nan,
                #
                "acc_id": row.get("acc_id", row.get("ACC_ID", "")), "source_tsv": row.get("_source_tsv", ""),
            })
            print(f"  OK {site_id}: PyTMs={rmsd:.3f}Å dih={pytms_dihedral:.0f}° | crystal={ref_dihedral:.0f}° | PF={pf_dihedral:.0f}° | oracle={oracle_angle:.0f}°" if all(v is not None for v in [rmsd, pytms_dihedral, ref_dihedral, pf_dihedral, oracle_angle]) else f"  OK {site_id}: PyTMs={rmsd:.3f}Å" if rmsd else f"  OK {site_id}")
        except Exception as exc:
            rows.append({"status": "FAILED", "message": str(exc), "restype_3": restype, "chain": chain, "position": pos, "reference_path": ref_path, "acc_id": row.get("acc_id", row.get("ACC_ID", "")), "source_tsv": row.get("_source_tsv", "")})
            print(f"  FAIL {site_id}: {exc}")
    detail = pd.DataFrame(rows)
    detail.to_csv(args.pytms_detail_output, sep="\t", index=False)
    return detail


def summarize_phosphofill_and_naive(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    work = df.copy()
    if "context_n_sites" in work.columns:
        work = work[work["context_n_sites"] == 1].copy()
    for restype in ["TPO", "SEP", "PTR", "Combined"]:
        sub = work if restype == "Combined" else work[work["restype_3"] == restype]
        if len(sub) == 0:
            continue
        metrics = [
            ("PhosphoFill Stage 0 seed", pd.to_numeric(sub["stage0_kabsch_rmsd"], errors="coerce")),
            ("PhosphoFill top1/S3", pd.to_numeric(sub["stage3_post_minimization_rmsd"], errors="coerce")),
        ]
        if "top3_best_rmsd" in sub.columns:
            metrics.append(("PhosphoFill best-of-3", pd.to_numeric(sub["top3_best_rmsd"], errors="coerce")))
        for method, vals in metrics:
            vals = vals.dropna()
            if len(vals) == 0:
                continue
            rows.append({"Residue type": restype, "Method": method, "n": len(vals), "median_RMSD_A": round(vals.median(), 3), "mean_RMSD_A": round(vals.mean(), 3), "le1_pct": round((vals <= 1.0).mean() * 100, 1), "le1p5_pct": round((vals <= 1.5).mean() * 100, 1)})
    return pd.DataFrame(rows)


def summarize_pytms(detail: pd.DataFrame) -> pd.DataFrame:
    ok = detail[detail["status"] == "OK"].copy()
    rows = []
    for restype in ["TPO", "SEP", "PTR", "Combined"]:
        sub = ok if restype == "Combined" else ok[ok["restype_3"] == restype]
        if len(sub) == 0:
            continue
        vals = pd.to_numeric(sub["pytms_phosphate_sym_rmsd"], errors="coerce").dropna()
        if len(vals) == 0:
            continue
        rows.append({"Residue type": restype, "Method": "PyTMs actual", "n": len(vals), "median_RMSD_A": round(vals.median(), 3), "mean_RMSD_A": round(vals.mean(), 3), "le1_pct": round((vals <= 1.0).mean() * 100, 1), "le1p5_pct": round((vals <= 1.5).mean() * 100, 1)})
    return pd.DataFrame(rows)


def summarize_dihedrals(detail: pd.DataFrame) -> None:
    """Print dihedral and frame selection summary."""
    ok = detail[detail["status"] == "OK"].copy()
    if "crystal_dihedral_deg" not in ok.columns:
        return

    print("\n" + "=" * 72)
    print("DIHEDRAL / FRAME SELECTION COMPARISON")
    print("=" * 72)

    for restype in ["TPO", "SEP", "PTR"]:
        sub = ok[ok["restype_3"] == restype]
        n = len(sub)
        if n == 0:
            continue

        print(f"\n--- {restype} (n={n}) ---")

        # Crystal dihedral distribution
        crystal_dihedrals = pd.to_numeric(sub["crystal_dihedral_deg"], errors="coerce").dropna()
        if len(crystal_dihedrals) > 0:
            print(f"  Crystal dihedral: median={crystal_dihedrals.median():.0f}°, mean={crystal_dihedrals.mean():.0f}°")

        # PyTMs dihedral
        pytms_dihedrals = pd.to_numeric(sub["pytms_dihedral_deg"], errors="coerce").dropna()
        if len(pytms_dihedrals) > 0:
            print(f"  PyTMs  dihedral: median={pytms_dihedrals.median():.0f}°, mean={pytms_dihedrals.mean():.0f}°")

        # PhosphoFill rotation
        pf_angles = pd.to_numeric(sub["pf_selected_angle"], errors="coerce").dropna()
        if len(pf_angles) > 0:
            print(f"  PF     selected: median={pf_angles.median():.0f}°, mean={pf_angles.mean():.0f}°")

        # Oracle
        oracle_angles = pd.to_numeric(sub["oracle_angle"], errors="coerce").dropna()
        if len(oracle_angles) > 0:
            print(f"  Oracle best:     median={oracle_angles.median():.0f}°, mean={oracle_angles.mean():.0f}°")

        # Frame match rates
        pytms_match = sub["pytms_matches_crystal_frame"].dropna()
        pf_oracle_match = sub["pf_matches_oracle_frame"].dropna()
        if len(pytms_match) > 0:
            print(f"  PyTMs picks correct crystal frame: {pytms_match.sum():.0f}/{len(pytms_match)} ({pytms_match.mean()*100:.0f}%)")
        if len(pf_oracle_match) > 0:
            print(f"  PF picks oracle-best frame:        {pf_oracle_match.sum():.0f}/{len(pf_oracle_match)} ({pf_oracle_match.mean()*100:.0f}%)")

        # Frame distribution for PyTMs
        pytms_frames = sub["pytms_frame"].dropna()
        if len(pytms_frames) > 0:
            frame_counts = pytms_frames.value_counts().sort_index()
            print(f"  PyTMs frame distribution:")
            for fi, cnt in frame_counts.items():
                angle = FRAME_ANGLES[int(fi)]
                print(f"    Frame {int(fi):2d} ({angle:5.0f}°): {cnt:4d} ({cnt/len(pytms_frames)*100:5.1f}%)")


def main() -> None:
    args = parse_args()
    print("=" * 72)
    print("PhosphoFill tool comparison benchmark")
    print("=" * 72)
    pf = load_phosphofill_results(args.phosphofill_tsv)
    print(f"Loaded {len(pf)} OK PhosphoFill rows")
    pf.to_csv(args.detail_output, sep="\t", index=False)
    summaries = [summarize_phosphofill_and_naive(pf)]
    if args.run_pytms:
        pytms_detail = run_pytms_benchmark(pf, args)
        pytms_summary = summarize_pytms(pytms_detail)
        if not pytms_summary.empty:
            summaries.append(pytms_summary)
        summarize_dihedrals(pytms_detail)
    final = pd.concat(summaries, ignore_index=True)
    final.to_csv(args.output, sep="\t", index=False)
    print("\nSummary:")
    print(final.to_string(index=False))
    print(f"\nSummary output: {args.output}")
    print(f"Detail output: {args.detail_output}")
    if args.run_pytms:
        print(f"PyTMs detail output: {args.pytms_detail_output}")


if __name__ == "__main__":
    main()
