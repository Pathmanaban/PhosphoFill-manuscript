#!/usr/bin/env python3
"""PhosphoFill tool comparison — benchmark against PTM-Psi.

Same strip-and-regraft protocol as the PyTMs benchmark:
1. Load crystal phosphostructure
2. Strip phosphate → parent residue (SER/THR/TYR)
3. Run PTM-Psi modify() on the stripped structure
4. Measure phosphate RMSD against crystal reference

PTM-Psi uses the NERF algorithm for coordinate placement with
canonical bond lengths, angles, and dihedrals. No rotamer search
or environment-aware scoring.

Installation:
    git clone https://github.com/pnnl/PTMPSI.git
    cd PTMPSI
    pip install -e .

Usage:
    python3 tool_comparison_ptmpsi.py \
        --phosphofill-tsv TPO_final_res.tsv SEP_final_res.tsv PTR_final_res.tsv \
        --pdb-cache-dir pdb_cache \
        --ptmpsi-source-dir PTMPSI \
        --ptmpsi-workdir ptmpsi_work \
        --ptmpsi-timeout-sec 300 \
        --output ptmpsi_comparison.tsv \
        --detail-output ptmpsi_detail.tsv

    # Without PTM-Psi installed (just PhosphoFill vs naive summary):
    python3 tool_comparison_ptmpsi.py \
        --phosphofill-tsv TPO_final_res.tsv SEP_final_res.tsv PTR_final_res.tsv \
        --skip-ptmpsi \
        --output ptmpsi_comparison.tsv
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import multiprocessing as mp
import os
import queue
import sys
import traceback
import types
import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

try:
    from Bio.PDB import MMCIFParser, PDBParser, PDBIO
    from Bio.PDB.Atom import Atom, DisorderedAtom
except ImportError:
    print("Required: pip install biopython")
    sys.exit(1)

PTMPSI_AVAILABLE = False
PTMPSI_IMPORT_ERROR = ""
PTMPSI_IMPORT_TRACEBACK = ""
PTMPSI_IMPORT_MODE = ""
PtmpsiProtein = None

CIF_PARSER = MMCIFParser(QUIET=True)
PDB_PARSER = PDBParser(QUIET=True)

RESTYPE_MAP = {"T": "TPO", "S": "SEP", "Y": "PTR", "TPO": "TPO", "SEP": "SEP", "PTR": "PTR"}
PARENT_MAP = {"TPO": "THR", "SEP": "SER", "PTR": "TYR"}
PHOSPHO_ATOMS = {"P", "O1P", "O2P", "O3P"}
PHOSPHO_ONLY_ATOMS = {"P", "O1P", "O2P", "O3P", "OP1", "OP2", "OP3", "H1P", "H2P", "H3P"}
PTMPSI_PHOSPHATE_ATOM_RENAMES = {"O1": "O1P", "O2": "O2P", "O3": "O3P"}

# Dihedral measurement
DIHEDRAL_ATOMS = {
    "TPO": ("CG2", "CB", "OG1", "P"),
    "SEP": ("CA", "CB", "OG", "P"),
    "PTR": ("CE1", "CZ", "OH", "P"),
}


def normalize_restype(x):
    val = str(x).strip().upper()
    return RESTYPE_MAP.get(val, val)


def normalize_pos(x):
    try:
        return int(float(x))
    except (ValueError, TypeError):
        return None


def norm_atom_name(name: str) -> str:
    n = str(name).strip().upper()
    return {
        "OP1": "O1P",
        "OP2": "O2P",
        "OP3": "O3P",
        "O1": "O1P",
        "O2": "O2P",
        "O3": "O3P",
    }.get(n, n)


def _ptmpsi_source_candidates(explicit: Optional[str] = None) -> List[Path]:
    """Return possible source roots containing a ptmpsi package."""
    candidates = []
    if explicit:
        candidates.append(Path(explicit))

    script_dir = Path(__file__).resolve().parent
    for base in [Path.cwd(), script_dir]:
        candidates.append(base / "PTMPSI")

    try:
        spec = importlib.util.find_spec("ptmpsi")
        if spec and spec.submodule_search_locations:
            candidates.append(Path(next(iter(spec.submodule_search_locations))).parent)
    except (ImportError, AttributeError, StopIteration, ValueError):
        pass

    valid = []
    seen = set()
    for candidate in candidates:
        try:
            resolved = candidate.expanduser().resolve()
        except OSError:
            continue
        if (resolved / "ptmpsi").is_dir() and resolved not in seen:
            valid.append(resolved)
            seen.add(resolved)
    return valid


def _clear_ptmpsi_modules() -> None:
    for name in list(sys.modules):
        if name == "ptmpsi" or name.startswith("ptmpsi."):
            del sys.modules[name]


def _install_ptmpsi_minimal_shims(source_dir: Path) -> None:
    """Load ptmpsi.protein without importing optional top-level modules.

    PTM-Psi's package __init__ imports broad optional stacks such as docking,
    AlphaFold, and Gromacs. The benchmark only needs Protein.modify() and
    Protein.write_pdb(), so we create a minimal package module and docking shim
    to avoid failing on optional dependencies that are unrelated here.
    """
    package_dir = source_dir / "ptmpsi"

    pkg = types.ModuleType("ptmpsi")
    pkg.__file__ = str(package_dir / "__init__.py")
    pkg.__path__ = [str(package_dir)]
    pkg.__package__ = "ptmpsi"
    pkg.__version__ = "0.1.0"
    sys.modules["ptmpsi"] = pkg

    docking = types.ModuleType("ptmpsi.docking")

    class Dock:
        def __init__(self):
            self.engine = None
            self.ligand = None
            self.receptor = None
            self.flexible = None
            self.output = None
            self.boxsize = None
            self.boxcenter = None
            self.exhaustiveness = None

    def dock_ligand(*_args, **_kwargs):
        raise RuntimeError("PTM-Psi docking is not available in benchmark-only mode")

    docking.Dock = Dock
    docking.dock_ligand = dock_ligand
    sys.modules["ptmpsi.docking"] = docking


class _DeepcopyCompat:
    def __call__(self, obj):
        import copy
        return copy.deepcopy(obj)

    def deepcopy(self, obj):
        import copy
        return copy.deepcopy(obj)


def _patch_ptmpsi_copy_bug() -> None:
    try:
        import ptmpsi.protein.mutate as mutate
        mutate.copy = _DeepcopyCompat()
    except Exception:
        pass


def load_ptmpsi(source_dir: Optional[str] = None) -> bool:
    """Import PTM-Psi and keep the real import error for diagnostics."""
    global PTMPSI_AVAILABLE, PTMPSI_IMPORT_ERROR, PTMPSI_IMPORT_TRACEBACK
    global PTMPSI_IMPORT_MODE, PtmpsiProtein

    PTMPSI_AVAILABLE = False
    PTMPSI_IMPORT_ERROR = ""
    PTMPSI_IMPORT_TRACEBACK = ""
    PTMPSI_IMPORT_MODE = ""
    PtmpsiProtein = None

    candidates = _ptmpsi_source_candidates(source_dir)
    for candidate in reversed(candidates):
        cstr = str(candidate)
        if cstr not in sys.path:
            sys.path.insert(0, cstr)

    try:
        from ptmpsi.protein import Protein as ProteinClass
        _patch_ptmpsi_copy_bug()
        PtmpsiProtein = ProteinClass
        PTMPSI_AVAILABLE = True
        PTMPSI_IMPORT_MODE = "normal import"
        return True
    except Exception as exc:
        PTMPSI_IMPORT_ERROR = repr(exc)
        PTMPSI_IMPORT_TRACEBACK = traceback.format_exc()

    for candidate in candidates:
        try:
            _clear_ptmpsi_modules()
            _install_ptmpsi_minimal_shims(candidate)
            from ptmpsi.protein import Protein as ProteinClass
            _patch_ptmpsi_copy_bug()
            PtmpsiProtein = ProteinClass
            PTMPSI_AVAILABLE = True
            PTMPSI_IMPORT_MODE = f"source import with optional-dependency shims ({candidate})"
            return True
        except Exception as exc:
            PTMPSI_IMPORT_ERROR = repr(exc)
            PTMPSI_IMPORT_TRACEBACK = traceback.format_exc()

    return False


def load_structure(path: str, sid: str = "s"):
    suffix = Path(path).suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return CIF_PARSER.get_structure(sid, path)
    return PDB_PARSER.get_structure(sid, path)


def download_pdb(pdb_id: str, cache_dir: str) -> Optional[str]:
    pdb_id = pdb_id.lower()
    cif_path = os.path.join(cache_dir, f"{pdb_id}.cif")
    if os.path.exists(cif_path):
        return cif_path
    os.makedirs(cache_dir, exist_ok=True)
    import urllib.request
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    try:
        urllib.request.urlretrieve(url, cif_path)
        return cif_path
    except Exception:
        return None


def find_residue(model, chain_id: str, resseq: int):
    if chain_id not in model:
        return None
    for hetfield in [" ", "H_TPO", "H_SEP", "H_PTR"]:
        try:
            return model[chain_id][(hetfield, int(resseq), " ")]
        except KeyError:
            pass
    for res in model[chain_id]:
        if res.id[1] == int(resseq):
            return res
    return None


def get_phosphate_coords(residue) -> Dict[str, np.ndarray]:
    coords = {}
    if residue is None:
        return coords
    for atom in residue.get_atoms():
        n = norm_atom_name(atom.get_name())
        if n in PHOSPHO_ATOMS:
            coords[n] = np.array(atom.coord, dtype=float).copy()
    return coords


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


def calc_dihedral(p1, p2, p3, p4):
    b1, b2, b3 = p2 - p1, p3 - p2, p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)
    if n1_norm < 1e-8 or n2_norm < 1e-8:
        return float("nan")
    n1 = n1 / n1_norm
    n2 = n2 / n2_norm
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    return float(np.degrees(np.arctan2(-np.dot(m1, n2), np.dot(n1, n2))))


def measure_phosphate_dihedral(residue, restype_3: str) -> Optional[float]:
    if residue is None or restype_3 not in DIHEDRAL_ATOMS:
        return None
    a1, a2, a3, a4 = DIHEDRAL_ATOMS[restype_3]
    coords = {}
    for atom in residue.get_atoms():
        n = norm_atom_name(atom.get_name())
        coords[n] = np.array(atom.coord, dtype=float)
    for name in [a1, a2, a3, a4]:
        if name not in coords:
            return None
    return calc_dihedral(coords[a1], coords[a2], coords[a3], coords[a4])


# =====================================================================
# Strip phosphate (same as PyTMs benchmark)
# =====================================================================

def normalize_chain_ids_for_pdbio(structure, target_chain_id: str) -> str:
    """Map chain IDs to one-character IDs required by the PDB format."""
    alphabet = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")
    mapping = {}
    used = set()

    for model in structure:
        for chain in model:
            cid = str(chain.id)
            if len(cid) == 1 and cid not in used:
                mapping[cid] = cid
                used.add(cid)

    for model in structure:
        for chain in model:
            cid = str(chain.id)
            if cid in mapping:
                continue
            choices = [c for c in alphabet if c not in used]
            if not choices:
                raise RuntimeError("too many chains to map into PDB one-character chain IDs")
            mapping[cid] = choices[0]
            used.add(choices[0])

    for model in structure:
        for chain in list(model):
            chain.id = mapping[str(chain.id)]

    return mapping.get(str(target_chain_id), str(target_chain_id))


def strip_phospho_to_parent(input_path: str, chain_id: str, resseq: int,
                            phospho_resname: str, output_pdb: str) -> Optional[str]:
    structure = load_structure(input_path, "strip")
    model = list(structure.get_models())[0]
    res = find_residue(model, chain_id, resseq)
    if res is None:
        return None
    parent = PARENT_MAP.get(phospho_resname)
    if parent is None:
        return None

    # Remove phosphate atoms
    atoms_to_remove = []
    for child in list(res.get_list()):
        child_name = norm_atom_name(child.get_name())
        if child_name in PHOSPHO_ONLY_ATOMS:
            atoms_to_remove.append(child.get_id())
    for aid in atoms_to_remove:
        res.detach_child(aid)

    res.resname = parent

    # Fix hetfield
    old_id = res.id
    if old_id[0] != " ":
        new_id = (" ", old_id[1], old_id[2])
        chain = res.get_parent()
        if old_id in chain.child_dict:
            del chain.child_dict[old_id]
        chain.child_dict[new_id] = res
        res._id = new_id

    # Handle DisorderedAtom
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

    ptmpsi_chain_id = normalize_chain_ids_for_pdbio(structure, chain_id)
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)
    return ptmpsi_chain_id if os.path.exists(output_pdb) else None


def ptmpsi_residue_index(pdb_path: str, chain_id: str, original_resseq: int) -> Optional[int]:
    """Return PTM-Psi's chain-local residue index for an original PDB residue.

    PTM-Psi selects residues by their position in the parsed chain. The PDB/mmCIF
    residue number is often different when structures have missing residues or
    non-1-based numbering.
    """
    structure = load_structure(pdb_path, "ptmpsi_index")
    model = list(structure.get_models())[0]
    if chain_id not in model:
        return None
    idx = 0
    for res in model[chain_id]:
        if res.id[0] != " ":
            continue
        idx += 1
        if int(res.id[1]) == int(original_resseq):
            return idx
    return None


def renumber_pdb_residues_for_ptmpsi(pdb_path: str) -> bool:
    """Rewrite ATOM/HETATM residue numbers to 1-based chain-local indices.

    PTM-Psi's parser can crash on perfectly valid PDB residue numbering such as
    a first residue numbered 0. It also selects by parsed chain-local index, so
    this makes the file numbering match the benchmark's PTM-Psi selection.
    """
    counters: Dict[str, int] = {}
    residue_ids: Dict[Tuple[str, str, str], int] = {}
    changed = False
    lines = []

    with open(pdb_path, "r") as handle:
        for line in handle:
            if line.startswith("ATOM  ") and len(line) >= 27:
                chain = line[21:22]
                original_key = (chain, line[22:26], line[26:27])
                if original_key not in residue_ids:
                    counters[chain] = counters.get(chain, 0) + 1
                    residue_ids[original_key] = counters[chain]
                new_resseq = residue_ids[original_key]
                new_field = f"{new_resseq:4d}"
                if line[22:26] != new_field:
                    line = f"{line[:22]}{new_field}{line[26:]}"
                    changed = True
            lines.append(line)

    if changed:
        with open(pdb_path, "w") as handle:
            handle.writelines(lines)
    return changed


def sanitize_ptmpsi_input_pdb(pdb_path: str) -> bool:
    """Drop records PTM-Psi ignores but Biopython/PTM-Psi may parse poorly.

    Some PDBs contain 4-character ligand names, which cannot be represented in
    standard PDB columns without shifting the chain/residue-number fields.
    PTM-Psi is run with delhet=True, so HETATM records are not part of the
    method comparison; removing them makes parsing deterministic.
    """
    changed = False
    kept = []
    drop_prefixes = ("HETATM", "ANISOU", "CONECT", "MASTER", "TER")
    with open(pdb_path, "r") as handle:
        for line in handle:
            if line.startswith(drop_prefixes):
                changed = True
                continue
            kept.append(line)

    if changed:
        with open(pdb_path, "w") as handle:
            handle.writelines(kept)
    return changed


def normalize_ptmpsi_output_pdb(pdb_path: str, chain_id: str, ptmpsi_resseq: int,
                                restype_3: str) -> bool:
    """Rename PTM-Psi's generic PTM residue to SEP/TPO/PTR in-place.

    PTM-Psi stores phosphorylated SER/THR/TYR as residue name PTM and uses
    O1/O2/O3 phosphate atom names. The benchmark normalizes those labels so
    saved PTM-Psi outputs can be inspected alongside crystal phosphoresidues.
    """
    changed = False
    target_seq = int(ptmpsi_resseq)
    lines = []
    with open(pdb_path, "r") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 26:
                line_chain = line[21:22].strip()
                try:
                    line_seq = int(line[22:26])
                except ValueError:
                    line_seq = None
                if line_chain == chain_id and line_seq == target_seq:
                    atom_name = line[12:16].strip().upper()
                    new_atom_name = PTMPSI_PHOSPHATE_ATOM_RENAMES.get(atom_name)
                    if new_atom_name:
                        line = f"{line[:12]}{new_atom_name:>4s}{line[16:]}"
                    if line[17:20].strip().upper() == "PTM":
                        line = f"{line[:17]}{restype_3:>3s}{line[20:]}"
                    changed = True
            lines.append(line)

    if changed:
        with open(pdb_path, "w") as handle:
            handle.writelines(lines)
    return changed


# =====================================================================
# PTM-Psi runner
# =====================================================================

def run_ptmpsi_single(input_pdb: str, chain: str, resseq: int,
                      parent_resname: str, output_pdb: str) -> Tuple[bool, str]:
    """Run PTM-Psi phosphorylation on a single site.

    PTM-Psi API: Protein.modify("A:SER154", "phosphorylation")
    The selection format is "chain:RESNAMEresseq"
    """
    if not PTMPSI_AVAILABLE:
        return False, "PTM-Psi not installed"

    try:
        prot = PtmpsiProtein(input_pdb)
        drop_empty_ptmpsi_residues(prot)
        selection = f"{chain}:{parent_resname}{resseq}"

        # PTM-Psi only codes this PTM as "phosphorylation".
        prot.modify(selection, "phosphorylation")

        prot.write_pdb(output_pdb)
        if os.path.exists(output_pdb):
            return True, "OK"
        return False, "output not created"
    except Exception as exc:
        return False, f"{exc}\n{traceback.format_exc()[-200:]}"


def drop_empty_ptmpsi_residues(prot) -> int:
    """Remove parser artifacts with no coordinates from a PTM-Psi Protein."""
    dropped = 0
    for chain in getattr(prot, "chains", []) or []:
        residues = []
        for residue in getattr(chain, "residues", []) or []:
            coords = getattr(residue, "coordinates", None)
            natoms = int(getattr(residue, "natoms", 0) or 0)
            if coords is None or natoms == 0 or len(coords) == 0:
                dropped += 1
                continue
            residues.append(residue)
        chain.residues = residues
    if dropped:
        prot.update()
    return dropped


def _ptmpsi_single_worker(result_queue, source_dir: Optional[str], input_pdb: str,
                          chain: str, resseq: int, parent_resname: str,
                          output_pdb: str) -> None:
    try:
        load_ptmpsi(source_dir)
        if not PTMPSI_AVAILABLE:
            result_queue.put((False, f"PTM-Psi import failed: {PTMPSI_IMPORT_ERROR}"))
            return
        result_queue.put(run_ptmpsi_single(input_pdb, chain, resseq, parent_resname, output_pdb))
    except BaseException as exc:
        result_queue.put((False, f"{exc}\n{traceback.format_exc()[-200:]}"))


def run_ptmpsi_single_with_timeout(input_pdb: str, chain: str, resseq: int,
                                   parent_resname: str, output_pdb: str,
                                   timeout_sec: Optional[float],
                                   source_dir: Optional[str]) -> Tuple[bool, str]:
    """Run one PTM-Psi site with a hard timeout."""
    if timeout_sec is None or timeout_sec <= 0:
        return run_ptmpsi_single(input_pdb, chain, resseq, parent_resname, output_pdb)

    ctx_name = "spawn" if os.name == "nt" else "fork"
    ctx = mp.get_context(ctx_name)
    result_queue = ctx.Queue()
    proc = ctx.Process(
        target=_ptmpsi_single_worker,
        args=(result_queue, source_dir, input_pdb, chain, resseq, parent_resname, output_pdb),
    )
    proc.daemon = True
    proc.start()
    proc.join(timeout_sec)

    if proc.is_alive():
        proc.terminate()
        proc.join(5)
        if proc.is_alive() and hasattr(proc, "kill"):
            proc.kill()
            proc.join(5)
        return False, f"PTM-Psi timed out after {timeout_sec:g} seconds"

    try:
        return result_queue.get_nowait()
    except queue.Empty:
        return False, f"PTM-Psi worker exited with code {proc.exitcode} without a result"


def write_detail_checkpoint(rows: List[Dict], detail_output: str) -> None:
    if rows:
        pd.DataFrame(rows).to_csv(detail_output, sep="\t", index=False)


def write_timeout_site_list(rows: List[Dict], output_path: Optional[str]) -> None:
    if not rows or not output_path:
        return
    df = pd.DataFrame(rows)
    if "site_id" not in df.columns or "message" not in df.columns:
        return
    messages = df["message"].fillna("").astype(str)
    timed_out = df[(df["status"] == "FAILED") & messages.str.contains("timed out", case=False, regex=False)]
    site_ids = timed_out["site_id"].dropna().drop_duplicates()
    out = Path(output_path)
    if out.parent != Path("."):
        out.parent.mkdir(parents=True, exist_ok=True)
    site_ids.to_csv(out, index=False, header=False)


def read_site_id_list(path: Optional[str]) -> Optional[set]:
    if not path:
        return None
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"site id list not found: {path}")
    try:
        df = pd.read_csv(p, sep=None, engine="python")
        if "site_id" in df.columns:
            return {str(x).strip() for x in df["site_id"].dropna() if str(x).strip()}
    except Exception:
        pass
    site_ids = set()
    with open(p, "r") as handle:
        for line in handle:
            val = line.strip().split()[0] if line.strip() else ""
            if val and val.lower() != "site_id":
                site_ids.add(val)
    return site_ids


def read_rerun_site_ids(detail_path: Optional[str], mode: str) -> Optional[set]:
    if not detail_path:
        return None
    df = pd.read_csv(detail_path, sep="\t")
    if "site_id" not in df.columns:
        raise ValueError(f"{detail_path} does not contain a site_id column")

    work = df.copy()
    if mode == "timeouts":
        if "message" not in work.columns:
            return set()
        messages = work["message"].fillna("").astype(str)
        work = work[(work["status"] == "FAILED") & messages.str.contains("timed out", case=False, regex=False)]
    elif mode == "failed":
        work = work[work["status"] == "FAILED"]
    elif mode == "non-ok":
        work = work[work["status"] != "OK"]
    else:
        raise ValueError(f"unsupported rerun mode: {mode}")

    return {str(x).strip() for x in work["site_id"].dropna() if str(x).strip()}


def build_site_id(row) -> str:
    restype = normalize_restype(row.get("restype_3", row.get("restype", "")))
    chain = str(row.get("chain", "")).strip()
    pos = normalize_pos(row.get("position"))
    acc_id = row.get("acc_id", row.get("ACC_ID", "NA"))
    pdb_id = row.get("pdb_id", row.get("PDBID_1", "NA"))
    return f"{acc_id}_{restype}_{pdb_id}_{chain}_{pos}"


# =====================================================================
# Load PhosphoFill results
# =====================================================================

def load_phosphofill_results(tsv_paths: List[str]) -> pd.DataFrame:
    frames = []
    for path in tsv_paths:
        df = pd.read_csv(path, sep="\t")
        if "status" in df.columns:
            df = df[df["status"] == "OK"].copy()
        if "restype" in df.columns:
            df["restype_3"] = df["restype"].apply(normalize_restype)

        # Filter to primary type from filename
        basename = Path(path).stem.upper()
        if "TPO" in basename:
            primary = "TPO"
        elif "SEP" in basename:
            primary = "SEP"
        elif "PTR" in basename:
            primary = "PTR"
        else:
            primary = None

        if primary and "restype_3" in df.columns:
            df = df[df["restype_3"] == primary]

        df["_source_tsv"] = path
        frames.append(df)

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# =====================================================================
# Resolve reference path
# =====================================================================

def resolve_reference_path(row, pdb_cache_dir: str) -> Optional[str]:
    """Get path to reference PDB/CIF file."""
    # Try ref_path column first
    for col in ["ref_path", "reference_path", "input_structure"]:
        val = str(row.get(col, "")).strip()
        if val and val != "nan" and os.path.exists(val):
            return val

    # Download from PDB
    for col in ["pdb_id", "PDBID_1", "PDB_ID"]:
        pdb_id = str(row.get(col, "")).strip()
        if pdb_id and pdb_id != "nan" and len(pdb_id) == 4:
            return download_pdb(pdb_id, pdb_cache_dir)

    return None


# =====================================================================
# Main benchmark
# =====================================================================

def run_ptmpsi_benchmark(pf_df: pd.DataFrame, args) -> pd.DataFrame:
    """Run PTM-Psi on single-site entries and compare to crystal."""
    os.makedirs(args.ptmpsi_workdir, exist_ok=True)

    work = pf_df.copy()
    if "context_n_sites" in work.columns:
        work = work[work["context_n_sites"] == 1].copy()

    work["_site_id"] = work.apply(build_site_id, axis=1)
    site_filters = []
    list_filter = read_site_id_list(args.site_id_list)
    if list_filter is not None:
        site_filters.append(("site list", list_filter, args.site_id_list))
    rerun_filter = read_rerun_site_ids(args.rerun_from_detail, args.rerun_mode)
    if rerun_filter is not None:
        site_filters.append((f"{args.rerun_mode} rows", rerun_filter, args.rerun_from_detail))

    site_filter = None
    if site_filters:
        site_filter = set().union(*(item[1] for item in site_filters))

    if site_filter is not None:
        work = work[work["_site_id"].isin(site_filter)].copy()
        missing = len(site_filter) - len(set(work["_site_id"]))
        sources = ", ".join(f"{label} from {path}" for label, _ids, path in site_filters)
        msg = f"Restricting to {len(work)} site IDs from {sources}"
        if missing:
            msg += f" ({missing} not found)"
        print(msg)

    if args.max_sites_per_restype and site_filter is None:
        rng = np.random.RandomState(42)
        subsets = []
        for rt, sub in work.groupby("restype_3"):
            if len(sub) > args.max_sites_per_restype:
                subsets.append(sub.sample(n=args.max_sites_per_restype, random_state=rng))
            else:
                subsets.append(sub)
        work = pd.concat(subsets, ignore_index=True)

    print(f"Running PTM-Psi on {len(work)} single-site entries...")
    if args.ptmpsi_timeout_sec and args.ptmpsi_timeout_sec > 0:
        print(f"Per-site PTM-Psi timeout: {args.ptmpsi_timeout_sec:g} seconds")

    rows = []
    n_ok = 0
    n_fail = 0

    for idx, row in work.iterrows():
        restype = normalize_restype(row.get("restype_3", row.get("restype", "")))
        if restype not in PARENT_MAP:
            continue
        parent = PARENT_MAP[restype]

        chain = str(row.get("chain", "")).strip()
        pos = normalize_pos(row.get("position"))
        if not chain or pos is None:
            continue

        ref_path = resolve_reference_path(row, args.pdb_cache_dir)
        acc_id = row.get("acc_id", row.get("ACC_ID", "NA"))
        pdb_id = row.get("pdb_id", row.get("PDBID_1", "NA"))
        site_id = row.get("_site_id", f"{acc_id}_{restype}_{pdb_id}_{chain}_{pos}")

        try:
            if ref_path is None or not os.path.exists(ref_path):
                raise RuntimeError("missing reference")

            # Load reference and get phosphate coords
            ref_structure = load_structure(ref_path, "ref")
            ref_model = list(ref_structure.get_models())[0]
            ref_res = find_residue(ref_model, chain, pos)
            ref_phos = get_phosphate_coords(ref_res)
            if "P" not in ref_phos:
                raise RuntimeError("reference phosphate not found")

            # Measure crystal dihedral
            ref_dihedral = measure_phosphate_dihedral(ref_res, restype)

            # Strip phosphate
            stripped_pdb = os.path.join(args.ptmpsi_workdir, f"{site_id}_stripped.pdb")
            ptmpsi_pdb = os.path.join(args.ptmpsi_workdir, f"{site_id}_ptmpsi.pdb")

            ptmpsi_chain = strip_phospho_to_parent(ref_path, chain, pos, restype, stripped_pdb)
            if not ptmpsi_chain:
                raise RuntimeError("failed to strip phospho residue")

            stripped_sanitized = sanitize_ptmpsi_input_pdb(stripped_pdb)
            ptmpsi_pos = ptmpsi_residue_index(stripped_pdb, ptmpsi_chain, pos)
            if ptmpsi_pos is None:
                raise RuntimeError(f"could not map {ptmpsi_chain}:{pos} to PTM-Psi residue index")
            stripped_renumbered = renumber_pdb_residues_for_ptmpsi(stripped_pdb)

            # Run PTM-Psi
            print(f"  RUN {site_id}: PTM-Psi residue index {ptmpsi_pos}", flush=True)
            ok_run, msg = run_ptmpsi_single_with_timeout(
                stripped_pdb,
                ptmpsi_chain,
                ptmpsi_pos,
                parent,
                ptmpsi_pdb,
                args.ptmpsi_timeout_sec,
                args.ptmpsi_source_dir,
            )
            if not ok_run:
                raise RuntimeError(msg)
            normalized_output = normalize_ptmpsi_output_pdb(ptmpsi_pdb, ptmpsi_chain, ptmpsi_pos, restype)

            # Load output and measure
            out_structure = load_structure(ptmpsi_pdb, "ptmpsi")
            out_model = list(out_structure.get_models())[0]
            out_res = find_residue(out_model, ptmpsi_chain, ptmpsi_pos)
            out_phos = get_phosphate_coords(out_res)

            if "P" not in out_phos:
                raise RuntimeError("PTM-Psi output missing phosphate atoms")

            rmsd = phosphate_sym_rmsd(ref_phos, out_phos)
            pdist = p_only_dist(ref_phos, out_phos)
            ptmpsi_dihedral = measure_phosphate_dihedral(out_res, restype)

            n_ok += 1
            rows.append({
                "status": "OK",
                "message": "",
                "site_id": site_id,
                "acc_id": acc_id,
                "pdb_id": pdb_id,
                "restype_3": restype,
                "chain": chain,
                "position": pos,
                "reference_path": ref_path,
                "stripped_input": stripped_pdb,
                "stripped_input_sanitized": stripped_sanitized,
                "stripped_input_renumbered": stripped_renumbered,
                "ptmpsi_output": ptmpsi_pdb,
                "ptmpsi_output_normalized": normalized_output,
                "ptmpsi_chain": ptmpsi_chain,
                "ptmpsi_residue_index": ptmpsi_pos,
                "ptmpsi_phosphate_sym_rmsd": round(rmsd, 4) if rmsd is not None else np.nan,
                "ptmpsi_p_only_dist": round(pdist, 4) if pdist is not None else np.nan,
                "crystal_dihedral_deg": round(ref_dihedral, 1) if ref_dihedral is not None else np.nan,
                "ptmpsi_dihedral_deg": round(ptmpsi_dihedral, 1) if ptmpsi_dihedral is not None else np.nan,
                "stage0_kabsch_rmsd": row.get("stage0_kabsch_rmsd", np.nan),
                "stage3_post_minimization_rmsd": row.get("stage3_post_minimization_rmsd", np.nan),
                "top3_best_rmsd": row.get("top3_best_rmsd", np.nan),
                "stage2_selected_angle": row.get("stage2_selected_angle", np.nan),
                "stage1_best_of_12_angle": row.get("stage1_best_of_12_angle", np.nan),
                "source_tsv": row.get("_source_tsv", ""),
            })
            write_detail_checkpoint(rows, args.detail_output)
            write_timeout_site_list(rows, args.timeout_sites_output)
            print(f"  OK {site_id}: PTM-Psi={rmsd:.3f}Å dih={ptmpsi_dihedral:.0f}° | crystal={ref_dihedral:.0f}°" if all(v is not None for v in [rmsd, ptmpsi_dihedral, ref_dihedral]) else f"  OK {site_id}: PTM-Psi={rmsd:.3f}Å" if rmsd else f"  OK {site_id}")

        except Exception as exc:
            n_fail += 1
            rows.append({
                "status": "FAILED",
                "message": str(exc),
                "site_id": site_id,
                "acc_id": acc_id,
                "pdb_id": pdb_id,
                "restype_3": restype,
                "chain": chain,
                "position": pos,
                "reference_path": ref_path if ref_path else "",
                "source_tsv": row.get("_source_tsv", ""),
            })
            write_detail_checkpoint(rows, args.detail_output)
            write_timeout_site_list(rows, args.timeout_sites_output)
            print(f"  FAIL {site_id}: {exc}", flush=True)

    print(f"\nPTM-Psi complete: {n_ok} OK, {n_fail} failed")
    detail = pd.DataFrame(rows)
    detail.to_csv(args.detail_output, sep="\t", index=False)
    write_timeout_site_list(rows, args.timeout_sites_output)
    return detail


def summarize_results(pf_df: pd.DataFrame, ptmpsi_detail: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Build combined summary table."""
    rows = []

    # PhosphoFill and naive from result TSVs
    work = pf_df.copy()
    if "context_n_sites" in work.columns:
        work = work[work["context_n_sites"] == 1].copy()

    for restype in ["TPO", "SEP", "PTR", "Combined"]:
        sub = work if restype == "Combined" else work[work["restype_3"] == restype]
        if len(sub) == 0:
            continue

        for method, col in [("PhosphoFill Stage 0 seed", "stage0_kabsch_rmsd"),
                            ("PhosphoFill top1/S3", "stage3_post_minimization_rmsd"),
                            ("PhosphoFill best-of-3", "top3_best_rmsd")]:
            vals = pd.to_numeric(sub[col], errors="coerce").dropna()
            if len(vals) == 0:
                continue
            rows.append({
                "Residue type": restype,
                "Method": method,
                "n": len(vals),
                "median_RMSD_A": round(vals.median(), 3),
                "mean_RMSD_A": round(vals.mean(), 3),
                "le1_pct": round((vals <= 1.0).mean() * 100, 1),
                "le1p5_pct": round((vals <= 1.5).mean() * 100, 1),
            })

    # PTM-Psi
    if ptmpsi_detail is not None:
        ok = ptmpsi_detail[ptmpsi_detail["status"] == "OK"]
        for restype in ["TPO", "SEP", "PTR", "Combined"]:
            sub = ok if restype == "Combined" else ok[ok["restype_3"] == restype]
            if len(sub) == 0:
                continue
            vals = pd.to_numeric(sub["ptmpsi_phosphate_sym_rmsd"], errors="coerce").dropna()
            if len(vals) == 0:
                continue
            rows.append({
                "Residue type": restype,
                "Method": "PTM-Psi (NERF)",
                "n": len(vals),
                "median_RMSD_A": round(vals.median(), 3),
                "mean_RMSD_A": round(vals.mean(), 3),
                "le1_pct": round((vals <= 1.0).mean() * 100, 1),
                "le1p5_pct": round((vals <= 1.5).mean() * 100, 1),
            })

    return pd.DataFrame(rows)


def print_dihedral_summary(detail: pd.DataFrame) -> None:
    """Print dihedral comparison summary."""
    ok = detail[detail["status"] == "OK"]
    if "crystal_dihedral_deg" not in ok.columns:
        return

    print(f"\n{'='*72}")
    print("DIHEDRAL COMPARISON: PTM-Psi vs Crystal")
    print(f"{'='*72}")

    for rt in ["TPO", "SEP", "PTR"]:
        sub = ok[ok["restype_3"] == rt]
        if len(sub) == 0:
            continue
        n = len(sub)
        crystal = pd.to_numeric(sub["crystal_dihedral_deg"], errors="coerce").dropna()
        ptmpsi = pd.to_numeric(sub["ptmpsi_dihedral_deg"], errors="coerce").dropna()
        pf_angle = pd.to_numeric(sub["stage2_selected_angle"], errors="coerce").dropna()

        print(f"\n--- {rt} (n={n}) ---")
        if len(crystal) > 0:
            print(f"  Crystal dihedral: median={crystal.median():.0f}°")
        if len(ptmpsi) > 0:
            print(f"  PTM-Psi dihedral: median={ptmpsi.median():.0f}°")
        if len(pf_angle) > 0:
            print(f"  PF selected angle: median={pf_angle.median():.0f}° (relative rotation)")


def parse_args():
    ap = argparse.ArgumentParser(description="PhosphoFill vs PTM-Psi benchmark")
    ap.add_argument("--phosphofill-tsv", nargs="+", required=True,
                    help="PhosphoFill benchmark result TSVs")
    ap.add_argument("--pdb-cache-dir", default="pdb_cache",
                    help="PDB download cache directory")
    ap.add_argument("--ptmpsi-source-dir", default=None,
                    help="Optional PTMPSI checkout/source directory containing ptmpsi/")
    ap.add_argument("--ptmpsi-workdir", default="ptmpsi_work",
                    help="Working directory for stripped/output PDBs (not the PTMPSI source dir)")
    ap.add_argument("--ptmpsi-timeout-sec", type=float, default=300.0,
                    help="Hard timeout per PTM-Psi site in seconds (<=0 disables timeout)")
    ap.add_argument("--output", default="ptmpsi_comparison.tsv",
                    help="Summary output TSV")
    ap.add_argument("--detail-output", default="ptmpsi_detail.tsv",
                    help="Per-site detail TSV")
    ap.add_argument("--timeout-sites-output", default="ptmpsi_timeout_sites.txt",
                    help="Write timed-out site IDs here as they occur")
    ap.add_argument("--site-id-list", default=None,
                    help="Only run site IDs from this text/CSV/TSV file")
    ap.add_argument("--rerun-from-detail", default=None,
                    help="Rerun sites selected from an existing detail TSV")
    ap.add_argument("--rerun-mode", choices=["timeouts", "failed", "non-ok"], default="timeouts",
                    help="Which rows to rerun from --rerun-from-detail")
    ap.add_argument("--max-sites-per-restype", type=int, default=None,
                    help="Max sites per residue type (None = all)")
    ap.add_argument("--skip-ptmpsi", action="store_true",
                    help="Skip PTM-Psi (just summarize PhosphoFill)")
    return ap.parse_args()


def main():
    args = parse_args()

    print("=" * 72)
    print("PhosphoFill vs PTM-Psi benchmark")
    print("=" * 72)

    if not args.skip_ptmpsi:
        load_ptmpsi(args.ptmpsi_source_dir)

    if not PTMPSI_AVAILABLE and not args.skip_ptmpsi:
        print("WARNING: PTM-Psi could not be imported. Install/check with:")
        print("  git clone https://github.com/pnnl/PTMPSI.git")
        print("  cd PTMPSI && pip install -e .")
        if PTMPSI_IMPORT_ERROR:
            print(f"\nImport error: {PTMPSI_IMPORT_ERROR}")
            last_line = PTMPSI_IMPORT_TRACEBACK.strip().splitlines()[-1] if PTMPSI_IMPORT_TRACEBACK else ""
            if last_line:
                print(f"Last traceback line: {last_line}")
        print("\nIf PTMPSI is checked out beside this script, use:")
        print("  --ptmpsi-source-dir PTMPSI --ptmpsi-workdir ptmpsi_work")
        print("Running in summary-only mode.\n")
        args.skip_ptmpsi = True
    elif PTMPSI_AVAILABLE:
        print(f"PTM-Psi loaded via {PTMPSI_IMPORT_MODE}")

    # Load PhosphoFill results
    pf = load_phosphofill_results(args.phosphofill_tsv)
    print(f"Loaded {len(pf)} OK PhosphoFill entries")

    ptmpsi_detail = None
    if not args.skip_ptmpsi:
        ptmpsi_detail = run_ptmpsi_benchmark(pf, args)
        print_dihedral_summary(ptmpsi_detail)

    # Summary
    summary = summarize_results(pf, ptmpsi_detail)
    summary.to_csv(args.output, sep="\t", index=False)

    print(f"\n{'='*72}")
    print("SUMMARY")
    print(f"{'='*72}")
    print(summary.to_string(index=False))
    print(f"\nSummary: {args.output}")
    if ptmpsi_detail is not None:
        print(f"Detail: {args.detail_output}")


if __name__ == "__main__":
    main()
