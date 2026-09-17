#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import itertools
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Superimposer
from Bio.PDB.mmcifio import MMCIFIO

PHOSPHO_MAP = {"SEP": "SER", "TPO": "THR", "PTR": "TYR"}
PHOSPHO_ATOMS = {"P", "O1P", "O2P", "O3P"}
PHOSPHO_OXYGENS = ("O1P", "O2P", "O3P")
PHOSPHO_HYDROGEN = {"H1P", "H2P", "H3P"}
BACKBONE_ATOMS = {"N", "CA", "C", "O"}

METAL_RESNAMES = {
    "MG", "MN", "ZN", "FE", "CA", "CO", "NI", "CU", "CD", "NA", "K",
    "FE2", "MN3", "CO2", "NI2", "CU1", "CU2", "ZN2",
}
WATER_RESNAMES = {"HOH", "WAT", "H2O", "DOD"}
BASIC_DONOR_ATOMS = {
    ("LYS", "NZ"),
    ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
    ("HIS", "ND1"), ("HIS", "NE2"), ("HID", "ND1"), ("HIE", "NE2"),
}
HBOND_DONOR_ATOMS = BASIC_DONOR_ATOMS | {
    ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"),
    ("ASN", "ND2"), ("GLN", "NE2"), ("TRP", "NE1"),
}
PDB_DOWNLOAD_URL = "https://files.rcsb.org/download/{}.cif"


def download_file(url: str, dest: str) -> str:
    if os.path.exists(dest):
        return dest
    os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
    print(f"  Downloading {url} ...", end=" ", flush=True)
    urllib.request.urlretrieve(url, dest)
    print("OK")
    return dest


def download_pdb(pdb_id: str, cache_dir: str) -> str:
    return download_file(PDB_DOWNLOAD_URL.format(pdb_id.lower()), os.path.join(cache_dir, f"{pdb_id.lower()}.cif"))


def parse_structure(path: str, name: str = "s"):
    suffix = Path(path).suffix.lower()
    parser = MMCIFParser(QUIET=True) if suffix in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    return parser.get_structure(name, path)


def get_residue(model, chain_id: str, resseq: int):
    if chain_id not in model:
        return None
    chain = model[chain_id]
    for hetflag in [" ", "H_SEP", "H_TPO", "H_PTR"]:
        rid = (hetflag, resseq, " ")
        if rid in chain:
            return chain[rid]
    for residue in chain.get_residues():
        if residue.id[1] == resseq:
            return residue
    return None


def strip_phosphate_from_structure(structure, chain_id: str, positions: List[int]) -> None:
    model = list(structure.get_models())[0]
    for resseq in positions:
        residue = get_residue(model, chain_id, resseq)
        if residue is None:
            continue
        resname = residue.get_resname().strip().upper()
        if resname not in PHOSPHO_MAP:
            continue
        atoms_to_remove = []
        for atom in residue.get_atoms():
            aname = atom.get_name().strip()
            if aname in PHOSPHO_ATOMS or aname in PHOSPHO_HYDROGEN:
                atoms_to_remove.append(atom.get_id())
        for aid in atoms_to_remove:
            residue.detach_child(aid)
        residue.resname = PHOSPHO_MAP[resname]
        # Fix hetfield: convert from HETATM (H_SEP/H_TPO/H_PTR) to standard
        # residue (' ') so the residue is found correctly after CIF round-trip.
        # IMPORTANT: we update the chain's internal dict key directly and set
        # residue._id to preserve residue ordering.  chain.detach_child +
        # chain.add would move the residue to the end of the chain.
        old_id = residue.id
        if old_id[0] != " ":
            new_id = (" ", old_id[1], old_id[2])
            chain = residue.get_parent()
            del chain.child_dict[old_id]
            chain.child_dict[new_id] = residue
            residue._id = new_id


def save_structure_file(structure, path: str) -> None:
    """Save structure in CIF or PDB format based on file extension.
    Defaults to CIF to avoid PDB format limitations (long ligand names,
    long chain IDs, >99999 atoms)."""
    suffix = Path(path).suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        io = MMCIFIO()
        io.set_structure(structure)
        io.save(path)
    else:
        io = PDBIO()
        io.set_structure(structure)
        io.save(path)
        standard_aa = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
            "TYR", "VAL",
        }
        lines = []
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                if line.startswith("HETATM") and line[17:20].strip() in standard_aa:
                    line = "ATOM  " + line[6:]
                lines.append(line)
        with open(path, "w", encoding="utf-8") as fh:
            fh.writelines(lines)

# Legacy alias
save_structure_pdb = save_structure_file


def get_phosphate_coords(residue) -> Optional[Dict[str, np.ndarray]]:
    coords = {}
    for atom in residue.get_atoms():
        aname = atom.get_name().strip()
        if aname in PHOSPHO_ATOMS:
            coords[aname] = atom.get_vector().get_array()
    return coords if "P" in coords else None


def phosphate_rmsd_named(ref_coords: Dict[str, np.ndarray], pred_coords: Dict[str, np.ndarray]) -> Optional[float]:
    shared = sorted(set(ref_coords) & set(pred_coords))
    if not shared:
        return None
    ref = np.array([ref_coords[a] for a in shared], dtype=float)
    pred = np.array([pred_coords[a] for a in shared], dtype=float)
    return float(np.sqrt(np.mean(np.sum((ref - pred) ** 2, axis=1))))


def _rmsd_from_arrays(ref: np.ndarray, pred: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((ref - pred) ** 2, axis=1))))


def phosphate_rmsd(ref_coords: Dict[str, np.ndarray], pred_coords: Dict[str, np.ndarray]) -> Optional[float]:
    if "P" not in ref_coords or "P" not in pred_coords:
        return phosphate_rmsd_named(ref_coords, pred_coords)
    ref_ox = [name for name in PHOSPHO_OXYGENS if name in ref_coords]
    pred_ox = [name for name in PHOSPHO_OXYGENS if name in pred_coords]
    if len(ref_ox) != 3 or len(pred_ox) != 3:
        return phosphate_rmsd_named(ref_coords, pred_coords)
    ref_p = np.array(ref_coords["P"], dtype=float)
    pred_p = np.array(pred_coords["P"], dtype=float)
    ref_ox_coords = [np.array(ref_coords[name], dtype=float) for name in PHOSPHO_OXYGENS]
    pred_ox_coords = [np.array(pred_coords[name], dtype=float) for name in PHOSPHO_OXYGENS]
    best = None
    for order in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
        ref_arr = np.array([ref_p] + ref_ox_coords, dtype=float)
        pred_arr = np.array([pred_p] + [pred_ox_coords[i] for i in order], dtype=float)
        rmsd = _rmsd_from_arrays(ref_arr, pred_arr)
        if best is None or rmsd < best:
            best = rmsd
    return best


def dihedral(p0, p1, p2, p3) -> Optional[float]:
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    p3 = np.asarray(p3, dtype=float)
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    n1 = np.cross(b0, b1)
    n2 = np.cross(b1, b2)
    n1n = np.linalg.norm(n1)
    n2n = np.linalg.norm(n2)
    b1n = np.linalg.norm(b1)
    if n1n < 1e-8 or n2n < 1e-8 or b1n < 1e-8:
        return None
    n1 /= n1n
    n2 /= n2n
    b1u = b1 / b1n
    m1 = np.cross(n1, b1u)
    x = float(np.dot(n1, n2))
    y = float(np.dot(m1, n2))
    return float(np.degrees(np.arctan2(y, x)))


def _canonical_residue_atom_coords(residue) -> Dict[str, np.ndarray]:
    atoms = {}
    for atom in residue.get_atoms():
        real_atom = atom.selected_child if atom.is_disordered() else atom
        atoms[real_atom.get_name().strip()] = np.asarray(real_atom.get_coord(), dtype=float)
    return atoms


def _bridging_anchor_names(resname: str) -> Optional[Tuple[str, str, str]]:
    resname = resname.strip().upper()
    if resname == "SEP":
        return ("CB", "OG", "P")
    if resname == "TPO":
        return ("CB", "OG1", "P")
    if resname == "PTR":
        return ("CZ", "OH", "P")
    return None


def best_permuted_torsion(residue) -> Optional[float]:
    if residue is None:
        return None
    atoms = _canonical_residue_atom_coords(residue)
    anchor = _bridging_anchor_names(residue.get_resname())
    if anchor is None:
        return None
    vals = []
    for oxy in PHOSPHO_OXYGENS:
        names = list(anchor) + [oxy]
        if all(name in atoms for name in names):
            angle = dihedral(atoms[names[0]], atoms[names[1]], atoms[names[2]], atoms[names[3]])
            if angle is not None:
                vals.append(angle)
    return vals[0] if vals else None


def torsion_error_label_invariant(ref_residue, pred_residue) -> Optional[float]:
    if ref_residue is None or pred_residue is None:
        return None

    def torsion_list(residue):
        resname = residue.get_resname().strip().upper()
        atoms = {a.get_name().strip(): a.get_vector().get_array() for a in residue.get_atoms()}
        if resname == "SEP":
            anchor = ["CB", "OG", "P"]
        elif resname == "TPO":
            anchor = ["CB", "OG1", "P"]
        elif resname == "PTR":
            anchor = ["CZ", "OH", "P"]
        else:
            return []
        vals = []
        for oxy in PHOSPHO_OXYGENS:
            if all(name in atoms for name in anchor + [oxy]):
                angle = dihedral(atoms[anchor[0]], atoms[anchor[1]], atoms[anchor[2]], atoms[oxy])
                if angle is not None:
                    vals.append(angle)
        return vals

    ref_vals = torsion_list(ref_residue)
    pred_vals = torsion_list(pred_residue)
    if not ref_vals or not pred_vals:
        return None
    best = None
    for a in ref_vals:
        for b in pred_vals:
            diff = abs(a - b)
            if diff > 180:
                diff = 360 - diff
            if best is None or diff < best:
                best = diff
    return best


def p_atom_displacement(ref_coords: Dict[str, np.ndarray], pred_coords: Dict[str, np.ndarray]) -> Optional[float]:
    if "P" not in ref_coords or "P" not in pred_coords:
        return None
    return float(np.linalg.norm(ref_coords["P"] - pred_coords["P"]))


def annotate_phosphosite_environment(model, chain_id: str, resseq: int,
                                     water_radius: float = 4.0, metal_radius: float = 5.0) -> Dict[str, object]:
    residue = get_residue(model, chain_id, resseq)
    if residue is None:
        return {"n_waters_4A": None, "metal_within_5A": None, "metal_names": "",
                "nearest_metal_dist": None, "nearest_water_dist": None,
                "site_secondary_structure": "unknown", "site_category": "unknown"}
    phos_coords = []
    for atom in residue.get_atoms():
        if atom.get_name().strip() in PHOSPHO_ATOMS:
            phos_coords.append(atom.get_vector().get_array())
    if not phos_coords:
        return {"n_waters_4A": None, "metal_within_5A": None, "metal_names": "",
                "nearest_metal_dist": None, "nearest_water_dist": None,
                "site_secondary_structure": "unknown", "site_category": "unknown"}
    phos_coords = np.array(phos_coords)

    n_waters = 0
    nearest_water_dist = None
    for chain in model.get_chains():
        for res in chain.get_residues():
            if res.get_resname().strip().upper() not in WATER_RESNAMES:
                continue
            counted = False
            for atom in res.get_atoms():
                wat_coord = atom.get_vector().get_array()
                for pc in phos_coords:
                    d = float(np.linalg.norm(wat_coord - pc))
                    if d <= water_radius and not counted:
                        n_waters += 1
                        counted = True
                    if nearest_water_dist is None or d < nearest_water_dist:
                        nearest_water_dist = d

    metal_within = False
    nearest_metal_dist = None
    metal_names_found = []
    for chain in model.get_chains():
        for res in chain.get_residues():
            resname = res.get_resname().strip().upper()
            if resname not in METAL_RESNAMES:
                continue
            for atom in res.get_atoms():
                metal_coord = atom.get_vector().get_array()
                for pc in phos_coords:
                    d = float(np.linalg.norm(metal_coord - pc))
                    if d <= metal_radius:
                        metal_within = True
                        metal_names_found.append(f"{resname}:{d:.1f}A")
                    if nearest_metal_dist is None or d < nearest_metal_dist:
                        nearest_metal_dist = d

    ss = _guess_secondary_structure(model, chain_id, resseq)
    category = _classify_site_category(n_waters, metal_within)
    return {
        "n_waters_4A": n_waters,
        "metal_within_5A": metal_within,
        "metal_names": ";".join(metal_names_found) if metal_names_found else "",
        "nearest_metal_dist": round(nearest_metal_dist, 2) if nearest_metal_dist is not None else None,
        "nearest_water_dist": round(nearest_water_dist, 2) if nearest_water_dist is not None else None,
        "site_secondary_structure": ss,
        "site_category": category,
    }


def _guess_secondary_structure(model, chain_id: str, resseq: int) -> str:
    if chain_id not in model:
        return "unknown"
    chain = model[chain_id]
    ca_coords = {}
    for res in chain.get_residues():
        rnum = res.id[1]
        if abs(rnum - resseq) <= 5:
            for atom in res.get_atoms():
                if atom.get_name().strip() == "CA":
                    ca_coords[rnum] = atom.get_vector().get_array()
                    break
    if resseq in ca_coords and (resseq + 4) in ca_coords:
        d_i4 = np.linalg.norm(ca_coords[resseq] - ca_coords[resseq + 4])
        if 5.0 < d_i4 < 7.0:
            return "helix"
    if resseq in ca_coords and (resseq + 2) in ca_coords:
        d_i2 = np.linalg.norm(ca_coords[resseq] - ca_coords[resseq + 2])
        if d_i2 > 6.5:
            return "strand"
    return "loop"


def _classify_site_category(n_waters: int, metal_within: bool) -> str:
    if metal_within:
        return "B_metal"
    if n_waters >= 3:
        return "B_water"
    return "A_protein"


def reference_pose_qc(model, chain_id: str, resseq: int) -> Dict[str, object]:
    residue = get_residue(model, chain_id, resseq)
    if residue is None:
        return {
            "ref_phosphate_b_mean": None,
            "ref_backbone_b_mean": None,
            "ref_b_ratio": None,
            "ref_phosphate_occupancy_mean": None,
            "ref_n_basic_contacts_4A": None,
            "ref_n_hbondlike_contacts_35A": None,
            "ref_nearest_basic_dist": None,
            "ref_nearest_any_contact": None,
            "ref_po_bond_mean": None,
            "ref_po_bond_std": None,
            "ref_pose_support_class": "unknown",
        }

    atoms = {}
    phosphate_b = []
    phosphate_occ = []
    backbone_b = []
    phosphate_coords = []
    p_coord = None
    for atom in residue.get_atoms():
        real = atom.selected_child if atom.is_disordered() else atom
        name = real.get_name().strip()
        atoms[name] = np.asarray(real.get_coord(), dtype=float)
        if name in PHOSPHO_ATOMS:
            phosphate_b.append(float(real.get_bfactor()))
            phosphate_occ.append(float(real.get_occupancy() if real.get_occupancy() is not None else 1.0))
            phosphate_coords.append(np.asarray(real.get_coord(), dtype=float))
            if name == "P":
                p_coord = np.asarray(real.get_coord(), dtype=float)
        if name in BACKBONE_ATOMS:
            backbone_b.append(float(real.get_bfactor()))

    if p_coord is None or not phosphate_coords:
        return {
            "ref_phosphate_b_mean": None,
            "ref_backbone_b_mean": None,
            "ref_b_ratio": None,
            "ref_phosphate_occupancy_mean": None,
            "ref_n_basic_contacts_4A": None,
            "ref_n_hbondlike_contacts_35A": None,
            "ref_nearest_basic_dist": None,
            "ref_nearest_any_contact": None,
            "ref_po_bond_mean": None,
            "ref_po_bond_std": None,
            "ref_pose_support_class": "unknown",
        }

    po_bonds = []
    for oxy in PHOSPHO_OXYGENS:
        if oxy in atoms:
            po_bonds.append(float(np.linalg.norm(atoms[oxy] - p_coord)))

    n_basic = 0
    n_hbondlike = 0
    nearest_basic = None
    nearest_any = None
    same_res_id = residue.get_full_id()
    for chain in model.get_chains():
        for other_res in chain.get_residues():
            if other_res.get_full_id() == same_res_id:
                continue
            resname = other_res.get_resname().strip().upper()
            for atom in other_res.get_atoms():
                real = atom.selected_child if atom.is_disordered() else atom
                aname = real.get_name().strip()
                coord = np.asarray(real.get_coord(), dtype=float)
                for pc in phosphate_coords:
                    d = float(np.linalg.norm(coord - pc))
                    if nearest_any is None or d < nearest_any:
                        nearest_any = d
                    if (resname, aname) in BASIC_DONOR_ATOMS and d <= 4.0:
                        n_basic += 1
                        if nearest_basic is None or d < nearest_basic:
                            nearest_basic = d
                    if (resname, aname) in HBOND_DONOR_ATOMS and d <= 3.5:
                        n_hbondlike += 1

    ref_phosphate_b_mean = float(np.mean(phosphate_b)) if phosphate_b else None
    ref_backbone_b_mean = float(np.mean(backbone_b)) if backbone_b else None
    ref_b_ratio = (ref_phosphate_b_mean / ref_backbone_b_mean) if (ref_phosphate_b_mean is not None and ref_backbone_b_mean not in (None, 0.0)) else None
    ref_occ_mean = float(np.mean(phosphate_occ)) if phosphate_occ else None
    po_mean = float(np.mean(po_bonds)) if po_bonds else None
    po_std = float(np.std(po_bonds)) if po_bonds else None

    support = "low"
    if (ref_b_ratio is not None and ref_b_ratio <= 1.5 and
            (n_basic >= 1 or n_hbondlike >= 2) and
            ref_occ_mean is not None and ref_occ_mean >= 0.95):
        support = "high"
    elif (ref_b_ratio is not None and ref_b_ratio <= 2.0 and
          ref_occ_mean is not None and ref_occ_mean >= 0.8):
        support = "medium"

    return {
        "ref_phosphate_b_mean": round(ref_phosphate_b_mean, 2) if ref_phosphate_b_mean is not None else None,
        "ref_backbone_b_mean": round(ref_backbone_b_mean, 2) if ref_backbone_b_mean is not None else None,
        "ref_b_ratio": round(ref_b_ratio, 3) if ref_b_ratio is not None else None,
        "ref_phosphate_occupancy_mean": round(ref_occ_mean, 3) if ref_occ_mean is not None else None,
        "ref_n_basic_contacts_4A": n_basic,
        "ref_n_hbondlike_contacts_35A": n_hbondlike,
        "ref_nearest_basic_dist": round(nearest_basic, 2) if nearest_basic is not None else None,
        "ref_nearest_any_contact": round(nearest_any, 2) if nearest_any is not None else None,
        "ref_po_bond_mean": round(po_mean, 3) if po_mean is not None else None,
        "ref_po_bond_std": round(po_std, 3) if po_std is not None else None,
        "ref_pose_support_class": support,
    }


def three_stage_rmsd(frames_data: Dict, ref_coords: Dict[str, np.ndarray]) -> Dict[str, object]:
    result = {
        "stage1_best_of_12_rmsd": None,
        "stage1_best_of_12_angle": None,
        "stage2_selected_rmsd": None,
        "stage2_selected_angle": None,
        "stage2_scoring_gap": None,
    }
    if isinstance(frames_data, dict) and "frames" in frames_data:
        frames = frames_data["frames"]
    elif isinstance(frames_data, list):
        frames = frames_data
    else:
        return result
    if not frames or "P" not in ref_coords:
        return result

    frame_rmsds = []
    best_score_idx = 0
    best_score = None
    selected_rank1_idx = None
    for i, frame in enumerate(frames):
        frame_coords = {}
        if "coords" in frame:
            for aname, xyz in frame["coords"].items():
                if aname in PHOSPHO_ATOMS:
                    frame_coords[aname] = np.array(xyz)
        if not frame_coords or "P" not in frame_coords:
            frame_rmsds.append(None)
            continue
        rmsd = phosphate_rmsd(ref_coords, frame_coords)
        frame_rmsds.append(rmsd)
        total_score = frame.get("total_score")
        if total_score is not None and (best_score is None or total_score < best_score):
            best_score = total_score
            best_score_idx = i
        if frame.get("selection_rank") == 1:
            selected_rank1_idx = i

    if selected_rank1_idx is not None:
        best_score_idx = selected_rank1_idx

    valid_rmsds = [(r, i) for i, r in enumerate(frame_rmsds) if r is not None]

    # Annotate each frame with ground-truth labels for supervised learning.
    # Mapping is airtight: frame_id + angle ensure unambiguous identification.
    best_frame_idx = min(valid_rmsds, key=lambda x: x[0])[1] if valid_rmsds else None
    for i, frame in enumerate(frames):
        rmsd = frame_rmsds[i]
        frame["frame_id"] = i
        frame["rmsd_to_ref"] = round(rmsd, 4) if rmsd is not None else None
        frame["is_best_rmsd"] = 1 if i == best_frame_idx else 0
        frame["is_good_1A"] = 1 if rmsd is not None and rmsd <= 1.0 else 0

    if valid_rmsds:
        best_rmsd, best_idx = min(valid_rmsds, key=lambda x: x[0])
        result["stage1_best_of_12_rmsd"] = round(best_rmsd, 4)
        result["stage1_best_of_12_angle"] = frames[best_idx].get("angle")
    if best_score_idx < len(frame_rmsds) and frame_rmsds[best_score_idx] is not None:
        selected_rmsd = frame_rmsds[best_score_idx]
        result["stage2_selected_rmsd"] = round(selected_rmsd, 4)
        result["stage2_selected_angle"] = frames[best_score_idx].get("angle")
        if result["stage1_best_of_12_rmsd"] is not None:
            result["stage2_scoring_gap"] = round(selected_rmsd - result["stage1_best_of_12_rmsd"], 4)
    return result


def benchmark_success_label(phosphate_rmsd_value: Optional[float], p_displacement_value: Optional[float],
                            clashes_after: Optional[float]) -> str:
    if phosphate_rmsd_value is None or p_displacement_value is None:
        return "unknown"
    clash_val = None if clashes_after in (None, "") else float(clashes_after)
    if phosphate_rmsd_value <= 1.5 and p_displacement_value <= 1.0 and (clash_val is None or clash_val <= 1):
        return "pass"
    if phosphate_rmsd_value <= 2.5 and p_displacement_value <= 2.0 and (clash_val is None or clash_val <= 2):
        return "acceptable"
    return "fail"


def _parse_optional_number(value):
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def pose_output_path(base_path: str, rank: int, n_poses: int) -> str:
    if n_poses <= 1:
        return base_path
    path = Path(base_path)
    return str(path.parent / f"{path.stem}_phosphoFill_{rank}{path.suffix}")


def run_phosphofill(graft_script: str, input_path: str, output_path: str, sites: List[str], cache_dir: str,
                    report_script: Optional[str] = None, extra_args: Optional[List[str]] = None,
                    graft_mode: str = "joint", n_poses: int = 1
                    ) -> Tuple[bool, str, Dict[str, Any]]:
    actual_n_poses = max(1, min(int(n_poses), 3))
    cmd = [sys.executable, graft_script, input_path, output_path, "--sites", *sites,
           "--cache-dir", cache_dir, "--n-poses", str(actual_n_poses)]
    if extra_args:
        cmd.extend(extra_args)
    pose_paths = [pose_output_path(output_path, rank, actual_n_poses)
                  for rank in range(1, actual_n_poses + 1)]
    graft_report_paths = [
        (output_path + ".report.tsv") if actual_n_poses == 1 else str(Path(path).with_suffix(".report.tsv"))
        for path in pose_paths
    ]
    empty_paths: Dict[str, Any] = {
        "graft_report": None, "site_report": None, "site_report_json": None,
        "prescan_frames_json": None, "output_structures": [], "graft_reports": [],
        "site_reports": [], "site_report_jsons": [],
    }
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if result.returncode != 0:
            return False, f"Exit code {result.returncode}: {(result.stderr or result.stdout)[:1000]}", empty_paths
        missing = [path for path in pose_paths + graft_report_paths if not os.path.exists(path)]
        if missing:
            return False, f"Expected multipose output not generated: {missing[0]}", empty_paths

        report_tsv_paths: List[Optional[str]] = []
        report_json_paths: List[Optional[str]] = []
        if report_script:
            for pose_path, graft_report_path in zip(pose_paths, graft_report_paths):
                report_prefix = pose_path + ".site_metrics"
                report_cmd = [sys.executable, report_script, input_path, pose_path, "--sites", *sites,
                              "--output-prefix", report_prefix, "--graft-report", graft_report_path,
                              "--graft-mode", graft_mode]
                report_result = subprocess.run(report_cmd, capture_output=True, text=True, timeout=600)
                candidate = report_prefix + ".tsv"
                report_tsv = candidate if report_result.returncode == 0 and os.path.exists(candidate) else None
                report_tsv_paths.append(report_tsv)
                report_json_paths.append(locate_sidecar_json(report_tsv))
        else:
            report_tsv_paths = [None] * actual_n_poses
            report_json_paths = [None] * actual_n_poses

        base_graft_report = output_path + ".report.tsv"
        return True, "OK", {
            "graft_report": graft_report_paths[0],
            "site_report": report_tsv_paths[0],
            "site_report_json": report_json_paths[0],
            "prescan_frames_json": locate_prescan_frames_path(base_graft_report),
            "output_structures": pose_paths,
            "graft_reports": graft_report_paths,
            "site_reports": report_tsv_paths,
            "site_report_jsons": report_json_paths,
        }
    except subprocess.TimeoutExpired:
        return False, "Timeout after 1800s", empty_paths
    except Exception as e:
        return False, str(e), empty_paths


def load_prescan_frames(report_tsv_path: str) -> Dict[str, Dict]:
    if report_tsv_path and str(report_tsv_path).endswith(".json") and os.path.exists(report_tsv_path):
        with open(report_tsv_path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
        return data if isinstance(data, dict) else {}
    candidates = [report_tsv_path + ".prescan_frames.json"]
    if report_tsv_path.endswith(".tsv"):
        candidates.append(report_tsv_path[:-4] + ".prescan_frames.json")
    try:
        candidates.append(str(Path(report_tsv_path).with_suffix(".prescan_frames.json")))
    except Exception:
        pass
    for path in candidates:
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
            if data:
                return data
    return {}


def locate_sidecar_json(tsv_path: Optional[str]) -> Optional[str]:
    if not tsv_path:
        return None
    candidates = []
    if tsv_path.endswith('.tsv'):
        candidates.append(tsv_path[:-4] + '.json')
    try:
        candidates.append(str(Path(tsv_path).with_suffix('.json')))
    except Exception:
        pass
    for path in candidates:
        if path and os.path.exists(path):
            return path
    return None


def locate_prescan_frames_path(report_tsv_path: Optional[str]) -> Optional[str]:
    if not report_tsv_path:
        return None
    candidates = [report_tsv_path + '.prescan_frames.json']
    if report_tsv_path.endswith('.tsv'):
        candidates.append(report_tsv_path[:-4] + '.prescan_frames.json')
    try:
        candidates.append(str(Path(report_tsv_path).with_suffix('.prescan_frames.json')))
    except Exception:
        pass
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def _json_default(obj: Any):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(f'Object of type {type(obj).__name__} is not JSON serializable')


def write_json(path: str, payload: Any) -> None:
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w', encoding='utf-8') as fh:
        json.dump(payload, fh, indent=2, default=_json_default)


def copy_if_exists(src: Optional[str], dst: str) -> Optional[str]:
    if not src or not os.path.exists(src):
        return None
    os.makedirs(os.path.dirname(dst) or '.', exist_ok=True)
    shutil.copy2(src, dst)
    return dst


def _round_coord_list(coord: Sequence[float], ndigits: int = 4) -> List[float]:
    return [round(float(x), ndigits) for x in coord]


def _round_coord_map(coord_map: Optional[Dict[str, Sequence[float]]], ndigits: int = 4) -> Optional[Dict[str, List[float]]]:
    if coord_map is None:
        return None
    return {str(k): _round_coord_list(v, ndigits=ndigits) for k, v in coord_map.items()}


def set_phosphate_coords_on_structure(structure, chain_id: str, resseq: int, coord_map: Dict[str, Sequence[float]]) -> bool:
    model = list(structure.get_models())[0]
    residue = get_residue(model, chain_id, resseq)
    if residue is None:
        return False
    changed = False
    for atom in residue.get_atoms():
        aname = atom.get_name().strip()
        if aname in coord_map:
            atom.coord = np.asarray(coord_map[aname], dtype=float)
            changed = True
    return changed


def write_stage_structure_from_base(base_structure_path: str, output_path: str, chain_id: str, per_site_coords: Dict[int, Dict[str, Sequence[float]]]) -> Optional[str]:
    if not per_site_coords or not os.path.exists(base_structure_path):
        return None
    structure = parse_structure(base_structure_path, 'stage_template')
    changed_any = False
    for pos, coord_map in per_site_coords.items():
        if coord_map:
            changed_any = set_phosphate_coords_on_structure(structure, chain_id, pos, coord_map) or changed_any
    if not changed_any:
        return None
    save_structure_pdb(structure, output_path)
    return output_path


def prescan_stage_details(frames_data: Dict, ref_coords: Dict[str, np.ndarray]) -> Dict[str, Any]:
    details: Dict[str, Any] = {
        'n_frames': 0,
        'stage1_best_index': None,
        'stage1_best_angle': None,
        'stage1_best_rmsd': None,
        'stage1_best_coords': None,
        'stage1_best_frame': None,
        'stage2_selected_index': None,
        'stage2_selected_angle': None,
        'stage2_selected_rmsd': None,
        'stage2_selected_coords': None,
        'stage2_selected_frame': None,
        'stage2_scoring_gap': None,
    }
    if isinstance(frames_data, dict) and 'frames' in frames_data:
        frames = frames_data['frames']
    elif isinstance(frames_data, list):
        frames = frames_data
    else:
        return details
    details['n_frames'] = len(frames)
    if not frames or 'P' not in ref_coords:
        return details

    frame_rmsds: List[Optional[float]] = []
    best_score_idx = 0
    best_score = None
    selected_rank1_idx = None
    best_rmsd_idx = None
    best_rmsd = None

    for i, frame in enumerate(frames):
        frame_coords = {}
        for aname, xyz in (frame.get('coords') or {}).items():
            if aname in PHOSPHO_ATOMS:
                frame_coords[aname] = np.array(xyz, dtype=float)
        if not frame_coords or 'P' not in frame_coords:
            frame_rmsds.append(None)
        else:
            rmsd = phosphate_rmsd(ref_coords, frame_coords)
            frame_rmsds.append(rmsd)
            if rmsd is not None and (best_rmsd is None or rmsd < best_rmsd):
                best_rmsd = rmsd
                best_rmsd_idx = i
        total_score = frame.get('total_score')
        if total_score is not None and (best_score is None or total_score < best_score):
            best_score = total_score
            best_score_idx = i
        if frame.get('selection_rank') == 1:
            selected_rank1_idx = i

    if selected_rank1_idx is not None:
        best_score_idx = selected_rank1_idx

    if best_rmsd_idx is not None:
        best_frame = frames[best_rmsd_idx]
        details['stage1_best_index'] = best_rmsd_idx
        details['stage1_best_angle'] = best_frame.get('angle')
        details['stage1_best_rmsd'] = round(float(best_rmsd), 4) if best_rmsd is not None else None
        details['stage1_best_coords'] = _round_coord_map(best_frame.get('coords'))
        details['stage1_best_frame'] = best_frame

    if best_score_idx < len(frames):
        sel_frame = frames[best_score_idx]
        sel_rmsd = frame_rmsds[best_score_idx]
        details['stage2_selected_index'] = best_score_idx
        details['stage2_selected_angle'] = sel_frame.get('angle')
        details['stage2_selected_rmsd'] = round(float(sel_rmsd), 4) if sel_rmsd is not None else None
        details['stage2_selected_coords'] = _round_coord_map(sel_frame.get('coords'))
        details['stage2_selected_frame'] = sel_frame
        if best_rmsd is not None and sel_rmsd is not None:
            details['stage2_scoring_gap'] = round(float(sel_rmsd - best_rmsd), 4)
    return details



def topk_prescan_combo_details(
    chain_id: str,
    target_positions: Sequence[int],
    prescan_frames: Dict[str, Dict],
    ref_phosphates: Dict[int, Dict[str, np.ndarray]],
    k: int = 3,
) -> Dict[str, Any]:
    """Return the top-k scored prescan frame combinations across target sites.

    For a single site this is just the top-k frames by total_score. For multiple
    sites, frame combinations are ranked by the sum of per-site total_score,
    assuming the independent prescan scores are additive.
    """
    per_site_ranked: Dict[int, List[Dict[str, Any]]] = {}
    for pos in target_positions:
        site_key = f"{chain_id}:{pos}"
        frames_data = prescan_frames.get(site_key, {})
        if isinstance(frames_data, dict) and 'frames' in frames_data:
            frames = frames_data['frames'] or []
        elif isinstance(frames_data, list):
            frames = frames_data
        else:
            frames = []
        ranked = []
        ref_coords = ref_phosphates.get(pos, {})
        for idx, frame in enumerate(frames):
            total_score = frame.get('total_score')
            if total_score is None:
                continue
            coords = frame.get('coords') or {}
            phosphate_coords = {aname: np.array(xyz, dtype=float) for aname, xyz in coords.items() if aname in PHOSPHO_ATOMS}
            rmsd = phosphate_rmsd(ref_coords, phosphate_coords) if phosphate_coords and 'P' in phosphate_coords and ref_coords else None
            ranked.append({
                'index': idx,
                'angle': frame.get('angle'),
                'score': float(total_score),
                'coords': _round_coord_map(coords),
                'frame': frame,
                'rmsd': round(float(rmsd), 4) if rmsd is not None else None,
            })
        ranked.sort(key=lambda x: (x['score'], x['index']))
        per_site_ranked[int(pos)] = ranked

    if any(len(v) == 0 for v in per_site_ranked.values()):
        return {'topk': [], 'per_site': {str(pos): {} for pos in target_positions}}

    combos = []
    ordered_positions = [int(p) for p in target_positions]
    frame_lists = [per_site_ranked[pos] for pos in ordered_positions]
    for combo in itertools.product(*frame_lists):
        total_score = sum(item['score'] for item in combo)
        combo_entry = {
            'combo_score': round(float(total_score), 6),
            'per_site': {},
        }
        for pos, item in zip(ordered_positions, combo):
            combo_entry['per_site'][str(pos)] = {
                'index': item['index'],
                'angle': item['angle'],
                'score': round(float(item['score']), 6),
                'coords': item['coords'],
                'rmsd': item['rmsd'],
            }
        combos.append(combo_entry)
    combos.sort(key=lambda x: x['combo_score'])
    topk = combos[:max(1, int(k))]

    per_site_summary: Dict[str, Dict[str, Any]] = {}
    for pos in ordered_positions:
        pos_str = str(pos)
        rmsd_vals = [entry['per_site'][pos_str].get('rmsd') for entry in topk]
        best_rmsd = min((v for v in rmsd_vals if v is not None), default=None)
        stage1_rmsd = None
        ranked = per_site_ranked[pos]
        if ranked:
            stage1_rmsd = min((item['rmsd'] for item in ranked if item['rmsd'] is not None), default=None)
        per_site_summary[pos_str] = {
            'top1_rmsd': rmsd_vals[0] if len(rmsd_vals) > 0 else None,
            'top2_rmsd': rmsd_vals[1] if len(rmsd_vals) > 1 else None,
            'top3_rmsd': rmsd_vals[2] if len(rmsd_vals) > 2 else None,
            'top3_best_rmsd': best_rmsd,
            'stage1_best_rmsd': stage1_rmsd,
            'top3_matches_stage1_best': (
                best_rmsd is not None and stage1_rmsd is not None and abs(best_rmsd - stage1_rmsd) <= 1e-4
            ),
        }
    return {'topk': topk, 'per_site': per_site_summary}


def stage_site_snapshot(structure_path: Optional[str], chain_id: str, positions: Sequence[int], ref_phosphates: Dict[int, Dict[str, np.ndarray]], ref_model) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    if not structure_path or not os.path.exists(structure_path):
        return out
    structure = parse_structure(structure_path, Path(structure_path).stem)
    model = list(structure.get_models())[0]
    for pos in positions:
        ref_coords = ref_phosphates.get(pos)
        pred_res = get_residue(model, chain_id, pos)
        ref_res = get_residue(ref_model, chain_id, pos)
        pred_coords = get_phosphate_coords(pred_res) if pred_res else None
        rmsd = phosphate_rmsd(ref_coords, pred_coords) if ref_coords is not None and pred_coords is not None else None
        rmsd_named = phosphate_rmsd_named(ref_coords, pred_coords) if ref_coords is not None and pred_coords is not None else None
        p_disp = p_atom_displacement(ref_coords, pred_coords) if ref_coords is not None and pred_coords is not None else None
        pred_torsion = best_permuted_torsion(pred_res) if pred_res is not None else None
        t_err = torsion_error_label_invariant(ref_res, pred_res) if ref_res is not None and pred_res is not None else None
        out[str(pos)] = {
            'coords': _round_coord_map(pred_coords),
            'phosphate_rmsd': round(float(rmsd), 4) if rmsd is not None else None,
            'phosphate_rmsd_named': round(float(rmsd_named), 4) if rmsd_named is not None else None,
            'p_displacement': round(float(p_disp), 4) if p_disp is not None else None,
            'pred_torsion': round(float(pred_torsion), 2) if pred_torsion is not None else None,
            'torsion_error_deg': round(float(t_err), 2) if t_err is not None else None,
        }
    return out


def case_archive_name(acc_id: str, pdb_id: str, chain_id: str, target_sig: str, context_sig: str, graft_mode: str) -> str:
    parts = [
        sanitize_label(acc_id or 'NA'),
        sanitize_label(pdb_id or 'pdb'),
        sanitize_label(chain_id or 'chain'),
        f'target-{sanitize_label(target_sig)}',
        f'ctx-{sanitize_label(context_sig)}',
        f'mode-{sanitize_label(graft_mode)}',
    ]
    return '__'.join(parts)


def archive_case_outputs(archive_case_dir: str, row: Dict[str, Any], run_meta: Dict[str, Any], chain_id: str,
                         positions: Sequence[int], target_positions: Sequence[int], ref_model, ref_phosphates: Dict[int, Dict[str, np.ndarray]],
                         stripped_path: str, pdb_path: str, stage0_output_path: str, stage0_reports: Dict[str, Optional[str]],
                         output_path: str, stage3_reports: Dict[str, Any], stage0_site_metrics: Dict[str, Dict[str, Any]],
                         stage3_site_metrics: Dict[str, Dict[str, Any]], prescan_frames: Dict[str, Dict],
                         ref_environment: Dict[int, Dict[str, Any]], ref_qc: Dict[int, Dict[str, Any]],
                         stage0_metrics: Dict[int, Dict[str, Any]], stage_rows: Sequence[Dict[str, Any]], save_stage_pdbs: bool = True) -> Dict[str, Any]:
    os.makedirs(archive_case_dir, exist_ok=True)
    for sub in ['input', 'ref', 'stage0', 'stage1', 'stage2', 'stage3', 'meta']:
        os.makedirs(os.path.join(archive_case_dir, sub), exist_ok=True)

    target_sig = context_signature(list(target_positions))
    context_sig = context_signature(list(positions))

    write_json(os.path.join(archive_case_dir, 'meta', 'benchmark_row.json'), row)
    write_json(os.path.join(archive_case_dir, 'meta', 'context.json'), {
        'acc_id': row.get('ACC_ID', ''),
        'pdb_id': row.get('PDBID_1', ''),
        'chain_id': chain_id,
        'target_positions': list(target_positions),
        'target_signature': target_sig,
        'context_positions': list(positions),
        'context_signature': context_sig,
        'context_n_sites': len(list(positions)),
        'graft_mode': run_meta.get('graft_mode', ''),
    })
    write_json(os.path.join(archive_case_dir, 'meta', 'run_meta.json'), run_meta)

    input_ref = copy_if_exists(pdb_path, os.path.join(archive_case_dir, 'input', 'reference_phospho.cif'))
    input_stripped = copy_if_exists(stripped_path, os.path.join(archive_case_dir, 'input', 'stripped_input.cif'))

    ref_payload = {}
    for pos in target_positions:
        ref_res = get_residue(ref_model, chain_id, pos)
        ref_payload[str(pos)] = {
            'coords': _round_coord_map(ref_phosphates.get(pos)),
            'ref_torsion': round(float(best_permuted_torsion(ref_res)), 2) if ref_res is not None and best_permuted_torsion(ref_res) is not None else None,
            'environment': ref_environment.get(pos, {}),
            'ref_qc': ref_qc.get(pos, {}),
        }
    write_json(os.path.join(archive_case_dir, 'ref', 'reference_sites.json'), ref_payload)

    stage0_paths = {
        'structure': copy_if_exists(stage0_output_path, os.path.join(archive_case_dir, 'stage0', 'structure_s0.cif')) if save_stage_pdbs else None,
        'graft_report_tsv': copy_if_exists(stage0_reports.get('graft_report'), os.path.join(archive_case_dir, 'stage0', 'graft_report.tsv')),
        'site_metrics_tsv': copy_if_exists(stage0_reports.get('site_report'), os.path.join(archive_case_dir, 'stage0', 'site_metrics.tsv')),
        'site_metrics_json': copy_if_exists(stage0_reports.get('site_report_json'), os.path.join(archive_case_dir, 'stage0', 'site_metrics.json')),
        'prescan_frames_json': copy_if_exists(stage0_reports.get('prescan_frames_json'), os.path.join(archive_case_dir, 'stage0', 'prescan_frames.json')),
    }
    write_json(os.path.join(archive_case_dir, 'stage0', 'per_site.json'), {str(pos): stage0_site_metrics.get(f'{chain_id}:{pos}', {}) for pos in target_positions})
    write_json(os.path.join(archive_case_dir, 'stage0', 'coords.json'), stage_site_snapshot(stage0_output_path, chain_id, target_positions, ref_phosphates, ref_model))
    write_json(os.path.join(archive_case_dir, 'stage0', 'metrics.json'), {str(pos): stage0_metrics.get(pos, {}) for pos in target_positions})

    stage3_paths = {
        'structure': copy_if_exists(output_path, os.path.join(archive_case_dir, 'stage3', 'structure_s3.cif')) if save_stage_pdbs else None,
        'graft_report_tsv': copy_if_exists(stage3_reports.get('graft_report'), os.path.join(archive_case_dir, 'stage3', 'graft_report.tsv')),
        'site_metrics_tsv': copy_if_exists(stage3_reports.get('site_report'), os.path.join(archive_case_dir, 'stage3', 'site_metrics.tsv')),
        'site_metrics_json': copy_if_exists(stage3_reports.get('site_report_json'), os.path.join(archive_case_dir, 'stage3', 'site_metrics.json')),
        'prescan_frames_json': copy_if_exists(stage3_reports.get('prescan_frames_json'), os.path.join(archive_case_dir, 'stage3', 'prescan_frames.json')),
    }
    # Write the in-memory prescan frames with ground-truth labels (rmsd_to_ref,
    # is_best_rmsd, is_good_1A) added by three_stage_rmsd().  The file above is
    # the raw graft-script output; this one has the benchmark annotations.
    labeled_frames = {}
    for pos in target_positions:
        site_key = f'{chain_id}:{pos}'
        if site_key in prescan_frames:
            labeled_frames[site_key] = prescan_frames[site_key]
    if labeled_frames:
        write_json(os.path.join(archive_case_dir, 'stage3', 'prescan_frames_labeled.json'), labeled_frames)
    write_json(os.path.join(archive_case_dir, 'stage3', 'per_site.json'), {str(pos): stage3_site_metrics.get(f'{chain_id}:{pos}', {}) for pos in target_positions})
    write_json(os.path.join(archive_case_dir, 'stage3', 'coords.json'), stage_site_snapshot(output_path, chain_id, target_positions, ref_phosphates, ref_model))

    stage1_per_site: Dict[str, Any] = {}
    stage2_per_site: Dict[str, Any] = {}
    stage1_coord_map: Dict[int, Dict[str, Sequence[float]]] = {}
    stage2_coord_map: Dict[int, Dict[str, Sequence[float]]] = {}
    for pos in target_positions:
        site_key = f'{chain_id}:{pos}'
        details = prescan_stage_details(prescan_frames.get(site_key, {}), ref_phosphates.get(pos, {}))
        stage1_per_site[str(pos)] = {
            'best_index': details.get('stage1_best_index'),
            'best_angle': details.get('stage1_best_angle'),
            'best_rmsd': details.get('stage1_best_rmsd'),
            'best_coords': details.get('stage1_best_coords'),
            'best_frame': details.get('stage1_best_frame'),
        }
        stage2_per_site[str(pos)] = {
            'selected_index': details.get('stage2_selected_index'),
            'selected_angle': details.get('stage2_selected_angle'),
            'selected_rmsd': details.get('stage2_selected_rmsd'),
            'selected_coords': details.get('stage2_selected_coords'),
            'selected_frame': details.get('stage2_selected_frame'),
            'scoring_gap': details.get('stage2_scoring_gap'),
        }
        if details.get('stage1_best_coords'):
            stage1_coord_map[int(pos)] = details['stage1_best_coords']
        if details.get('stage2_selected_coords'):
            stage2_coord_map[int(pos)] = details['stage2_selected_coords']
    write_json(os.path.join(archive_case_dir, 'stage1', 'per_site.json'), stage1_per_site)
    write_json(os.path.join(archive_case_dir, 'stage2', 'per_site.json'), stage2_per_site)
    if save_stage_pdbs and stage0_output_path and os.path.exists(stage0_output_path):
        stage1_structure_path = write_stage_structure_from_base(stage0_output_path, os.path.join(archive_case_dir, 'stage1', 'structure_s1.cif'), chain_id, stage1_coord_map)
        stage2_structure_path = write_stage_structure_from_base(stage0_output_path, os.path.join(archive_case_dir, 'stage2', 'structure_s2.cif'), chain_id, stage2_coord_map)
    else:
        stage1_structure_path = None
        stage2_structure_path = None
    write_json(os.path.join(archive_case_dir, 'stage1', 'coords.json'), stage_site_snapshot(stage1_structure_path, chain_id, target_positions, ref_phosphates, ref_model))
    write_json(os.path.join(archive_case_dir, 'stage2', 'coords.json'), stage_site_snapshot(stage2_structure_path, chain_id, target_positions, ref_phosphates, ref_model))

    topk_paths: Dict[str, Optional[str]] = {}
    topk_meta: Dict[str, Any] = {}
    pose_sources = list(stage3_reports.get('output_structures', []))
    pose_graft_reports = list(stage3_reports.get('graft_reports', []))
    pose_site_reports = list(stage3_reports.get('site_reports', []))
    pose_site_jsons = list(stage3_reports.get('site_report_jsons', []))
    for rank, source_path in enumerate(pose_sources, start=1):
        rank_path = None
        if save_stage_pdbs:
            rank_path = copy_if_exists(
                source_path,
                os.path.join(archive_case_dir, 'stage3', f'phosphoFill{rank}.cif'),
            )
        if rank_path:
            topk_paths[f'phosphoFill{rank}'] = rank_path
        topk_meta[f'phosphoFill{rank}'] = {
            'rank': rank,
            'post_minimization': True,
            'coords_snapshot': stage_site_snapshot(rank_path or source_path, chain_id, target_positions, ref_phosphates, ref_model),
            'graft_report': copy_if_exists(
                pose_graft_reports[rank - 1] if rank - 1 < len(pose_graft_reports) else None,
                os.path.join(archive_case_dir, 'stage3', f'graft_report_rank{rank}.tsv'),
            ),
            'site_metrics': copy_if_exists(
                pose_site_reports[rank - 1] if rank - 1 < len(pose_site_reports) else None,
                os.path.join(archive_case_dir, 'stage3', f'site_metrics_rank{rank}.tsv'),
            ),
            'site_metrics_json': copy_if_exists(
                pose_site_jsons[rank - 1] if rank - 1 < len(pose_site_jsons) else None,
                os.path.join(archive_case_dir, 'stage3', f'site_metrics_rank{rank}.json'),
            ),
        }
    write_json(os.path.join(archive_case_dir, 'stage3', 'top3_ranked_poses.json'), {
        'definition': 'three independently minimized rank-matched prescan candidates',
        'paths': topk_paths,
        'structures': topk_meta,
    })

    summary_rows = {str(r.get('target_position')): r for r in stage_rows}
    write_json(os.path.join(archive_case_dir, 'meta', 'result_rows.json'), summary_rows)

    write_json(os.path.join(archive_case_dir, 'meta', 'paths.json'), {
        'input_reference_phospho': input_ref,
        'input_stripped': input_stripped,
        'stage0': stage0_paths,
        'stage1': {'structure': stage1_structure_path},
        'stage2': {'structure': stage2_structure_path},
        'stage3': {**stage3_paths, **topk_paths},
    })

    return {
        'archive_case_dir': archive_case_dir,
        'archive_stage0_structure': stage0_paths.get('structure') or '',
        'archive_stage1_structure': stage1_structure_path or '',
        'archive_stage2_structure': stage2_structure_path or '',
        'archive_stage3_structure': stage3_paths.get('structure') or '',
        'archive_stage3_top1_structure': topk_paths.get('phosphoFill1') or '',
        'archive_stage3_top2_structure': topk_paths.get('phosphoFill2') or '',
        'archive_stage3_top3_structure': topk_paths.get('phosphoFill3') or '',
        'archive_stage3_prescan_json': stage3_paths.get('prescan_frames_json') or '',
        'archive_stage3_site_metrics_json': stage3_paths.get('site_metrics_json') or '',
    }


def load_report_metrics(report_tsv_path: str) -> Dict[str, Dict]:
    out = {}
    with open(report_tsv_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = f"{row.get('chain_id', '')}:{row.get('resseq', '')}"
            out[key] = row
    return out


def parse_positions(text: str) -> List[int]:
    out = []
    for chunk in str(text).split(","):
        chunk = chunk.strip()
        if chunk and chunk.lstrip("-").isdigit():
            out.append(int(chunk))
    return out


def sanitize_label(text: str) -> str:
    keep = []
    for ch in str(text):
        if ch.isalnum() or ch in {"-", "_", "."}:
            keep.append(ch)
        else:
            keep.append("_")
    return "".join(keep)


def context_signature(positions: List[int]) -> str:
    return ",".join(str(p) for p in sorted(positions))


def run_mode_label(graft_mode: str) -> str:
    mode = str(graft_mode or "joint").strip().lower()
    return mode if mode in {"joint", "single"} else "joint"


def parse_restypes(text: str) -> List[str]:
    return [x.strip() for x in str(text).split(",") if x.strip()]


def _run_stage0(graft_script: str, report_script: Optional[str], input_path: str, stage0_output_path: str,
                sites: List[str], cache_dir: str, ref_coords: Dict[int, Dict[str, np.ndarray]],
                chain_id: str, positions: List[int]) -> Tuple[Dict[int, Dict[str, Optional[float]]], Dict[str, Optional[str]], Dict[str, Dict[str, Any]]]:
    stage0_args = ["--no-prescan-rotamer", "--relax-mode", "none"]
    success, _, report_paths = run_phosphofill(graft_script, input_path, stage0_output_path, sites, cache_dir,
                                               report_script=report_script, extra_args=stage0_args, graft_mode="joint")
    metrics: Dict[int, Dict[str, Optional[float]]] = {}
    report_metrics = load_report_metrics(report_paths.get('site_report') or report_paths.get('graft_report')) if success and (report_paths.get('site_report') or report_paths.get('graft_report')) else {}
    if not success or not os.path.exists(stage0_output_path):
        for pos in positions:
            metrics[pos] = {"stage0_kabsch_rmsd": None, "stage0_p_displacement": None}
        return metrics, report_paths, report_metrics
    stage0_structure = parse_structure(stage0_output_path, "stage0")
    stage0_model = list(stage0_structure.get_models())[0]
    for pos in positions:
        pred_res = get_residue(stage0_model, chain_id, pos)
        pred_coords = get_phosphate_coords(pred_res) if pred_res else None
        refc = ref_coords.get(pos)
        if refc is None or pred_coords is None:
            metrics[pos] = {"stage0_kabsch_rmsd": None, "stage0_p_displacement": None}
        else:
            s0_rmsd = phosphate_rmsd(refc, pred_coords)
            s0_p = p_atom_displacement(refc, pred_coords)
            metrics[pos] = {
                "stage0_kabsch_rmsd": round(s0_rmsd, 4) if s0_rmsd is not None else None,
                "stage0_p_displacement": round(s0_p, 4) if s0_p is not None else None,
            }
    return metrics, report_paths, report_metrics


def validate_tier1(row: Dict, graft_script: str, cache_dir: str, work_dir: str,
                   report_script: Optional[str] = None, extra_args: Optional[List[str]] = None,
                   graft_mode: str = "joint", archive_root: Optional[str] = None,
                   archive_run_tag: Optional[str] = None, save_stage_pdbs: bool = True) -> List[Dict]:
    pdb_id = row["PDBID_1"].strip().lower()
    chain_id = row["CHAINID_1"].strip()
    positions = parse_positions(row["newpos_1"])
    restypes = parse_restypes(row["MODRES_1"])
    if not positions:
        return [{"tier": 1, "status": "SKIP", "message": "no positions"}]

    graft_mode = run_mode_label(graft_mode)
    context_sig = context_signature(positions)
    context_tag = sanitize_label(f"{pdb_id}_{chain_id}_{context_sig}")

    pdb_path = download_pdb(pdb_id, cache_dir)
    ref_structure = parse_structure(pdb_path, "ref")
    ref_model = list(ref_structure.get_models())[0]

    ref_phosphates, ref_torsions, ref_environment, ref_qc = {}, {}, {}, {}
    for pos in positions:
        ref_res = get_residue(ref_model, chain_id, pos)
        if ref_res is not None:
            ref_phosphates[pos] = get_phosphate_coords(ref_res)
            ref_torsions[pos] = best_permuted_torsion(ref_res)
            ref_environment[pos] = annotate_phosphosite_environment(ref_model, chain_id, pos)
            ref_qc[pos] = reference_pose_qc(ref_model, chain_id, pos)

    stripped_structure = parse_structure(pdb_path, "stripped")
    strip_phosphate_from_structure(stripped_structure, chain_id, positions)
    stripped_path = os.path.join(work_dir, f"{context_tag}_stripped.cif")
    save_structure_file(stripped_structure, stripped_path)

    if graft_mode == "single":
        target_runs = [[pos] for pos in positions]
    else:
        target_runs = [positions]

    results: List[Dict] = []
    for target_positions in target_runs:
        sites = [f"{chain_id}:{p}" for p in target_positions]
        target_sig = context_signature(target_positions)
        target_tag = sanitize_label(f"{context_tag}_{graft_mode}_{target_sig}")

        stage0_output_path = os.path.join(work_dir, f"{target_tag}_stage0_naive.cif")
        stage0_metrics, stage0_report_paths, stage0_report_metrics = _run_stage0(graft_script, report_script, stripped_path, stage0_output_path, sites,
                                         cache_dir, ref_phosphates, chain_id, target_positions)

        output_path = os.path.join(work_dir, f"{target_tag}_regrafted.cif")
        success, message, report_paths = run_phosphofill(
            graft_script, stripped_path, output_path, sites, cache_dir,
            report_script=report_script,
            extra_args=(["--prescan-rotamer", "--relax-mode", "openmm_local"] + (extra_args or [])),
            graft_mode=graft_mode,
            n_poses=3,
        )
        if not success:
            for pos in target_positions:
                results.append({
                    "tier": 1, "acc_id": row.get("ACC_ID", ""), "pdb_id": pdb_id,
                    "chain": chain_id, "position": pos, "target_position": pos,
                    "target_position_signature": str(pos),
                    "context_positions": row.get("newpos_1", ""),
                    "context_signature": context_sig,
                    "graft_mode": graft_mode,
                    "context_n_sites": len(positions),
                    "n_sites": len(target_positions),
                    "status": "FAILED", "message": message,
                })
            continue

        pose_paths = list(report_paths.get("output_structures", []))
        pose_structures = [parse_structure(path, f"pred_rank_{rank}")
                           for rank, path in enumerate(pose_paths, start=1)]
        pose_models = [list(structure.get_models())[0] for structure in pose_structures]
        site_report_paths = list(report_paths.get("site_reports", []))
        graft_report_paths = list(report_paths.get("graft_reports", []))
        report_metrics_by_rank: List[Dict[str, Dict]] = []
        for rank in range(len(pose_paths)):
            site_report_path = site_report_paths[rank] if rank < len(site_report_paths) else None
            graft_report_path = graft_report_paths[rank] if rank < len(graft_report_paths) else None
            metrics_path = site_report_path or graft_report_path
            report_metrics_by_rank.append(load_report_metrics(metrics_path) if metrics_path else {})
        report_metrics = report_metrics_by_rank[0] if report_metrics_by_rank else {}
        frames_path = report_paths.get("prescan_frames_json")
        prescan_frames = load_prescan_frames(frames_path) if frames_path else {}

        stage_rows: List[Dict[str, Any]] = []
        for pos in target_positions:
            idx = positions.index(pos) if pos in positions else -1
            ref_coords = ref_phosphates.get(pos)
            ref_res = get_residue(ref_model, chain_id, pos)
            pred_residues = [get_residue(model, chain_id, pos) for model in pose_models]
            pred_coords_by_rank = [get_phosphate_coords(residue) if residue else None
                                   for residue in pred_residues]
            pred_res = pred_residues[0] if pred_residues else None
            pred_coords = pred_coords_by_rank[0] if pred_coords_by_rank else None
            site_key = f"{chain_id}:{pos}"
            site_metric_row = report_metrics.get(site_key, {})
            site_status = str(site_metric_row.get("status", "")).strip().upper()
            site_message = str(site_metric_row.get("message", "")).strip()
            if ref_coords is None or pred_coords is None:
                stage_rows.append({
                    "tier": 1, "acc_id": row.get("ACC_ID", ""), "pdb_id": pdb_id, "chain": chain_id,
                    "position": pos, "target_position": pos, "target_position_signature": str(pos),
                    "context_positions": row.get("newpos_1", ""),
                    "context_signature": context_sig,
                    "graft_mode": graft_mode,
                    "context_n_sites": len(positions),
                    "n_sites": len(target_positions),
                    "restype": restypes[idx] if 0 <= idx < len(restypes) else "?",
                    "status": ("FAILED" if site_status == "FAILED" else "NO_COORDS"),
                    "message": (site_message if site_status == "FAILED" and site_message else f"ref={ref_coords is not None} pred={pred_coords is not None}"),
                })
                continue
            rank_rmsds = [phosphate_rmsd(ref_coords, coords) if coords is not None else None
                          for coords in pred_coords_by_rank]
            rank_named_rmsds = [phosphate_rmsd_named(ref_coords, coords) if coords is not None else None
                                 for coords in pred_coords_by_rank]
            rank_p_displacements = [p_atom_displacement(ref_coords, coords) if coords is not None else None
                                    for coords in pred_coords_by_rank]
            rmsd = rank_rmsds[0] if rank_rmsds else None
            rmsd_named = rank_named_rmsds[0] if rank_named_rmsds else None
            p_disp = rank_p_displacements[0] if rank_p_displacements else None
            valid_rank_rmsds = [(rank + 1, value) for rank, value in enumerate(rank_rmsds) if value is not None]
            top3_best_rank, top3_best_rmsd = min(valid_rank_rmsds, key=lambda item: item[1]) if valid_rank_rmsds else (None, None)
            pred_torsion = best_permuted_torsion(pred_res) if pred_res else None
            t_err = torsion_error_label_invariant(ref_res, pred_res) if ref_res is not None and pred_res is not None else None
            site_key = f"{chain_id}:{pos}"
            stages = three_stage_rmsd(prescan_frames.get(site_key, {}), ref_coords)
            stage3_minimization_delta = None
            if rmsd is not None and stages["stage2_selected_rmsd"] is not None:
                stage3_minimization_delta = round(rmsd - stages["stage2_selected_rmsd"], 4)
            rmetrics = report_metrics.get(site_key, {})
            s0 = stage0_metrics.get(pos, {"stage0_kabsch_rmsd": None, "stage0_p_displacement": None})
            stage_rows.append({
                "tier": 1, "acc_id": row.get("ACC_ID", ""), "pdb_id": pdb_id, "chain": chain_id, "position": pos,
                "target_position": pos,
                "target_position_signature": str(pos),
                "context_positions": row.get("newpos_1", ""),
                "context_signature": context_sig,
                "graft_mode": graft_mode,
                "context_n_sites": len(positions),
                "restype": restypes[idx] if 0 <= idx < len(restypes) else "?", "n_sites": len(target_positions),
                "resolution": row.get("Resolution_1", ""), "status": "OK",
                "stage0_kabsch_rmsd": s0.get("stage0_kabsch_rmsd"),
                "stage0_seed_rmsd": s0.get("stage0_kabsch_rmsd"),
                "stage0_p_displacement": s0.get("stage0_p_displacement"),
                "phosphate_rmsd": round(rmsd, 4) if rmsd is not None else None,
                "phosphate_rmsd_named": round(rmsd_named, 4) if rmsd_named is not None else None,
                "p_displacement": round(p_disp, 4) if p_disp is not None else None,
                "torsion_error_deg": round(t_err, 2) if t_err is not None else None,
                "ref_torsion": round(ref_torsions.get(pos, 0), 2) if ref_torsions.get(pos) is not None else None,
                "pred_torsion": round(pred_torsion, 2) if pred_torsion is not None else None,
                **stages,
                "rank1_postmin_rmsd": round(rank_rmsds[0], 4) if len(rank_rmsds) > 0 and rank_rmsds[0] is not None else None,
                "rank2_postmin_rmsd": round(rank_rmsds[1], 4) if len(rank_rmsds) > 1 and rank_rmsds[1] is not None else None,
                "rank3_postmin_rmsd": round(rank_rmsds[2], 4) if len(rank_rmsds) > 2 and rank_rmsds[2] is not None else None,
                "top3_best_postmin_rmsd": round(top3_best_rmsd, 4) if top3_best_rmsd is not None else None,
                "top3_best_postmin_rank": top3_best_rank,
                # Backward-compatible aliases now explicitly contain
                # post-minimization values.
                "top3_rank1_rmsd": round(rank_rmsds[0], 4) if len(rank_rmsds) > 0 and rank_rmsds[0] is not None else None,
                "top3_rank2_rmsd": round(rank_rmsds[1], 4) if len(rank_rmsds) > 1 and rank_rmsds[1] is not None else None,
                "top3_rank3_rmsd": round(rank_rmsds[2], 4) if len(rank_rmsds) > 2 and rank_rmsds[2] is not None else None,
                "top3_best_rmsd": round(top3_best_rmsd, 4) if top3_best_rmsd is not None else None,
                "top3_matches_stage1_best": "",
                "stage3_post_minimization_rmsd": round(rmsd, 4) if rmsd is not None else None,
                "stage3_minimization_delta": stage3_minimization_delta,
                "prescan_n_zero_clash": rmetrics.get("prescan_n_zero_clash", ""),
                "prescan_zero_clash_arc": rmetrics.get("prescan_zero_clash_arc", ""),
                "pf_confidence": rmetrics.get("phosphofill_confidence", ""),
                "quality_label": rmetrics.get("phospho_quality_label", ""),
                "clashes_after": rmetrics.get("clashes_after_relax", ""),
                "benchmark_label": benchmark_success_label(rmsd, p_disp, _parse_optional_number(rmetrics.get("clashes_after_relax", ""))),
                **{k: v for k, v in ref_environment.get(pos, {}).items()},
                **{k: v for k, v in ref_qc.get(pos, {}).items()},
            })

        archive_meta = {}
        if archive_root:
            run_tag = sanitize_label(archive_run_tag or 'run')
            archive_case_dir = os.path.join(archive_root, run_tag, case_archive_name(row.get('ACC_ID', ''), pdb_id, chain_id, target_sig, context_sig, graft_mode))
            archive_meta = archive_case_outputs(
                archive_case_dir=archive_case_dir,
                row=row,
                run_meta={
                    'graft_mode': graft_mode,
                    'archive_run_tag': archive_run_tag or '',
                    'extra_args': list(extra_args or []),
                },
                chain_id=chain_id,
                positions=positions,
                target_positions=target_positions,
                ref_model=ref_model,
                ref_phosphates=ref_phosphates,
                stripped_path=stripped_path,
                pdb_path=pdb_path,
                stage0_output_path=stage0_output_path,
                stage0_reports=stage0_report_paths,
                output_path=pose_paths[0],
                stage3_reports=report_paths,
                stage0_site_metrics=stage0_report_metrics,
                stage3_site_metrics=report_metrics,
                prescan_frames=prescan_frames,
                ref_environment=ref_environment,
                ref_qc=ref_qc,
                stage0_metrics=stage0_metrics,
                stage_rows=stage_rows,
                save_stage_pdbs=save_stage_pdbs,
            )
        for r in stage_rows:
            r.update(archive_meta)
            results.append(r)
    return results


def compute_summary(results: List[Dict]) -> Dict:
    ok_results = [r for r in results if r.get("status") == "OK"]
    def numeric_values(column: str) -> List[float]:
        values = [_optional_float(row.get(column)) for row in ok_results]
        return [value for value in values if value is not None]

    stage0s = numeric_values("stage0_seed_rmsd")
    stage1s = numeric_values("stage1_best_of_12_rmsd")
    stage2s = numeric_values("stage2_selected_rmsd")
    rmsds = numeric_values("rank1_postmin_rmsd")
    top3_rmsds = numeric_values("top3_best_postmin_rmsd")
    scoring_gaps = numeric_values("stage2_scoring_gap")

    def stats(values):
        if not values:
            return {"n": 0}
        arr = np.array(values)
        return {"n": len(arr), "mean": round(float(np.mean(arr)), 4), "median": round(float(np.median(arr)), 4),
                "std": round(float(np.std(arr)), 4), "min": round(float(np.min(arr)), 4), "max": round(float(np.max(arr)), 4)}

    return {
        "total_rows": len(results),
        "ok": len(ok_results),
        "stage0_seed_rmsd": stats(stage0s),
        "stage1_best_of_12_rmsd": stats(stage1s),
        "stage2_selected_rmsd": stats(stage2s),
        "stage3_post_minimization_rmsd": stats(rmsds),
        "top3_best_postmin_rmsd": stats(top3_rmsds),
        "stage2_scoring_gap": stats(scoring_gaps),
    }


def parse_weight_sweep(text_value: Optional[str]) -> List[Optional[float]]:
    if text_value is None or str(text_value).strip() == "":
        return [None]
    return [float(chunk.strip()) for chunk in str(text_value).split(",") if chunk.strip()] or [None]


def _optional_float(value: Any) -> Optional[float]:
    if value is None or str(value).strip() == "":
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if np.isfinite(result) else None


def _same_optional_number(value: Any, expected: Optional[float]) -> bool:
    parsed = _optional_float(value)
    if expected is None:
        return parsed is None
    return parsed is not None and abs(parsed - expected) <= 1e-12


def _preferred_output_columns() -> List[str]:
    return [
        "tier", "benchmark_row_index", "benchmark_case_key", "run_graft_mode", "graft_mode",
        "prescan_contact_weight", "prescan_pack_weight", "acc_id", "pdb_id", "chain", "position",
        "target_position", "target_position_signature", "context_positions", "context_signature", "context_n_sites",
        "restype", "n_sites", "resolution", "status", "message",
        "stage0_seed_rmsd", "stage0_kabsch_rmsd", "stage0_p_displacement",
        "stage1_best_of_12_rmsd", "stage1_best_of_12_angle",
        "stage2_selected_rmsd", "stage2_selected_angle", "stage2_scoring_gap",
        "rank1_postmin_rmsd", "rank2_postmin_rmsd", "rank3_postmin_rmsd",
        "top3_best_postmin_rmsd", "top3_best_postmin_rank",
        "top3_rank1_rmsd", "top3_rank2_rmsd", "top3_rank3_rmsd", "top3_best_rmsd", "top3_matches_stage1_best",
        "stage3_post_minimization_rmsd", "stage3_minimization_delta",
        "phosphate_rmsd", "phosphate_rmsd_named", "p_displacement", "torsion_error_deg",
        "ref_torsion", "pred_torsion", "pf_confidence", "quality_label", "clashes_after", "benchmark_label",
        "n_waters_4A", "metal_within_5A", "metal_names", "nearest_metal_dist", "nearest_water_dist",
        "site_secondary_structure", "site_category",
        "ref_phosphate_b_mean", "ref_backbone_b_mean", "ref_b_ratio", "ref_phosphate_occupancy_mean",
        "ref_n_basic_contacts_4A", "ref_n_hbondlike_contacts_35A", "ref_nearest_basic_dist",
        "ref_nearest_any_contact", "ref_po_bond_mean", "ref_po_bond_std", "ref_pose_support_class",
        "archive_case_dir", "archive_stage0_structure", "archive_stage1_structure", "archive_stage2_structure", "archive_stage3_structure",
        "archive_stage3_top1_structure", "archive_stage3_top2_structure", "archive_stage3_top3_structure",
        "archive_stage3_prescan_json", "archive_stage3_site_metrics_json",
    ]


def write_results_tsv(path: str, results: Sequence[Dict[str, Any]]) -> None:
    all_keys: List[str] = []
    seen = set()
    for key in _preferred_output_columns():
        if any(key in result for result in results):
            all_keys.append(key)
            seen.add(key)
    for result in results:
        for key in result:
            if key not in seen:
                all_keys.append(key)
                seen.add(key)
    if not all_keys:
        return
    temporary_path = f"{path}.tmp"
    with open(temporary_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_keys, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    os.replace(temporary_path, path)


def main():
    ap = argparse.ArgumentParser(description="Tier 1 benchmark with Stage 0 and crystal-pose QC.")
    ap.add_argument("benchmark_tsv", help="Tab-separated table with paired phospho/unmod PDBs")
    ap.add_argument("--graft-script", required=True, help="Path to graft_phospho_openmm_v4.py")
    ap.add_argument("--report-script", default=None, help="Optional path to phosphofill_report_v4.py")
    ap.add_argument("--output-prefix", default="benchmark_results_refqc")
    ap.add_argument("--cache-dir", default="benchmark_cache")
    ap.add_argument("--work-dir", default=None)
    ap.add_argument("--max-rows", type=int, default=None)
    ap.add_argument("--prescan-contact-weight", type=float, default=None)
    ap.add_argument("--prescan-weight-sweep", default=None)
    ap.add_argument("--prescan-basic-weight", type=float, default=None)
    ap.add_argument("--prescan-hbond-weight", type=float, default=None)
    ap.add_argument("--prescan-polar-clash-scale", type=float, default=None)
    ap.add_argument("--prescan-pack-weight", type=float, default=None)
    ap.add_argument("--prescan-pack-weight-sweep", default=None)
    ap.add_argument("--prescan-geom-prior-weight", type=float, default=None)
    ap.add_argument("--prescan-tpo-torsion-sigma", type=float, default=None)
    ap.add_argument("--prescan-keep-kabsch-if-good", action="store_true")
    ap.add_argument("--graft-mode", choices=["joint", "single"], default="joint", help="Benchmark mode: graft all context sites together or one target site at a time")
    ap.add_argument("--archive-root", default=None, help="Optional root directory for per-case archived stage outputs")
    ap.add_argument("--no-save-stage-pdbs", action="store_true", help="Archive JSON/metrics only and skip writing stage0/1/2/3 PDB files")
    ap.add_argument("--resume", action="store_true", help="Continue from an existing output TSV, skipping checkpointed benchmark cases")
    ap.add_argument("--checkpoint-every", type=int, default=1, help="Rewrite the output TSV after this many newly processed benchmark rows (default: 1)")
    args = ap.parse_args()

    if args.checkpoint_every < 1:
        raise ValueError("--checkpoint-every must be at least 1")

    os.makedirs(args.cache_dir, exist_ok=True)
    with open(args.benchmark_tsv, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if args.max_rows:
        rows = rows[:args.max_rows]

    if args.prescan_weight_sweep and args.prescan_contact_weight is not None:
        raise ValueError("Use either --prescan-contact-weight or --prescan-weight-sweep, not both")
    if args.prescan_pack_weight_sweep and args.prescan_pack_weight is not None:
        raise ValueError("Use either --prescan-pack-weight or --prescan-pack-weight-sweep, not both")

    weights = parse_weight_sweep(args.prescan_weight_sweep)
    if args.prescan_contact_weight is not None:
        weights = [args.prescan_contact_weight]
    pack_weights = parse_weight_sweep(args.prescan_pack_weight_sweep)
    if args.prescan_pack_weight is not None:
        pack_weights = [args.prescan_pack_weight]

    print(f"Loaded {len(rows)} benchmark pairs")
    print(f"Prescan contact weights: {weights}")
    print(f"Prescan pack weights: {pack_weights}")

    results_path = f"{args.output_prefix}.tsv"
    os.makedirs(os.path.dirname(os.path.abspath(results_path)), exist_ok=True)
    all_results: List[Dict] = []
    if args.resume and os.path.exists(results_path):
        with open(results_path, newline="", encoding="utf-8") as handle:
            all_results = list(csv.DictReader(handle, delimiter="\t"))
        print(f"Loaded {len(all_results)} checkpointed result rows from {results_path}")
    completed_case_keys = {
        str(result.get("benchmark_case_key", "")).strip()
        for result in all_results
        if str(result.get("benchmark_case_key", "")).strip()
    }
    processed_since_checkpoint = 0
    work_dir = args.work_dir or tempfile.mkdtemp(prefix="phosphofill_refqc_")
    os.makedirs(work_dir, exist_ok=True)
    save_stage_pdbs = not args.no_save_stage_pdbs
    if args.archive_root:
        os.makedirs(args.archive_root, exist_ok=True)
        write_json(os.path.join(args.archive_root, "manifest.json"), {
            "benchmark_tsv": args.benchmark_tsv,
            "graft_script": args.graft_script,
            "report_script": args.report_script,
            "output_prefix": args.output_prefix,
            "cache_dir": args.cache_dir,
            "work_dir": args.work_dir or work_dir,
            "graft_mode": args.graft_mode,
            "weights": weights,
            "pack_weights": pack_weights,
            "prescan_basic_weight": args.prescan_basic_weight,
            "prescan_hbond_weight": args.prescan_hbond_weight,
            "prescan_polar_clash_scale": args.prescan_polar_clash_scale,
            "prescan_geom_prior_weight": args.prescan_geom_prior_weight,
            "prescan_tpo_torsion_sigma": args.prescan_tpo_torsion_sigma,
            "prescan_keep_kabsch_if_good": args.prescan_keep_kabsch_if_good,
            "save_stage_pdbs": save_stage_pdbs,
        })

    try:
        for weight in weights:
            for pack_weight in pack_weights:
                extra_args: List[str] = []
                weight_label = "default" if weight is None else str(weight)
                pack_label = "default" if pack_weight is None else str(pack_weight)
                if weight is not None:
                    extra_args.extend(["--prescan-contact-weight", str(weight)])
                if pack_weight is not None:
                    extra_args.extend(["--prescan-pack-weight", str(pack_weight)])
                if args.prescan_basic_weight is not None:
                    extra_args.extend(["--prescan-basic-weight", str(args.prescan_basic_weight)])
                if args.prescan_hbond_weight is not None:
                    extra_args.extend(["--prescan-hbond-weight", str(args.prescan_hbond_weight)])
                if args.prescan_polar_clash_scale is not None:
                    extra_args.extend(["--prescan-polar-clash-scale", str(args.prescan_polar_clash_scale)])
                if args.prescan_geom_prior_weight is not None:
                    extra_args.extend(["--prescan-geom-prior-weight", str(args.prescan_geom_prior_weight)])
                if args.prescan_tpo_torsion_sigma is not None:
                    extra_args.extend(["--prescan-tpo-torsion-sigma", str(args.prescan_tpo_torsion_sigma)])
                if args.prescan_keep_kabsch_if_good:
                    extra_args.append("--prescan-keep-kabsch-if-good")

                print("\n" + "#" * 70)
                print(f"Running Tier 1 with graft_mode={args.graft_mode}, prescan_contact_weight={weight_label}, prescan_pack_weight={pack_label}")
                print("#" * 70)
                archive_run_tag = sanitize_label(f"mode-{args.graft_mode}__w-{weight_label}__pack-{pack_label}")

                for i, row in enumerate(rows):
                    acc = row.get("ACC_ID", "?")
                    pdb1 = row.get("PDBID_1", "?")
                    case_key = f"row-{i + 1}__mode-{args.graft_mode}__w-{weight_label}__pack-{pack_label}"
                    if case_key in completed_case_keys:
                        print(f"\n[{i+1}/{len(rows)}] {acc} | phospho={pdb1} | checkpointed; skipping")
                        continue
                    print(f"\n[{i+1}/{len(rows)}] {acc} | phospho={pdb1} | weight={weight_label} | pack={pack_label}")
                    try:
                        results = validate_tier1(row, args.graft_script, args.cache_dir, work_dir,
                                                 report_script=args.report_script, extra_args=extra_args,
                                                 graft_mode=args.graft_mode, archive_root=args.archive_root,
                                                 archive_run_tag=archive_run_tag, save_stage_pdbs=save_stage_pdbs)
                        for r in results:
                            r["prescan_contact_weight"] = weight if weight is not None else ""
                            r["prescan_pack_weight"] = pack_weight if pack_weight is not None else ""
                            r["run_graft_mode"] = args.graft_mode
                            r["benchmark_row_index"] = i + 1
                            r["benchmark_case_key"] = case_key
                        all_results.extend(results)
                        for r in results:
                            if r.get("status") == "OK":
                                print(
                                    f"    {r.get('chain','')}:{r.get('position','')} "
                                    f"S0={r.get('stage0_kabsch_rmsd','?')}Å "
                                    f"S1={r.get('stage1_best_of_12_rmsd','?')}Å "
                                    f"S2={r.get('stage2_selected_rmsd','?')}Å "
                                    f"S3={r.get('stage3_post_minimization_rmsd','?')}Å "
                                    f"refQC={r.get('ref_pose_support_class','?')}",
                                    flush=True,
                                )
                            else:
                                print(
                                    f"    {r.get('chain','')}:{r.get('position','?')} "
                                    f"status={r.get('status','?')} "
                                    f"message={r.get('message','')}",
                                    flush=True,
                                )
                    except Exception as e:
                        print(f"    ERROR: {e}")
                        all_results.append({"tier": 1, "pdb_id": pdb1, "status": "ERROR",
                                            "message": str(e), "prescan_contact_weight": weight if weight is not None else "",
                                            "prescan_pack_weight": pack_weight if pack_weight is not None else "", "run_graft_mode": args.graft_mode,
                                            "benchmark_row_index": i + 1, "benchmark_case_key": case_key})
                    completed_case_keys.add(case_key)
                    processed_since_checkpoint += 1
                    if processed_since_checkpoint >= args.checkpoint_every:
                        write_results_tsv(results_path, all_results)
                        processed_since_checkpoint = 0
    finally:
        if args.work_dir is None:
            shutil.rmtree(work_dir, ignore_errors=True)

    write_results_tsv(results_path, all_results)
    print(f"\nWrote {len(all_results)} results to {results_path}")

    summary = {}
    for weight in weights:
        weight_key = "default" if weight is None else str(weight)
        for pack_weight in pack_weights:
            pack_key = "default" if pack_weight is None else str(pack_weight)
            weight_results = [
                r for r in all_results
                if (_same_optional_number(r.get("prescan_contact_weight", ""), weight)
                    and _same_optional_number(r.get("prescan_pack_weight", ""), pack_weight))
            ]
            summary[f"weight_{weight_key}__pack_{pack_key}"] = compute_summary(weight_results)

    summary_path = f"{args.output_prefix}_summary.json"
    with open(summary_path, "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"Wrote summary to {summary_path}")


if __name__ == "__main__":
    main()
