#!/usr/bin/env python3
"""Compare unmodified, phospho-crystal, and predicted site environments."""
from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
from PIL import Image, ImageDraw, ImageFont

import case_study_tool_analysis as base


ROOT = base.ROOT
OUT_DIR = base.OUT_DIR
FIG_DIR = base.FIG_DIR
PYMOL_DIR = OUT_DIR / "pymol"

CONTACT_CUTOFF_A = 5.0
BASIC_CUTOFF_A = 4.0
PHOSPHATE_NAMES = ("P", "O1P", "O2P", "O3P")
PHOSPHATE_OXYGEN_NAMES = ("O1P", "O2P", "O3P")
PHOSPHATE_SELECTION = "P+O1P+O2P+O3P+OP1+OP2+OP3"
PHOSPHATE_OXYGEN_SELECTION = "O1P+O2P+O3P+OP1+OP2+OP3"
PHOSPHO_LINKAGE_ATOMS = {
    "SEP": "OG",
    "TPO": "OG1",
    "PTR": "OH",
}
UNMODIFIED_MARKER_ATOMS = {
    "TPO": ("OG1",),
    "PTR": ("OH",),
}
UNMODIFIED_SITE_MAP = {
    ("CDK2", "TPO", 160): {
        "pdb_id": "1b38",
        "chain": "A",
        "position": 160,
        "requested_position": 160,
        "note": "CDK2 1B38 A:160 is the unmodified Thr160 site.",
    },
    ("ERK2", "TPO", 185): {
        "pdb_id": "5umo",
        "chain": "A",
        "position": 183,
        "requested_position": 185,
        "note": "5UMO author numbering places the unmodified ERK TEY Thr at A:183; A:185 is Tyr.",
    },
    ("ERK2", "PTR", 187): {
        "pdb_id": "5umo",
        "chain": "A",
        "position": 185,
        "requested_position": 187,
        "note": "5UMO author numbering places the unmodified ERK TEY Tyr at A:185; A:187 is Ala.",
    },
}


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]+", "_", str(text)).strip("_")


def font(size: int, bold: bool = False):
    candidates = [
        "arialbd.ttf" if bold else "arial.ttf",
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
        "C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size=size)
        except Exception:
            pass
    return ImageFont.load_default()


def text(draw: ImageDraw.ImageDraw, xy: Tuple[float, float], value, size=12, fill="#222222", anchor="mm", bold=False):
    draw.text(xy, str(value), font=font(size, bold), fill=fill, anchor=anchor)


def has_backbone_rmsd_column(df: pd.DataFrame) -> bool:
    return "af2_crystal_local_backbone_rmsd" in df.columns


def backbone_rmsd_value(row) -> float:
    return base.safe_float(row.get("af2_crystal_local_backbone_rmsd", ""))


def backbone_rmsd_label(row) -> str:
    value = backbone_rmsd_value(row)
    return "NA" if math.isnan(value) else f"{value:.1f} A"


def backbone_rmsd_fill(row) -> str:
    value = backbone_rmsd_value(row)
    if math.isnan(value):
        return "#f2f2f2"
    if value <= 2.0:
        return "#cfead1"
    if value <= 3.0:
        return "#fff1bd"
    return "#f4c7c3"


def marker_coords_for_unmodified(residue, residue_type: str) -> Dict[str, object]:
    coords = base.atom_coords(residue)
    out = {}
    for name in UNMODIFIED_MARKER_ATOMS.get(residue_type, ()):
        if name in coords:
            out[name] = coords[name]
    return out


def marker_coords_for_phospho(residue) -> Dict[str, object]:
    return base.phosphate_coords(residue)


def is_heavy_atom(atom) -> bool:
    name = atom.get_name().strip().upper()
    element = getattr(atom, "element", "").strip().upper()
    return not (name.startswith("H") or element == "H")


def is_contact_residue(residue) -> bool:
    resname = residue.get_resname().strip().upper()
    if resname in {"HOH", "WAT"}:
        return False
    if residue.id[0].startswith("W"):
        return False
    return True


def residue_contacts(
    model,
    marker_coords: Dict[str, object],
    exclude: Optional[Tuple[str, int]],
    cutoff: float = CONTACT_CUTOFF_A,
) -> List[Dict]:
    if not marker_coords:
        return []
    rows = []
    for chain in model:
        for residue in chain:
            if not is_contact_residue(residue):
                continue
            resseq = int(residue.id[1])
            if exclude and chain.id == exclude[0] and resseq == int(exclude[1]):
                continue
            best = None
            best_atom = ""
            best_marker = ""
            for atom in residue.get_atoms():
                if not is_heavy_atom(atom):
                    continue
                atom_coord = atom.coord
                for marker_name, marker_coord in marker_coords.items():
                    dist = float(math.dist(atom_coord, marker_coord))
                    if best is None or dist < best:
                        best = dist
                        best_atom = atom.get_name().strip()
                        best_marker = marker_name
            if best is not None and best <= cutoff:
                resname = residue.get_resname().strip().upper()
                rows.append({
                    "chain": chain.id,
                    "resname": resname,
                    "resseq": resseq,
                    "residue": base.residue_label(chain.id, resname, resseq),
                    "contact_atom": best_atom,
                    "marker_atom": best_marker,
                    "distance": best,
                })
    rows.sort(key=lambda row: (row["distance"], row["chain"], row["resseq"]))
    return rows


def basic_donor_contacts(
    model,
    marker_coords: Dict[str, object],
    exclude: Optional[Tuple[str, int]],
    cutoff: float = BASIC_CUTOFF_A,
) -> List[Dict]:
    if not marker_coords:
        return []
    rows = []
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip().upper()
            donors = base.BASIC_DONORS.get(resname)
            if not donors:
                continue
            resseq = int(residue.id[1])
            if exclude and chain.id == exclude[0] and resseq == int(exclude[1]):
                continue
            best = None
            best_atom = ""
            best_marker = ""
            for atom in residue.get_atoms():
                atom_name = atom.get_name().strip().upper()
                if atom_name not in donors:
                    continue
                for marker_name, marker_coord in marker_coords.items():
                    dist = float(math.dist(atom.coord, marker_coord))
                    if best is None or dist < best:
                        best = dist
                        best_atom = atom_name
                        best_marker = marker_name
            if best is not None and best <= cutoff:
                rows.append({
                    "chain": chain.id,
                    "resname": resname,
                    "resseq": resseq,
                    "residue": base.residue_label(chain.id, resname, resseq),
                    "donor_atom": best_atom,
                    "marker_atom": best_marker,
                    "distance": best,
                })
    rows.sort(key=lambda row: (row["distance"], row["chain"], row["resseq"]))
    return rows


def basic_atomic_contacts(
    model,
    marker_coords: Dict[str, object],
    exclude: Optional[Tuple[str, int]],
    cutoff: float = BASIC_CUTOFF_A,
) -> List[Dict]:
    oxygen_coords = {name: marker_coords[name] for name in PHOSPHATE_OXYGEN_NAMES if name in marker_coords}
    if not oxygen_coords:
        return []
    rows = []
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip().upper()
            donors = base.BASIC_DONORS.get(resname)
            if not donors:
                continue
            resseq = int(residue.id[1])
            if exclude and chain.id == exclude[0] and resseq == int(exclude[1]):
                continue
            for atom in residue.get_atoms():
                atom_name = atom.get_name().strip().upper()
                if atom_name not in donors:
                    continue
                distances = []
                for marker_name, marker_coord in oxygen_coords.items():
                    dist = float(math.dist(atom.coord, marker_coord))
                    distances.append((marker_name, dist))
                within = [(name, dist) for name, dist in distances if dist <= cutoff]
                if within:
                    nearest_name, nearest_dist = min(within, key=lambda item: item[1])
                    residue_label = base.residue_label(chain.id, resname, resseq)
                    rows.append({
                        "chain": chain.id,
                        "resname": resname,
                        "resseq": resseq,
                        "residue": residue_label,
                        "donor_atom": atom_name,
                        "marker_atom": "Ophos",
                        "nearest_marker_atom": nearest_name,
                        "distance": nearest_dist,
                        "oxygen_contacts": ";".join(f"{name}:{dist:.2f}" for name, dist in sorted(within)),
                        "atomic_contact": f"{residue_label}:{atom_name}-Ophos",
                    })
    rows.sort(key=lambda row: (row["distance"], row["chain"], row["resseq"], row["donor_atom"]))
    return rows


def rows_to_tuples(rows: Sequence[Dict]) -> List[Tuple[str, str, int]]:
    tuples = []
    seen = set()
    for row in rows:
        key = (row["chain"], row["resname"], int(row["resseq"]))
        if key not in seen:
            seen.add(key)
            tuples.append(key)
    return tuples


def labels_from_tuples(rows: Sequence[Tuple[str, str, int]]) -> List[str]:
    return [base.residue_label(*row) for row in rows]


def labels_from_rows(rows: Sequence[Dict]) -> List[str]:
    return [row["residue"] for row in rows]


def atomic_labels_from_rows(rows: Sequence[Dict]) -> List[str]:
    return [row["atomic_contact"] for row in rows]


def mapped_atomic_labels_from_rows(rows: Sequence[Dict], reference_model, tool: str) -> List[str]:
    if tool != "PTM-Psi":
        return atomic_labels_from_rows(rows)
    labels = []
    for row in rows:
        local = base.chain_local_index(reference_model, row["chain"], int(row["resseq"]))
        if local is None:
            continue
        residue_label = base.residue_label(row["chain"], row["resname"], local)
        labels.append(f"{residue_label}:{row['donor_atom']}-{row['marker_atom']}")
    return labels


def join_labels(labels: Iterable[str]) -> str:
    return ";".join(str(label) for label in labels)


def residue_selection(contact_tuples: Sequence[Tuple[str, str, int]], object_name: str) -> str:
    by_chain: Dict[str, List[int]] = {}
    for chain_id, _resname, resseq in contact_tuples:
        by_chain.setdefault(chain_id, []).append(int(resseq))
    parts = []
    for chain_id, positions in sorted(by_chain.items()):
        resi = "+".join(str(pos) for pos in sorted(set(positions)))
        parts.append(f"({object_name} and chain {chain_id} and resi {resi})")
    return " or ".join(parts) if parts else "none"


def contact_label_tuples(value) -> List[Tuple[str, str, int]]:
    tuples = []
    for label in contact_items(value):
        match = re.match(r"^([^:]+):([A-Za-z]{3})(-?\d+)", label)
        if not match:
            continue
        tuples.append((match.group(1), match.group(2).upper(), int(match.group(3))))
    return tuples


def basic_donor_selection_from_labels(value, object_name: str) -> str:
    tuples = contact_label_tuples(value)
    by_chain: Dict[str, List[int]] = {}
    for chain_id, _resname, resseq in tuples:
        by_chain.setdefault(chain_id, []).append(int(resseq))
    parts = []
    for chain_id, positions in sorted(by_chain.items()):
        resi = "+".join(str(pos) for pos in sorted(set(positions)))
        parts.append(
            f"({object_name} and chain {chain_id} and resi {resi} and "
            "polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2)"
        )
    return " or ".join(parts) if parts else "none"


def atomic_contact_tuples(value) -> List[Tuple[str, str, int, str, str]]:
    tuples = []
    for label in contact_items(value):
        match = re.match(r"^([^:]+):([A-Za-z]{3})(-?\d+):([A-Za-z0-9']+)-([A-Za-z0-9']+)$", label)
        if not match:
            continue
        tuples.append((
            match.group(1),
            match.group(2).upper(),
            int(match.group(3)),
            match.group(4).upper(),
            match.group(5).upper(),
        ))
    return tuples


def phosphate_atom_selection_name(atom_name: str) -> str:
    atom_name = str(atom_name).strip().upper()
    return {
        "O1P": "O1P+OP1",
        "O2P": "O2P+OP2",
        "O3P": "O3P+OP3",
        "OPHOS": PHOSPHATE_OXYGEN_SELECTION,
    }.get(atom_name, atom_name)


def phosphosite_bond_cleanup_pml(site_selection: str, residue_type: str) -> List[str]:
    linkage_atom = PHOSPHO_LINKAGE_ATOMS.get(str(residue_type).upper())
    lines = [
        f"unbond ({site_selection} and name {PHOSPHATE_SELECTION}), ({site_selection})",
    ]
    if linkage_atom:
        lines.append(f"bond ({site_selection} and name {linkage_atom}), ({site_selection} and name P)")
    for oxygen_names in ("O1P+OP1", "O2P+OP2", "O3P+OP3"):
        lines.append(f"bond ({site_selection} and name P), ({site_selection} and name {oxygen_names})")
    return lines


def colored_basic_contact_pml(
    object_name: str,
    phosphate_selection: str,
    maintained_labels,
    new_labels,
    lost_labels,
    lost_cutoff: float = 12.0,
) -> List[str]:
    lines = []
    classes = [
        ("maintained", maintained_labels, "contact_maintained_gray", 4.0),
        ("new", new_labels, "contact_new_blue", 4.0),
        ("lost", lost_labels, "contact_lost_red", lost_cutoff),
    ]
    for suffix, labels, color_name, cutoff in classes:
        donor_selection = basic_donor_selection_from_labels(labels, object_name)
        if donor_selection == "none":
            lines.append(f"select {object_name}_{suffix}_basic_residues, none")
            lines.append(f"select {object_name}_{suffix}_basic_contacts, none")
            continue
        dist_name = f"{object_name}_{suffix}_basic_contacts"
        lines.extend([
            f"select {object_name}_{suffix}_basic_residues, byres ({donor_selection})",
            f"show sticks, {object_name}_{suffix}_basic_residues",
            f"color {color_name}, {object_name}_{suffix}_basic_residues",
            f"distance {dist_name}, {phosphate_selection}, ({donor_selection}), {cutoff:.1f}",
            f"color {color_name}, {dist_name}",
            f"set dash_color, {color_name}, {dist_name}",
            f"set dash_width, 2.2, {dist_name}",
            f"hide labels, {dist_name}",
        ])
    return lines


def colored_basic_atomic_contact_pml(
    object_name: str,
    phosphate_selection: str,
    maintained_labels,
    new_labels,
    lost_labels,
    lost_cutoff: float = 12.0,
) -> List[str]:
    lines = []
    classes = [
        ("maintained", maintained_labels, "contact_maintained_gray", 4.0),
        ("new", new_labels, "contact_new_blue", 4.0),
        ("lost", lost_labels, "contact_lost_red", lost_cutoff),
    ]
    for suffix, labels, color_name, cutoff in classes:
        tuples = atomic_contact_tuples(labels)
        residue_parts = []
        distance_members = []
        if not tuples:
            lines.extend([
                f"select {object_name}_{suffix}_basic_atomic_residues, none",
                f"select {object_name}_{suffix}_basic_atomic_atoms, none",
                f"select {object_name}_{suffix}_basic_atomic_contacts, none",
            ])
            continue
        for idx, (chain_id, _resname, resseq, donor_atom, phosphate_atom) in enumerate(tuples, start=1):
            residue_sel = (
                f"({object_name} and chain {chain_id} and resi {resseq} and polymer.protein "
                f"and resn ARG+LYS+HIS and name {donor_atom})"
            )
            phosphate_sel = f"({phosphate_selection} and name {phosphate_atom_selection_name(phosphate_atom)})"
            donor_sel_name = f"{object_name}_{suffix}_basic_atomic_donor_{idx}"
            phosphate_sel_name = f"{object_name}_{suffix}_basic_atomic_phos_{idx}"
            distance_name = f"{object_name}_{suffix}_basic_atomic_contact_{idx}"
            residue_parts.append(residue_sel)
            distance_members.append(distance_name)
            lines.extend([
                f"select {donor_sel_name}, {residue_sel}",
                f"select {phosphate_sel_name}, {phosphate_sel}",
                f"distance {distance_name}, {phosphate_sel_name}, {donor_sel_name}, {cutoff:.1f}",
                f"color {color_name}, {distance_name}",
                f"set dash_color, {color_name}, {distance_name}",
                f"set dash_width, 2.0, {distance_name}",
                f"hide labels, {distance_name}",
                f"hide dashes, {distance_name}",
            ])
        residue_selection_text = " or ".join(residue_parts)
        lines.extend([
            f"select {object_name}_{suffix}_basic_atomic_residues, byres ({residue_selection_text})",
            f"select {object_name}_{suffix}_basic_atomic_atoms, ({residue_selection_text}) or ({phosphate_selection})",
            f"group {object_name}_{suffix}_basic_atomic_contacts, {' '.join(distance_members)}",
        ])
    return lines


def find_site_residue(model, chain_hint: str, position_hint: int, residue_type: str):
    exact = base.find_residue(model, chain_hint, position_hint)
    if exact is not None and "P" in marker_coords_for_phospho(exact):
        return exact, chain_hint, int(position_hint)
    candidates = []
    for chain in model:
        for residue in chain:
            if residue.get_resname().strip().upper() == residue_type and "P" in marker_coords_for_phospho(residue):
                candidates.append((residue, chain.id, int(residue.id[1])))
    if not candidates:
        return exact, chain_hint, int(position_hint)
    return min(candidates, key=lambda item: abs(item[2] - int(position_hint)))


def compare_sets(predicted: Sequence[str], crystal: Sequence[str]) -> Dict:
    pred_set = set(predicted)
    crystal_set = set(crystal)
    recovered = sorted(crystal_set & pred_set)
    missing = sorted(crystal_set - pred_set)
    extra = sorted(pred_set - crystal_set)
    recall = len(recovered) / len(crystal_set) if crystal_set else float("nan")
    precision = len(recovered) / len(pred_set) if pred_set else float("nan")
    return {
        "recovered_crystal_contact_count_5A": len(recovered),
        "missing_crystal_contact_count_5A": len(missing),
        "extra_predicted_contact_count_5A": len(extra),
        "contact_recall_5A": recall,
        "contact_precision_5A": precision,
        "recovered_crystal_contacts_5A": join_labels(recovered),
        "missing_crystal_contacts_5A": join_labels(missing),
        "extra_predicted_contacts_5A": join_labels(extra),
    }


def compare_basic_atomic_sets(predicted: Sequence[str], crystal: Sequence[str]) -> Dict:
    comparison = compare_sets(predicted, crystal)
    return {
        "recovered_crystal_basic_atomic_contact_count_4A": comparison["recovered_crystal_contact_count_5A"],
        "missing_crystal_basic_atomic_contact_count_4A": comparison["missing_crystal_contact_count_5A"],
        "extra_predicted_basic_atomic_contact_count_4A": comparison["extra_predicted_contact_count_5A"],
        "basic_atomic_contact_recall_4A": comparison["contact_recall_5A"],
        "basic_atomic_contact_precision_4A": comparison["contact_precision_5A"],
        "recovered_crystal_basic_atomic_contacts_4A": comparison["recovered_crystal_contacts_5A"],
        "missing_crystal_basic_atomic_contacts_4A": comparison["missing_crystal_contacts_5A"],
        "extra_predicted_basic_atomic_contacts_4A": comparison["extra_predicted_contacts_5A"],
    }


def summarize_environment(model, residue, chain_id: str, position: int, marker_coords: Dict[str, object]) -> Dict:
    contacts = residue_contacts(model, marker_coords, (chain_id, position), CONTACT_CUTOFF_A)
    basic = basic_donor_contacts(model, marker_coords, (chain_id, position), BASIC_CUTOFF_A)
    basic_atomic = basic_atomic_contacts(model, marker_coords, (chain_id, position), BASIC_CUTOFF_A)
    contact_tuples = rows_to_tuples(contacts)
    basic_tuples = rows_to_tuples(basic)
    return {
        "site_resname": residue.get_resname().strip().upper() if residue is not None else "",
        "contact_rows": contacts,
        "basic_rows": basic,
        "basic_atomic_rows": basic_atomic,
        "contact_tuples": contact_tuples,
        "basic_tuples": basic_tuples,
        "contact_labels": labels_from_tuples(contact_tuples),
        "basic_labels": labels_from_tuples(basic_tuples),
        "basic_atomic_labels": atomic_labels_from_rows(basic_atomic),
        "contact_count": len(contact_tuples),
        "basic_contact_count": len(basic_tuples),
        "basic_atomic_contact_count": len(basic_atomic),
        "contact_set": join_labels(labels_from_tuples(contact_tuples)),
        "basic_contact_set": join_labels(labels_from_tuples(basic_tuples)),
        "basic_atomic_contact_set": join_labels(atomic_labels_from_rows(basic_atomic)),
    }


def add_contact_detail_rows(out_rows: List[Dict], case: Dict, pdb_id: str, env_type: str, structure_label: str,
                            variant: str, rows: Sequence[Dict], structure_path: str):
    for row in rows:
        out_rows.append({
            "protein": case["protein"],
            "site_label": case["site_label"],
            "case_id": base.case_key(case, pdb_id),
            "pdb_id": pdb_id,
            "environment_type": env_type,
            "structure_label": structure_label,
            "variant": variant,
            "contact_cutoff_A": CONTACT_CUTOFF_A,
            "contact_residue": row["residue"],
            "contact_chain": row["chain"],
            "contact_resname": row["resname"],
            "contact_resseq": row["resseq"],
            "contact_atom": row["contact_atom"],
            "marker_atom": row["marker_atom"],
            "distance": row["distance"],
            "structure_path": structure_path,
        })


def add_basic_atomic_detail_rows(out_rows: List[Dict], case: Dict, pdb_id: str, env_type: str, structure_label: str,
                                 variant: str, rows: Sequence[Dict], structure_path: str):
    for row in rows:
        out_rows.append({
            "protein": case["protein"],
            "site_label": case["site_label"],
            "case_id": base.case_key(case, pdb_id),
            "pdb_id": pdb_id,
            "environment_type": env_type,
            "structure_label": structure_label,
            "variant": variant,
            "cutoff_A": BASIC_CUTOFF_A,
            "contact_residue": row["residue"],
            "contact_chain": row["chain"],
            "contact_resname": row["resname"],
            "contact_resseq": row["resseq"],
            "donor_atom": row["donor_atom"],
            "phosphate_atom": row["marker_atom"],
            "nearest_phosphate_oxygen": row.get("nearest_marker_atom", row["marker_atom"]),
            "all_phosphate_oxygen_contacts": row.get("oxygen_contacts", ""),
            "atomic_contact": row["atomic_contact"],
            "distance": row["distance"],
            "structure_path": structure_path,
        })


def build_environment_tables(metrics: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict, Dict]:
    unmodified_by_site = {}
    phospho_by_case = {}
    summary_rows = []
    contact_detail_rows = []
    basic_atomic_detail_rows = []
    unmodified_rows = []

    for case in base.CASES:
        site_key = (case["protein"], case["residue_type"], int(case["position"]))
        unmod = UNMODIFIED_SITE_MAP[site_key]
        unmod_path = base.PDB_CACHE / f"{unmod['pdb_id']}.cif"
        unmod_model = base.first_model(unmod_path)
        unmod_res = base.find_residue(unmod_model, unmod["chain"], unmod["position"])
        unmod_marker = marker_coords_for_unmodified(unmod_res, case["residue_type"])
        unmod_summary = summarize_environment(unmod_model, unmod_res, unmod["chain"], unmod["position"], unmod_marker)
        unmod_summary.update({
            "pdb_id": unmod["pdb_id"],
            "chain": unmod["chain"],
            "position": unmod["position"],
            "requested_position": unmod["requested_position"],
            "path": str(unmod_path),
            "note": unmod["note"],
        })
        unmodified_by_site[site_key] = unmod_summary
        unmodified_rows.append({
            "protein": case["protein"],
            "site_label": case["site_label"],
            "residue_type": case["residue_type"],
            "unmodified_pdb_id": unmod["pdb_id"],
            "unmodified_chain": unmod["chain"],
            "unmodified_author_position_used": unmod["position"],
            "requested_unmodified_position": unmod["requested_position"],
            "unmodified_site_resname": unmod_summary["site_resname"],
            "unmodified_contact_count_5A": unmod_summary["contact_count"],
            "unmodified_basic_contact_count_4A": unmod_summary["basic_contact_count"],
            "unmodified_basic_atomic_contact_count_4A": unmod_summary["basic_atomic_contact_count"],
            "unmodified_contact_set_5A": unmod_summary["contact_set"],
            "unmodified_basic_contact_set_4A": unmod_summary["basic_contact_set"],
            "unmodified_basic_atomic_contact_set_4A": unmod_summary["basic_atomic_contact_set"],
            "unmodified_structure_path": str(unmod_path),
            "numbering_note": unmod["note"],
        })

    for case in base.CASES:
        site_key = (case["protein"], case["residue_type"], int(case["position"]))
        unmod_summary = unmodified_by_site[site_key]
        for pdb_id in case["pdb_ids"]:
            ref_path = base.PDB_CACHE / f"{pdb_id}.cif"
            ref_model = base.first_model(ref_path)
            ref_res = base.find_residue(ref_model, case["chain"], case["position"])
            ref_marker = marker_coords_for_phospho(ref_res)
            phospho_summary = summarize_environment(ref_model, ref_res, case["chain"], case["position"], ref_marker)
            phospho_summary.update({"path": str(ref_path), "model": ref_model})
            phospho_by_case[(case["site_label"], pdb_id)] = phospho_summary
            add_contact_detail_rows(contact_detail_rows, case, pdb_id, "phospho_crystal",
                                    f"{pdb_id}_phospho_crystal", "Reference", phospho_summary["contact_rows"], str(ref_path))
            add_basic_atomic_detail_rows(basic_atomic_detail_rows, case, pdb_id, "phospho_crystal",
                                         f"{pdb_id}_phospho_crystal", "Reference",
                                         phospho_summary["basic_atomic_rows"], str(ref_path))

            add_contact_detail_rows(contact_detail_rows, case, pdb_id, "unmodified_crystal",
                                    f"{unmod_summary['pdb_id']}_unmodified", "Unmodified", unmod_summary["contact_rows"], unmod_summary["path"])
            add_basic_atomic_detail_rows(basic_atomic_detail_rows, case, pdb_id, "unmodified_crystal",
                                         f"{unmod_summary['pdb_id']}_unmodified", "Unmodified",
                                         unmod_summary["basic_atomic_rows"], unmod_summary["path"])

            sub = metrics[(metrics["site_label"] == case["site_label"]) & (metrics["pdb_id"].str.lower() == pdb_id.lower())].copy()
            for _, metric in sub.iterrows():
                structure_path = base.resolve_path(metric.get("structure_path", ""))
                pred_contacts = []
                pred_basic = []
                pred_basic_atomic = []
                pred_site_chain = str(metric.get("chain", case["chain"]))
                pred_site_pos = int(metric.get("position", case["position"]))
                notes = []
                if structure_path is None or not structure_path.exists():
                    notes.append("missing_prediction_structure")
                else:
                    pred_model = base.first_model(structure_path)
                    pred_res, pred_site_chain, pred_site_pos = find_site_residue(
                        pred_model, pred_site_chain, pred_site_pos, case["residue_type"]
                    )
                    pred_marker = marker_coords_for_phospho(pred_res)
                    pred_summary = summarize_environment(pred_model, pred_res, pred_site_chain, pred_site_pos, pred_marker)
                    pred_contacts = pred_summary["contact_rows"]
                    pred_basic = pred_summary["basic_rows"]
                    pred_basic_atomic = pred_summary["basic_atomic_rows"]
                    add_contact_detail_rows(contact_detail_rows, case, pdb_id, "prediction",
                                            str(metric["variant"]), str(metric["variant"]), pred_contacts, str(structure_path))
                    add_basic_atomic_detail_rows(basic_atomic_detail_rows, case, pdb_id, "prediction",
                                                 str(metric["variant"]), str(metric["variant"]),
                                                 pred_basic_atomic, str(structure_path))

                pred_labels = labels_from_rows(pred_contacts)
                mapped_crystal_tuples = base.map_contacts_for_output(
                    phospho_summary["contact_tuples"], ref_model, str(metric["tool"])
                )
                mapped_crystal_labels = labels_from_tuples(mapped_crystal_tuples)
                comparison = compare_sets(pred_labels, mapped_crystal_labels)
                pred_basic_labels = labels_from_rows(pred_basic)
                mapped_crystal_basic_tuples = base.map_contacts_for_output(
                    phospho_summary["basic_tuples"], ref_model, str(metric["tool"])
                )
                mapped_crystal_basic_labels = labels_from_tuples(mapped_crystal_basic_tuples)
                basic_comparison = compare_sets(pred_basic_labels, mapped_crystal_basic_labels)
                pred_basic_atomic_labels = atomic_labels_from_rows(pred_basic_atomic)
                mapped_crystal_basic_atomic_labels = mapped_atomic_labels_from_rows(
                    phospho_summary["basic_atomic_rows"], ref_model, str(metric["tool"])
                )
                basic_atomic_comparison = compare_basic_atomic_sets(
                    pred_basic_atomic_labels, mapped_crystal_basic_atomic_labels
                )
                crystal_basic_atomic_count = phospho_summary["basic_atomic_contact_count"]
                pred_basic_atomic_count = len(pred_basic_atomic)
                summary_rows.append({
                    "protein": case["protein"],
                    "site_label": case["site_label"],
                    "case_id": base.case_key(case, pdb_id),
                    "pdb_id": pdb_id,
                    "chain": case["chain"],
                    "position": case["position"],
                    "residue_type": case["residue_type"],
                    "tool": metric["tool"],
                    "variant": metric["variant"],
                    "unmodified_pdb_id": unmod_summary["pdb_id"],
                    "unmodified_author_position_used": unmod_summary["position"],
                    "requested_unmodified_position": unmod_summary["requested_position"],
                    "unmodified_site_resname": unmod_summary["site_resname"],
                    "unmodified_contact_count_5A": unmod_summary["contact_count"],
                    "unmodified_basic_contact_count_4A": unmod_summary["basic_contact_count"],
                    "unmodified_basic_atomic_contact_count_4A": unmod_summary["basic_atomic_contact_count"],
                    "phospho_crystal_contact_count_5A": phospho_summary["contact_count"],
                    "phospho_crystal_basic_contact_count_4A": phospho_summary["basic_contact_count"],
                    "phospho_crystal_basic_atomic_contact_count_4A": crystal_basic_atomic_count,
                    "predicted_contact_count_5A": len(rows_to_tuples(pred_contacts)),
                    "predicted_basic_contact_count_4A": len(rows_to_tuples(pred_basic)),
                    "predicted_basic_atomic_contact_count_4A": pred_basic_atomic_count,
                    "basic_atomic_contact_delta_4A": pred_basic_atomic_count - crystal_basic_atomic_count,
                    "predicted_site_chain": pred_site_chain,
                    "predicted_site_position": pred_site_pos,
                    "unmodified_contact_set_5A": unmod_summary["contact_set"],
                    "phospho_crystal_contact_set_5A": join_labels(mapped_crystal_labels),
                    "predicted_contact_set_5A": join_labels(pred_labels),
                    "unmodified_basic_contact_set_4A": unmod_summary["basic_contact_set"],
                    "unmodified_basic_atomic_contact_set_4A": unmod_summary["basic_atomic_contact_set"],
                    "phospho_crystal_basic_contact_set_4A": phospho_summary["basic_contact_set"],
                    "phospho_crystal_basic_contact_set_4A_mapped": join_labels(mapped_crystal_basic_labels),
                    "phospho_crystal_basic_atomic_contact_set_4A": phospho_summary["basic_atomic_contact_set"],
                    "phospho_crystal_basic_atomic_contact_set_4A_mapped": join_labels(mapped_crystal_basic_atomic_labels),
                    "predicted_basic_contact_set_4A": join_labels(pred_basic_labels),
                    "predicted_basic_atomic_contact_set_4A": join_labels(pred_basic_atomic_labels),
                    "recovered_crystal_basic_contact_count_4A": basic_comparison["recovered_crystal_contact_count_5A"],
                    "missing_crystal_basic_contact_count_4A": basic_comparison["missing_crystal_contact_count_5A"],
                    "extra_predicted_basic_contact_count_4A": basic_comparison["extra_predicted_contact_count_5A"],
                    "basic_contact_recall_4A": basic_comparison["contact_recall_5A"],
                    "basic_contact_precision_4A": basic_comparison["contact_precision_5A"],
                    "recovered_crystal_basic_contacts_4A": basic_comparison["recovered_crystal_contacts_5A"],
                    "missing_crystal_basic_contacts_4A": basic_comparison["missing_crystal_contacts_5A"],
                    "extra_predicted_basic_contacts_4A": basic_comparison["extra_predicted_contacts_5A"],
                    **basic_atomic_comparison,
                    "structure_path": metric.get("structure_path", ""),
                    "unmodified_structure_path": unmod_summary["path"],
                    "phospho_crystal_structure_path": str(ref_path),
                    "numbering_note": unmod_summary["note"],
                    "notes": ";".join(notes),
                    **comparison,
                })

    return (
        pd.DataFrame(summary_rows),
        pd.DataFrame(contact_detail_rows),
        pd.DataFrame(basic_atomic_detail_rows),
        pd.DataFrame(unmodified_rows),
        unmodified_by_site,
        phospho_by_case,
    )


def draw_grouped_bar_png(df: pd.DataFrame, path: Path, variants: Sequence[str], title: str):
    work = df[df["variant"].isin(variants)].copy()
    case_ids = list(dict.fromkeys(work["case_id"]))
    max_y = max(0.5, float(pd.to_numeric(work["rmsd"], errors="coerce").max()) * 1.15)
    width = max(1100, 150 + len(case_ids) * 145)
    height = 640
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    left, top, plot_h = 90, 70, 390
    plot_w = width - 170
    group_w = plot_w / max(1, len(case_ids))
    bar_w = min(20, group_w / (len(variants) + 2))
    text(draw, (width / 2, 32), title, size=22, bold=True)
    draw.line((left, top, left, top + plot_h), fill="#333333", width=2)
    draw.line((left, top + plot_h, left + plot_w, top + plot_h), fill="#333333", width=2)
    for tick in [max_y * i / 5 for i in range(6)]:
        y = top + plot_h - (tick / max_y) * plot_h
        draw.line((left - 5, y, left + plot_w, y), fill="#e7e7e7")
        text(draw, (left - 12, y), f"{tick:.1f}", size=12, anchor="rm")
    for i, cid in enumerate(case_ids):
        cx = left + group_w * i + group_w / 2
        sub = work[work["case_id"] == cid]
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            if row.empty:
                continue
            value = base.safe_float(row.iloc[0]["rmsd"])
            if math.isnan(value):
                continue
            x = cx - (len(variants) * bar_w) / 2 + j * bar_w
            h = (value / max_y) * plot_h
            y = top + plot_h - h
            draw.rectangle((x, y, x + bar_w - 3, top + plot_h), fill=base.COLORS[variant])
        for k, part in enumerate(cid.split(" ")):
            text(draw, (cx, top + plot_h + 24 + k * 14), part, size=11)
    for i, variant in enumerate(variants):
        x = left + 10 + i * 210
        y = height - 92
        draw.rectangle((x, y, x + 16, y + 16), fill=base.COLORS[variant])
        text(draw, (x + 24, y + 8), variant, size=13, anchor="lm")
    text(draw, (27, top + plot_h / 2), "RMSD (A)", size=14)
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def draw_torsion_png(df: pd.DataFrame, path: Path, variants: Optional[Sequence[str]] = None,
                     title: str = "Predicted phosphate torsion vs crystal"):
    variants = list(variants or base.PRIMARY_VARIANTS)
    work = df[df["variant"].isin(variants)].copy()
    case_ids = list(dict.fromkeys(work["case_id"]))
    width = max(1100, 150 + len(case_ids) * 145)
    height = 640
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    left, top, plot_h = 90, 70, 390
    plot_w = width - 170
    group_w = plot_w / max(1, len(case_ids))

    def y_for(angle):
        return top + plot_h - ((float(angle) + 180.0) / 360.0) * plot_h

    text(draw, (width / 2, 32), title, size=22, bold=True)
    draw.line((left, top, left, top + plot_h), fill="#333333", width=2)
    draw.line((left, top + plot_h, left + plot_w, top + plot_h), fill="#333333", width=2)
    for tick in [-180, -120, -60, 0, 60, 120, 180]:
        y = y_for(tick)
        draw.line((left - 5, y, left + plot_w, y), fill="#e7e7e7")
        text(draw, (left - 12, y), str(tick), size=12, anchor="rm")
    for i, cid in enumerate(case_ids):
        cx = left + group_w * i + group_w / 2
        sub = work[work["case_id"] == cid]
        ref = base.safe_float(sub.iloc[0]["torsion_ref"]) if not sub.empty else float("nan")
        if not math.isnan(ref):
            y = y_for(ref)
            draw.line((cx - group_w * 0.34, y, cx + group_w * 0.34, y), fill="#222222", width=3)
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            if row.empty:
                continue
            value = base.safe_float(row.iloc[0]["torsion_pred"])
            if math.isnan(value):
                continue
            x = cx - 42 + j * 28
            y = y_for(value)
            draw.ellipse((x - 6, y - 6, x + 6, y + 6), fill=base.COLORS[variant])
        for k, part in enumerate(cid.split(" ")):
            text(draw, (cx, top + plot_h + 24 + k * 14), part, size=11)
    legend = ["Reference"] + list(variants)
    for i, name in enumerate(legend):
        x = left + 10 + i * 180
        y = height - 92
        if name == "Reference":
            draw.line((x, y + 8, x + 18, y + 8), fill="#222222", width=3)
        else:
            draw.ellipse((x + 3, y + 3, x + 15, y + 15), fill=base.COLORS[name])
        text(draw, (x + 24, y + 8), name, size=13, anchor="lm")
    text(draw, (34, top + plot_h / 2), "Torsion angle", size=14)
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def _rgb(hex_color: str) -> Tuple[int, int, int]:
    color = str(hex_color).lstrip("#")
    return tuple(int(color[i:i + 2], 16) for i in (0, 2, 4))


def _mix(hex_color: str, other: str = "#ffffff", weight: float = 0.55) -> Tuple[int, int, int]:
    a = _rgb(hex_color)
    b = _rgb(other)
    return tuple(int(a[i] * (1.0 - weight) + b[i] * weight) for i in range(3))


def _angle_point(cx: float, cy: float, radius: float, angle: float) -> Tuple[float, float]:
    theta = math.radians(float(angle) - 90.0)
    return cx + radius * math.cos(theta), cy + radius * math.sin(theta)


def _signed_circular_delta(target: float, reference: float) -> float:
    return ((float(target) - float(reference) + 180.0) % 360.0) - 180.0


def _arc_points(cx: float, cy: float, radius: float, reference: float, target: float) -> List[Tuple[float, float]]:
    delta = _signed_circular_delta(target, reference)
    steps = max(8, int(abs(delta) / 6.0) + 1)
    return [_angle_point(cx, cy, radius, float(reference) + delta * i / steps) for i in range(steps + 1)]


def _torsion_case_label(case_id: str) -> List[str]:
    parts = str(case_id).split()
    if "AF2" in parts:
        idx = parts.index("AF2")
        return [" ".join(parts[:idx]), " ".join(parts[idx:])]
    if len(parts) >= 4:
        return [" ".join(parts[:-1]), parts[-1]]
    return [str(case_id)]


def draw_torsion_polar_png(df: pd.DataFrame, path: Path, variants: Optional[Sequence[str]] = None,
                           title: str = "Circular torsion mode comparison"):
    variants = list(variants or base.PRIMARY_VARIANTS)
    work = df[df["variant"].isin(variants)].copy()
    case_ids = list(dict.fromkeys(work["case_id"]))
    if not case_ids:
        return

    cols = min(3, len(case_ids))
    rows = int(math.ceil(len(case_ids) / cols))
    panel_w, panel_h = 360, 305
    left, top = 28, 86
    width = max(1100, left * 2 + cols * panel_w)
    legend_h = 105 if len(variants) <= 5 else 130
    height = top + rows * panel_h + legend_h
    scale = 4
    image = Image.new("RGB", (width * scale, height * scale), "white")
    raw_draw = ImageDraw.Draw(image)

    class ScaledDraw:
        def __init__(self, inner: ImageDraw.ImageDraw, factor: int):
            self.inner = inner
            self.factor = factor

        def _coords(self, coords):
            if isinstance(coords, list):
                return [(x * self.factor, y * self.factor) for x, y in coords]
            return tuple(value * self.factor for value in coords)

        def line(self, coords, fill=None, width=1, **kwargs):
            self.inner.line(self._coords(coords), fill=fill, width=max(1, int(round(width * self.factor))), **kwargs)

        def ellipse(self, coords, fill=None, outline=None, width=1, **kwargs):
            self.inner.ellipse(
                self._coords(coords),
                fill=fill,
                outline=outline,
                width=max(1, int(round(width * self.factor))),
                **kwargs,
            )

    draw = ScaledDraw(raw_draw, scale)

    def draw_label(xy: Tuple[float, float], value, size=12, fill="#222222", anchor="mm", bold=False):
        raw_draw.text(
            (xy[0] * scale, xy[1] * scale),
            str(value),
            font=font(int(size * scale), bold),
            fill=fill,
            anchor=anchor,
        )

    draw_label((width / 2, 30), title, size=22, bold=True)
    draw_label((width / 2, 58), "Angle is circular: 0 at top, +/-180 at bottom; radial offsets only separate tools.",
               size=12, fill="#555555")

    for idx, cid in enumerate(case_ids):
        col = idx % cols
        row_idx = idx // cols
        x0 = left + col * panel_w
        y0 = top + row_idx * panel_h
        cx = x0 + panel_w / 2
        cy = y0 + 122
        outer_r = 92
        inner_r = 42
        point_outer_r = 84

        for radius in (inner_r, 62, point_outer_r, outer_r):
            draw.ellipse((cx - radius, cy - radius, cx + radius, cy + radius), outline="#e4e4e4", width=1)
        for tick in (-180, -90, 0, 90, 180):
            x, y = _angle_point(cx, cy, outer_r, tick)
            draw.line((cx, cy, x, y), fill="#ededed", width=1)
            tx, ty = _angle_point(cx, cy, outer_r + 18, tick)
            tick_label = "+/-180" if tick in (-180, 180) else str(tick)
            if tick == -180:
                continue
            draw_label((tx, ty), tick_label, size=10, fill="#555555")
        draw.ellipse((cx - outer_r, cy - outer_r, cx + outer_r, cy + outer_r), outline="#333333", width=2)

        sub = work[work["case_id"] == cid]
        ref = base.safe_float(sub.iloc[0]["torsion_ref"]) if not sub.empty else float("nan")
        if not math.isnan(ref):
            rx, ry = _angle_point(cx, cy, outer_r + 4, ref)
            draw.line((cx, cy, rx, ry), fill=base.COLORS["Reference"], width=3)
            draw.ellipse((rx - 4, ry - 4, rx + 4, ry + 4), fill=base.COLORS["Reference"])

        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            if row.empty:
                continue
            value = base.safe_float(row.iloc[0]["torsion_pred"])
            if math.isnan(value):
                continue
            radius = point_outer_r if len(variants) == 1 else inner_r + j * (point_outer_r - inner_r) / (len(variants) - 1)
            color = base.COLORS.get(variant, "#777777")
            if not math.isnan(ref):
                points = _arc_points(cx, cy, radius, ref, value)
                if len(points) > 1:
                    draw.line(points, fill=_mix(color, "#ffffff", 0.50), width=3)
            px, py = _angle_point(cx, cy, radius, value)
            draw.ellipse((px - 6, py - 6, px + 6, py + 6), fill=color, outline="#222222", width=1)

        for line_idx, case_label in enumerate(_torsion_case_label(cid)):
            draw_label((cx, y0 + 242 + line_idx * 17), case_label, size=12, bold=(line_idx == 0))

    legend_items = ["Reference"] + variants
    per_row = max(1, min(len(legend_items), int((width - 70) // 165)))
    legend_top = height - legend_h + 28
    for i, name in enumerate(legend_items):
        lx = 45 + (i % per_row) * 165
        ly = legend_top + (i // per_row) * 28
        if name == "Reference":
            draw.line((lx, ly + 8, lx + 25, ly + 8), fill=base.COLORS["Reference"], width=3)
            draw.ellipse((lx + 21, ly + 4, lx + 29, ly + 12), fill=base.COLORS["Reference"])
        else:
            color = base.COLORS.get(name, "#777777")
            draw.line((lx, ly + 8, lx + 25, ly + 8), fill=_mix(color, "#ffffff", 0.50), width=3)
            draw.ellipse((lx + 8, ly + 2, lx + 20, ly + 14), fill=color, outline="#222222")
        draw_label((lx + 34, ly + 8), name, size=12, anchor="lm")

    path.parent.mkdir(parents=True, exist_ok=True)
    image = image.resize((width, height), Image.Resampling.LANCZOS)
    image.save(path)


def draw_status_table_png(df: pd.DataFrame, path: Path, value_col: str, dist_col: str, title: str):
    variants = base.PRIMARY_VARIANTS
    cases = list(dict.fromkeys(df["case_id"]))
    row_h = 38
    col_w = 170
    width = 360 + len(variants) * col_w
    height = 105 + len(cases) * row_h
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    text(draw, (width / 2, 32), title, size=22, bold=True)
    x0, y0 = 24, 70
    text(draw, (x0, y0), "Case", size=14, anchor="lm", bold=True)
    for j, variant in enumerate(variants):
        text(draw, (315 + j * col_w, y0), variant, size=12, bold=True)
    for i, cid in enumerate(cases):
        y = y0 + 30 + i * row_h
        text(draw, (x0, y), cid, size=12, anchor="lm")
        sub = df[df["case_id"] == cid]
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            x = 245 + j * col_w
            label = "NA"
            fill = "#f2f2f2"
            if not row.empty:
                val = str(row.iloc[0][value_col]).lower()
                ok = val == "true" or val == "1" or val == "1.0"
                dist = base.safe_float(row.iloc[0].get(dist_col, ""))
                label = ("yes" if ok else "no") + (f" ({dist:.2f})" if not math.isnan(dist) else "")
                fill = "#cfead1" if ok else "#f4c7c3"
            draw.rectangle((x, y - 16, x + col_w - 14, y + 12), fill=fill, outline="#dddddd")
            text(draw, (x + (col_w - 14) / 2, y - 2), label, size=12)
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def draw_environment_png(env: pd.DataFrame, path: Path, variants: Optional[Sequence[str]] = None,
                         title: str = "Crystal environment contact recovery (5 A)"):
    variants = list(variants or base.PRIMARY_VARIANTS)
    work = env[env["variant"].isin(variants)].copy()
    cases = list(dict.fromkeys(work["case_id"]))
    has_bb = has_backbone_rmsd_column(work)
    bb_w = 110 if has_bb else 0
    row_h = 38
    col_w = 170
    width = 520 + bb_w + len(variants) * col_w
    height = 120 + len(cases) * row_h
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    text(draw, (width / 2, 32), title, size=22, bold=True)
    x0, y0 = 24, 72
    headers = [(x0, "Case"), (260 + bb_w, "Unmod"), (360 + bb_w, "Phospho")]
    if has_bb:
        headers.insert(1, (232, "Local BB"))
    for x, header in headers:
        text(draw, (x, y0), header, size=13, anchor="lm", bold=True)
    for j, variant in enumerate(variants):
        text(draw, (500 + bb_w + j * col_w, y0), variant, size=12, bold=True)
    for i, cid in enumerate(cases):
        y = y0 + 30 + i * row_h
        sub = work[work["case_id"] == cid]
        first = sub.iloc[0]
        text(draw, (x0, y), cid, size=12, anchor="lm")
        if has_bb:
            draw.rectangle((225, y - 16, 310, y + 12), fill=backbone_rmsd_fill(first), outline="#dddddd")
            text(draw, (267, y - 2), backbone_rmsd_label(first), size=12)
        text(draw, (278 + bb_w, y), int(first["unmodified_contact_count_5A"]), size=12)
        text(draw, (384 + bb_w, y), int(first["phospho_crystal_contact_count_5A"]), size=12)
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            x = 430 + bb_w + j * col_w
            if row.empty:
                label = "NA"
                fill = "#f2f2f2"
            else:
                r = row.iloc[0]
                recall = base.safe_float(r["contact_recall_5A"])
                recovered = int(r["recovered_crystal_contact_count_5A"])
                total = int(r["phospho_crystal_contact_count_5A"])
                label = f"{recovered}/{total}" + (f" ({recall:.0%})" if not math.isnan(recall) else "")
                fill = "#cfead1" if recall >= 0.75 else "#fff1bd" if recall >= 0.4 else "#f4c7c3"
            draw.rectangle((x, y - 16, x + col_w - 14, y + 12), fill=fill, outline="#dddddd")
            text(draw, (x + (col_w - 14) / 2, y - 2), label, size=12)
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def contact_items(value) -> List[str]:
    if value is None:
        return []
    try:
        if pd.isna(value):
            return []
    except Exception:
        pass
    text_value = str(value).strip()
    if not text_value or text_value.lower() == "nan":
        return []
    return [item.strip() for item in text_value.split(";") if item.strip() and item.strip().lower() != "nan"]


def contact_line(prefix: str, value) -> str:
    items = contact_items(value)
    return f"{prefix}: {', '.join(items) if items else '-'}"


def short_atomic_contact(label: str) -> str:
    match = re.match(r"^([^:]+):([A-Z]{3})(-?\d+):(.+)$", str(label).strip())
    if not match:
        return str(label).strip()
    chain, resname, resseq, atom_pair = match.groups()
    residue_short = {"ARG": "R", "LYS": "K", "HIS": "H"}.get(resname, resname)
    return f"{chain}:{residue_short}{resseq}:{atom_pair}"


def compact_contact_line(prefix: str, value, max_items: int = 3, atomic: bool = False) -> str:
    items = contact_items(value)
    if atomic:
        items = [short_atomic_contact(item) for item in items]
    visible = items[:max_items]
    if len(items) > max_items:
        visible.append(f"+{len(items) - max_items} more")
    return f"{prefix}: {', '.join(visible) if visible else '-'}"


def draw_multiline(draw: ImageDraw.ImageDraw, x: float, y: float, lines: Sequence[str], size: int = 11,
                   fill: str = "#222222", bold: bool = False, line_h: int = 15):
    for i, line in enumerate(lines):
        text(draw, (x, y + i * line_h), line, size=size, fill=fill, anchor="lm", bold=bold)


def draw_basic_recovery_png(env: pd.DataFrame, path: Path, variants: Optional[Sequence[str]] = None,
                            title: str = "Basic residue recovery from crystal contacts (4 A)"):
    variants = list(variants or base.PRIMARY_VARIANTS)
    work = env[env["variant"].isin(variants)].copy()
    cases = list(dict.fromkeys(work["case_id"]))
    has_bb = has_backbone_rmsd_column(work)
    bb_w = 95 if has_bb else 0
    row_h = 96
    case_w = 220
    crystal_w = 250
    col_w = 255
    width = 36 + case_w + bb_w + crystal_w + len(variants) * col_w + 28
    height = 122 + len(cases) * row_h
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    text(draw, (width / 2, 32), title, size=22, bold=True)
    x0, y0 = 24, 76
    text(draw, (x0, y0), "Case", size=13, anchor="lm", bold=True)
    bb_x = x0 + case_w
    crystal_x = x0 + case_w + bb_w
    if has_bb:
        text(draw, (bb_x, y0), "Local BB", size=12, anchor="lm", bold=True)
    text(draw, (crystal_x, y0), "Crystal contacts", size=13, anchor="lm", bold=True)
    for j, variant in enumerate(variants):
        text(draw, (crystal_x + crystal_w + j * col_w + 8, y0), variant, size=12, anchor="lm", bold=True)
    for i, cid in enumerate(cases):
        y = y0 + 34 + i * row_h
        sub = work[work["case_id"] == cid]
        first = sub.iloc[0]
        if i % 2:
            draw.rectangle((10, y - 24, width - 10, y + row_h - 28), fill="#fafafa")
        draw_multiline(draw, x0, y - 8, cid.split(" "), size=12, bold=False, line_h=15)
        if has_bb:
            draw.rectangle((bb_x, y - 22, bb_x + bb_w - 14, y + 16), fill=backbone_rmsd_fill(first), outline="#dddddd")
            text(draw, (bb_x + (bb_w - 14) / 2, y - 3), backbone_rmsd_label(first), size=11)
        total = int(first["phospho_crystal_basic_contact_count_4A"])
        crystal_lines = [f"{total} total", contact_line("C", first["phospho_crystal_basic_contact_set_4A"])]
        draw_multiline(draw, crystal_x, y - 8, crystal_lines, size=11, line_h=15)
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            x = crystal_x + crystal_w + j * col_w
            if row.empty:
                lines = ["NA", "M: -", "X: -"]
                fill = "#f2f2f2"
            else:
                r = row.iloc[0]
                recall = base.safe_float(r["basic_contact_recall_4A"])
                recovered = int(r["recovered_crystal_basic_contact_count_4A"])
                total = int(r["phospho_crystal_basic_contact_count_4A"])
                lines = [
                    f"{recovered}/{total}" + (f" ({recall:.0%})" if not math.isnan(recall) else ""),
                    contact_line("M", r["missing_crystal_basic_contacts_4A"]),
                    contact_line("X", r["extra_predicted_basic_contacts_4A"]),
                ]
                fill = "#cfead1" if recall >= 0.75 else "#fff1bd" if recall >= 0.4 else "#f4c7c3"
            draw.rectangle((x, y - 22, x + col_w - 12, y + 54), fill=fill, outline="#dddddd")
            draw_multiline(draw, x + 8, y - 11, lines, size=10, line_h=18)
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def draw_basic_recovery_simple_png(env: pd.DataFrame, path: Path, variants: Optional[Sequence[str]] = None,
                                   title: str = "Basic residue recovery from crystal contacts (4 A)"):
    variants = list(variants or base.PRIMARY_VARIANTS)
    work = env[env["variant"].isin(variants)].copy()
    cases = list(dict.fromkeys(work["case_id"]))
    has_bb = has_backbone_rmsd_column(work)
    bb_w = 95 if has_bb else 0
    row_h = 38
    col_w = 170
    width = 560 + bb_w + len(variants) * col_w
    height = 120 + len(cases) * row_h
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    text(draw, (width / 2, 32), title, size=22, bold=True)
    x0, y0 = 24, 72
    text(draw, (x0, y0), "Case", size=13, anchor="lm", bold=True)
    if has_bb:
        text(draw, (270, y0), "Local BB", size=12, bold=True)
    text(draw, (300 + bb_w, y0), "Crystal Basic", size=13, bold=True)
    for j, variant in enumerate(variants):
        text(draw, (540 + bb_w + j * col_w, y0), variant, size=12, bold=True)
    for i, cid in enumerate(cases):
        y = y0 + 30 + i * row_h
        sub = work[work["case_id"] == cid]
        first = sub.iloc[0]
        text(draw, (x0, y), cid, size=12, anchor="lm")
        if has_bb:
            draw.rectangle((245, y - 16, 330, y + 12), fill=backbone_rmsd_fill(first), outline="#dddddd")
            text(draw, (287, y - 2), backbone_rmsd_label(first), size=12)
        total = int(first["phospho_crystal_basic_contact_count_4A"])
        text(draw, (300 + bb_w, y), str(total), size=12)
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            x = 470 + bb_w + j * col_w
            if row.empty:
                label = "NA"
                fill = "#f2f2f2"
            else:
                r = row.iloc[0]
                recall = base.safe_float(r["basic_contact_recall_4A"])
                recovered = int(r["recovered_crystal_basic_contact_count_4A"])
                total = int(r["phospho_crystal_basic_contact_count_4A"])
                label = f"{recovered}/{total}" + (f" ({recall:.0%})" if not math.isnan(recall) else "")
                fill = "#cfead1" if recall >= 0.75 else "#fff1bd" if recall >= 0.4 else "#f4c7c3"
            draw.rectangle((x, y - 16, x + col_w - 14, y + 12), fill=fill, outline="#dddddd")
            text(draw, (x + (col_w - 14) / 2, y - 2), label, size=12)
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def draw_basic_atomic_png(env: pd.DataFrame, path: Path, variants: Optional[Sequence[str]] = None,
                          title: str = "Basic residue recovery and donor atom-phosphate oxygen contacts (4 A)"):
    variants = list(variants or base.PRIMARY_VARIANTS)
    work = env[env["variant"].isin(variants)].copy()
    cases = list(dict.fromkeys(work["case_id"]))
    has_bb = has_backbone_rmsd_column(work)
    bb_w = 95 if has_bb else 0
    row_h = 42
    col_w = 185
    width = 640 + bb_w + len(variants) * col_w
    height = 128 + len(cases) * row_h
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    text(draw, (width / 2, 32), title, size=22, bold=True)
    text(draw, (width / 2, 56), "Cell format: recovered residues / crystal residues; predicted donor atoms (delta vs crystal)", size=12)
    x0, y0 = 24, 88
    text(draw, (x0, y0), "Case", size=13, anchor="lm", bold=True)
    if has_bb:
        text(draw, (260, y0), "Local BB", size=12, bold=True)
    text(draw, (315 + bb_w, y0), "Crystal residues", size=12, bold=True)
    text(draw, (450 + bb_w, y0), "Crystal donor atoms", size=12, bold=True)
    for j, variant in enumerate(variants):
        text(draw, (620 + bb_w + j * col_w, y0), variant, size=12, bold=True)
    for i, cid in enumerate(cases):
        y = y0 + 32 + i * row_h
        sub = work[work["case_id"] == cid]
        first = sub.iloc[0]
        text(draw, (x0, y), cid, size=12, anchor="lm")
        if has_bb:
            draw.rectangle((232, y - 16, 318, y + 12), fill=backbone_rmsd_fill(first), outline="#dddddd")
            text(draw, (275, y - 2), backbone_rmsd_label(first), size=12)
        crystal_res = int(first["phospho_crystal_basic_contact_count_4A"])
        crystal_atom = int(first.get("phospho_crystal_basic_atomic_contact_count_4A", 0))
        text(draw, (315 + bb_w, y), str(crystal_res), size=12)
        text(draw, (450 + bb_w, y), str(crystal_atom), size=12)
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            x = 535 + bb_w + j * col_w
            if row.empty:
                label = "NA"
                fill = "#f2f2f2"
            else:
                r = row.iloc[0]
                recall = base.safe_float(r["basic_contact_recall_4A"])
                recovered = int(r["recovered_crystal_basic_contact_count_4A"])
                total = int(r["phospho_crystal_basic_contact_count_4A"])
                pred_atom = int(r.get("predicted_basic_atomic_contact_count_4A", 0))
                delta = int(r.get("basic_atomic_contact_delta_4A", 0))
                label = f"res {recovered}/{total}; donor {pred_atom} ({delta:+d})"
                fill = "#cfead1" if recall >= 0.75 else "#fff1bd" if recall >= 0.4 else "#f4c7c3"
            draw.rectangle((x, y - 16, x + col_w - 12, y + 12), fill=fill, outline="#dddddd")
            text(draw, (x + (col_w - 12) / 2, y - 2), label, size=11)
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def draw_basic_atomic_detail_png(env: pd.DataFrame, path: Path, variants: Optional[Sequence[str]] = None,
                                 title: str = "Detailed donor atom-phosphate oxygen recovery (4 A)"):
    variants = list(variants or base.PRIMARY_VARIANTS)
    work = env[env["variant"].isin(variants)].copy()
    cases = list(dict.fromkeys(work["case_id"]))
    has_bb = has_backbone_rmsd_column(work)
    bb_w = 95 if has_bb else 0
    row_h = 108
    case_w = 220
    crystal_w = 360
    col_w = 340
    width = 36 + case_w + bb_w + crystal_w + len(variants) * col_w + 28
    height = 136 + len(cases) * row_h
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    text(draw, (width / 2, 30), title, size=22, bold=True)
    text(
        draw,
        (width / 2, 55),
        "Donor atoms are compared against any phosphate oxygen: M = missing crystal donors, X = extra predicted donors",
        size=12,
    )
    x0, y0 = 24, 88
    bb_x = x0 + case_w
    crystal_x = x0 + case_w + bb_w
    text(draw, (x0, y0), "Case", size=13, anchor="lm", bold=True)
    if has_bb:
        text(draw, (bb_x, y0), "Local BB", size=12, anchor="lm", bold=True)
    text(draw, (crystal_x, y0), "Crystal donor atoms", size=13, anchor="lm", bold=True)
    for j, variant in enumerate(variants):
        text(draw, (crystal_x + crystal_w + j * col_w + 8, y0), variant, size=12, anchor="lm", bold=True)
    for i, cid in enumerate(cases):
        y = y0 + 36 + i * row_h
        sub = work[work["case_id"] == cid]
        first = sub.iloc[0]
        if i % 2:
            draw.rectangle((10, y - 26, width - 10, y + row_h - 30), fill="#fafafa")
        draw_multiline(draw, x0, y - 8, cid.split(" "), size=12, line_h=15)
        if has_bb:
            draw.rectangle((bb_x, y - 24, bb_x + bb_w - 14, y + 14), fill=backbone_rmsd_fill(first), outline="#dddddd")
            text(draw, (bb_x + (bb_w - 14) / 2, y - 5), backbone_rmsd_label(first), size=11)
        crystal_total = int(first.get("phospho_crystal_basic_atomic_contact_count_4A", 0))
        crystal_set = first.get(
            "phospho_crystal_basic_atomic_contact_set_4A_mapped",
            first.get("phospho_crystal_basic_atomic_contact_set_4A", ""),
        )
        crystal_lines = [
            f"{crystal_total} donor atoms",
            compact_contact_line("C", crystal_set, max_items=4, atomic=True),
        ]
        draw_multiline(draw, crystal_x, y - 10, crystal_lines, size=10, line_h=17)
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            x = crystal_x + crystal_w + j * col_w
            if row.empty:
                lines = ["NA", "M: -", "X: -"]
                fill = "#f2f2f2"
            else:
                r = row.iloc[0]
                recall = base.safe_float(r.get("basic_atomic_contact_recall_4A", float("nan")))
                recovered = int(r.get("recovered_crystal_basic_atomic_contact_count_4A", 0))
                total = int(r.get("phospho_crystal_basic_atomic_contact_count_4A", 0))
                lines = [
                    f"{recovered}/{total}" + (f" ({recall:.0%})" if not math.isnan(recall) else ""),
                    compact_contact_line("M", r.get("missing_crystal_basic_atomic_contacts_4A", ""), max_items=3, atomic=True),
                    compact_contact_line("X", r.get("extra_predicted_basic_atomic_contacts_4A", ""), max_items=3, atomic=True),
                ]
                fill = "#cfead1" if recall >= 0.75 else "#fff1bd" if recall >= 0.4 else "#f4c7c3"
            draw.rectangle((x, y - 24, x + col_w - 12, y + 56), fill=fill, outline="#dddddd")
            draw_multiline(draw, x + 8, y - 12, lines, size=9, line_h=19)
    path.parent.mkdir(parents=True, exist_ok=True)
    image.save(path)


def write_environment_svg(env: pd.DataFrame, path: Path, variants: Optional[Sequence[str]] = None,
                          title: str = "Crystal environment contact recovery (5 A)"):
    variants = list(variants or base.PRIMARY_VARIANTS)
    work = env[env["variant"].isin(variants)].copy()
    cases = list(dict.fromkeys(work["case_id"]))
    has_bb = has_backbone_rmsd_column(work)
    bb_w = 95 if has_bb else 0
    row_h = 34
    col_w = 150
    width = 500 + bb_w + len(variants) * col_w
    height = 90 + len(cases) * row_h
    body = [base.svg_text(width / 2, 28, title, size=18, weight="bold")]
    x0, y0 = 20, 58
    body.append(base.svg_text(x0, y0, "Case", anchor="start", weight="bold"))
    if has_bb:
        body.append(base.svg_text(228, y0, "Local BB", weight="bold", size=11))
    body.append(base.svg_text(265 + bb_w, y0, "Unmod", weight="bold", size=11))
    body.append(base.svg_text(350 + bb_w, y0, "Phospho", weight="bold", size=11))
    for j, variant in enumerate(variants):
        body.append(base.svg_text(500 + bb_w + j * col_w, y0, variant, weight="bold", size=11))
    for i, cid in enumerate(cases):
        y = y0 + 26 + i * row_h
        sub = work[work["case_id"] == cid]
        first = sub.iloc[0]
        body.append(base.svg_text(x0, y, cid, anchor="start", size=11))
        if has_bb:
            body.append(f'<rect x="215" y="{y-18}" width="76" height="24" fill="{backbone_rmsd_fill(first)}" stroke="#dddddd"/>')
            body.append(base.svg_text(253, y - 2, backbone_rmsd_label(first), size=10))
        body.append(base.svg_text(265 + bb_w, y, int(first["unmodified_contact_count_5A"]), size=11))
        body.append(base.svg_text(350 + bb_w, y, int(first["phospho_crystal_contact_count_5A"]), size=11))
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            x = 430 + bb_w + j * col_w
            label = "NA"
            fill = "#f2f2f2"
            if not row.empty:
                r = row.iloc[0]
                recall = base.safe_float(r["contact_recall_5A"])
                recovered = int(r["recovered_crystal_contact_count_5A"])
                total = int(r["phospho_crystal_contact_count_5A"])
                label = f"{recovered}/{total}" + (f" ({recall:.0%})" if not math.isnan(recall) else "")
                fill = "#cfead1" if recall >= 0.75 else "#fff1bd" if recall >= 0.4 else "#f4c7c3"
            body.append(f'<rect x="{x}" y="{y-18}" width="{col_w-10}" height="24" fill="{fill}" stroke="#dddddd"/>')
            body.append(base.svg_text(x + (col_w - 10) / 2, y - 2, label, size=11))
    base.write_svg(path, width, height, body)


def write_basic_recovery_svg(env: pd.DataFrame, path: Path, variants: Optional[Sequence[str]] = None,
                             title: str = "Basic residue recovery from crystal contacts (4 A)"):
    variants = list(variants or base.PRIMARY_VARIANTS)
    work = env[env["variant"].isin(variants)].copy()
    cases = list(dict.fromkeys(work["case_id"]))
    has_bb = has_backbone_rmsd_column(work)
    bb_w = 95 if has_bb else 0
    row_h = 92
    case_w = 220
    crystal_w = 250
    col_w = 255
    width = 30 + case_w + bb_w + crystal_w + len(variants) * col_w + 30
    height = 110 + len(cases) * row_h
    body = [base.svg_text(width / 2, 28, title, size=18, weight="bold")]
    x0, y0 = 20, 62
    bb_x = x0 + case_w
    crystal_x = x0 + case_w + bb_w
    body.append(base.svg_text(x0, y0, "Case", anchor="start", weight="bold"))
    if has_bb:
        body.append(base.svg_text(bb_x, y0, "Local BB", anchor="start", weight="bold", size=11))
    body.append(base.svg_text(crystal_x, y0, "Crystal contacts", anchor="start", weight="bold", size=11))
    for j, variant in enumerate(variants):
        body.append(base.svg_text(crystal_x + crystal_w + j * col_w + 8, y0, variant, anchor="start", weight="bold", size=11))
    for i, cid in enumerate(cases):
        y = y0 + 32 + i * row_h
        sub = work[work["case_id"] == cid]
        first = sub.iloc[0]
        if i % 2:
            body.append(f'<rect x="10" y="{y-24}" width="{width-20}" height="{row_h-8}" fill="#fafafa"/>')
        for k, part in enumerate(cid.split(" ")):
            body.append(base.svg_text(x0, y - 6 + k * 14, part, anchor="start", size=10))
        if has_bb:
            body.append(f'<rect x="{bb_x}" y="{y-22}" width="{bb_w-14}" height="38" fill="{backbone_rmsd_fill(first)}" stroke="#dddddd"/>')
            body.append(base.svg_text(bb_x + (bb_w - 14) / 2, y - 3, backbone_rmsd_label(first), size=10))
        total = int(first["phospho_crystal_basic_contact_count_4A"])
        crystal_lines = [f"{total} total", contact_line("C", first["phospho_crystal_basic_contact_set_4A"])]
        for k, line in enumerate(crystal_lines):
            body.append(base.svg_text(crystal_x, y - 6 + k * 15, line, anchor="start", size=10))
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            x = crystal_x + crystal_w + j * col_w
            lines = ["NA", "M: -", "X: -"]
            fill = "#f2f2f2"
            if not row.empty:
                r = row.iloc[0]
                recall = base.safe_float(r["basic_contact_recall_4A"])
                recovered = int(r["recovered_crystal_basic_contact_count_4A"])
                total = int(r["phospho_crystal_basic_contact_count_4A"])
                lines = [
                    f"{recovered}/{total}" + (f" ({recall:.0%})" if not math.isnan(recall) else ""),
                    contact_line("M", r["missing_crystal_basic_contacts_4A"]),
                    contact_line("X", r["extra_predicted_basic_contacts_4A"]),
                ]
                fill = "#cfead1" if recall >= 0.75 else "#fff1bd" if recall >= 0.4 else "#f4c7c3"
            body.append(f'<rect x="{x}" y="{y-22}" width="{col_w-12}" height="76" fill="{fill}" stroke="#dddddd"/>')
            for k, line in enumerate(lines):
                body.append(base.svg_text(x + 8, y - 8 + k * 18, line, anchor="start", size=10))
    base.write_svg(path, width, height, body)


def write_basic_recovery_simple_svg(env: pd.DataFrame, path: Path, variants: Optional[Sequence[str]] = None,
                                    title: str = "Basic residue recovery from crystal contacts (4 A)"):
    variants = list(variants or base.PRIMARY_VARIANTS)
    work = env[env["variant"].isin(variants)].copy()
    cases = list(dict.fromkeys(work["case_id"]))
    has_bb = has_backbone_rmsd_column(work)
    bb_w = 95 if has_bb else 0
    row_h = 34
    col_w = 150
    width = 535 + bb_w + len(variants) * col_w
    height = 90 + len(cases) * row_h
    body = [base.svg_text(width / 2, 28, title, size=18, weight="bold")]
    x0, y0 = 20, 58
    body.append(base.svg_text(x0, y0, "Case", anchor="start", weight="bold"))
    if has_bb:
        body.append(base.svg_text(270, y0, "Local BB", weight="bold", size=11))
    body.append(base.svg_text(300 + bb_w, y0, "Crystal Basic", weight="bold", size=11))
    for j, variant in enumerate(variants):
        body.append(base.svg_text(535 + bb_w + j * col_w, y0, variant, weight="bold", size=11))
    for i, cid in enumerate(cases):
        y = y0 + 26 + i * row_h
        sub = work[work["case_id"] == cid]
        first = sub.iloc[0]
        body.append(base.svg_text(x0, y, cid, anchor="start", size=11))
        if has_bb:
            body.append(f'<rect x="245" y="{y-18}" width="78" height="24" fill="{backbone_rmsd_fill(first)}" stroke="#dddddd"/>')
            body.append(base.svg_text(284, y - 2, backbone_rmsd_label(first), size=10))
        body.append(base.svg_text(300 + bb_w, y, int(first["phospho_crystal_basic_contact_count_4A"]), size=11))
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            x = 465 + bb_w + j * col_w
            label = "NA"
            fill = "#f2f2f2"
            if not row.empty:
                r = row.iloc[0]
                recall = base.safe_float(r["basic_contact_recall_4A"])
                recovered = int(r["recovered_crystal_basic_contact_count_4A"])
                total = int(r["phospho_crystal_basic_contact_count_4A"])
                label = f"{recovered}/{total}" + (f" ({recall:.0%})" if not math.isnan(recall) else "")
                fill = "#cfead1" if recall >= 0.75 else "#fff1bd" if recall >= 0.4 else "#f4c7c3"
            body.append(f'<rect x="{x}" y="{y-18}" width="{col_w-10}" height="24" fill="{fill}" stroke="#dddddd"/>')
            body.append(base.svg_text(x + (col_w - 10) / 2, y - 2, label, size=11))
    base.write_svg(path, width, height, body)


def write_basic_atomic_svg(env: pd.DataFrame, path: Path, variants: Optional[Sequence[str]] = None,
                           title: str = "Basic residue recovery and donor atom-phosphate oxygen contacts (4 A)"):
    variants = list(variants or base.PRIMARY_VARIANTS)
    work = env[env["variant"].isin(variants)].copy()
    cases = list(dict.fromkeys(work["case_id"]))
    has_bb = has_backbone_rmsd_column(work)
    bb_w = 95 if has_bb else 0
    row_h = 36
    col_w = 170
    width = 595 + bb_w + len(variants) * col_w
    height = 106 + len(cases) * row_h
    body = [
        base.svg_text(width / 2, 26, title, size=18, weight="bold"),
        base.svg_text(width / 2, 48, "Cell format: recovered residues / crystal residues; predicted donor atoms (delta vs crystal)", size=10),
    ]
    x0, y0 = 20, 76
    body.append(base.svg_text(x0, y0, "Case", anchor="start", weight="bold", size=11))
    if has_bb:
        body.append(base.svg_text(260, y0, "Local BB", weight="bold", size=11))
    body.append(base.svg_text(315 + bb_w, y0, "Crystal res", weight="bold", size=11))
    body.append(base.svg_text(430 + bb_w, y0, "Crystal donor", weight="bold", size=11))
    for j, variant in enumerate(variants):
        body.append(base.svg_text(585 + bb_w + j * col_w, y0, variant, weight="bold", size=10))
    for i, cid in enumerate(cases):
        y = y0 + 26 + i * row_h
        sub = work[work["case_id"] == cid]
        first = sub.iloc[0]
        body.append(base.svg_text(x0, y, cid, anchor="start", size=10))
        if has_bb:
            body.append(f'<rect x="230" y="{y-18}" width="82" height="24" fill="{backbone_rmsd_fill(first)}" stroke="#dddddd"/>')
            body.append(base.svg_text(271, y - 2, backbone_rmsd_label(first), size=10))
        body.append(base.svg_text(315 + bb_w, y, int(first["phospho_crystal_basic_contact_count_4A"]), size=10))
        body.append(base.svg_text(430 + bb_w, y, int(first.get("phospho_crystal_basic_atomic_contact_count_4A", 0)), size=10))
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            x = 500 + bb_w + j * col_w
            label = "NA"
            fill = "#f2f2f2"
            if not row.empty:
                r = row.iloc[0]
                recall = base.safe_float(r["basic_contact_recall_4A"])
                recovered = int(r["recovered_crystal_basic_contact_count_4A"])
                total = int(r["phospho_crystal_basic_contact_count_4A"])
                pred_atom = int(r.get("predicted_basic_atomic_contact_count_4A", 0))
                delta = int(r.get("basic_atomic_contact_delta_4A", 0))
                label = f"res {recovered}/{total}; donor {pred_atom} ({delta:+d})"
                fill = "#cfead1" if recall >= 0.75 else "#fff1bd" if recall >= 0.4 else "#f4c7c3"
            body.append(f'<rect x="{x}" y="{y-18}" width="{col_w-10}" height="24" fill="{fill}" stroke="#dddddd"/>')
            body.append(base.svg_text(x + (col_w - 10) / 2, y - 2, label, size=9))
    base.write_svg(path, width, height, body)


def write_basic_atomic_detail_svg(env: pd.DataFrame, path: Path, variants: Optional[Sequence[str]] = None,
                                  title: str = "Detailed donor atom-phosphate oxygen recovery (4 A)"):
    variants = list(variants or base.PRIMARY_VARIANTS)
    work = env[env["variant"].isin(variants)].copy()
    cases = list(dict.fromkeys(work["case_id"]))
    has_bb = has_backbone_rmsd_column(work)
    bb_w = 95 if has_bb else 0
    row_h = 104
    case_w = 220
    crystal_w = 360
    col_w = 340
    width = 30 + case_w + bb_w + crystal_w + len(variants) * col_w + 30
    height = 126 + len(cases) * row_h
    body = [
        base.svg_text(width / 2, 26, title, size=18, weight="bold"),
        base.svg_text(
            width / 2,
            48,
            "Donor atoms are compared against any phosphate oxygen: M = missing crystal donors, X = extra predicted donors",
            size=10,
        ),
    ]
    x0, y0 = 20, 78
    bb_x = x0 + case_w
    crystal_x = x0 + case_w + bb_w
    body.append(base.svg_text(x0, y0, "Case", anchor="start", weight="bold", size=11))
    if has_bb:
        body.append(base.svg_text(bb_x, y0, "Local BB", anchor="start", weight="bold", size=11))
    body.append(base.svg_text(crystal_x, y0, "Crystal donor atoms", anchor="start", weight="bold", size=11))
    for j, variant in enumerate(variants):
        body.append(base.svg_text(crystal_x + crystal_w + j * col_w + 8, y0, variant, anchor="start", weight="bold", size=11))
    for i, cid in enumerate(cases):
        y = y0 + 34 + i * row_h
        sub = work[work["case_id"] == cid]
        first = sub.iloc[0]
        if i % 2:
            body.append(f'<rect x="10" y="{y-26}" width="{width-20}" height="{row_h-8}" fill="#fafafa"/>')
        for k, part in enumerate(cid.split(" ")):
            body.append(base.svg_text(x0, y - 6 + k * 14, part, anchor="start", size=10))
        if has_bb:
            body.append(f'<rect x="{bb_x}" y="{y-24}" width="{bb_w-14}" height="38" fill="{backbone_rmsd_fill(first)}" stroke="#dddddd"/>')
            body.append(base.svg_text(bb_x + (bb_w - 14) / 2, y - 5, backbone_rmsd_label(first), size=10))
        crystal_total = int(first.get("phospho_crystal_basic_atomic_contact_count_4A", 0))
        crystal_set = first.get(
            "phospho_crystal_basic_atomic_contact_set_4A_mapped",
            first.get("phospho_crystal_basic_atomic_contact_set_4A", ""),
        )
        crystal_lines = [
            f"{crystal_total} donor atoms",
            compact_contact_line("C", crystal_set, max_items=4, atomic=True),
        ]
        for k, line in enumerate(crystal_lines):
            body.append(base.svg_text(crystal_x, y - 6 + k * 16, line, anchor="start", size=9))
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            x = crystal_x + crystal_w + j * col_w
            lines = ["NA", "M: -", "X: -"]
            fill = "#f2f2f2"
            if not row.empty:
                r = row.iloc[0]
                recall = base.safe_float(r.get("basic_atomic_contact_recall_4A", float("nan")))
                recovered = int(r.get("recovered_crystal_basic_atomic_contact_count_4A", 0))
                total = int(r.get("phospho_crystal_basic_atomic_contact_count_4A", 0))
                lines = [
                    f"{recovered}/{total}" + (f" ({recall:.0%})" if not math.isnan(recall) else ""),
                    compact_contact_line("M", r.get("missing_crystal_basic_atomic_contacts_4A", ""), max_items=3, atomic=True),
                    compact_contact_line("X", r.get("extra_predicted_basic_atomic_contacts_4A", ""), max_items=3, atomic=True),
                ]
                fill = "#cfead1" if recall >= 0.75 else "#fff1bd" if recall >= 0.4 else "#f4c7c3"
            body.append(f'<rect x="{x}" y="{y-24}" width="{col_w-12}" height="78" fill="{fill}" stroke="#dddddd"/>')
            for k, line in enumerate(lines):
                body.append(base.svg_text(x + 8, y - 8 + k * 18, line, anchor="start", size=9))
    base.write_svg(path, width, height, body)


def pml_quote(path: str | Path) -> str:
    return str(path).replace("\\", "/").replace('"', '\\"')


PYMOL_COLORS = [
    "set_color unmodified_gray80, [0.80, 0.80, 0.80]",
    "set_color crystal_wheat_tint, [0.93, 0.82, 0.62]",
    "set_color phosphofill_green_1, [0.00, 0.48, 0.18]",
    "set_color phosphofill_green_2, [0.20, 0.65, 0.30]",
    "set_color phosphofill_green_3, [0.50, 0.78, 0.45]",
    "set_color pytms_blue_default, [0.25, 0.58, 0.95]",
    "set_color pytms_blue_optimized, [0.02, 0.25, 0.75]",
    "set_color ptmpsi_gold, [0.95, 0.62, 0.05]",
    "set_color contact_crystal_green, [0.00, 0.48, 0.18]",
    "set_color contact_maintained_gray, [0.45, 0.45, 0.45]",
    "set_color contact_new_blue, [0.05, 0.35, 0.95]",
    "set_color contact_lost_red, [0.88, 0.12, 0.10]",
]


def phosphofill_visual_structures(case: Dict, pdb_id: str) -> List[Dict]:
    pf_df = base.read_tsv(case["pf_file"])
    pf_row = base.find_case_row(pf_df, case, pdb_id, case["residue_type"])
    if pf_row is None:
        return []
    structures = []
    for rank, color_name in [
        (1, "phosphofill_green_1"),
        (2, "phosphofill_green_2"),
        (3, "phosphofill_green_3"),
    ]:
        path_value = str(pf_row.get(f"archive_stage3_top{rank}_structure", "")).strip()
        resolved = base.resolve_path(path_value)
        if not path_value or resolved is None or not resolved.exists():
            continue
        structures.append({
            "object_name": f"PhosphoFill{rank}",
            "path": resolved,
            "color": color_name,
            "tool": "PhosphoFill",
            "site_chain": case["chain"],
            "site_position": case["position"],
            "label": f"PhosphoFill top-{rank}",
        })
    return structures


def write_pymol_scripts(metrics: pd.DataFrame, env: pd.DataFrame, unmodified_by_site: Dict, phospho_by_case: Dict):
    PYMOL_DIR.mkdir(parents=True, exist_ok=True)
    manifest = []
    for case in base.CASES:
        site_key = (case["protein"], case["residue_type"], int(case["position"]))
        unmod = unmodified_by_site[site_key]
        for pdb_id in case["pdb_ids"]:
            case_id = base.case_key(case, pdb_id)
            safe_case = safe_name(case_id)
            script_path = PYMOL_DIR / f"{safe_case}.pml"
            ref_path = base.PDB_CACHE / f"{pdb_id}.cif"
            phospho_summary = phospho_by_case[(case["site_label"], pdb_id)]
            visible_chain = case["chain"]
            distance_alias_members = {
                "PhosphoFill": {"maintained": [], "new": [], "lost": []},
                "PyTMs": {"maintained": [], "new": [], "lost": []},
            }
            atomic_distance_alias_members = {
                "PhosphoFill": {"maintained": [], "new": [], "lost": []},
                "PyTMs": {"maintained": [], "new": [], "lost": []},
            }
            lines = [
                "reinitialize",
                "bg_color white",
                "set cartoon_transparency, 0.65",
                "set stick_radius, 0.16",
                "set dash_radius, 0.06",
                "set label_size, 16",
                *PYMOL_COLORS,
                f"# {case_id}",
                f"# {unmod['note']}",
                "# Contact colors per prediction: maintained/recovered grey, new/extra blue, lost/missing red.",
                f'load "{pml_quote(unmod["path"])}", unmodified',
                f'load "{pml_quote(ref_path)}", phospho_crystal',
                f"align unmodified and chain {unmod['chain']}, phospho_crystal and chain {visible_chain}",
                "hide everything, all",
                f"show cartoon, unmodified and chain {unmod['chain']}",
                f"show cartoon, phospho_crystal and chain {visible_chain}",
                f"hide everything, unmodified and not chain {unmod['chain']}",
                f"hide everything, phospho_crystal and not chain {visible_chain}",
                f"color unmodified_gray80, unmodified and chain {unmod['chain']}",
                f"color crystal_wheat_tint, phospho_crystal and chain {visible_chain}",
                f"select unmodified_site, unmodified and chain {unmod['chain']} and resi {unmod['position']}",
                f"select unmodified_marker, unmodified_site and name {'+'.join(UNMODIFIED_MARKER_ATOMS[case['residue_type']])}",
                f"select phospho_site, phospho_crystal and chain {visible_chain} and resi {case['position']}",
                f"select phospho_phosphate, phospho_site and name {PHOSPHATE_SELECTION}",
                f"select phospho_phosphate_oxygen, phospho_site and name {PHOSPHATE_OXYGEN_SELECTION}",
                *phosphosite_bond_cleanup_pml("phospho_site", case["residue_type"]),
                "show sticks, unmodified_site or phospho_site",
                "color crystal_wheat_tint, phospho_phosphate",
                f"select phospho_contacts_5A, byres (((phospho_crystal and chain {visible_chain}) within 5 of phospho_phosphate) and polymer.protein and not phospho_site)",
                "show sticks, phospho_contacts_5A",
                f"select phospho_basic_residues_4A, byres (((phospho_crystal and chain {visible_chain} and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2) within 4 of phospho_phosphate_oxygen))",
                "show sticks, phospho_basic_residues_4A",
                "color contact_crystal_green, phospho_basic_residues_4A",
                f"distance phospho_basic_4A, phospho_phosphate_oxygen, ((phospho_crystal and chain {visible_chain} and polymer.protein and resn ARG+LYS+HIS and name NE+NH1+NH2+NZ+ND1+NE2) within 4 of phospho_phosphate_oxygen)",
                "color contact_crystal_green, phospho_basic_4A",
                "set dash_color, contact_crystal_green, phospho_basic_4A",
                "set dash_width, 2.4, phospho_basic_4A",
                "hide labels, phospho_basic_4A",
            ]
            crystal_selection = residue_selection(phospho_summary["contact_tuples"], "phospho_crystal")
            lines.append(f"select phospho_crystal_contact_residues, {crystal_selection}")

            def computed_contact_classes(structure_path: Path, tool: str, site_chain: str, site_pos: int) -> Dict[str, Dict[str, str]]:
                pred_model = base.first_model(structure_path)
                pred_res, pred_site_chain, pred_site_pos = find_site_residue(
                    pred_model, site_chain, int(site_pos), case["residue_type"]
                )
                pred_summary = summarize_environment(
                    pred_model, pred_res, pred_site_chain, pred_site_pos, marker_coords_for_phospho(pred_res)
                )
                mapped_crystal_basic = base.map_contacts_for_output(
                    phospho_summary["basic_tuples"], phospho_summary["model"], tool
                )
                comparison = compare_sets(pred_summary["basic_labels"], labels_from_tuples(mapped_crystal_basic))
                mapped_crystal_atomic = mapped_atomic_labels_from_rows(
                    phospho_summary["basic_atomic_rows"], phospho_summary["model"], tool
                )
                atomic_comparison = compare_basic_atomic_sets(
                    pred_summary["basic_atomic_labels"], mapped_crystal_atomic
                )
                return {
                    "basic": {
                        "maintained": comparison["recovered_crystal_contacts_5A"],
                        "new": comparison["extra_predicted_contacts_5A"],
                        "lost": comparison["missing_crystal_contacts_5A"],
                    },
                    "atomic": {
                        "maintained": atomic_comparison["recovered_crystal_basic_atomic_contacts_4A"],
                        "new": atomic_comparison["extra_predicted_basic_atomic_contacts_4A"],
                        "lost": atomic_comparison["missing_crystal_basic_atomic_contacts_4A"],
                    },
                }

            def add_prediction_object(
                obj: str,
                structure_path: Path,
                color_name: str,
                tool: str,
                site_chain: str,
                site_pos: int,
                contact_classes: Optional[Dict[str, str]] = None,
                atomic_contact_classes: Optional[Dict[str, str]] = None,
            ):
                computed_classes = None
                if contact_classes is None or atomic_contact_classes is None:
                    computed_classes = computed_contact_classes(structure_path, tool, site_chain, site_pos)
                if contact_classes is None:
                    contact_classes = computed_classes["basic"]
                if atomic_contact_classes is None:
                    atomic_contact_classes = computed_classes["atomic"]
                for suffix in ("maintained", "new", "lost"):
                    if contact_items(contact_classes.get(suffix, "")):
                        distance_name = f"{obj}_{suffix}_basic_contacts"
                        if obj.startswith("PhosphoFill"):
                            distance_alias_members["PhosphoFill"][suffix].append(distance_name)
                        elif obj.startswith("PyTMs_"):
                            distance_alias_members["PyTMs"][suffix].append(distance_name)
                    if contact_items(atomic_contact_classes.get(suffix, "")):
                        atomic_distance_name = f"{obj}_{suffix}_basic_atomic_contacts"
                        if obj.startswith("PhosphoFill"):
                            atomic_distance_alias_members["PhosphoFill"][suffix].append(atomic_distance_name)
                        elif obj.startswith("PyTMs_"):
                            atomic_distance_alias_members["PyTMs"][suffix].append(atomic_distance_name)
                lines.extend([
                    f'load "{pml_quote(structure_path)}", {obj}',
                    f"align {obj} and chain {site_chain}, phospho_crystal and chain {visible_chain}",
                    f"hide everything, {obj}",
                    f"show cartoon, {obj} and chain {site_chain}",
                    f"hide everything, {obj} and not chain {site_chain}",
                    f"color {color_name}, {obj} and chain {site_chain}",
                    f"select {obj}_site, {obj} and chain {site_chain} and resi {site_pos}",
                    f"select {obj}_phosphate, {obj}_site and name {PHOSPHATE_SELECTION}",
                    *phosphosite_bond_cleanup_pml(f"{obj}_site", case["residue_type"]),
                    f"show sticks, {obj}_site",
                    f"select {obj}_contacts_5A, byres ((({obj} and chain {site_chain}) within 5 of {obj}_phosphate) and polymer.protein and not {obj}_site)",
                    f"show sticks, {obj}_contacts_5A",
                ])
                lines.extend(colored_basic_contact_pml(
                    obj,
                    f"{obj}_phosphate",
                    contact_classes.get("maintained", ""),
                    contact_classes.get("new", ""),
                    contact_classes.get("lost", ""),
                ))
                lines.extend(colored_basic_atomic_contact_pml(
                    obj,
                    f"{obj}_phosphate",
                    atomic_contact_classes.get("maintained", ""),
                    atomic_contact_classes.get("new", ""),
                    atomic_contact_classes.get("lost", ""),
                ))
                mapped = base.map_contacts_for_output(phospho_summary["contact_tuples"], phospho_summary["model"], tool)
                selection = residue_selection(mapped, obj)
                lines.append(f"select {obj}_crystal_contact_residues, {selection}")

            for pf_visual in phosphofill_visual_structures(case, pdb_id):
                add_prediction_object(
                    pf_visual["object_name"],
                    pf_visual["path"],
                    pf_visual["color"],
                    pf_visual["tool"],
                    pf_visual["site_chain"],
                    int(pf_visual["site_position"]),
                )

            case_metrics = metrics[(metrics["case_id"] == case_id)].copy()
            for _, row in case_metrics.iterrows():
                if str(row["tool"]) == "PhosphoFill":
                    continue
                obj = visual_object_name(row["variant"])
                structure_path = base.resolve_path(row["structure_path"])
                if structure_path is None or not structure_path.exists():
                    continue
                env_row = env[(env["case_id"] == case_id) & (env["variant"] == row["variant"])].iloc[0]
                site_chain = env_row["predicted_site_chain"]
                site_pos = int(env_row["predicted_site_position"])
                contact_classes = {
                    "maintained": str(env_row.get("recovered_crystal_basic_contacts_4A", "")),
                    "new": str(env_row.get("extra_predicted_basic_contacts_4A", "")),
                    "lost": str(env_row.get("missing_crystal_basic_contacts_4A", "")),
                }
                atomic_contact_classes = {
                    "maintained": str(env_row.get("recovered_crystal_basic_atomic_contacts_4A", "")),
                    "new": str(env_row.get("extra_predicted_basic_atomic_contacts_4A", "")),
                    "lost": str(env_row.get("missing_crystal_basic_atomic_contacts_4A", "")),
                }
                add_prediction_object(
                    obj,
                    structure_path,
                    variant_color_name(row["variant"]),
                    str(row["tool"]),
                    str(site_chain),
                    site_pos,
                    contact_classes,
                    atomic_contact_classes,
                )
            lines.extend([
                "# Convenience groups: show dashes, PyTMs_maintained_basic_contacts / PyTMs_new_basic_contacts / PyTMs_lost_basic_contacts",
                "# Donor-atom groups: show dashes, PyTMs_maintained_basic_atomic_contacts / PyTMs_new_basic_atomic_contacts / PyTMs_lost_basic_atomic_contacts",
            ])
            for alias_root, suffix_map in distance_alias_members.items():
                for suffix, members in suffix_map.items():
                    alias = f"{alias_root}_{suffix}_basic_contacts"
                    if members:
                        lines.append(f"group {alias}, {' '.join(members)}")
                    else:
                        lines.append(f"select {alias}, none")
            for alias_root, suffix_map in atomic_distance_alias_members.items():
                for suffix, members in suffix_map.items():
                    alias = f"{alias_root}_{suffix}_basic_atomic_contacts"
                    if members:
                        lines.append(f"group {alias}, {' '.join(members)}")
                    else:
                        lines.append(f"select {alias}, none")
            lines.extend([
                "orient phospho_site",
                "zoom phospho_site or phospho_contacts_5A, 12",
                "set_view (\\",
                "    0.7, 0.1, 0.7,\\",
                "    0.0, 1.0, -0.2,\\",
                "    -0.7, 0.2, 0.7,\\",
                "    0.0, 0.0, -120.0,\\",
                "    0.0, 0.0, 0.0,\\",
                "    80.0, 160.0, -20.0 )",
            ])
            script_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            manifest.append({
                "case_id": case_id,
                "pdb_id": pdb_id,
                "site_label": case["site_label"],
                "pymol_script": str(script_path),
            })
    pd.DataFrame(manifest).to_csv(OUT_DIR / "case_study_pymol_manifest.tsv", sep="\t", index=False)


def variant_color_name(variant: str) -> str:
    return {
        "PyTMs default": "pytms_blue_default",
        "PyTMs optimized": "pytms_blue_optimized",
        "PTM-Psi": "ptmpsi_gold",
    }.get(str(variant), "white")


def visual_object_name(variant: str) -> str:
    return {
        "PyTMs default": "PyTMs_default",
        "PyTMs optimized": "PyTMs_optimized",
        "PTM-Psi": "PTM_Psi",
    }.get(str(variant), safe_name(variant)[:40])


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    PYMOL_DIR.mkdir(parents=True, exist_ok=True)

    metrics_path = OUT_DIR / "case_study_metrics_long.tsv"
    if not metrics_path.exists():
        raise FileNotFoundError(f"Run case_study_tool_analysis.py first: {metrics_path}")
    metrics = pd.read_csv(metrics_path, sep="\t")

    env, contact_details, basic_atomic_details, unmodified, unmodified_by_site, phospho_by_case = build_environment_tables(metrics)
    variant_order = {variant: i for i, variant in enumerate(base.SUPPLEMENTARY_VARIANTS)}
    if not env.empty:
        env["_variant_order"] = env["variant"].map(variant_order).fillna(99)
        env.sort_values(["protein", "position", "pdb_id", "_variant_order", "variant"], inplace=True)
        env.drop(columns=["_variant_order"], inplace=True)
    env.to_csv(OUT_DIR / "case_study_environment_recovery.tsv", sep="\t", index=False)
    env[env["variant"].isin(base.PRIMARY_VARIANTS)].to_csv(
        OUT_DIR / "case_study_environment_recovery_primary.tsv", sep="\t", index=False
    )
    env[[
        "protein", "site_label", "case_id", "pdb_id", "tool", "variant",
        "phospho_crystal_basic_contact_count_4A", "predicted_basic_contact_count_4A",
        "phospho_crystal_basic_atomic_contact_count_4A", "predicted_basic_atomic_contact_count_4A",
        "basic_atomic_contact_delta_4A",
        "recovered_crystal_basic_contact_count_4A", "missing_crystal_basic_contact_count_4A",
        "extra_predicted_basic_contact_count_4A", "basic_contact_recall_4A",
        "basic_contact_precision_4A", "phospho_crystal_basic_contact_set_4A",
        "phospho_crystal_basic_contact_set_4A_mapped", "predicted_basic_contact_set_4A",
        "recovered_crystal_basic_atomic_contact_count_4A", "missing_crystal_basic_atomic_contact_count_4A",
        "extra_predicted_basic_atomic_contact_count_4A", "basic_atomic_contact_recall_4A",
        "basic_atomic_contact_precision_4A",
        "phospho_crystal_basic_atomic_contact_set_4A", "phospho_crystal_basic_atomic_contact_set_4A_mapped",
        "predicted_basic_atomic_contact_set_4A", "recovered_crystal_basic_atomic_contacts_4A",
        "missing_crystal_basic_atomic_contacts_4A", "extra_predicted_basic_atomic_contacts_4A",
        "recovered_crystal_basic_contacts_4A", "missing_crystal_basic_contacts_4A",
        "extra_predicted_basic_contacts_4A", "structure_path",
    ]].to_csv(OUT_DIR / "case_study_crystal_basic_contact_recovery.tsv", sep="\t", index=False)
    contact_details.to_csv(OUT_DIR / "case_study_environment_contacts.tsv", sep="\t", index=False)
    basic_atomic_details.to_csv(OUT_DIR / "case_study_basic_atomic_contacts.tsv", sep="\t", index=False)
    unmodified.to_csv(OUT_DIR / "case_study_unmodified_environment.tsv", sep="\t", index=False)

    draw_grouped_bar_png(metrics, FIG_DIR / "rmsd_grouped_bar_with_pytms_default.png", base.SUPPLEMENTARY_VARIANTS,
                         "RMSD comparison including PyTMs default")
    draw_torsion_png(metrics, FIG_DIR / "torsion_mode_comparison_with_pytms_default.png", base.SUPPLEMENTARY_VARIANTS,
                     "Predicted phosphate torsion vs crystal including PyTMs default")
    draw_torsion_polar_png(metrics, FIG_DIR / "torsion_mode_polar_with_pytms_default.png",
                           base.SUPPLEMENTARY_VARIANTS,
                           "Circular torsion mode comparison including PyTMs default")
    draw_environment_png(env, FIG_DIR / "environment_recovery_table_with_pytms_default.png",
                         base.SUPPLEMENTARY_VARIANTS,
                         "Crystal environment contact recovery including PyTMs default (5 A)")
    write_environment_svg(env, FIG_DIR / "environment_recovery_table_with_pytms_default.svg",
                          base.SUPPLEMENTARY_VARIANTS,
                          "Crystal environment contact recovery including PyTMs default (5 A)")
    draw_basic_recovery_png(env, FIG_DIR / "salt_bridge_recovery_table_with_pytms_default.png",
                            base.SUPPLEMENTARY_VARIANTS,
                            "Basic residue recovery including PyTMs default (4 A)")
    write_basic_recovery_svg(env, FIG_DIR / "salt_bridge_recovery_table_with_pytms_default.svg",
                             base.SUPPLEMENTARY_VARIANTS,
                             "Basic residue recovery including PyTMs default (4 A)")
    draw_basic_recovery_simple_png(env, FIG_DIR / "salt_bridge_recovery_simple_table_with_pytms_default.png",
                                   base.SUPPLEMENTARY_VARIANTS,
                                   "Basic residue recovery including PyTMs default (4 A)")
    write_basic_recovery_simple_svg(env, FIG_DIR / "salt_bridge_recovery_simple_table_with_pytms_default.svg",
                                    base.SUPPLEMENTARY_VARIANTS,
                                    "Basic residue recovery including PyTMs default (4 A)")
    draw_basic_atomic_png(env, FIG_DIR / "salt_bridge_atomic_contact_table_with_pytms_default.png",
                          base.SUPPLEMENTARY_VARIANTS,
                          "Basic residue recovery and donor atom contacts including PyTMs default (4 A)")
    write_basic_atomic_svg(env, FIG_DIR / "salt_bridge_atomic_contact_table_with_pytms_default.svg",
                           base.SUPPLEMENTARY_VARIANTS,
                           "Basic residue recovery and donor atom contacts including PyTMs default (4 A)")
    draw_basic_atomic_detail_png(env, FIG_DIR / "salt_bridge_atomic_contact_detail_table_with_pytms_default.png",
                                 base.SUPPLEMENTARY_VARIANTS,
                                 "Detailed donor atom contact recovery including PyTMs default (4 A)")
    write_basic_atomic_detail_svg(env, FIG_DIR / "salt_bridge_atomic_contact_detail_table_with_pytms_default.svg",
                                  base.SUPPLEMENTARY_VARIANTS,
                                  "Detailed donor atom contact recovery including PyTMs default (4 A)")

    write_pymol_scripts(metrics, env, unmodified_by_site, phospho_by_case)

    print(f"Wrote {len(env)} environment recovery rows")
    print(f"Wrote {len(contact_details)} detailed contact rows")
    print(f"Wrote PyMOL scripts to {PYMOL_DIR}")
    print(f"Wrote PNG figures to {FIG_DIR}")


if __name__ == "__main__":
    main()
