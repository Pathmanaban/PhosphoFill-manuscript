#!/usr/bin/env python3
"""Consolidate CDK2/ERK2 strip-and-regraft case-study metrics."""
from __future__ import annotations

import html
import math
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

ROOT = Path(__file__).resolve().parent
LOCAL_BIO_DEPS = ROOT / "case_study_outputs" / "python_deps_bioonly"
if LOCAL_BIO_DEPS.exists():
    sys.path.insert(0, str(LOCAL_BIO_DEPS))

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser

OUT_DIR = ROOT / "case_study_outputs"
FIG_DIR = OUT_DIR / "figures"
PDB_CACHE = ROOT / "pdb_cache"

CIF_PARSER = MMCIFParser(QUIET=True)
PDB_PARSER = PDBParser(QUIET=True)

DIHEDRAL_ATOMS = {
    "TPO": ("CG2", "CB", "OG1", "P"),
    "SEP": ("CA", "CB", "OG", "P"),
    "PTR": ("CE1", "CZ", "OH", "P"),
}
PHOSPHATE_ATOMS = {"P", "O1P", "O2P", "O3P"}
BASIC_DONORS = {
    "ARG": {"NE", "NH1", "NH2"},
    "LYS": {"NZ"},
    "HIS": {"ND1", "NE2"},
}
RESTYPE_MAP = {"T": "TPO", "THR": "TPO", "TPO": "TPO", "Y": "PTR", "TYR": "PTR", "PTR": "PTR"}

PRIMARY_VARIANTS = ["PhosphoFill top1", "PhosphoFill best-of-3", "PyTMs optimized", "PTM-Psi"]
SUPPLEMENTARY_VARIANTS = [
    "PhosphoFill top1",
    "PhosphoFill best-of-3",
    "PyTMs default",
    "PyTMs optimized",
    "PTM-Psi",
]
COLORS = {
    "PhosphoFill top1": "#2f6fbb",
    "PhosphoFill best-of-3": "#7aa6d8",
    "PyTMs optimized": "#d97932",
    "PyTMs default": "#e8b07a",
    "PTM-Psi": "#4c9a5f",
    "Reference": "#222222",
}

CASES = [
    {
        "protein": "CDK2",
        "site_label": "CDK2 pThr160",
        "acc_id": "P24941",
        "pdb_ids": ["2cch", "4eoj"],
        "chain": "A",
        "position": 160,
        "residue_type": "TPO",
        "pf_file": "TPO_final_res.tsv",
        "expected_contacts": [("A", "ARG", 150)],
    },
    {
        "protein": "ERK2",
        "site_label": "ERK2 pThr185",
        "acc_id": "P28482",
        "pdb_ids": ["4iza", "5v62"],
        "chain": "A",
        "position": 185,
        "residue_type": "TPO",
        "pf_file": "TPO_final_res.tsv",
        "expected_contacts": [("A", "ARG", 65), ("A", "ARG", 68), ("A", "ARG", 146)],
    },
    {
        "protein": "ERK2",
        "site_label": "ERK2 pTyr187",
        "acc_id": "P28482",
        "pdb_ids": ["4iza", "5v62"],
        "chain": "A",
        "position": 187,
        "residue_type": "PTR",
        "pf_file": "case_study_inputs/PTR_P28482_ERK2_Y187.tsv",
        "expected_contacts": [("A", "ARG", 189)],
        "include_observed_arg_contacts": True,
    },
]


def read_tsv(path: str | Path) -> pd.DataFrame:
    p = ROOT / path if not Path(path).is_absolute() else Path(path)
    if not p.exists():
        return pd.DataFrame()
    return pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False)


def safe_float(value) -> float:
    try:
        if value is None or str(value).strip() == "":
            return float("nan")
        return float(value)
    except Exception:
        return float("nan")


def norm_restype(value: str) -> str:
    return RESTYPE_MAP.get(str(value).strip().upper(), str(value).strip().upper())


def rel_or_abs(path_value: str | Path | None) -> str:
    if path_value is None:
        return ""
    text = str(path_value).strip()
    if not text:
        return ""
    return text


def resolve_path(path_value: str | Path | None) -> Optional[Path]:
    text = rel_or_abs(path_value)
    if not text:
        return None
    p = Path(text)
    if not p.is_absolute():
        p = ROOT / p
    return p


def pdb_from_reference_path(path_value: str) -> str:
    return Path(str(path_value).replace("\\", os.sep)).stem.lower()


def load_structure(path: str | Path, sid: str = "structure"):
    p = Path(path)
    if p.suffix.lower() in {".cif", ".mmcif"}:
        return CIF_PARSER.get_structure(sid, str(p))
    return PDB_PARSER.get_structure(sid, str(p))


def first_model(path: str | Path):
    structure = load_structure(path)
    return list(structure.get_models())[0]


def find_residue(model, chain_id: str, resseq: int):
    if chain_id not in model:
        return None
    chain = model[chain_id]
    for hetfield in [" ", "H_TPO", "H_SEP", "H_PTR", "H_PTM"]:
        try:
            return chain[(hetfield, int(resseq), " ")]
        except KeyError:
            pass
    for residue in chain:
        if int(residue.id[1]) == int(resseq):
            return residue
    return None


def norm_atom_name(name: str) -> str:
    name = str(name).strip().upper()
    return {"OP1": "O1P", "OP2": "O2P", "OP3": "O3P"}.get(name, name)


def atom_coords(residue) -> Dict[str, np.ndarray]:
    coords = {}
    if residue is None:
        return coords
    for atom in residue.get_atoms():
        coords[norm_atom_name(atom.get_name())] = np.array(atom.coord, dtype=float)
    return coords


def phosphate_coords(residue) -> Dict[str, np.ndarray]:
    coords = {}
    for name, coord in atom_coords(residue).items():
        if name in PHOSPHATE_ATOMS:
            coords[name] = coord
    return coords


def calc_dihedral(p1, p2, p3, p4) -> Optional[float]:
    b1, b2, b3 = p2 - p1, p3 - p2, p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    if np.linalg.norm(n1) < 1e-8 or np.linalg.norm(n2) < 1e-8 or np.linalg.norm(b2) < 1e-8:
        return None
    n1 = n1 / np.linalg.norm(n1)
    n2 = n2 / np.linalg.norm(n2)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    return float(np.degrees(np.arctan2(-np.dot(m1, n2), np.dot(n1, n2))))


def measure_torsion(residue, restype: str) -> Optional[float]:
    names = DIHEDRAL_ATOMS.get(restype)
    if residue is None or names is None:
        return None
    coords = atom_coords(residue)
    if any(name not in coords for name in names):
        return None
    return calc_dihedral(*(coords[name] for name in names))


def circular_delta(a: Optional[float], b: Optional[float]) -> float:
    if a is None or b is None or math.isnan(float(a)) or math.isnan(float(b)):
        return float("nan")
    return abs(((float(a) - float(b) + 180.0) % 360.0) - 180.0)


def frame_index(angle: Optional[float]) -> Optional[int]:
    if angle is None or math.isnan(float(angle)):
        return None
    normalized = float(angle) % 360.0
    return int(round(normalized / 30.0)) % 12


def phosphate_oxygens(phos: Dict[str, np.ndarray]) -> List[np.ndarray]:
    return [phos[name] for name in ("O1P", "O2P", "O3P") if name in phos]


def residue_label(chain_id: str, resname: str, resseq: int) -> str:
    return f"{chain_id}:{resname}{int(resseq)}"


def basic_donor_distances(model, phos: Dict[str, np.ndarray], exclude: Tuple[str, int] | None = None) -> List[Dict]:
    oxy = phosphate_oxygens(phos)
    if not oxy:
        return []
    rows = []
    for chain in model:
        for residue in chain:
            resname = residue.get_resname().strip().upper()
            donors = BASIC_DONORS.get(resname)
            if not donors:
                continue
            resseq = int(residue.id[1])
            if exclude and chain.id == exclude[0] and resseq == int(exclude[1]):
                continue
            for atom in residue.get_atoms():
                atom_name = atom.get_name().strip().upper()
                if atom_name not in donors:
                    continue
                coord = np.array(atom.coord, dtype=float)
                for oxygen_name, oxygen_coord in [(k, phos[k]) for k in ("O1P", "O2P", "O3P") if k in phos]:
                    rows.append({
                        "chain": chain.id,
                        "resname": resname,
                        "resseq": resseq,
                        "residue": residue_label(chain.id, resname, resseq),
                        "donor_atom": atom_name,
                        "phosphate_atom": oxygen_name,
                        "distance": float(np.linalg.norm(coord - oxygen_coord)),
                    })
    return rows


def nearest_basic_distance(model, phos: Dict[str, np.ndarray], exclude: Tuple[str, int] | None = None) -> Tuple[float, str]:
    rows = basic_donor_distances(model, phos, exclude=exclude)
    if not rows:
        return float("nan"), ""
    best = min(rows, key=lambda r: r["distance"])
    return best["distance"], f"{best['residue']}:{best['donor_atom']}-{best['phosphate_atom']}"


def reference_contacts(model, phos: Dict[str, np.ndarray], site_chain: str, site_pos: int, cutoff: float = 4.0) -> List[Tuple[str, str, int]]:
    contacts = basic_donor_distances(model, phos, exclude=(site_chain, site_pos))
    residues = []
    seen = set()
    for row in sorted(contacts, key=lambda r: r["distance"]):
        if row["distance"] <= cutoff:
            key = (row["chain"], row["resname"], int(row["resseq"]))
            if key not in seen:
                seen.add(key)
                residues.append(key)
    return residues


def chain_local_index(model, chain_id: str, original_resseq: int) -> Optional[int]:
    if chain_id not in model:
        return None
    idx = 0
    for residue in model[chain_id]:
        if residue.id[0] != " ":
            continue
        idx += 1
        if int(residue.id[1]) == int(original_resseq):
            return idx
    return None


def map_contacts_for_output(contacts: Sequence[Tuple[str, str, int]], reference_model, tool: str) -> List[Tuple[str, str, int]]:
    if tool != "PTM-Psi":
        return list(contacts)
    mapped = []
    for chain_id, resname, resseq in contacts:
        local = chain_local_index(reference_model, chain_id, resseq)
        if local is not None:
            mapped.append((chain_id, resname, local))
    return mapped


def target_contact_distance(model, phos: Dict[str, np.ndarray], contacts: Sequence[Tuple[str, str, int]]) -> Tuple[float, str, List[str]]:
    oxy = phosphate_oxygens(phos)
    if not oxy:
        return float("nan"), "", []
    hits = []
    missing = []
    best_distance = float("nan")
    best_detail = ""
    for chain_id, resname, resseq in contacts:
        residue = find_residue(model, chain_id, resseq)
        if residue is None:
            missing.append(residue_label(chain_id, resname, resseq))
            continue
        actual_resname = residue.get_resname().strip().upper()
        donors = BASIC_DONORS.get(actual_resname)
        if not donors:
            missing.append(f"{residue_label(chain_id, actual_resname, resseq)}(expected_{resname};nonbasic)")
            continue
        for atom in residue.get_atoms():
            atom_name = atom.get_name().strip().upper()
            if atom_name not in donors:
                continue
            coord = np.array(atom.coord, dtype=float)
            for oxygen_name, oxygen_coord in [(k, phos[k]) for k in ("O1P", "O2P", "O3P") if k in phos]:
                dist = float(np.linalg.norm(coord - oxygen_coord))
                if math.isnan(best_distance) or dist < best_distance:
                    best_distance = dist
                    best_detail = f"{residue_label(chain_id, residue.get_resname().strip().upper(), resseq)}:{atom_name}-{oxygen_name}"
                if dist <= 4.0:
                    hits.append(f"{residue_label(chain_id, residue.get_resname().strip().upper(), resseq)}:{atom_name}-{oxygen_name}:{dist:.2f}")
    return best_distance, best_detail, hits or missing


def find_case_row(df: pd.DataFrame, case: Dict, pdb_id: str, restype: str) -> Optional[pd.Series]:
    if df.empty:
        return None
    work = df.copy()
    for col in ["acc_id", "pdb_id", "chain", "position", "status"]:
        if col not in work.columns:
            return None
    mask = (
        (work["acc_id"].str.upper() == case["acc_id"].upper())
        & (work["pdb_id"].str.lower() == pdb_id.lower())
        & (work["chain"].astype(str) == case["chain"])
        & (work["position"].astype(str) == str(case["position"]))
        & (work["status"].astype(str).str.upper() == "OK")
    )
    rows = work[mask].copy()
    if rows.empty:
        return None
    rows["_restype_3"] = rows.get("restype", "").apply(norm_restype)
    if restype == "PTR":
        # The corrected mini TSV has Y/PTR. The original PTR_final_res rows for ERK2 Y187 say T.
        if any(rows["_restype_3"] == "PTR"):
            rows = rows[rows["_restype_3"] == "PTR"]
    else:
        rows = rows[rows["_restype_3"] == restype]
    return None if rows.empty else rows.iloc[0]


def find_ptmpsi_row(df: pd.DataFrame, case: Dict, pdb_id: str) -> Optional[pd.Series]:
    if df.empty or "pdb_id" not in df.columns:
        return None
    mask = (
        (df["status"].astype(str).str.upper() == "OK")
        & (df["acc_id"].str.upper() == case["acc_id"].upper())
        & (df["pdb_id"].str.lower() == pdb_id.lower())
        & (df["position"].astype(str) == str(case["position"]))
        & (df["restype_3"].str.upper() == case["residue_type"])
    )
    rows = df[mask]
    return None if rows.empty else rows.iloc[0]


def find_pytms_row(df: pd.DataFrame, case: Dict, pdb_id: str) -> Optional[pd.Series]:
    if df.empty or "reference_path" not in df.columns:
        return None
    work = df.copy()
    work["_pdb_id"] = work["reference_path"].apply(pdb_from_reference_path)
    mask = (
        (work["status"].astype(str).str.upper() == "OK")
        & (work["acc_id"].str.upper() == case["acc_id"].upper())
        & (work["_pdb_id"].str.lower() == pdb_id.lower())
        & (work["position"].astype(str) == str(case["position"]))
        & (work["restype_3"].str.upper() == case["residue_type"])
    )
    rows = work[mask]
    return None if rows.empty else rows.iloc[0]


def best_of_three_path(row: pd.Series) -> Tuple[float, str]:
    ranks = [
        (safe_float(row.get("top3_rank1_rmsd")), row.get("archive_stage3_top1_structure", "")),
        (safe_float(row.get("top3_rank2_rmsd")), row.get("archive_stage3_top2_structure", "")),
        (safe_float(row.get("top3_rank3_rmsd")), row.get("archive_stage3_top3_structure", "")),
    ]
    valid = [(rmsd, path) for rmsd, path in ranks if not math.isnan(rmsd) and str(path).strip()]
    if not valid:
        return safe_float(row.get("top3_best_rmsd")), row.get("archive_stage3_top1_structure", "")
    return min(valid, key=lambda item: item[0])


def case_key(case: Dict, pdb_id: str) -> str:
    return f"{case['protein']} {case['residue_type']} {case['position']} {pdb_id.upper()}"


def make_record(
    case: Dict,
    pdb_id: str,
    tool: str,
    variant: str,
    rmsd: float,
    p_only_dist: float,
    structure_path: str,
    source_file: str,
    reference_model,
    reference_torsion: Optional[float],
    reference_nearest: float,
    reference_nearest_detail: str,
    target_contacts: Sequence[Tuple[str, str, int]],
    expected_contacts: Sequence[Tuple[str, str, int]],
    site_query: Tuple[str, int],
) -> Dict:
    resolved = resolve_path(structure_path)
    exists = resolved.exists() if resolved is not None else False
    torsion_pred = None
    nearest_dist = float("nan")
    nearest_detail = ""
    target_min = float("nan")
    target_detail = ""
    target_hits: List[str] = []
    expected_min = float("nan")
    expected_detail = ""
    expected_hits: List[str] = []
    notes = []

    if not exists:
        notes.append("missing_structure_path")
    else:
        try:
            model = first_model(resolved)
            site_res = find_residue(model, site_query[0], site_query[1])
            phos = phosphate_coords(site_res)
            if "P" not in phos:
                notes.append("site_phosphate_not_found")
            torsion_pred = measure_torsion(site_res, case["residue_type"])
            nearest_dist, nearest_detail = nearest_basic_distance(model, phos, exclude=site_query)
            mapped_target = map_contacts_for_output(target_contacts, reference_model, tool)
            mapped_expected = map_contacts_for_output(expected_contacts, reference_model, tool)
            target_min, target_detail, target_hits = target_contact_distance(model, phos, mapped_target)
            expected_min, expected_detail, expected_hits = target_contact_distance(model, phos, mapped_expected)
        except Exception as exc:
            notes.append(f"geometry_error:{exc}")

    torsion_error = circular_delta(torsion_pred, reference_torsion)
    return {
        "protein": case["protein"],
        "site_label": case["site_label"],
        "case_id": case_key(case, pdb_id),
        "pdb_id": pdb_id.lower(),
        "chain": case["chain"],
        "position": case["position"],
        "residue_type": case["residue_type"],
        "tool": tool,
        "variant": variant,
        "rmsd": rmsd,
        "p_only_dist": p_only_dist,
        "torsion_ref": reference_torsion,
        "torsion_pred": torsion_pred,
        "torsion_error": torsion_error,
        "torsion_ref_frame": frame_index(reference_torsion),
        "torsion_pred_frame": frame_index(torsion_pred),
        "torsion_frame_match": frame_index(reference_torsion) == frame_index(torsion_pred)
        if frame_index(reference_torsion) is not None and frame_index(torsion_pred) is not None
        else "",
        "nearest_basic_dist": nearest_dist,
        "nearest_basic_detail": nearest_detail,
        "reference_nearest_basic_dist": reference_nearest,
        "reference_nearest_basic_detail": reference_nearest_detail,
        "nearest_basic_delta": nearest_dist - reference_nearest
        if not math.isnan(nearest_dist) and not math.isnan(reference_nearest)
        else float("nan"),
        "salt_bridge_recovered": target_min <= 4.0 if not math.isnan(target_min) else (False if target_contacts else ""),
        "salt_bridge_min_dist": target_min,
        "salt_bridge_detail": target_detail,
        "salt_bridge_contacts_or_missing": ";".join(target_hits),
        "expected_contact_recovered": expected_min <= 4.0 if not math.isnan(expected_min) else (False if expected_contacts else ""),
        "expected_contact_min_dist": expected_min,
        "expected_contact_detail": expected_detail,
        "expected_contacts_or_missing": ";".join(expected_hits),
        "target_contact_set": ";".join(residue_label(*c) for c in target_contacts),
        "expected_contact_set": ";".join(residue_label(*c) for c in expected_contacts),
        "structure_path": rel_or_abs(structure_path),
        "structure_exists": exists,
        "reference_path": str(PDB_CACHE / f"{pdb_id.lower()}.cif"),
        "source_file": source_file,
        "notes": ";".join(notes),
    }


def svg_text(x, y, text, size=12, anchor="middle", weight="normal", fill="#222", rotate=None):
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" font-family="Arial, sans-serif" font-size="{size}" '
        f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}"{transform}>{html.escape(str(text))}</text>'
    )


def write_svg(path: Path, width: int, height: int, body: Iterable[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    content = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        *body,
        "</svg>",
    ]
    path.write_text("\n".join(content), encoding="utf-8")


def write_grouped_bar_svg(df: pd.DataFrame, path: Path, variants: Sequence[str], title: str):
    work = df[df["variant"].isin(variants)].copy()
    case_ids = list(dict.fromkeys(work["case_id"]))
    max_y = max(0.5, float(pd.to_numeric(work["rmsd"], errors="coerce").max()) * 1.15)
    width = max(980, 120 + len(case_ids) * 130)
    height = 560
    left, top, plot_h = 80, 60, 350
    plot_w = width - 140
    group_w = plot_w / max(1, len(case_ids))
    bar_w = min(18, group_w / (len(variants) + 2))
    body = [svg_text(width / 2, 28, title, size=18, weight="bold")]
    body.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333"/>')
    body.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333"/>')
    for tick in np.linspace(0, max_y, 6):
        y = top + plot_h - (tick / max_y) * plot_h
        body.append(f'<line x1="{left-5}" y1="{y:.1f}" x2="{left+plot_w}" y2="{y:.1f}" stroke="#e8e8e8"/>')
        body.append(svg_text(left - 10, y + 4, f"{tick:.1f}", anchor="end", size=11))
    for i, cid in enumerate(case_ids):
        cx = left + group_w * i + group_w / 2
        sub = work[work["case_id"] == cid]
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            if row.empty:
                continue
            value = safe_float(row.iloc[0]["rmsd"])
            if math.isnan(value):
                continue
            x = cx - (len(variants) * bar_w) / 2 + j * bar_w
            h = (value / max_y) * plot_h
            y = top + plot_h - h
            body.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w-2:.1f}" height="{h:.1f}" fill="{COLORS[variant]}"/>')
        label = cid.replace(" ", "\n")
        parts = label.split("\n")
        for k, part in enumerate(parts):
            body.append(svg_text(cx, top + plot_h + 24 + k * 14, part, size=10))
    legend_x = left + 10
    legend_y = height - 80
    for i, variant in enumerate(variants):
        x = legend_x + i * 190
        body.append(f'<rect x="{x}" y="{legend_y}" width="14" height="14" fill="{COLORS[variant]}"/>')
        body.append(svg_text(x + 20, legend_y + 12, variant, size=12, anchor="start"))
    body.append(svg_text(22, top + plot_h / 2, "RMSD (A)", size=13, rotate=-90))
    write_svg(path, width, height, body)


def write_torsion_svg(df: pd.DataFrame, path: Path):
    variants = PRIMARY_VARIANTS
    work = df[df["variant"].isin(variants)].copy()
    case_ids = list(dict.fromkeys(work["case_id"]))
    width = max(1000, 120 + len(case_ids) * 140)
    height = 560
    left, top, plot_h = 80, 60, 350
    plot_w = width - 140
    group_w = plot_w / max(1, len(case_ids))

    def y_for(angle):
        return top + plot_h - ((float(angle) + 180.0) / 360.0) * plot_h

    body = [svg_text(width / 2, 28, "Predicted phosphate torsion vs crystal", size=18, weight="bold")]
    body.append(f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333"/>')
    body.append(f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333"/>')
    for tick in [-180, -120, -60, 0, 60, 120, 180]:
        y = y_for(tick)
        body.append(f'<line x1="{left-5}" y1="{y:.1f}" x2="{left+plot_w}" y2="{y:.1f}" stroke="#e8e8e8"/>')
        body.append(svg_text(left - 10, y + 4, str(tick), anchor="end", size=11))
    for i, cid in enumerate(case_ids):
        cx = left + group_w * i + group_w / 2
        sub = work[work["case_id"] == cid]
        ref = safe_float(sub.iloc[0]["torsion_ref"]) if not sub.empty else float("nan")
        if not math.isnan(ref):
            y = y_for(ref)
            body.append(f'<line x1="{cx-group_w*0.34:.1f}" y1="{y:.1f}" x2="{cx+group_w*0.34:.1f}" y2="{y:.1f}" stroke="{COLORS["Reference"]}" stroke-width="2"/>')
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            if row.empty:
                continue
            value = safe_float(row.iloc[0]["torsion_pred"])
            if math.isnan(value):
                continue
            x = cx - 36 + j * 24
            y = y_for(value)
            body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" fill="{COLORS[variant]}"/>')
        for k, part in enumerate(cid.replace(" ", "\n").split("\n")):
            body.append(svg_text(cx, top + plot_h + 24 + k * 14, part, size=10))
    legend_y = height - 82
    legend = ["Reference"] + list(variants)
    for i, name in enumerate(legend):
        x = left + i * 160
        if name == "Reference":
            body.append(f'<line x1="{x}" y1="{legend_y+7}" x2="{x+18}" y2="{legend_y+7}" stroke="{COLORS[name]}" stroke-width="2"/>')
        else:
            body.append(f'<circle cx="{x+8}" cy="{legend_y+7}" r="6" fill="{COLORS[name]}"/>')
        body.append(svg_text(x + 24, legend_y + 12, name, size=12, anchor="start"))
    body.append(svg_text(22, top + plot_h / 2, "Torsion angle (deg)", size=13, rotate=-90))
    write_svg(path, width, height, body)


def write_salt_table_svg(df: pd.DataFrame, path: Path):
    variants = PRIMARY_VARIANTS
    cases = list(dict.fromkeys(df["case_id"]))
    row_h = 34
    col_w = 150
    width = 330 + len(variants) * col_w
    height = 90 + len(cases) * row_h
    body = [svg_text(width / 2, 28, "Salt bridge recovery (target contact <= 4.0 A)", size=18, weight="bold")]
    x0, y0 = 20, 58
    body.append(svg_text(x0, y0, "Case", anchor="start", weight="bold"))
    for j, variant in enumerate(variants):
        body.append(svg_text(300 + j * col_w, y0, variant, weight="bold", size=11))
    for i, cid in enumerate(cases):
        y = y0 + 26 + i * row_h
        body.append(svg_text(x0, y, cid, anchor="start", size=11))
        sub = df[df["case_id"] == cid]
        for j, variant in enumerate(variants):
            row = sub[sub["variant"] == variant]
            x = 235 + j * col_w
            recovered = False
            label = "NA"
            fill = "#f2f2f2"
            if not row.empty:
                val = str(row.iloc[0]["salt_bridge_recovered"]).lower()
                recovered = val == "true"
                dist = safe_float(row.iloc[0]["salt_bridge_min_dist"])
                label = ("yes" if recovered else "no") + (f" ({dist:.2f})" if not math.isnan(dist) else "")
                fill = "#cfead1" if recovered else "#f4c7c3"
            body.append(f'<rect x="{x}" y="{y-18}" width="{col_w-10}" height="24" fill="{fill}" stroke="#dddddd"/>')
            body.append(svg_text(x + (col_w - 10) / 2, y - 2, label, size=11))
    write_svg(path, width, height, body)


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    tpo_pf = read_tsv("TPO_final_res.tsv")
    ptr_pf = read_tsv("case_study_inputs/PTR_P28482_ERK2_Y187.tsv")
    ptmpsi_full = read_tsv("ptmpsi_detail.tsv")
    ptmpsi_ptr = read_tsv("case_study_outputs/ptmpsi_detail_ERK2_Y187.tsv")
    pytms_default_full = read_tsv("pytms_detail.tsv")
    pytms_default_ptr = read_tsv("case_study_outputs/pytms_detail_ERK2_Y187.tsv")
    pytms_opt_full = read_tsv("pytms_detail_optimized.tsv")
    pytms_opt_ptr = read_tsv("case_study_outputs/pytms_rotation_detail_ERK2_Y187.tsv")

    records = []
    reference_rows = []

    for case in CASES:
        pf_df = ptr_pf if case["residue_type"] == "PTR" else tpo_pf
        for pdb_id in case["pdb_ids"]:
            ref_path = PDB_CACHE / f"{pdb_id}.cif"
            ref_model = first_model(ref_path)
            ref_res = find_residue(ref_model, case["chain"], case["position"])
            ref_phos = phosphate_coords(ref_res)
            ref_torsion = measure_torsion(ref_res, case["residue_type"])
            ref_nearest, ref_nearest_detail = nearest_basic_distance(ref_model, ref_phos, exclude=(case["chain"], case["position"]))
            observed_contacts = reference_contacts(ref_model, ref_phos, case["chain"], case["position"])
            expected_contacts = list(case["expected_contacts"])
            target_contacts = []
            seen = set()
            observed_arg_contacts = [contact for contact in observed_contacts if contact[1] == "ARG"]
            contacts_to_add = expected_contacts + (observed_arg_contacts if case.get("include_observed_arg_contacts") else [])
            for contact in contacts_to_add:
                if contact not in seen:
                    seen.add(contact)
                    target_contacts.append(contact)
            reference_rows.append({
                "protein": case["protein"],
                "site_label": case["site_label"],
                "case_id": case_key(case, pdb_id),
                "pdb_id": pdb_id,
                "chain": case["chain"],
                "position": case["position"],
                "residue_type": case["residue_type"],
                "reference_path": str(ref_path),
                "reference_torsion": ref_torsion,
                "reference_nearest_basic_dist": ref_nearest,
                "reference_nearest_basic_detail": ref_nearest_detail,
                "expected_contact_set": ";".join(residue_label(*c) for c in expected_contacts),
                "observed_reference_contact_set": ";".join(residue_label(*c) for c in observed_contacts),
                "observed_reference_arg_contact_set": ";".join(residue_label(*c) for c in observed_arg_contacts),
                "target_contact_set": ";".join(residue_label(*c) for c in target_contacts),
            })

            pf_row = find_case_row(pf_df, case, pdb_id, case["residue_type"])
            if pf_row is not None:
                records.append(make_record(
                    case, pdb_id, "PhosphoFill", "PhosphoFill top1",
                    safe_float(pf_row.get("stage3_post_minimization_rmsd")),
                    safe_float(pf_row.get("p_displacement")),
                    pf_row.get("archive_stage3_top1_structure", pf_row.get("archive_stage3_structure", "")),
                    case["pf_file"], ref_model, ref_torsion, ref_nearest, ref_nearest_detail,
                    target_contacts, expected_contacts, (case["chain"], case["position"]),
                ))
                best_rmsd, best_path = best_of_three_path(pf_row)
                records.append(make_record(
                    case, pdb_id, "PhosphoFill", "PhosphoFill best-of-3",
                    best_rmsd,
                    float("nan"),
                    best_path,
                    case["pf_file"], ref_model, ref_torsion, ref_nearest, ref_nearest_detail,
                    target_contacts, expected_contacts, (case["chain"], case["position"]),
                ))

            ptmpsi_df = ptmpsi_ptr if case["residue_type"] == "PTR" else ptmpsi_full
            ptmpsi_row = find_ptmpsi_row(ptmpsi_df, case, pdb_id)
            if ptmpsi_row is not None:
                site_chain = str(ptmpsi_row.get("ptmpsi_chain", case["chain"])) or case["chain"]
                site_pos = int(float(ptmpsi_row.get("ptmpsi_residue_index", case["position"])))
                records.append(make_record(
                    case, pdb_id, "PTM-Psi", "PTM-Psi",
                    safe_float(ptmpsi_row.get("ptmpsi_phosphate_sym_rmsd")),
                    safe_float(ptmpsi_row.get("ptmpsi_p_only_dist")),
                    ptmpsi_row.get("ptmpsi_output", ""),
                    "case_study_outputs/ptmpsi_detail_ERK2_Y187.tsv" if case["residue_type"] == "PTR" else "ptmpsi_detail.tsv",
                    ref_model, ref_torsion, ref_nearest, ref_nearest_detail,
                    target_contacts, expected_contacts, (site_chain, site_pos),
                ))

            for variant, source_df, source_file in [
                ("PyTMs optimized", pytms_opt_ptr if case["residue_type"] == "PTR" else pytms_opt_full,
                 "case_study_outputs/pytms_rotation_detail_ERK2_Y187.tsv" if case["residue_type"] == "PTR" else "pytms_detail_optimized.tsv"),
                ("PyTMs default", pytms_default_ptr if case["residue_type"] == "PTR" else pytms_default_full,
                 "case_study_outputs/pytms_detail_ERK2_Y187.tsv" if case["residue_type"] == "PTR" else "pytms_detail.tsv"),
            ]:
                pytms_row = find_pytms_row(source_df, case, pdb_id)
                if pytms_row is None:
                    continue
                records.append(make_record(
                    case, pdb_id, "PyTMs", variant,
                    safe_float(pytms_row.get("pytms_phosphate_sym_rmsd")),
                    safe_float(pytms_row.get("pytms_p_only_dist")),
                    pytms_row.get("pytms_output", ""),
                    source_file,
                    ref_model, ref_torsion, ref_nearest, ref_nearest_detail,
                    target_contacts, expected_contacts, (case["chain"], case["position"]),
                ))

    metrics = pd.DataFrame(records)
    references = pd.DataFrame(reference_rows)
    variant_order = {variant: i for i, variant in enumerate(SUPPLEMENTARY_VARIANTS)}
    metrics["_variant_order"] = metrics["variant"].map(variant_order).fillna(99)
    metrics.sort_values(["protein", "position", "pdb_id", "_variant_order", "variant"], inplace=True)
    metrics.drop(columns=["_variant_order"], inplace=True)

    metrics.to_csv(OUT_DIR / "case_study_metrics_long.tsv", sep="\t", index=False)
    references.to_csv(OUT_DIR / "case_study_reference_contacts.tsv", sep="\t", index=False)
    metrics[[
        "protein", "site_label", "pdb_id", "residue_type", "tool", "variant",
        "salt_bridge_recovered", "salt_bridge_min_dist", "salt_bridge_detail",
        "salt_bridge_contacts_or_missing", "expected_contact_recovered", "expected_contact_min_dist",
        "expected_contacts_or_missing", "target_contact_set", "expected_contact_set", "structure_path",
    ]].to_csv(OUT_DIR / "case_study_salt_bridge_recovery.tsv", sep="\t", index=False)
    metrics[[
        "protein", "site_label", "pdb_id", "residue_type", "tool", "variant",
        "reference_nearest_basic_dist", "nearest_basic_dist", "nearest_basic_delta",
        "reference_nearest_basic_detail", "nearest_basic_detail", "structure_path",
    ]].to_csv(OUT_DIR / "case_study_nearest_basic_distance.tsv", sep="\t", index=False)
    metrics[[
        "protein", "site_label", "pdb_id", "residue_type", "tool", "variant",
        "torsion_ref", "torsion_pred", "torsion_error", "torsion_ref_frame", "torsion_pred_frame",
        "torsion_frame_match", "structure_path",
    ]].to_csv(OUT_DIR / "case_study_torsion_comparison.tsv", sep="\t", index=False)
    metrics[metrics["variant"].isin(PRIMARY_VARIANTS)].to_csv(OUT_DIR / "case_study_metrics_primary.tsv", sep="\t", index=False)
    metrics[metrics["variant"].isin(SUPPLEMENTARY_VARIANTS)].to_csv(OUT_DIR / "case_study_metrics_with_pytms_default.tsv", sep="\t", index=False)

    manifest_cols = [
        "protein", "site_label", "case_id", "pdb_id", "chain", "position", "residue_type",
        "tool", "variant", "structure_path", "reference_path", "structure_exists",
    ]
    metrics[manifest_cols].to_csv(OUT_DIR / "case_study_visual_structure_manifest.tsv", sep="\t", index=False)

    write_grouped_bar_svg(
        metrics[metrics["variant"].isin(PRIMARY_VARIANTS)],
        FIG_DIR / "rmsd_grouped_bar_primary.svg",
        PRIMARY_VARIANTS,
        "Case-study phosphate RMSD by tool",
    )
    write_grouped_bar_svg(
        metrics[metrics["variant"].isin(SUPPLEMENTARY_VARIANTS)],
        FIG_DIR / "rmsd_grouped_bar_with_pytms_default.svg",
        SUPPLEMENTARY_VARIANTS,
        "Case-study phosphate RMSD by tool, including PyTMs default",
    )
    write_torsion_svg(metrics, FIG_DIR / "torsion_mode_comparison.svg")
    write_salt_table_svg(metrics[metrics["variant"].isin(PRIMARY_VARIANTS)], FIG_DIR / "salt_bridge_recovery_table.svg")

    primary = metrics[metrics["variant"].isin(PRIMARY_VARIANTS)]
    print(f"Wrote {len(metrics)} total metric rows")
    print(f"Wrote {len(primary)} primary metric rows")
    print(f"All primary structures exist: {bool(primary['structure_exists'].all())}")
    print(f"Outputs: {OUT_DIR}")


if __name__ == "__main__":
    main()
