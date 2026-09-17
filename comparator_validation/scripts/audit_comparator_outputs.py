#!/usr/bin/env python3
"""Audit stored PyTMs and PTM-Psi strip-and-regraft benchmark outputs.

The audit deliberately ignores the PhosphoFill metric columns embedded in the
older comparator TSVs.  It checks only comparator provenance and coordinates,
then matches the comparator site keys to the frozen corrected PhosphoFill run.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser


PARENTS = {"TPO": "THR", "SEP": "SER", "PTR": "TYR"}
PHOSPHATE = {"P", "O1P", "O2P", "O3P"}
DIHEDRAL_ATOMS = {
    "TPO": ("CG2", "CB", "OG1", "P"),
    "SEP": ("CA", "CB", "OG", "P"),
    "PTR": ("CE1", "CZ", "OH", "P"),
}
REFERENCE_CACHE: dict[tuple[str, str, int], tuple[dict[str, np.ndarray], str]] = {}


class PDBAtom:
    """Minimal atom adapter for large PyMOL-written PDB files.

    Some archived comparator structures contain hybrid-36 identifiers elsewhere
    in the file.  Bio.PDB aborts while parsing those unrelated records, so the
    audit reads only the requested residue from fixed-width PDB columns.
    """

    def __init__(self, name: str, coord: np.ndarray):
        self._name = name
        self.coord = coord

    def get_name(self) -> str:
        return self._name


class PDBResidue:
    def __init__(self, resname: str, atoms: list[PDBAtom]):
        self.resname = resname
        self._atoms = atoms

    def get_atoms(self):
        return iter(self._atoms)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--phosphofill-tsv", nargs=3, required=True)
    p.add_argument("--pytms-default", required=True)
    p.add_argument("--pytms-optimised", required=True)
    p.add_argument("--ptmpsi", required=True)
    p.add_argument(
        "--comparator-root",
        help="Directory containing pdb_cache, pytms_work, pytms_work_optimized, and ptmpsi_work",
    )
    p.add_argument("--output-json", required=True)
    p.add_argument("--output-tsv", required=True)
    return p.parse_args()


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def norm_restype(value: object) -> str:
    return {"T": "TPO", "S": "SEP", "Y": "PTR"}.get(
        str(value).strip().upper(), str(value).strip().upper()
    )


def site_key(acc_id: object, pdb_id: object, restype: object, chain: object, position: object) -> str:
    return "|".join(
        [
            str(acc_id).strip(),
            str(pdb_id).strip().lower(),
            norm_restype(restype),
            str(chain).strip(),
            str(int(float(position))),
        ]
    )


def load_phosphofill(paths: list[str]) -> pd.DataFrame:
    frames = []
    for raw in paths:
        path = Path(raw)
        primary = next(rt for rt in ("TPO", "SEP", "PTR") if rt in path.stem.upper())
        df = pd.read_csv(path, sep="\t")
        df = df[(df["status"] == "OK") & (df["context_n_sites"] == 1)].copy()
        df["restype_3"] = df["restype"].map(norm_restype)
        df = df[df["restype_3"] == primary].copy()
        df["site_key"] = [
            site_key(a, p, r, c, n)
            for a, p, r, c, n in zip(
                df["acc_id"], df["pdb_id"], df["restype_3"], df["chain"], df["position"]
            )
        ]
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def load_comparator(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    pdb_ids = df["reference_path"].astype(str).map(lambda value: Path(value).stem.lower())
    df["site_key"] = [
        site_key(a, p, r, c, n)
        for a, p, r, c, n in zip(
            df["acc_id"], pdb_ids, df["restype_3"], df["chain"], df["position"]
        )
    ]
    return df


def resolve_path(raw: object, detail_path: Path, comparator_root: Path | None = None) -> Path | None:
    text = str(raw).strip()
    if not text or text.lower() == "nan":
        return None
    candidate = Path(text)
    choices = [candidate] if candidate.is_absolute() else [
        detail_path.parent / candidate,
        *((comparator_root / candidate,) if comparator_root is not None else ()),
        Path.cwd() / candidate,
    ]
    for choice in choices:
        if choice.exists():
            return choice.resolve()
    return choices[0].resolve()


def parser_for(path: Path):
    if path.suffix.lower() in {".cif", ".mmcif"}:
        return MMCIFParser(QUIET=True)
    return PDBParser(QUIET=True)


def find_residue(path: Path, chain_id: str, position: int):
    if path.suffix.lower() not in {".cif", ".mmcif"}:
        atoms = []
        resname = None
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith(("ATOM  ", "HETATM")) or len(line) < 54:
                    continue
                if line[21].strip() != chain_id.strip():
                    continue
                try:
                    resseq = int(line[22:26].strip())
                except ValueError:
                    continue
                if resseq != int(position):
                    continue
                altloc = line[16].strip()
                if altloc not in {"", "A"}:
                    continue
                try:
                    coord = np.array(
                        [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                        dtype=float,
                    )
                except ValueError:
                    continue
                resname = line[17:20].strip().upper()
                atoms.append(PDBAtom(line[12:16].strip(), coord))
        return PDBResidue(resname, atoms) if atoms and resname else None

    structure = parser_for(path).get_structure("x", str(path))
    model = next(structure.get_models())
    if chain_id not in model:
        return None
    matches = [res for res in model[chain_id] if int(res.id[1]) == int(position)]
    if not matches:
        return None
    for residue in matches:
        if residue.resname.strip().upper() in set(PARENTS) | set(PARENTS.values()):
            return residue
    return matches[0]


def cached_reference(path: Path, chain_id: str, position: int) -> tuple[dict[str, np.ndarray], str]:
    key = (str(path), chain_id, int(position))
    if key not in REFERENCE_CACHE:
        residue = find_residue(path, chain_id, position)
        if residue is None:
            raise ValueError("target residue was not found in the reference coordinate file")
        REFERENCE_CACHE[key] = (phosphate_coords(residue), residue.resname.strip().upper())
    return REFERENCE_CACHE[key]


def atom_coord(residue, name: str) -> np.ndarray | None:
    aliases = {"O1P": ("O1P", "OP1", "O1"), "O2P": ("O2P", "OP2", "O2"), "O3P": ("O3P", "OP3", "O3")}
    names = aliases.get(name, (name,))
    for atom in residue.get_atoms():
        if atom.get_name().strip().upper() in names:
            return np.asarray(atom.coord, dtype=float)
    return None


def phosphate_coords(residue) -> dict[str, np.ndarray]:
    if residue is None:
        return {}
    result = {}
    for name in ("P", "O1P", "O2P", "O3P"):
        coord = atom_coord(residue, name)
        if coord is not None:
            result[name] = coord
    return result


def phosphate_rmsd(reference: dict[str, np.ndarray], predicted: dict[str, np.ndarray]) -> float | None:
    if "P" not in reference or "P" not in predicted:
        return None
    ref_ox = [name for name in ("O1P", "O2P", "O3P") if name in reference]
    pred_ox = [name for name in ("O1P", "O2P", "O3P") if name in predicted]
    n_ox = min(len(ref_ox), len(pred_ox))
    if n_ox == 0:
        return float(np.linalg.norm(reference["P"] - predicted["P"]))
    best = np.inf
    for perm in itertools.permutations(pred_ox, n_ox):
        values = [np.sum((reference["P"] - predicted["P"]) ** 2)]
        values.extend(
            np.sum((reference[r] - predicted[p]) ** 2)
            for r, p in zip(ref_ox[:n_ox], perm)
        )
        best = min(best, float(np.sqrt(np.mean(values))))
    return best


def dihedral_degrees(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> float:
    b0 = -(b - a)
    b1 = c - b
    b2 = d - c
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def residue_dihedral(residue, restype: str) -> float | None:
    coords = [atom_coord(residue, name) for name in DIHEDRAL_ATOMS[restype]]
    if any(coord is None for coord in coords):
        return None
    return dihedral_degrees(*coords)


def audit_method(
    name: str,
    detail_path: Path,
    df: pd.DataFrame,
    comparator_root: Path | None = None,
) -> tuple[dict, list[dict]]:
    metric_col = "ptmpsi_phosphate_sym_rmsd" if name == "PTM-Psi" else "pytms_phosphate_sym_rmsd"
    output_col = "ptmpsi_output" if name == "PTM-Psi" else "pytms_output"
    records = []
    default_dihedrals: dict[str, list[float]] = {rt: [] for rt in PARENTS}

    ok_df = df[df["status"] == "OK"]
    for row_number, (_, row) in enumerate(ok_df.iterrows(), start=1):
        reference_path = resolve_path(row["reference_path"], detail_path, comparator_root)
        stripped_path = resolve_path(row["stripped_input"], detail_path, comparator_root)
        output_path = resolve_path(row[output_col], detail_path, comparator_root)
        ref_chain = str(row["chain"]).strip()
        ref_position = int(float(row["position"]))
        restype = norm_restype(row["restype_3"])
        out_chain = str(row.get("ptmpsi_chain", ref_chain)).strip() if name == "PTM-Psi" else ref_chain
        out_position = int(float(row.get("ptmpsi_residue_index", ref_position))) if name == "PTM-Psi" else ref_position

        checks = {
            "reference_exists": bool(reference_path and reference_path.exists()),
            "stripped_exists": bool(stripped_path and stripped_path.exists()),
            "output_exists": bool(output_path and output_path.exists()),
        }
        recomputed = None
        difference = None
        stripped_has_p = None
        output_has_p = None
        output_resname = None
        message = ""
        try:
            if not all(checks.values()):
                raise FileNotFoundError("one or more archived coordinate files are missing")
            reference, _reference_resname = cached_reference(reference_path, ref_chain, ref_position)
            stripped_residue = find_residue(stripped_path, out_chain, out_position)
            output_residue = find_residue(output_path, out_chain, out_position)
            if stripped_residue is None or output_residue is None:
                raise ValueError("target residue was not found in one or more coordinate files")
            stripped_has_p = "P" in phosphate_coords(stripped_residue)
            predicted = phosphate_coords(output_residue)
            output_has_p = "P" in predicted
            output_resname = output_residue.resname.strip().upper()
            recomputed = phosphate_rmsd(reference, predicted)
            recorded = float(row[metric_col])
            difference = None if recomputed is None else abs(recomputed - recorded)
            if name == "PyTMs default":
                dih = residue_dihedral(output_residue, restype)
                if dih is not None:
                    default_dihedrals[restype].append(dih)
        except Exception as exc:
            message = str(exc)

        passed = (
            all(checks.values())
            and stripped_has_p is False
            and output_has_p is True
            and recomputed is not None
            and difference is not None
            and difference <= 0.002
        )
        records.append(
            {
                "method": name,
                "site_key": row["site_key"],
                "restype": restype,
                **checks,
                "stripped_has_phosphorus": stripped_has_p,
                "output_has_phosphorus": output_has_p,
                "output_resname": output_resname,
                "recorded_rmsd_A": row[metric_col],
                "recomputed_rmsd_A": recomputed,
                "absolute_difference_A": difference,
                "passed": passed,
                "message": message,
            }
        )
        if row_number % 100 == 0 or row_number == len(ok_df):
            print(f"{name}: checked {row_number}/{len(ok_df)} archived outputs", flush=True)

    record_df = pd.DataFrame(records)
    dihedral_summary = {}
    for restype, values in default_dihedrals.items():
        arr = np.asarray(values, dtype=float)
        dihedral_summary[restype] = {
            "n": int(arr.size),
            "median_deg": float(np.median(arr)) if arr.size else None,
            "circular_resultant_length": float(np.hypot(np.mean(np.cos(np.radians(arr))), np.mean(np.sin(np.radians(arr))))) if arr.size else None,
        }
    summary = {
        "method": name,
        "requested_rows": int(len(df)),
        "unique_requested_sites": int(df["site_key"].nunique()),
        "ok_rows": int((df["status"] == "OK").sum()),
        "failed_rows": int((df["status"] != "OK").sum()),
        "coordinate_rows_checked": int(len(record_df)),
        "coordinate_rows_passed": int(record_df["passed"].sum()),
        "maximum_rmsd_recalculation_difference_A": float(record_df["absolute_difference_A"].dropna().max()),
        "default_dihedral_summary": dihedral_summary if name == "PyTMs default" else None,
        "failed_messages": df.loc[df["status"] != "OK", ["site_key", "message"]].to_dict("records"),
        "input_sha256": sha256(detail_path),
    }
    return summary, records


def main() -> int:
    args = parse_args()
    comparator_root = Path(args.comparator_root).resolve() if args.comparator_root else None
    pf = load_phosphofill(args.phosphofill_tsv)
    method_specs = [
        ("PyTMs default", Path(args.pytms_default)),
        ("PyTMs optimised", Path(args.pytms_optimised)),
        ("PTM-Psi", Path(args.ptmpsi)),
    ]
    loaded = [(name, path.resolve(), load_comparator(path)) for name, path in method_specs]

    requested_sets = {name: set(df["site_key"]) for name, _, df in loaded}
    pf_set = set(pf["site_key"])
    intersection = set(pf_set)
    for values in requested_sets.values():
        intersection &= values
    ok_intersection = set(pf_set)
    for name, _, df in loaded:
        ok_intersection &= set(df.loc[df["status"] == "OK", "site_key"])

    summaries = []
    records = []
    for name, path, df in loaded:
        summary, detail = audit_method(name, path, df, comparator_root)
        summary["requested_keys_match_phosphofill"] = requested_sets[name] == pf_set
        summaries.append(summary)
        records.extend(detail)

    output_json = Path(args.output_json)
    output_tsv = Path(args.output_tsv)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_csv(output_tsv, sep="\t", index=False)
    report = {
        "phosphofill_single_site_n": int(len(pf)),
        "phosphofill_unique_site_keys": int(pf["site_key"].nunique()),
        "all_requested_site_intersection_n": len(intersection),
        "all_successful_complete_case_n": len(ok_intersection),
        "methods": summaries,
        "conclusion": (
            "Comparator coordinates are independent of the corrected PhosphoFill rank-2/rank-3 "
            "minimisation. A comparator rerun is unnecessary if every requested key matches and "
            "every archived successful coordinate row passes recalculation."
        ),
    }
    output_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"complete_case_n": len(ok_intersection), "methods": summaries}, indent=2))
    return 0 if all(item["coordinate_rows_checked"] == item["coordinate_rows_passed"] for item in summaries) else 1


if __name__ == "__main__":
    raise SystemExit(main())
