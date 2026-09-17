#!/usr/bin/env python3
"""Independently compare extractor torsions with Bio.PDB on source CIF atoms."""

from __future__ import annotations

import argparse
from math import degrees
from pathlib import Path

import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.vectors import calc_dihedral


SPECS = (
    ("mod", "TPO", ("N", "CA", "CB", "OG1"), "torsion_N_CA_CB_OG1"),
    ("mod", "TPO", ("CG2", "CB", "OG1", "P"), "torsion_CG2_CB_OG1_P"),
    ("mod", "SEP", ("CA", "CB", "OG", "P"), "torsion_CA_CB_OG_P"),
    ("mod", "PTR", ("CE1", "CZ", "OH", "P"), "torsion_CE1_CZ_OH_P"),
    ("unmod", "THR", ("N", "CA", "CB", "OG1"), "torsion_N_CA_CB_OG1"),
    ("unmod", "SER", ("N", "CA", "CB", "OG"), "torsion_N_CA_CB_OG"),
)


def wrap(angle: float) -> float:
    return ((angle + 180.0) % 360.0) - 180.0


def main() -> None:
    benchmarkv2 = Path(__file__).resolve().parents[4]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phospho-root", type=Path, default=benchmarkv2.parent)
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).resolve().parents[1] / "outputs"
                        / "torsion_convention_examples.tsv")
    args = parser.parse_args()
    source = args.phospho_root / "Validation-allPDB"
    parser_cif = MMCIFParser(QUIET=True)
    rows = []
    for kind, residue, atoms, column in SPECS:
        path = source / f"{kind}_geometry_prior_all" / f"{residue}_geometry_prior.tsv"
        data = pd.read_csv(path, sep="\t", low_memory=False)
        data = data.loc[data["status"].eq("ok") & data["resname"].eq(residue)]
        successes = 0
        for record in data.itertuples(index=False):
            if successes >= 5:
                break
            cif = source / "pdb_cache" / f"{record.pdb_id.lower()}.cif"
            if not cif.exists():
                continue
            try:
                structure = parser_cif.get_structure("check", str(cif))
                chain = structure[0][record.chain_id]
                target = next(r for r in chain if r.id[1] == int(record.target_position)
                              and r.get_resname() == residue)
                bio = degrees(calc_dihedral(*[target[name].get_vector() for name in atoms]))
            except (KeyError, StopIteration, ValueError):
                continue
            extracted = float(getattr(record, column))
            predicted_bio = wrap(180.0 - extracted)
            error = wrap(bio - predicted_bio)
            if abs(error) > 1e-4:
                raise ValueError(f"Torsion convention mismatch for {record.pdb_id} {residue}: {error}")
            rows.append({"residue": residue, "pdb_id": record.pdb_id,
                         "chain": record.chain_id, "position": int(record.target_position),
                         "torsion": column, "extractor_deg": extracted,
                         "biopython_deg": bio, "wrap_180_minus_extractor_deg": predicted_bio,
                         "circular_error_deg": error})
            successes += 1
        if successes != 5:
            raise ValueError(f"Only {successes} CIF records could be checked for {residue}/{column}")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.out, sep="\t", index=False)
    print(f"Verified {len(rows)} direct-coordinate torsions: {args.out}")


if __name__ == "__main__":
    main()
