#!/usr/bin/env python3
"""
05_add_dssp.py
==============
Merges Structure_AF.txt (DSSP SS8, ASA, RSA) into
proteome_results_enriched.tsv on ACC_ID + UniProt_pos.

Usage
-----
    python3 05_add_dssp.py \
        --results  proteome_results_enriched.tsv \
        --dssp     Structure_AF.txt \
        --output   proteome_results_final.tsv
"""

import argparse
import pandas as pd

# SS8 → simplified SS3 mapping
SS8_TO_SS3 = {
    "H": "helix",   # Alpha helix
    "G": "helix",   # 3-10 helix
    "I": "helix",   # Pi helix
    "E": "sheet",   # Beta strand
    "B": "sheet",   # Beta bridge
    "T": "loop",    # Turn
    "S": "loop",    # Bend
    "C": "loop",    # Coil/loop
    "-": "disorder",
}


def main(results_path: str, dssp_path: str, output_path: str) -> None:

    print(f"Loading results: {results_path}")
    df = pd.read_csv(results_path, sep="\t")
    print(f"  {len(df):,} rows, {df.acc_id.nunique():,} proteins")

    print(f"Loading DSSP: {dssp_path}")
    dssp = pd.read_csv(dssp_path, sep="\t")
    print(f"  {len(dssp):,} rows total")

    # Filter to phosphorylation only
    dssp = dssp[dssp["Position_Label"] == "Phosphorylation"].copy()
    print(f"  {len(dssp):,} phosphorylation rows")

    # Rename for merge
    dssp = dssp.rename(columns={
        "ACC_ID":      "acc_id",
        "UniProt_pos": "resseq",
        "SS8":         "ss8",
        "ASA":         "asa",
        "RSA":         "rsa",
        "Residue":     "residue_aa",
    })[["acc_id", "resseq", "ss8", "asa", "rsa", "residue_aa"]]

    dssp["ss3"] = dssp["ss8"].map(SS8_TO_SS3).fillna("loop")
    dssp["surface_exposed"] = dssp["rsa"] >= 0.25  # standard threshold

    # Merge
    before = len(df)
    df = df.merge(dssp, on=["acc_id", "resseq"], how="left")
    print(f"\nMerge: {before:,} → {len(df):,} rows")

    # Check coverage
    df1 = df[df["pose_rank"] == 1]
    n_ss = df1["ss3"].notna().sum()
    print(f"SS coverage (top-1): {n_ss:,}/{len(df1):,} ({100*n_ss/len(df1):.1f}%)")

    print("\nSS3 distribution (top-1):")
    print(df1["ss3"].value_counts())
    print("\nSurface exposed (top-1):")
    print(df1["surface_exposed"].value_counts())

    df.to_csv(output_path, sep="\t", index=False)
    print(f"\nSaved: {output_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--dssp",    required=True)
    ap.add_argument("--output",  default="proteome_results_final.tsv")
    args = ap.parse_args()
    main(args.results, args.dssp, args.output)
