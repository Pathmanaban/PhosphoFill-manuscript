#!/usr/bin/env python3
"""
04_collect_rich_reports.py
==========================
Walks pf_runs/ and extracts rich metrics from each
*_PF{rank}_site_report.tsv, then merges into proteome_results.tsv.

Usage
-----
    python3 04_collect_rich_reports.py \
        --results  proteome_results.tsv \
        --pf-runs  pf_runs \
        --output   proteome_results_enriched.tsv
"""

import argparse
from pathlib import Path
import pandas as pd

# Key columns from the rich report TSV
RICH_COLS = [
    # Torsion
    "torsion_CA_CB_OG_P",       # SEP key torsion
    "torsion_CG2_CB_OG1_P",     # TPO key torsion
    "torsion_CE1_CZ_OH_P",      # PTR key torsion
    # Contacts
    "nearby_basic_count_after",
    "nearby_acidic_count_after",
    "salt_bridge_count_after",
    "salt_bridge_details_after",
    "salt_bridge_count_before",
    "polar_contact_count_after",
    "hydrogen_bond_count_after",
    "electrostatic_score_weighted_after",
    # Confidence
    "phosphofill_confidence",
    "confidence_tier",
    "phospho_quality_label",
    "intrinsic_site_fit",
    "contextual_site_fit",
    # Structure
    "plddt_site",
    "plddt_window_mean",
    "plddt_class",
    "sasa_unmodified",
    "sasa_modified",
    "sasa_delta",
    "exposure_class_unmodified",
    "exposure_class_modified",
    "local_ca_rmsd",
    # Energy
    "site_repulsion_energy_before",
    "site_repulsion_energy_after",
    "site_repulsion_energy_delta",
    # Geometry
    "d_anchor_p_modified",
    "angle_base_anchor_p_modified",
    # Interpretation
    "overall_interpretation",
    "interpretation_flags",
    "cb_neighborhood_after",
]


def find_tsv(pf_runs: Path, acc_id: str, rank: int) -> Path | None:
    """Find the rich report TSV for a given protein/rank."""
    direct = pf_runs / acc_id / f"{acc_id}_PF{rank}_site_report.tsv"
    if direct.is_file():
        return direct
    matches = sorted((pf_runs / acc_id).glob(f"**/*_PF{rank}_site_report.tsv"))
    if not matches:
        # Backward compatibility with archived pre-release runs.
        matches = sorted((pf_runs / acc_id).glob(f"**/*_pf{rank}_site_report.tsv"))
    return matches[0] if matches else None


def main(results_path: str, pf_runs_dir: str, output_path: str) -> None:
    pf_runs = Path(pf_runs_dir)

    print(f"Loading {results_path}...")
    df = pd.read_csv(results_path, sep="\t")
    print(f"  {len(df):,} rows, {df.acc_id.nunique():,} proteins")

    # Add empty rich columns
    for col in RICH_COLS:
        df[col] = None

    # Build index: (acc_id, resseq, rank) -> row index in df
    idx_map = {}
    for idx, row in df.iterrows():
        key = (row["acc_id"], row["resseq"], row["pose_rank"])
        idx_map[key] = idx

    print("Collecting rich report TSVs...")
    n_found = n_missing = n_rows_merged = 0

    for acc_id in df["acc_id"].unique():
        for rank in df[df["acc_id"] == acc_id]["pose_rank"].unique():
            tsv = find_tsv(pf_runs, acc_id, int(rank))
            if tsv is None:
                n_missing += 1
                continue
            try:
                rich = pd.read_csv(tsv, sep="\t")
            except Exception:
                n_missing += 1
                continue
            n_found += 1
            for _, row in rich.iterrows():
                resseq = int(row.get("residue_number", 0))
                key = (acc_id, resseq, rank)
                if key in idx_map:
                    idx = idx_map[key]
                    for col in RICH_COLS:
                        if col in row:
                            df.at[idx, col] = row[col]
                    n_rows_merged += 1

        if n_found % 1000 == 0 and n_found > 0:
            print(f"  {n_found:,} TSVs processed...", end="\r")

    print(f"\n  Found: {n_found:,}  Missing: {n_missing:,}  Rows merged: {n_rows_merged:,}")

    # Add convenience torsion column: pick the relevant one per residue type
    def get_torsion(row):
        rn = row.get("new_resname", "")
        if rn == "SEP":
            return row.get("torsion_CA_CB_OG_P")
        elif rn == "TPO":
            return row.get("torsion_CG2_CB_OG1_P")
        elif rn == "PTR":
            return row.get("torsion_CE1_CZ_OH_P")
        return None

    df["torsion_key"] = df.apply(get_torsion, axis=1)

    # Summary
    print("\nRich column population (top-1):")
    df1 = df[df["pose_rank"] == 1]
    for col in ["phosphofill_confidence", "salt_bridge_count_after",
                "torsion_key", "plddt_site", "exposure_class_modified",
                "phospho_quality_label"]:
        n = df1[col].notna().sum()
        print(f"  {col}: {n:,} / {len(df1):,} ({100*n/len(df1):.1f}%)")


    df.to_csv(output_path, sep="\t", index=False)
    print(f"\nSaved: {output_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results",  required=True)
    ap.add_argument("--pf-runs",  default="pf_runs")
    ap.add_argument("--output",   default="proteome_results_enriched.tsv")
    args = ap.parse_args()
    main(args.results, args.pf_runs, args.output)
