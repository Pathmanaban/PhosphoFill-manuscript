#!/usr/bin/env python3
"""Summarise corrected PhosphoFill benchmark outputs and recreate Figure 3.

This script only analyses measured values written by
``phosphofill_benchmark_production.py``.  In particular, it requires three
independently minimised pose RMSDs and never substitutes the old pre-minimised
top-k scan values.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


RESTYPE_MAP = {"T": "TPO", "S": "SEP", "Y": "PTR", "TPO": "TPO", "SEP": "SEP", "PTR": "PTR"}
RESTYPES = ("TPO", "SEP", "PTR")
TYPE_COLORS = {"TPO": "#1565c0", "SEP": "#2e7d32", "PTR": "#6a1b9a"}

REQUIRED_COLUMNS = (
    "status",
    "restype",
    "context_n_sites",
    "stage0_seed_rmsd",
    "stage1_best_of_12_rmsd",
    "stage2_selected_rmsd",
    "rank1_postmin_rmsd",
    "rank2_postmin_rmsd",
    "rank3_postmin_rmsd",
    "top3_best_postmin_rmsd",
    "top3_best_postmin_rank",
)

METRIC_COLUMNS = (
    "stage0_seed_rmsd",
    "stage1_best_of_12_rmsd",
    "stage2_selected_rmsd",
    "rank1_postmin_rmsd",
    "rank2_postmin_rmsd",
    "rank3_postmin_rmsd",
    "top3_best_postmin_rmsd",
    "top3_best_postmin_rank",
)

STAGES = (
    ("S0\nSeed", "stage0_seed_rmsd"),
    ("S1\nOracle\n(best scan)", "stage1_best_of_12_rmsd"),
    ("S2\nSelected\n(pre-min.)", "stage2_selected_rmsd"),
    ("S3\nTop-1\n(post-min.)", "rank1_postmin_rmsd"),
    ("Best of 3\n(post-min.)", "top3_best_postmin_rmsd"),
)

RC = {
    "axes.titlesize": 13,
    "axes.titleweight": "bold",
    "axes.labelsize": 14,
    "axes.labelweight": "bold",
    "xtick.labelsize": 11,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "legend.frameon": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def norm_restype(value: object) -> str:
    key = str(value).strip().upper()
    return RESTYPE_MAP.get(key, key)


def infer_primary_restype(path: str) -> str:
    name = Path(path).stem.upper()
    tokens = set(re.split(r"[^A-Z0-9]+", name))
    for restype in RESTYPES:
        if restype in tokens or name.startswith(restype):
            return restype
    raise ValueError(
        f"Cannot infer the primary residue type from {path!r}; "
        "include TPO, SEP, or PTR in each input filename."
    )


def _site_key_columns(df: pd.DataFrame) -> List[str]:
    candidates = ["acc_id", "pdb_id", "chain", "position", "target_position"]
    return [column for column in candidates if column in df.columns]


def load_outputs(paths: Sequence[str]) -> pd.DataFrame:
    frames: List[pd.DataFrame] = []
    for path in paths:
        df = pd.read_csv(path, sep="\t", low_memory=False)
        missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
        if missing:
            raise ValueError(
                f"{path} is not a corrected three-pose benchmark output; "
                f"missing columns: {', '.join(missing)}"
            )

        primary = infer_primary_restype(path)
        df = df.loc[df["status"].astype(str).str.upper().eq("OK")].copy()
        df["restype_3"] = df["restype"].map(norm_restype)
        df = df.loc[df["restype_3"].eq(primary)].copy()
        df["source_tsv"] = str(Path(path).resolve())
        df["primary_restype"] = primary
        for column in METRIC_COLUMNS:
            df[column] = pd.to_numeric(df[column], errors="coerce")
        df["context_n_sites"] = pd.to_numeric(df["context_n_sites"], errors="coerce")

        incomplete = df[list(METRIC_COLUMNS)].isna().any(axis=1)
        if incomplete.any():
            keys = _site_key_columns(df)
            examples = df.loc[incomplete, keys].head(5).to_dict("records") if keys else []
            raise ValueError(
                f"{path} contains {int(incomplete.sum())} OK rows without all three "
                f"post-minimisation RMSDs. Example rows: {examples}"
            )
        if df["context_n_sites"].isna().any():
            raise ValueError(f"{path} contains OK rows with missing context_n_sites")
        frames.append(df)

    if not frames:
        raise ValueError("No benchmark files were supplied")
    combined = pd.concat(frames, ignore_index=True)
    if combined.empty:
        raise ValueError("No primary-residue OK rows were found in the supplied outputs")

    key_columns = _site_key_columns(combined)
    duplicate_subset = ["primary_restype", *key_columns]
    if key_columns and combined.duplicated(duplicate_subset, keep=False).any():
        examples = combined.loc[
            combined.duplicated(duplicate_subset, keep=False), duplicate_subset
        ].head(10)
        raise ValueError(
            "Duplicate benchmark sites were found across the supplied TSVs. "
            f"Do not mix overlapping full runs and shards. Examples:\n{examples.to_string(index=False)}"
        )
    return combined


def pct(values: pd.Series, threshold: float) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna()
    return float((numeric <= threshold).mean() * 100.0) if len(numeric) else math.nan


def cohort_row(label: str, df: pd.DataFrame) -> Dict[str, object]:
    rank1 = pd.to_numeric(df["rank1_postmin_rmsd"], errors="coerce")
    top3 = pd.to_numeric(df["top3_best_postmin_rmsd"], errors="coerce")
    return {
        "type": label,
        "n": int(len(df)),
        "stage0_seed_median": float(df["stage0_seed_rmsd"].median()),
        "stage1_oracle_median": float(df["stage1_best_of_12_rmsd"].median()),
        "stage2_selected_median": float(df["stage2_selected_rmsd"].median()),
        "top1_postmin_median": float(rank1.median()),
        "top3_best_postmin_median": float(top3.median()),
        "top1_le_1A_pct": pct(rank1, 1.0),
        "top3_le_1A_pct": pct(top3, 1.0),
        "top1_le_1_5A_pct": pct(rank1, 1.5),
        "top3_le_1_5A_pct": pct(top3, 1.5),
        "top3_improved_n": int((top3 < rank1 - 1e-9).sum()),
        "top3_improved_pct": float((top3 < rank1 - 1e-9).mean() * 100.0),
    }


def make_single_summary(single: pd.DataFrame) -> pd.DataFrame:
    rows = [cohort_row(restype, single.loc[single["restype_3"].eq(restype)]) for restype in RESTYPES]
    rows = [row for row in rows if row["n"]]
    rows.append(cohort_row("Combined", single))
    return pd.DataFrame(rows)


def make_context_summary(all_sites: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, object]] = []
    for restype in RESTYPES:
        type_df = all_sites.loc[all_sites["restype_3"].eq(restype)]
        cohorts = (
            ("single", type_df.loc[type_df["context_n_sites"].eq(1)]),
            ("2 sites", type_df.loc[type_df["context_n_sites"].eq(2)]),
            ("3+ sites", type_df.loc[type_df["context_n_sites"].ge(3)]),
        )
        for context, cohort in cohorts:
            if cohort.empty:
                continue
            row = cohort_row(restype, cohort)
            row["context"] = context
            rows.append(row)
    columns = ["type", "context", *[column for column in cohort_row("x", all_sites).keys() if column != "type"]]
    return pd.DataFrame(rows)[columns]


def fmt_number(value: object, digits: int = 2) -> str:
    if pd.isna(value):
        return "--"
    return f"{float(value):.{digits}f}"


def fmt_pct(value: object, digits: int = 1) -> str:
    if pd.isna(value):
        return "--"
    return f"{float(value):.{digits}f}\\%"


def write_table2_tex(summary: pd.DataFrame, path: Path) -> None:
    lines = [
        r"\begin{table}[h!]",
        r"\caption{Strip-and-regraft validation for single-site phosphosites. RMSD values (\AA) are permutation-invariant phosphate RMSD against experimental coordinates. Top-1 is the highest-scored pose after minimisation. Top-3 best reports the minimum RMSD across three independently minimised ranked output poses.}",
        r"\label{tab:validation}",
        r"\centering",
        r"\small",
        r"\begin{tabularx}{\textwidth}{l c c c c c c c c}",
        r"\toprule",
        r"& & \multicolumn{3}{c}{\textbf{Median RMSD (\AA)}} & \multicolumn{4}{c}{\textbf{\% sites within threshold}} \\",
        r"\cmidrule(lr){3-5} \cmidrule(lr){6-9}",
        r"\textbf{Type} & \textbf{$n$} & \textbf{Top-1} & \textbf{Top-3} & \textbf{S1} & \textbf{Top-1 $\leq$1.0} & \textbf{Top-3 $\leq$1.0} & \textbf{Top-1 $\leq$1.5} & \textbf{Top-3 $\leq$1.5} \\",
        r"\midrule",
    ]
    for _, row in summary.iterrows():
        if row["type"] == "Combined":
            lines.append(r"\midrule")
        lines.append(
            f"{row['type']} & {int(row['n'])} & "
            f"{fmt_number(row['top1_postmin_median'])} & {fmt_number(row['top3_best_postmin_median'])} & "
            f"{fmt_number(row['stage1_oracle_median'])} & {fmt_pct(row['top1_le_1A_pct'])} & "
            f"{fmt_pct(row['top3_le_1A_pct'])} & {fmt_pct(row['top1_le_1_5A_pct'])} & "
            f"{fmt_pct(row['top3_le_1_5A_pct'])} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabularx}", r"\end{table}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_table3_tex(summary: pd.DataFrame, path: Path) -> None:
    lines = [
        r"\begin{table}[h!]",
        r"\caption{Strip-and-regraft validation stratified by the number of phosphosites in each structure. Top-3 best is the minimum RMSD across three independently minimised ranked output poses.}",
        r"\label{tab:multisite}",
        r"\centering",
        r"\small",
        r"\begin{tabularx}{\textwidth}{l l c c c c c}",
        r"\toprule",
        r"\textbf{Type} & \textbf{Context} & \textbf{$n$} & \textbf{Top-1 median} & \textbf{Top-3 median} & \textbf{Top-1 $\leq$1.0} & \textbf{Top-3 $\leq$1.5} \\",
        r"\midrule",
    ]
    previous = None
    for _, row in summary.iterrows():
        if previous is not None and row["type"] != previous:
            lines.append(r"\addlinespace")
        lines.append(
            f"{row['type']} & {row['context']} & {int(row['n'])} & "
            f"{fmt_number(row['top1_postmin_median'])} & {fmt_number(row['top3_best_postmin_median'])} & "
            f"{fmt_pct(row['top1_le_1A_pct'])} & {fmt_pct(row['top3_le_1_5A_pct'])} \\\\"
        )
        previous = row["type"]
    lines.extend([r"\bottomrule", r"\end{tabularx}", r"\end{table}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def plot_figure3(single: pd.DataFrame, output_dir: Path, ymax: float) -> None:
    plt.rcParams.update(RC)
    fig, axes = plt.subplots(1, 3, figsize=(17.2, 5.6), sharey=True)
    positions = np.arange(1, len(STAGES) + 1)

    for index, restype in enumerate(RESTYPES):
        ax = axes[index]
        cohort = single.loc[single["restype_3"].eq(restype)]
        data = [pd.to_numeric(cohort[column], errors="coerce").dropna().to_numpy() for _, column in STAGES]

        if data and all(len(values) >= 2 and np.ptp(values) > 0 for values in data):
            violin = ax.violinplot(
                data,
                positions=positions,
                showmeans=False,
                showmedians=False,
                showextrema=False,
                widths=0.72,
                bw_method=0.25,
            )
            for body in violin["bodies"]:
                body.set_facecolor(TYPE_COLORS[restype])
                body.set_edgecolor(TYPE_COLORS[restype])
                body.set_alpha(0.25)

        if data and all(len(values) for values in data):
            box = ax.boxplot(
                data,
                positions=positions,
                widths=0.20,
                patch_artist=True,
                showfliers=False,
                medianprops={"color": "black", "linewidth": 2.0},
                whiskerprops={"linewidth": 1.0},
                capprops={"linewidth": 1.0},
            )
            for patch in box["boxes"]:
                patch.set_facecolor(TYPE_COLORS[restype])
                patch.set_edgecolor(TYPE_COLORS[restype])
                patch.set_alpha(0.75)

        ax.set_xticks(positions)
        ax.set_xticklabels([label for label, _ in STAGES])
        ax.set_title(f"{restype} (n={len(cohort)})")
        ax.set_ylabel("Phosphate RMSD (Å)" if index == 0 else "")
        ax.axhline(1.0, color="#c62828", alpha=0.45, linestyle="--", linewidth=1.1)
        ax.axhline(1.5, color="#ef6c00", alpha=0.45, linestyle="--", linewidth=1.1)
        ax.set_ylim(0, ymax)
        ax.grid(axis="y", alpha=0.16, linestyle="--")
        ax.set_axisbelow(True)

        for pos, values in zip(positions, data):
            if not len(values):
                continue
            median = float(np.median(values))
            q75 = float(np.percentile(values, 75))
            q25 = float(np.percentile(values, 25))
            upper = min(q75 + 1.5 * (q75 - q25), float(np.max(values)), ymax - 0.28)
            label_y = min(max(upper + 0.12, median + 0.16), ymax - 0.12)
            ax.text(pos, label_y, f"{median:.2f}", ha="center", va="bottom", fontsize=9, fontweight="bold")

        ax.text(-0.10, 1.035, chr(ord("a") + index), transform=ax.transAxes, fontsize=18, fontweight="bold")

    fig.tight_layout(w_pad=1.5)
    for extension in ("pdf", "png"):
        fig.savefig(output_dir / f"fig3_validation_stages_postmin.{extension}", bbox_inches="tight", dpi=300)
    plt.close(fig)


def safe_stat(value: float) -> object:
    return None if not np.isfinite(value) else float(value)


def statistical_summary(single: pd.DataFrame) -> Dict[str, object]:
    result: Dict[str, object] = {}
    groups = [
        single.loc[single["restype_3"].eq(restype), "rank1_postmin_rmsd"].dropna().to_numpy()
        for restype in RESTYPES
    ]
    if all(len(group) for group in groups):
        test = stats.kruskal(*groups)
        result["kruskal_wallis_rank1_by_residue"] = {
            "H": safe_stat(float(test.statistic)),
            "p": safe_stat(float(test.pvalue)),
        }

    if "stage2_scoring_gap" in single.columns:
        gap = pd.to_numeric(single["stage2_scoring_gap"], errors="coerce")
        rmsd = pd.to_numeric(single["rank1_postmin_rmsd"], errors="coerce")
        valid = gap.notna() & rmsd.notna()
        if int(valid.sum()) >= 3:
            test = stats.spearmanr(gap.loc[valid], rmsd.loc[valid])
            result["spearman_scoring_gap_vs_rank1_postmin"] = {
                "n": int(valid.sum()),
                "rho": safe_stat(float(test.statistic)),
                "p": safe_stat(float(test.pvalue)),
            }

    ranks = single["top3_best_postmin_rank"].round().astype(int).value_counts().sort_index()
    result["best_postmin_rank_counts"] = {str(int(rank)): int(count) for rank, count in ranks.items()}
    result["all_three_poses_independently_measured"] = True
    return result


def write_caption(single_summary: pd.DataFrame, path: Path) -> None:
    counts = {row["type"]: int(row["n"]) for _, row in single_summary.iterrows()}
    combined = single_summary.loc[single_summary["type"].eq("Combined")].iloc[0]
    caption = (
        r"\caption{\textbf{Validation results and RMSD decomposition.} "
        r"Per-site RMSD distributions for internal-coordinate seed placement (S0), the oracle-best of 12 scanned orientations (S1), "
        r"the orientation selected by scoring before minimisation (S2), the independently minimised top-ranked pose (S3), and the best "
        rf"of three independently minimised ranked poses for TPO ($n={counts.get('TPO', 0)}$), SEP ($n={counts.get('SEP', 0)}$), "
        rf"and PTR ($n={counts.get('PTR', 0)}$) single-site phosphosites. Median values are annotated above each box. Dashed lines "
        rf"indicate 1.0~\AA\ and 1.5~\AA\ thresholds. Across all residue types, {combined['top3_le_1_5A_pct']:.1f}\% of sites "
        r"are recovered within 1.5~\AA\ by at least one of the three minimised outputs.}"
    )
    path.write_text(caption + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create corrected validation summaries, manuscript tables, and Figure 3 from independently minimised top-3 benchmark outputs."
    )
    parser.add_argument("--phosphofill-tsv", nargs="+", required=True, help="Corrected TPO/SEP/PTR benchmark TSVs or non-overlapping shards")
    parser.add_argument("--outdir", default="benchmark_postmin_analysis")
    parser.add_argument("--ymax", type=float, default=4.5, help="Figure 3 y-axis maximum in Å")
    args = parser.parse_args()

    output_dir = Path(args.outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    all_sites = load_outputs(args.phosphofill_tsv)
    single = all_sites.loc[all_sites["context_n_sites"].eq(1)].copy()
    multi = all_sites.loc[all_sites["context_n_sites"].gt(1)].copy()
    if single.empty:
        raise ValueError("No single-site rows were found; Figure 3 and Table 2 cannot be generated")

    all_sites.to_csv(output_dir / "benchmark_primary_sites.tsv", sep="\t", index=False)
    single.to_csv(output_dir / "benchmark_single_site.tsv", sep="\t", index=False)
    multi.to_csv(output_dir / "benchmark_multi_site.tsv", sep="\t", index=False)

    single_summary = make_single_summary(single)
    context_summary = make_context_summary(all_sites)
    single_summary.to_csv(output_dir / "table2_validation_values.tsv", sep="\t", index=False, float_format="%.6g")
    context_summary.to_csv(output_dir / "table3_multisite_values.tsv", sep="\t", index=False, float_format="%.6g")
    write_table2_tex(single_summary, output_dir / "table2_validation.tex")
    write_table3_tex(context_summary, output_dir / "table3_multisite.tex")
    write_caption(single_summary, output_dir / "figure3_caption.tex")
    plot_figure3(single, output_dir, args.ymax)

    payload = {
        "input_files": [str(Path(path).resolve()) for path in args.phosphofill_tsv],
        "filter": "status=OK, primary residue type per input file",
        "single_site_n": int(len(single)),
        "multi_site_n": int(len(multi)),
        "statistics": statistical_summary(single),
        "single_site_summary": single_summary.to_dict("records"),
        "context_summary": context_summary.to_dict("records"),
    }
    (output_dir / "benchmark_postmin_summary.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    print(f"Loaded {len(all_sites)} primary-type OK sites ({len(single)} single-site; {len(multi)} multi-site)")
    print(f"Wrote corrected Figure 3, Tables 2/3, and summaries to {output_dir.resolve()}")


if __name__ == "__main__":
    main()
