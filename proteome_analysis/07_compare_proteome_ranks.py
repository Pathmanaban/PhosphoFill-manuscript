#!/usr/bin/env python3
"""Compare independently minimized PhosphoFill pose ranks proteome-wide."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.ndimage import gaussian_filter1d


SITE_KEY = ["acc_id", "chain_id", "resseq", "new_resname"]
MOD_ORDER = ["TPO", "SEP", "PTR"]
RANK_ORDER = [1, 2, 3]
RANK_COLORS = {1: "#173f5f", 2: "#e07a24", 3: "#7a7a7a"}
MOD_COLORS = {"TPO": "#1565c0", "SEP": "#2e7d32", "PTR": "#6a1b9a"}
PLDDT_COLORS = {"high": "#427f49", "low": "#b84b4b"}

TPO_MODE = -58.9
SEP_SHARP_MODE = 66.3
SEP_BROAD = (-60.0, 50.0)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def circular_distance_deg(a: np.ndarray, b: np.ndarray | float) -> np.ndarray:
    return np.abs(((np.asarray(a, dtype=float) - np.asarray(b, dtype=float) + 180.0) % 360.0) - 180.0)


def load_and_validate(path: Path) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    required = {
        *SITE_KEY,
        "pose_rank",
        "status",
        "phosphofill_confidence",
        "torsion_key",
        "salt_bridge_count_after",
        "nearby_basic_count_after",
        "plddt_site",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df["pose_rank"] = pd.to_numeric(df["pose_rank"], errors="raise").astype(int)
    duplicate_rows = int(df.duplicated(SITE_KEY + ["pose_rank"]).sum())
    rank_sets = df.groupby(SITE_KEY, dropna=False)["pose_rank"].agg(lambda values: tuple(sorted(values)))
    invalid_rank_sets = int((rank_sets != (1, 2, 3)).sum())
    non_ok_rows = int((df["status"].astype(str).str.upper() != "OK").sum())
    if duplicate_rows or invalid_rank_sets or non_ok_rows:
        raise ValueError(
            "Rank comparison requires a complete successful paired cohort: "
            f"duplicates={duplicate_rows}, invalid rank sets={invalid_rank_sets}, non-OK rows={non_ok_rows}"
        )

    df["mod_type"] = df["new_resname"].astype(str)
    df["confidence"] = pd.to_numeric(df["phosphofill_confidence"], errors="coerce")
    df["torsion"] = pd.to_numeric(df["torsion_key"], errors="coerce")
    df["salt_bridge"] = pd.to_numeric(df["salt_bridge_count_after"], errors="coerce").fillna(0) > 0
    df["n_basic"] = pd.to_numeric(df["nearby_basic_count_after"], errors="coerce")
    df["plddt"] = pd.to_numeric(df["plddt_site"], errors="coerce")
    df["plddt_group"] = np.where(df["plddt"] >= 70.0, "high", "low")

    missing_metrics = {
        metric: int(df[metric].isna().sum()) for metric in ["confidence", "torsion", "n_basic", "plddt"]
    }
    if any(missing_metrics.values()):
        raise ValueError(f"Missing rank-comparison metrics: {missing_metrics}")

    audit = {
        "pose_rows": int(len(df)),
        "sites": int(rank_sets.shape[0]),
        "proteins": int(df["acc_id"].nunique()),
        "rows_per_rank": {str(rank): int((df["pose_rank"] == rank).sum()) for rank in RANK_ORDER},
        "duplicate_site_rank_rows": duplicate_rows,
        "invalid_rank_sets": invalid_rank_sets,
        "non_ok_rows": non_ok_rows,
        "missing_metrics": missing_metrics,
    }
    return df, audit


def paired_rank_table(df: pd.DataFrame, value: str) -> pd.DataFrame:
    return df.pivot(index=SITE_KEY, columns="pose_rank", values=value).sort_index()


def calculate_statistics(df: pd.DataFrame, source: Path, audit: dict) -> tuple[dict, pd.DataFrame]:
    summary_rows: list[dict] = []

    for mod in MOD_ORDER:
        for rank in RANK_ORDER:
            sub = df[(df["mod_type"] == mod) & (df["pose_rank"] == rank)]
            summary_rows.append(
                {
                    "section": "rank_by_residue",
                    "residue_type": mod,
                    "pose_rank": rank,
                    "plddt_group": "all",
                    "metric": "confidence",
                    "n_sites": int(len(sub)),
                    "median": float(sub["confidence"].median()),
                    "q1": float(sub["confidence"].quantile(0.25)),
                    "q3": float(sub["confidence"].quantile(0.75)),
                    "percentage": np.nan,
                }
            )
            summary_rows.append(
                {
                    "section": "rank_by_residue",
                    "residue_type": mod,
                    "pose_rank": rank,
                    "plddt_group": "all",
                    "metric": "predicted_salt_bridge",
                    "n_sites": int(len(sub)),
                    "median": np.nan,
                    "q1": np.nan,
                    "q3": np.nan,
                    "percentage": float(sub["salt_bridge"].mean() * 100.0),
                }
            )

    for rank in RANK_ORDER:
        for group in ["high", "low"]:
            sub = df[(df["pose_rank"] == rank) & (df["plddt_group"] == group)]
            summary_rows.append(
                {
                    "section": "rank_by_plddt",
                    "residue_type": "ALL",
                    "pose_rank": rank,
                    "plddt_group": group,
                    "metric": "predicted_salt_bridge",
                    "n_sites": int(len(sub)),
                    "median": np.nan,
                    "q1": np.nan,
                    "q3": np.nan,
                    "percentage": float(sub["salt_bridge"].mean() * 100.0),
                }
            )

    torsion_wide = paired_rank_table(df, "torsion")
    torsion_mod = torsion_wide.reset_index()["new_resname"]
    separation_stats = {}
    for mod in MOD_ORDER:
        mask = (torsion_mod == mod).to_numpy()
        sub = torsion_wide.iloc[np.flatnonzero(mask)]
        separation_stats[mod] = {}
        for rank in [2, 3]:
            delta = circular_distance_deg(sub[rank].to_numpy(), sub[1].to_numpy())
            separation_stats[mod][str(rank)] = {
                "n": int(len(delta)),
                "median_deg": float(np.median(delta)),
                "q1_deg": float(np.quantile(delta, 0.25)),
                "q3_deg": float(np.quantile(delta, 0.75)),
                "pct_at_least_15_deg": float(np.mean(delta >= 15.0) * 100.0),
                "pct_at_least_30_deg": float(np.mean(delta >= 30.0) * 100.0),
                "pct_at_least_60_deg": float(np.mean(delta >= 60.0) * 100.0),
            }
            summary_rows.append(
                {
                    "section": "paired_torsion_separation",
                    "residue_type": mod,
                    "pose_rank": rank,
                    "plddt_group": "all",
                    "metric": "absolute_circular_difference_from_rank1_deg",
                    "n_sites": int(len(delta)),
                    "median": float(np.median(delta)),
                    "q1": float(np.quantile(delta, 0.25)),
                    "q3": float(np.quantile(delta, 0.75)),
                    "percentage": np.nan,
                }
            )

    salt_wide = paired_rank_table(df, "salt_bridge").astype(bool)
    salt_mod = salt_wide.reset_index()["new_resname"]
    salt_rescue = {}
    for mod in [*MOD_ORDER, "ALL"]:
        mask = np.ones(len(salt_wide), dtype=bool) if mod == "ALL" else (salt_mod == mod).to_numpy()
        sub = salt_wide.iloc[np.flatnonzero(mask)]
        rank1 = sub[1]
        any_rank = sub.any(axis=1)
        alt_only = (~rank1) & (sub[2] | sub[3])
        salt_rescue[mod] = {
            "n_sites": int(len(sub)),
            "rank1_n": int(rank1.sum()),
            "rank1_pct": float(rank1.mean() * 100.0),
            "any_rank_n": int(any_rank.sum()),
            "any_rank_pct": float(any_rank.mean() * 100.0),
            "rank2_or_rank3_only_n": int(alt_only.sum()),
            "rank2_or_rank3_only_pct": float(alt_only.mean() * 100.0),
        }

    rank_distinctness = {}
    for mod in MOD_ORDER:
        mask = (torsion_mod == mod).to_numpy()
        sub = torsion_wide.iloc[np.flatnonzero(mask)]
        d12 = circular_distance_deg(sub[1], sub[2])
        d13 = circular_distance_deg(sub[1], sub[3])
        d23 = circular_distance_deg(sub[2], sub[3])
        rank_distinctness[mod] = {
            "all_three_identical_within_0.1_deg_n": int(((d12 < 0.1) & (d13 < 0.1) & (d23 < 0.1)).sum()),
            "all_three_identical_within_0.1_deg_pct": float(np.mean((d12 < 0.1) & (d13 < 0.1) & (d23 < 0.1)) * 100.0),
            "at_least_one_alternative_30_deg_from_rank1_pct": float(np.mean((d12 >= 30.0) | (d13 >= 30.0)) * 100.0),
        }

    statistics = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": str(source.resolve()),
        "source_sha256": sha256_file(source),
        "analysis_scope": "paired rank comparison on complete canonical single-F1 AlphaFold Database v6 cohort",
        "interpretation": "Alternative-pose diversity and model-derived annotations; not experimental placement accuracy",
        "validation": audit,
        "paired_torsion_separation": separation_stats,
        "predicted_salt_bridge_coverage": salt_rescue,
        "rank_distinctness": rank_distinctness,
    }
    return statistics, pd.DataFrame(summary_rows)


def style_axis(ax: plt.Axes) -> None:
    sns.despine(ax=ax)
    ax.grid(axis="y", color="#e8e8e8", linewidth=0.7)
    ax.set_axisbelow(True)


def circular_density(values: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    bins = np.arange(-180.0, 180.0001, 5.0)
    centres = (bins[:-1] + bins[1:]) / 2.0
    density, _ = np.histogram(values.to_numpy(dtype=float), bins=bins, density=True)
    return centres, gaussian_filter1d(density, sigma=1.6, mode="wrap")


def make_figure(df: pd.DataFrame, statistics: dict, output_pdf: Path, output_png: Path) -> None:
    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8,
        }
    )

    fig, axes = plt.subplots(2, 3, figsize=(15.2, 9.1))
    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.085, top=0.965, wspace=0.29, hspace=0.35)
    ax_a, ax_b, ax_c, ax_d, ax_e, ax_f = axes.ravel()

    for ax, mod, label in zip([ax_a, ax_b, ax_c], MOD_ORDER, ["(a)", "(b)", "(c)"]):
        for rank in RANK_ORDER:
            values = df.loc[(df["mod_type"] == mod) & (df["pose_rank"] == rank), "torsion"]
            centres, density = circular_density(values)
            ax.plot(centres, density, color=RANK_COLORS[rank], linewidth=2.0, label=f"Rank {rank}")
        if mod == "TPO":
            ax.axvline(TPO_MODE, color="#222222", linestyle="--", linewidth=1.2, label="-58.9° prior")
        elif mod == "SEP":
            ax.axvspan(SEP_BROAD[0], SEP_BROAD[1], color=MOD_COLORS["SEP"], alpha=0.08, label="broad basin")
            ax.axvline(SEP_SHARP_MODE, color="#222222", linestyle="--", linewidth=1.2, label="+66.3° prior")
        else:
            ax.text(
                0.02,
                0.95,
                "No preferred orientation prior",
                transform=ax.transAxes,
                ha="left",
                va="top",
                color=MOD_COLORS["PTR"],
                fontsize=8.5,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.5},
            )
        ax.set_xlim(-180, 180)
        ax.set_xticks([-180, -90, 0, 90, 180])
        ax.set_xlabel("Final torsion (°)")
        ax.set_ylabel("Circular density")
        ax.set_title(f"{mod} torsion by pose rank")
        ax.text(-0.16, 1.06, label, transform=ax.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")
        ax.legend(frameon=False, loc="upper left" if mod == "SEP" else "upper right", ncol=1)
        style_axis(ax)

    offsets = {"TPO": -0.23, "SEP": 0.0, "PTR": 0.23}
    for mod in MOD_ORDER:
        medians, q1s, q3s = [], [], []
        for rank in RANK_ORDER:
            values = df.loc[(df["mod_type"] == mod) & (df["pose_rank"] == rank), "confidence"]
            medians.append(float(values.median()))
            q1s.append(float(values.quantile(0.25)))
            q3s.append(float(values.quantile(0.75)))
        x = np.arange(3) + offsets[mod]
        medians_arr = np.asarray(medians)
        ax_d.errorbar(
            x,
            medians_arr,
            yerr=np.vstack([medians_arr - np.asarray(q1s), np.asarray(q3s) - medians_arr]),
            fmt="o-",
            color=MOD_COLORS[mod],
            capsize=3,
            linewidth=1.5,
            markersize=5,
            label=mod,
        )
    ax_d.set_xticks(range(3), ["Rank 1", "Rank 2", "Rank 3"])
    ax_d.set_ylabel("PhosphoFill confidence score")
    ax_d.set_title("Median confidence and interquartile range")
    ax_d.legend(frameon=False, ncol=3, loc="lower center")
    ax_d.text(-0.16, 1.06, "(d)", transform=ax_d.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")
    style_axis(ax_d)

    x = np.arange(3)
    width = 0.34
    for offset, group, color, label in [
        (-width / 2, "high", PLDDT_COLORS["high"], "pLDDT ≥70"),
        (width / 2, "low", PLDDT_COLORS["low"], "pLDDT <70"),
    ]:
        values = []
        for rank in RANK_ORDER:
            sub = df[(df["pose_rank"] == rank) & (df["plddt_group"] == group)]
            values.append(float(sub["salt_bridge"].mean() * 100.0))
        bars = ax_e.bar(x + offset, values, width, color=color, label=label)
        for bar, value in zip(bars, values):
            ax_e.text(bar.get_x() + bar.get_width() / 2, value + 0.45, f"{value:.1f}%", ha="center", va="bottom", fontsize=7.5)
    ax_e.set_xticks(x, ["Rank 1", "Rank 2", "Rank 3"])
    ax_e.set_ylabel("Sites with a predicted salt bridge (%)")
    ax_e.set_title("Salt-bridge frequency by rank and pLDDT")
    ax_e.legend(frameon=False)
    ax_e.text(-0.16, 1.06, "(e)", transform=ax_e.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")
    style_axis(ax_e)

    torsion_wide = paired_rank_table(df, "torsion")
    torsion_mod = torsion_wide.reset_index()["new_resname"]
    positions, datasets, colors = [], [], []
    for i, mod in enumerate(MOD_ORDER):
        mask = (torsion_mod == mod).to_numpy()
        sub = torsion_wide.iloc[np.flatnonzero(mask)]
        for rank, offset in [(2, -0.16), (3, 0.16)]:
            positions.append(i + offset)
            datasets.append(circular_distance_deg(sub[rank].to_numpy(), sub[1].to_numpy()))
            colors.append(RANK_COLORS[rank])
    vp = ax_f.violinplot(datasets, positions=positions, widths=0.28, showmedians=False, showextrema=False)
    for body, color in zip(vp["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("none")
        body.set_alpha(0.58)
    for position, values, color in zip(positions, datasets, colors):
        median = float(np.median(values))
        ax_f.hlines(median, position - 0.08, position + 0.08, color=color, linewidth=2.0)
    ax_f.plot([], [], color=RANK_COLORS[2], linewidth=6, alpha=0.58, label="Rank 2 vs rank 1")
    ax_f.plot([], [], color=RANK_COLORS[3], linewidth=6, alpha=0.58, label="Rank 3 vs rank 1")
    ax_f.set_xticks(range(3), MOD_ORDER)
    ax_f.set_ylim(0, 180)
    ax_f.set_yticks([0, 30, 60, 90, 120, 150, 180])
    ax_f.axhline(30, color="#777777", linestyle=":", linewidth=1.0)
    ax_f.set_ylabel("Absolute circular torsion difference (°)")
    ax_f.set_title("Alternative-pose separation from rank 1")
    ax_f.legend(frameon=False, loc="upper left")
    ax_f.text(-0.16, 1.06, "(f)", transform=ax_f.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")
    style_axis(ax_f)

    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    fig.savefig(output_png, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_caption(statistics: dict, path: Path) -> None:
    rescue = statistics["predicted_salt_bridge_coverage"]["ALL"]
    text = (
        "\\caption{\\textbf{Proteome-scale comparison of independently minimised PhosphoFill pose ranks.} "
        "All panels use the same 102,300 successfully processed phosphosites represented by complete canonical "
        "single-F1 AlphaFold Database v6 models. \\textbf{(a--c)} Circular distributions of final torsions for "
        "ranks~1--3 for TPO, SEP and PTR. The dashed lines and shaded SEP interval indicate the empirical criteria "
        "used during placement; no preferred-orientation torsion prior is applied to PTR. \\textbf{(d)} Median "
        "model-derived PhosphoFill confidence with interquartile ranges by residue type and pose rank. "
        "\\textbf{(e)} Predicted salt-bridge frequency by pose rank and AlphaFold pLDDT class. "
        "\\textbf{(f)} Paired absolute circular torsion difference between rank~1 and each alternative rank. "
        f"Across all sites, rank~1 contains a predicted salt bridge at {rescue['rank1_pct']:.1f}\\%, whereas at least "
        f"one of the three ranks contains one at {rescue['any_rank_pct']:.1f}\\%; ranks~2 or~3 add a contact absent "
        f"from rank~1 at {rescue['rank2_or_rank3_only_pct']:.1f}\\% of sites. These analyses quantify diversity among "
        "the predicted alternatives and do not compare pose ranks with experimental phosphate coordinates.}"
    )
    path.write_text(text + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    df, audit = load_and_validate(args.results)
    statistics, summary = calculate_statistics(df, args.results, audit)

    (args.outdir / "proteome_rank_comparison_statistics.json").write_text(
        json.dumps(statistics, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    summary.round(4).to_csv(args.outdir / "proteome_rank_comparison_statistics.tsv", sep="\t", index=False, na_rep="n.a.")
    write_caption(statistics, args.outdir / "proteome_rank_comparison_caption.tex")
    make_figure(
        df,
        statistics,
        args.outdir / "proteome_rank_comparison.pdf",
        args.outdir / "proteome_rank_comparison.png",
    )
    print(json.dumps({"sites": audit["sites"], "pose_rows": audit["pose_rows"], "outdir": str(args.outdir)}, indent=2))


if __name__ == "__main__":
    main()
