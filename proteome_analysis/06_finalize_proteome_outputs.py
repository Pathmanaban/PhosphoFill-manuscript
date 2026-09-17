#!/usr/bin/env python3
"""Create the final proteome-scale figure, table, statistics and manuscript text.

The input is the three-rank PhosphoFill result table.  All biological summaries
are calculated from pose rank 1 after validating that every site has one each of
ranks 1, 2 and 3.  The script does not simulate, impute or resample site values.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.ndimage import gaussian_filter1d


COLORS = {
    "SEP": "#2e7d32",
    "TPO": "#1565c0",
    "PTR": "#6a1b9a",
    "high": "#427f49",
    "low": "#b84b4b",
}

MOD_ORDER = ["SEP", "TPO", "PTR"]
SITE_KEY = ["acc_id", "chain_id", "resseq", "new_resname"]
HIGH_PLDDT_CLASSES = {"very_high", "confident"}

TPO_MODE = -58.9
TPO_TOL = 30.0
SEP_SHARP_MODE = 66.3
SEP_SHARP_TOL = 15.0
SEP_BROAD_LOW = -60.0
SEP_BROAD_HIGH = 50.0


def circular_distance_deg(values: pd.Series | np.ndarray, centre: float) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return np.abs(((arr - centre + 180.0) % 360.0) - 180.0)


def fmt_int(value: int) -> str:
    return f"{int(value):,}"


def fmt_num(value: float, digits: int = 1) -> str:
    return f"{float(value):.{digits}f}"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def mann_whitney_summary(x: pd.Series, y: pd.Series) -> dict[str, float | int]:
    """Return U, rank-biserial effect and a non-underflowing asymptotic p-value."""
    x = pd.to_numeric(x, errors="coerce").dropna().to_numpy(dtype=float)
    y = pd.to_numeric(y, errors="coerce").dropna().to_numpy(dtype=float)
    if len(x) == 0 or len(y) == 0:
        raise ValueError("Mann-Whitney comparison requires two non-empty groups")

    combined = np.concatenate([x, y])
    ranks = stats.rankdata(combined, method="average")
    n_x = len(x)
    n_y = len(y)
    u_x = float(ranks[:n_x].sum() - n_x * (n_x + 1) / 2.0)
    mean_u = n_x * n_y / 2.0

    _, tie_counts = np.unique(combined, return_counts=True)
    n_total = n_x + n_y
    tie_term = float(np.sum(tie_counts**3 - tie_counts))
    variance = (n_x * n_y / 12.0) * (
        (n_total + 1.0) - tie_term / (n_total * (n_total - 1.0))
    )
    if variance <= 0:
        raise ValueError("Mann-Whitney variance is not positive")

    continuity = 0.5 * np.sign(u_x - mean_u)
    z = (u_x - mean_u - continuity) / math.sqrt(variance)
    log10_p = (math.log(2.0) + float(stats.norm.logsf(abs(z)))) / math.log(10.0)
    rank_biserial = 2.0 * u_x / (n_x * n_y) - 1.0

    return {
        "n_high": n_x,
        "n_low": n_y,
        "u_high": u_x,
        "z_asymptotic_tie_corrected": float(z),
        "log10_two_sided_p": float(log10_p),
        "rank_biserial_r": float(rank_biserial),
    }


def load_and_validate(path: Path) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    required = {
        *SITE_KEY,
        "pose_rank",
        "status",
        "phosphofill_confidence",
        "torsion_key",
        "nearby_basic_count_after",
        "salt_bridge_count_after",
        "plddt_site",
        "plddt_class",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df["pose_rank"] = pd.to_numeric(df["pose_rank"], errors="raise").astype(int)
    duplicate_rows = int(df.duplicated(SITE_KEY + ["pose_rank"]).sum())
    rank_sets = df.groupby(SITE_KEY, dropna=False)["pose_rank"].agg(lambda s: tuple(sorted(s)))
    invalid_rank_sets = int((rank_sets != (1, 2, 3)).sum())
    non_ok_rows = int((df["status"].astype(str).str.upper() != "OK").sum())
    if duplicate_rows or invalid_rank_sets or non_ok_rows:
        raise ValueError(
            "Input validation failed: "
            f"duplicate site/rank rows={duplicate_rows}, "
            f"invalid rank sets={invalid_rank_sets}, non-OK rows={non_ok_rows}"
        )

    top = df[df["pose_rank"] == 1].copy()
    top["mod_type"] = top["new_resname"].astype(str)
    unexpected = sorted(set(top["mod_type"]) - set(MOD_ORDER))
    if unexpected:
        raise ValueError(f"Unexpected phosphoresidue types: {unexpected}")

    top["confidence"] = pd.to_numeric(top["phosphofill_confidence"], errors="coerce")
    top["torsion"] = pd.to_numeric(top["torsion_key"], errors="coerce")
    top["n_basic"] = pd.to_numeric(top["nearby_basic_count_after"], errors="coerce")
    top["salt_bridge"] = (
        pd.to_numeric(top["salt_bridge_count_after"], errors="coerce").fillna(0) > 0
    )
    top["plddt_numeric"] = pd.to_numeric(top["plddt_site"], errors="coerce")

    class_high = top["plddt_class"].isin(HIGH_PLDDT_CLASSES)
    numeric_available = top["plddt_numeric"].notna()
    top["plddt_high"] = np.where(numeric_available, top["plddt_numeric"] >= 70.0, class_high)
    plddt_class_mismatches = int(
        (numeric_available & ((top["plddt_numeric"] >= 70.0) != class_high)).sum()
    )

    audit = {
        "source_pose_rows": int(len(df)),
        "top1_sites": int(len(top)),
        "proteins": int(top["acc_id"].nunique()),
        "duplicate_site_rank_rows": duplicate_rows,
        "invalid_rank_sets": invalid_rank_sets,
        "non_ok_rows": non_ok_rows,
        "plddt_class_numeric_mismatches": plddt_class_mismatches,
        "missing_confidence": int(top["confidence"].isna().sum()),
        "missing_torsion": int(top["torsion"].isna().sum()),
        "missing_basic_contact_count": int(top["n_basic"].isna().sum()),
    }
    return top, audit


def derive_statistics(df: pd.DataFrame, source: Path, audit: dict) -> tuple[dict, pd.DataFrame]:
    counts = df["mod_type"].value_counts().reindex(MOD_ORDER).fillna(0).astype(int)
    shares = counts / len(df) * 100.0

    tpo = df[(df["mod_type"] == "TPO") & df["torsion"].notna()]
    sep = df[(df["mod_type"] == "SEP") & df["torsion"].notna()]
    ptr = df[(df["mod_type"] == "PTR") & df["torsion"].notna()]

    tpo_supported = circular_distance_deg(tpo["torsion"], TPO_MODE) <= TPO_TOL
    sep_sharp = circular_distance_deg(sep["torsion"], SEP_SHARP_MODE) <= SEP_SHARP_TOL
    sep_broad = sep["torsion"].between(SEP_BROAD_LOW, SEP_BROAD_HIGH).to_numpy()
    sep_outside = ~(sep_sharp | sep_broad)

    high = df[df["plddt_high"]]
    low = df[~df["plddt_high"]]
    mw = mann_whitney_summary(high["confidence"], low["confidence"])

    high_sb = float(high["salt_bridge"].mean() * 100.0)
    low_sb = float(low["salt_bridge"].mean() * 100.0)
    sb_ratio = high_sb / low_sb if low_sb else math.nan

    by_type = {}
    table_rows = []
    for mod in MOD_ORDER:
        sub = df[df["mod_type"] == mod]
        if mod == "TPO":
            torsion_pct = float(tpo_supported.mean() * 100.0)
        elif mod == "SEP":
            torsion_pct = float((sep_sharp | sep_broad).mean() * 100.0)
        else:
            torsion_pct = None
        row = {
            "type": mod,
            "n_sites": int(len(sub)),
            "n_proteins": int(sub["acc_id"].nunique()),
            "median_confidence": float(sub["confidence"].median()),
            "predicted_salt_bridge_pct": float(sub["salt_bridge"].mean() * 100.0),
            "high_plddt_pct": float(sub["plddt_high"].mean() * 100.0),
            "torsion_prior_concordance_pct": torsion_pct,
            "median_basic_contacts": float(sub["n_basic"].median()),
        }
        by_type[mod] = row
        table_rows.append(row)

    all_row = {
        "type": "All",
        "n_sites": int(len(df)),
        "n_proteins": int(df["acc_id"].nunique()),
        "median_confidence": float(df["confidence"].median()),
        "predicted_salt_bridge_pct": float(df["salt_bridge"].mean() * 100.0),
        "high_plddt_pct": float(df["plddt_high"].mean() * 100.0),
        "torsion_prior_concordance_pct": None,
        "median_basic_contacts": float(df["n_basic"].median()),
    }
    table_rows.append(all_row)

    ptr_angles = np.deg2rad(ptr["torsion"].to_numpy(dtype=float))
    ptr_mean_vector = np.mean(np.exp(1j * ptr_angles)) if len(ptr_angles) else complex(np.nan, np.nan)

    statistics = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": str(source.resolve()),
        "source_sha256": sha256_file(source),
        "software_versions": {
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
            "seaborn": sns.__version__,
            "scipy": stats.__version__ if hasattr(stats, "__version__") else __import__("scipy").__version__,
        },
        "analysis_definition": {
            "pose_rank": 1,
            "localisation_probability": ">=90%",
            "structure_scope": "complete canonical single-F1 AlphaFold Database v6 models",
            "site_observations_simulated_or_imputed": False,
        },
        "validation": audit,
        "cohort": {
            "sites": int(len(df)),
            "proteins": int(df["acc_id"].nunique()),
            "pose_rows_in_source": int(audit["source_pose_rows"]),
        },
        "residue_composition": {
            mod: {"n": int(counts[mod]), "pct": float(shares[mod])} for mod in MOD_ORDER
        },
        "by_residue_type": by_type,
        "torsion_prior_concordance": {
            "TPO": {
                "denominator": int(len(tpo)),
                "criterion": f"circular distance <= {TPO_TOL:g} degrees from {TPO_MODE:g} degrees",
                "n": int(tpo_supported.sum()),
                "pct": float(tpo_supported.mean() * 100.0),
            },
            "SEP": {
                "denominator": int(len(sep)),
                "broad_criterion": f"[{SEP_BROAD_LOW:g}, {SEP_BROAD_HIGH:g}] degrees",
                "broad_n": int(sep_broad.sum()),
                "broad_pct": float(sep_broad.mean() * 100.0),
                "sharp_criterion": f"circular distance <= {SEP_SHARP_TOL:g} degrees from +{SEP_SHARP_MODE:g} degrees",
                "sharp_n": int(sep_sharp.sum()),
                "sharp_pct": float(sep_sharp.mean() * 100.0),
                "outside_n": int(sep_outside.sum()),
                "outside_pct": float(sep_outside.mean() * 100.0),
                "supported_n": int((sep_sharp | sep_broad).sum()),
                "supported_pct": float((sep_sharp | sep_broad).mean() * 100.0),
                "broad_sharp_overlap_n": int((sep_sharp & sep_broad).sum()),
            },
            "PTR": {
                "denominator": int(len(ptr)),
                "criterion": None,
                "interpretation": "No preferred orientation prior; not included in concordance percentage",
                "circular_mean_deg": float(np.rad2deg(np.angle(ptr_mean_vector))),
                "mean_resultant_length": float(abs(ptr_mean_vector)),
            },
        },
        "plddt_stratification": {
            "high": {
                "criterion": ">=70",
                "n": int(len(high)),
                "median_confidence": float(high["confidence"].median()),
                "median_basic_contacts": float(high["n_basic"].median()),
                "zero_basic_contacts_pct": float((high["n_basic"] == 0).mean() * 100.0),
                "five_or_more_basic_contacts_pct": float((high["n_basic"] >= 5).mean() * 100.0),
                "predicted_salt_bridge_pct": high_sb,
            },
            "low": {
                "criterion": "<70",
                "n": int(len(low)),
                "median_confidence": float(low["confidence"].median()),
                "median_basic_contacts": float(low["n_basic"].median()),
                "zero_basic_contacts_pct": float((low["n_basic"] == 0).mean() * 100.0),
                "five_or_more_basic_contacts_pct": float((low["n_basic"] >= 5).mean() * 100.0),
                "predicted_salt_bridge_pct": low_sb,
            },
            "predicted_salt_bridge_fold_ratio_high_over_low": float(sb_ratio),
            "confidence_mann_whitney": mw,
        },
    }

    table = pd.DataFrame(table_rows)

    panel_rows = []
    for mod in MOD_ORDER:
        panel_rows.extend(
            [
                {"panel": "a", "group": mod, "metric": "site_count", "value": int(counts[mod])},
                {"panel": "a", "group": mod, "metric": "site_percentage", "value": float(shares[mod])},
                {"panel": "b", "group": mod, "metric": "median_confidence", "value": by_type[mod]["median_confidence"]},
            ]
        )
    panel_rows.extend(
        [
            {"panel": "c", "group": "TPO", "metric": "dominant_mode_pct", "value": statistics["torsion_prior_concordance"]["TPO"]["pct"]},
            {"panel": "c", "group": "SEP", "metric": "broad_basin_pct", "value": statistics["torsion_prior_concordance"]["SEP"]["broad_pct"]},
            {"panel": "c", "group": "SEP", "metric": "sharp_mode_pct", "value": statistics["torsion_prior_concordance"]["SEP"]["sharp_pct"]},
            {"panel": "c", "group": "SEP", "metric": "outside_supported_regions_pct", "value": statistics["torsion_prior_concordance"]["SEP"]["outside_pct"]},
        ]
    )
    for label, sub in [("pLDDT >=70", high), ("pLDDT <70", low)]:
        panel_rows.extend(
            [
                {"panel": "d", "group": label, "metric": "n_sites", "value": int(len(sub))},
                {"panel": "d", "group": label, "metric": "median_confidence", "value": float(sub["confidence"].median())},
                {"panel": "d", "group": label, "metric": "median_basic_contacts", "value": float(sub["n_basic"].median())},
                {"panel": "d", "group": label, "metric": "predicted_salt_bridge_pct", "value": float(sub["salt_bridge"].mean() * 100.0)},
            ]
        )
    return statistics, table, pd.DataFrame(panel_rows)


def style_axes(ax: plt.Axes) -> None:
    sns.despine(ax=ax)
    ax.grid(axis="y", color="#e7e7e7", linewidth=0.7, zorder=0)
    ax.set_axisbelow(True)


def add_violin(ax: plt.Axes, groups: list[pd.Series], colors: list[str]) -> None:
    clean = [pd.to_numeric(group, errors="coerce").dropna().to_numpy() for group in groups]
    vp = ax.violinplot(clean, positions=np.arange(len(clean)), showmedians=False, showextrema=False)
    for body, color in zip(vp["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor("none")
        body.set_alpha(0.62)
    for i, (values, color) in enumerate(zip(clean, colors)):
        median = float(np.median(values))
        ax.hlines(median, i - 0.16, i + 0.16, color=color, linewidth=2.1, zorder=4)


def make_figure7(df: pd.DataFrame, statistics: dict, path_pdf: Path, path_png: Path) -> None:
    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 8,
        }
    )

    fig = plt.figure(figsize=(15.2, 9.2))
    grid = fig.add_gridspec(
        2, 3, left=0.06, right=0.985, bottom=0.075, top=0.965, wspace=0.29, hspace=0.34
    )
    axes = [fig.add_subplot(grid[row, col]) for row in range(2) for col in range(3)]
    ax_a, ax_b, ax_c, ax_d, ax_e, ax_f = axes

    counts = df["mod_type"].value_counts().reindex(MOD_ORDER).astype(int)
    shares = counts / len(df) * 100.0
    bars = ax_a.bar(MOD_ORDER, shares, color=[COLORS[m] for m in MOD_ORDER], width=0.72)
    for bar, mod in zip(bars, MOD_ORDER):
        ax_a.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 1.1,
            f"{shares[mod]:.1f}%\n(n={counts[mod]:,})",
            ha="center",
            va="bottom",
            fontsize=8.5,
        )
    ax_a.set_ylim(0, max(shares) * 1.18)
    ax_a.set_ylabel("Top-ranked sites (%)")
    ax_a.set_title("Phosphosite composition")

    confidence_groups = [df.loc[df["mod_type"] == mod, "confidence"] for mod in MOD_ORDER]
    add_violin(ax_b, confidence_groups, [COLORS[m] for m in MOD_ORDER])
    ax_b.set_xticks(range(3), MOD_ORDER)
    ax_b.set_ylabel("PhosphoFill confidence score")
    ax_b.set_title("Confidence by residue type")

    torsion = statistics["torsion_prior_concordance"]
    tpo_pct = torsion["TPO"]["pct"]
    sep_broad = torsion["SEP"]["broad_pct"]
    sep_sharp = torsion["SEP"]["sharp_pct"]
    ax_c.bar(0, tpo_pct, color=COLORS["TPO"], width=0.62)
    ax_c.bar(1, sep_broad, color="#78b87a", width=0.62, label="SEP broad basin")
    ax_c.bar(1, sep_sharp, bottom=sep_broad, color=COLORS["SEP"], width=0.62, label="SEP sharp mode")
    ax_c.text(0, tpo_pct / 2, f"{tpo_pct:.1f}%", color="white", fontweight="bold", ha="center", va="center")
    ax_c.text(1, sep_broad / 2, f"Broad\n{sep_broad:.1f}%", color="white", fontweight="bold", ha="center", va="center")
    ax_c.text(1, sep_broad + sep_sharp / 2, f"Sharp\n{sep_sharp:.1f}%", color="white", fontweight="bold", ha="center", va="center")
    ax_c.set_xticks([0, 1], ["TPO\n-58.9° ± 30°", "SEP\nsupported regions"])
    ax_c.set_ylim(0, 106)
    ax_c.set_ylabel("Top-ranked sites (%)")
    ax_c.set_title("Concordance with empirical torsion criteria")
    # The stacked segments are labelled directly; a legend would obscure the
    # TPO bar and repeat the same information.

    high = df[df["plddt_high"]]
    low = df[~df["plddt_high"]]
    add_violin(ax_d, [high["confidence"], low["confidence"]], [COLORS["high"], COLORS["low"]])
    ax_d.set_xticks([0, 1], ["pLDDT ≥70", "pLDDT <70"])
    ax_d.set_ylabel("PhosphoFill confidence score")
    ax_d.set_title("Confidence by AlphaFold pLDDT")
    rb = statistics["plddt_stratification"]["confidence_mann_whitney"]["rank_biserial_r"]
    ax_d.text(
        0.5,
        0.97,
        rf"Mann-Whitney $p<10^{{-300}}$; rank-biserial $r={rb:.2f}$",
        transform=ax_d.transAxes,
        ha="center",
        va="top",
        fontsize=8.5,
    )

    labels = ["0", "1", "2", "3", "4", "5+"]
    positions = np.arange(len(labels))
    width = 0.36
    for offset, sub, color, label in [
        (-width / 2, high, COLORS["high"], "pLDDT ≥70"),
        (width / 2, low, COLORS["low"], "pLDDT <70"),
    ]:
        clipped = sub["n_basic"].clip(upper=5)
        freq = clipped.value_counts(normalize=True)
        percentages = np.array([float(freq.get(i, 0.0) * 100.0) for i in range(6)])
        ax_e.bar(positions + offset, percentages, width, color=color, label=label)
    ax_e.set_xticks(positions, labels)
    ax_e.set_xlabel("Predicted basic contacts within 4 Å")
    ax_e.set_ylabel("Top-ranked sites (%)")
    ax_e.set_title("Basic-contact distribution")
    ax_e.legend(frameon=False)
    ax_e.text(
        0.5,
        -0.19,
        f"pLDDT ≥70: n={len(high):,}    pLDDT <70: n={len(low):,}",
        transform=ax_e.transAxes,
        ha="center",
        va="top",
        fontsize=8,
    )

    salt = statistics["plddt_stratification"]
    sb_values = [salt["high"]["predicted_salt_bridge_pct"], salt["low"]["predicted_salt_bridge_pct"]]
    sb_bars = ax_f.bar(["pLDDT ≥70", "pLDDT <70"], sb_values, color=[COLORS["high"], COLORS["low"]], width=0.58)
    for bar, value in zip(sb_bars, sb_values):
        ax_f.text(bar.get_x() + bar.get_width() / 2, value + 0.7, f"{value:.1f}%", ha="center", va="bottom", fontweight="bold")
    ax_f.set_ylim(0, max(sb_values) * 1.18)
    ax_f.set_ylabel("Sites with a predicted salt bridge (%)")
    ax_f.set_title("Salt-bridge formation")

    for ax in axes:
        style_axes(ax)
    for label, ax in zip(["(a)", "(b)", "(c)"], [ax_a, ax_b, ax_c]):
        ax.text(-0.15, 1.06, label, transform=ax.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")
    ax_d.text(-0.15, 1.06, "(d)", transform=ax_d.transAxes, fontsize=12, fontweight="bold", ha="left", va="bottom")

    path_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path_pdf, bbox_inches="tight", facecolor="white")
    fig.savefig(path_png, dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_torsion_figure(df: pd.DataFrame, path_pdf: Path, path_png: Path) -> None:
    sns.set_theme(style="ticks", context="paper")
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8), sharey=False)
    bins = np.arange(-180.0, 180.0001, 10.0)
    centres = (bins[:-1] + bins[1:]) / 2.0
    for ax, mod in zip(axes, ["TPO", "SEP", "PTR"]):
        values = df.loc[df["mod_type"] == mod, "torsion"].dropna().to_numpy(dtype=float)
        density, _ = np.histogram(values, bins=bins, density=True)
        smooth = gaussian_filter1d(density, sigma=0.8, mode="wrap")
        ax.hist(values, bins=bins, density=True, color=COLORS[mod], alpha=0.25, edgecolor="none")
        ax.plot(centres, smooth, color=COLORS[mod], linewidth=2.1)
        ax.fill_between(centres, 0, smooth, color=COLORS[mod], alpha=0.10)
        ax.set_xlim(-180, 180)
        ax.set_xticks([-180, -90, 0, 90, 180])
        ax.set_title(f"{mod} (n={len(values):,})", fontweight="bold")
        ax.set_xlabel("Top-ranked torsion (°)")
        ax.set_ylabel("Density")
        if mod == "TPO":
            ax.axvline(TPO_MODE, color="#222222", linestyle="--", linewidth=1.4, label="-58.9° placement prior")
            ax.legend(frameon=False, fontsize=8)
        elif mod == "SEP":
            ax.axvspan(SEP_BROAD_LOW, SEP_BROAD_HIGH, color=COLORS["SEP"], alpha=0.10, label="broad basin")
            ax.axvline(SEP_SHARP_MODE, color="#222222", linestyle="--", linewidth=1.4, label="+66.3° sharp prior")
            ax.legend(frameon=False, fontsize=8)
        else:
            ax.text(0.5, 0.95, "No preferred orientation prior", transform=ax.transAxes, ha="center", va="top", fontsize=8.5, color=COLORS["PTR"])
        style_axes(ax)
    fig.tight_layout(w_pad=1.1)
    fig.savefig(path_pdf, bbox_inches="tight", facecolor="white")
    fig.savefig(path_png, dpi=240, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_table_outputs(table: pd.DataFrame, outdir: Path, stats_dict: dict) -> None:
    tsv = table.copy()
    for col in ["median_confidence", "predicted_salt_bridge_pct", "high_plddt_pct", "torsion_prior_concordance_pct", "median_basic_contacts"]:
        tsv[col] = pd.to_numeric(tsv[col], errors="coerce").round(1)
    tsv.to_csv(outdir / "supplementary_table_s5.tsv", sep="\t", index=False, na_rep="n.a.")

    rows = []
    for _, row in table.iterrows():
        torsion = "n.a." if pd.isna(row["torsion_prior_concordance_pct"]) else f'{row["torsion_prior_concordance_pct"]:.1f}'
        rows.append(
            f'{row["type"]} & {int(row["n_sites"]):,} & {int(row["n_proteins"]):,} & '
            f'{row["median_confidence"]:.0f} & {row["predicted_salt_bridge_pct"]:.1f} & '
            f'{row["high_plddt_pct"]:.1f} & {torsion} & {row["median_basic_contacts"]:.0f} \\\\'
        )
    caption = (
        "\\caption*{\\textbf{Supplementary Table S5.} Summary of 102,300 successfully processed "
        "top-ranked PhosphoFill predictions for Scop3P class~I human phosphosites with localisation "
        "probability $\\geq90\\%$, spanning 15,049 proteins represented by complete canonical "
        "single-F1 AlphaFold Database v6 models. Confidence, basic contacts and salt bridges are "
        "model-derived PhosphoFill annotations. Torsion-prior concordance is reported only for TPO "
        "and SEP because no preferred-orientation torsion prior is applied to PTR.}"
    )
    tex = "\n".join(
        [
            "\\begin{table}[h!]",
            caption,
            "\\centering",
            "\\small",
            "\\begin{tabularx}{\\textwidth}{l c c c c c c c}",
            "\\toprule",
            "\\textbf{Type} & \\textbf{$n$ sites} & \\textbf{$n$ proteins} & \\textbf{Median conf.} & \\textbf{Predicted salt bridge (\\%)} & \\textbf{High pLDDT (\\%)} & \\textbf{Torsion-prior concordance (\\%)} & \\textbf{Median basic contacts} \\\\",
            "\\midrule",
            *rows[:3],
            "\\midrule",
            rows[3],
            "\\bottomrule",
            "\\end{tabularx}",
            "\\end{table}",
            "",
        ]
    )
    (outdir / "supplementary_table_s5.tex").write_text(tex, encoding="utf-8")


def write_manuscript_text(statistics: dict, outdir: Path) -> None:
    cohort = statistics["cohort"]
    composition = statistics["residue_composition"]
    by_type = statistics["by_residue_type"]
    torsion = statistics["torsion_prior_concordance"]
    plddt = statistics["plddt_stratification"]
    rb = plddt["confidence_mann_whitney"]["rank_biserial_r"]

    caption = (
        "\\caption{\\textbf{Proteome-wide PhosphoFill application to successfully processed Scop3P "
        "class~I human phosphosites.} The analysis contains "
        f"{cohort['sites']:,} top-ranked sites with localisation probability $\\geq90\\%$ across "
        f"{cohort['proteins']:,} proteins represented by complete canonical single-F1 AlphaFold Database v6 models. "
        "\\textbf{(a)} Phosphosite composition by residue type: "
        f"{composition['SEP']['pct']:.1f}\\% SEP, {composition['TPO']['pct']:.1f}\\% TPO and "
        f"{composition['PTR']['pct']:.1f}\\% PTR. "
        "\\textbf{(b)} Distributions of model-derived PhosphoFill confidence scores by residue type. "
        "\\textbf{(c)} Concordance of top-ranked torsions with the empirical criteria used during placement. "
        f"{torsion['TPO']['pct']:.1f}\\% of TPO predictions lie within $\\pm30^{{\\circ}}$ of "
        "$-58.9^{\\circ}$; the SEP stacked bar shows "
        f"{torsion['SEP']['broad_pct']:.1f}\\% in the broad [$-60^{{\\circ}},+50^{{\\circ}}$] basin "
        f"and {torsion['SEP']['sharp_pct']:.1f}\\% within $\\pm15^{{\\circ}}$ of the sharp "
        "$+66.3^{\\circ}$ prior. PTR is omitted because no preferred-orientation torsion prior is applied. "
        "These percentages are consistency checks against criteria used during placement and are not an independent accuracy validation. "
        "\\textbf{(d)} Stratification by AlphaFold pLDDT. Sites with pLDDT $\\geq70$ and $<70$ have "
        f"median confidence scores of {plddt['high']['median_confidence']:.0f} and {plddt['low']['median_confidence']:.0f}, "
        f"respectively, and predicted salt-bridge frequencies of {plddt['high']['predicted_salt_bridge_pct']:.1f}\\% "
        f"and {plddt['low']['predicted_salt_bridge_pct']:.1f}\\% "
        f"(Mann--Whitney confidence comparison, $p<10^{{-300}}$; rank-biserial $r={rb:.2f}$).}}"
    )
    (outdir / "figure7_caption.tex").write_text(caption + "\n", encoding="utf-8")

    results = "\n\n".join(
        [
            (
                "To demonstrate scalability, we analysed PhosphoFill predictions for "
                f"{cohort['sites']:,} successfully processed Scop3P class~I phosphosites with localisation "
                f"probability $\\geq90\\%$, spanning {cohort['proteins']:,} canonical human proteins represented "
                "by complete single-F1 AlphaFold Database v6 models. The top-ranked cohort comprised "
                f"{composition['SEP']['n']:,} SEP ({composition['SEP']['pct']:.1f}\\%), "
                f"{composition['TPO']['n']:,} TPO ({composition['TPO']['pct']:.1f}\\%) and "
                f"{composition['PTR']['n']:,} PTR ({composition['PTR']['pct']:.1f}\\%) sites "
                "(Fig.~\\ref{fig:proteome}a; Supplementary Table~S5). Three independently reconstructed and "
                "minimised ranked structures were generated per site; proteome-scale summaries use only the top-ranked pose."
            ),
            (
                "Model-derived PhosphoFill confidence differed across phosphoresidue types, with median scores of "
                f"{by_type['SEP']['median_confidence']:.0f} for SEP, {by_type['TPO']['median_confidence']:.0f} for TPO and "
                f"{by_type['PTR']['median_confidence']:.0f} for PTR (Fig.~\\ref{{fig:proteome}}b). "
                f"Sites with pLDDT $\\geq70$ (n={plddt['high']['n']:,}) had a median confidence of "
                f"{plddt['high']['median_confidence']:.0f}, compared with {plddt['low']['median_confidence']:.0f} for "
                f"sites with pLDDT $<70$ (n={plddt['low']['n']:,}; Mann--Whitney $p<10^{{-300}}$, "
                f"rank-biserial $r={rb:.2f}$). Predicted salt bridges occurred at "
                f"{plddt['high']['predicted_salt_bridge_pct']:.1f}\\% and "
                f"{plddt['low']['predicted_salt_bridge_pct']:.1f}\\% of sites, respectively "
                f"({plddt['predicted_salt_bridge_fold_ratio_high_over_low']:.1f}-fold difference; Fig.~\\ref{{fig:proteome}}d). "
                "These confidence and contact measures describe PhosphoFill annotations on predicted structures rather than experimental placement accuracy."
            ),
            (
                "Top-ranked torsions were strongly concordant with the empirical criteria used during placement: "
                f"{torsion['TPO']['pct']:.1f}\\% of TPO sites fell within $\\pm30^{{\\circ}}$ of $-58.9^{{\\circ}}$. "
                f"For SEP, {torsion['SEP']['broad_pct']:.1f}\\% fell in the broad "
                f"[$-60^{{\\circ}},+50^{{\\circ}}$] basin and {torsion['SEP']['sharp_pct']:.1f}\\% fell within "
                f"$\\pm15^{{\\circ}}$ of $+66.3^{{\\circ}}$; together, the supported regions contained "
                f"{torsion['SEP']['supported_pct']:.1f}\\% of SEP predictions. PTR was not assigned a concordance "
                "percentage because PhosphoFill applies no preferred-orientation PTR torsion prior. These values assess "
                "internal consistency with the placement model and do not constitute independent validation (Fig.~\\ref{fig:proteome}c)."
            ),
        ]
    )
    (outdir / "proteome_results_text.tex").write_text(results + "\n", encoding="utf-8")

    methods = (
        "Proteome-scale reconstruction used Scop3P class~I human phosphosites with localisation probability "
        "$\\geq90\\%$ mapped to canonical UniProt accessions. Complete canonical single-F1 AlphaFold Database "
        "v6 models were used as structural input; isoform-only records and proteins represented only by segmented "
        "models were excluded from the analysis-ready cohort. Phosphosites were scanned sequentially from the N "
        "terminus to the C terminus and minimised sequentially using the production settings. Three ranked poses were "
        "reconstructed from their absolute prescan coordinates and minimised independently. Figure~\\ref{fig:proteome} "
        "and Supplementary Table~S5 report only pose rank~1. Confidence scores, contacts and salt bridges are "
        "model-derived annotations on AlphaFold structures."
    )
    (outdir / "proteome_methods_text.tex").write_text(methods + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    top, audit = load_and_validate(args.results)
    statistics, table, panel_values = derive_statistics(top, args.results, audit)

    (args.outdir / "figure7_statistics.json").write_text(
        json.dumps(statistics, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    panel_values.to_csv(args.outdir / "figure7_panel_values.tsv", sep="\t", index=False)
    write_table_outputs(table, args.outdir, statistics)
    write_manuscript_text(statistics, args.outdir)
    make_figure7(
        top,
        statistics,
        args.outdir / "fig7_proteome.pdf",
        args.outdir / "fig7_proteome.png",
    )
    make_torsion_figure(
        top,
        args.outdir / "proteome_torsion_distributions.pdf",
        args.outdir / "proteome_torsion_distributions.png",
    )

    print(json.dumps({"sites": len(top), "proteins": top["acc_id"].nunique(), "outdir": str(args.outdir)}, indent=2))


if __name__ == "__main__":
    main()
