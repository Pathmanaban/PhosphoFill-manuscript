#!/usr/bin/env python3
"""
03_analyse_proteome.py
======================
Analyses proteome_results.tsv across six stratifications and produces
figures for the supplementary section of the PhosphoFill manuscript.

Stratifications
---------------
  1. pLDDT ≥70 vs <70        — confidence vs salt bridge formation
  2. TPO / SEP / PTR          — residue type distribution, torsion modes
  3. Secondary structure       — where phosphosites sit, PhosphoFill confidence
  4. Salt bridge yes/no        — fraction by type and pLDDT class
  5. Surface exposed vs buried — DSSP RSA or report-derived SASA class
  6. Multi-site vs single-site — crowding effect on confidence

Usage
-----
    python3 03_analyse_proteome.py \
        --results  proteome_run/proteome_results.tsv \
        --sites    proteome_run/proteome_sites.tsv \
        --outdir   proteome_run/analysis

Outputs
-------
    proteome_run/analysis/
        fig_s_plddt_stratification.png
        fig_s_residue_type_distribution.png
        fig_s_torsion_mode_frequencies.png
        fig_s_secondary_structure.png
        fig_s_salt_bridge_fraction.png
        fig_s_multisite_crowding.png
        proteome_summary_table.tsv
        torsion_mode_frequencies.tsv
        salt_bridge_fraction_by_class.tsv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

# ── Palette ────────────────────────────────────────────────────────────────
COLORS = {
    "SEP": "#2e7d32",
    "TPO": "#1565c0",
    "PTR": "#6a1b9a",
    "high_plddt": "#1b5e20",
    "low_plddt":  "#b71c1c",
    "helix":      "#e53935",
    "sheet":      "#1e88e5",
    "loop":       "#fb8c00",
    "disorder":   "#757575",
}

PLDDT_THRESH = 70.0

TPO_MODE = -58.9
TPO_MODE_TOL = 30.0
SEP_SHARP_MODE = 66.3
SEP_SHARP_TOL = 15.0
SEP_BROAD_RANGE = (-60.0, 50.0)


def circular_distance_deg(value: float, centre: float) -> float:
    return abs(((value - centre + 180.0) % 360.0) - 180.0)


def load_data(results_path: str, sites_path: str = None) -> pd.DataFrame:
    df = pd.read_csv(results_path, sep="\t")


    # Top-1 only for most analyses
    df1 = df[df["pose_rank"] == 1].copy() if "pose_rank" in df.columns else df.copy()

    # Column already has mod_type from new_resname
    df1["mod_type"] = df1["new_resname"]

    # Derived columns — use actual column names from enriched TSV
    df1["plddt_high"] = df1["plddt_class"].isin(["very_high", "confident"])
    df1["plddt_class_binary"] = df1["plddt_high"].map({True: "high", False: "low"})
    df1["salt_bridge"] = pd.to_numeric(df1["salt_bridge_count_after"], errors="coerce").fillna(0) > 0
    df1["n_basic"] = pd.to_numeric(df1["nearby_basic_count_after"], errors="coerce").fillna(0)
    df1["confidence"] = pd.to_numeric(df1["phosphofill_confidence"], errors="coerce")
    df1["torsion"] = pd.to_numeric(df1["torsion_key"], errors="coerce")
    # Use plddt_class_binary as plddt_class for stratification
    df1["plddt_class"] = df1["plddt_class_binary"]

    # Secondary structure is reported only when DSSP-derived SS3 data exist.
    # pLDDT is a confidence metric and must not be used as a structural-class
    # surrogate.
    if "ss3" in df1.columns and df1["ss3"].notna().sum() > 100:
        df1["ss_class"] = df1["ss3"]
    else:
        df1["ss_class"] = pd.NA

    # Surface exposure from RSA
    if "rsa" in df1.columns:
        df1["surface_exposed"] = pd.to_numeric(df1["rsa"], errors="coerce") >= 0.25
    elif "exposure_class_modified" in df1.columns:
        df1["surface_exposed"] = df1["exposure_class_modified"].isin(
            ["exposed", "partially_exposed"])

    # Concordance with the empirical torsion criteria used during placement.
    # PTR deliberately has no torsion prior and is therefore left unclassified.
    def torsion_class(row):
        t = row["torsion"]
        m = row["mod_type"]
        if pd.isna(t) or pd.isna(m):
            return None
        if m == "TPO":
            return "dominant mode" if circular_distance_deg(t, TPO_MODE) <= TPO_MODE_TOL else "outside"
        if m == "SEP":
            if circular_distance_deg(t, SEP_SHARP_MODE) <= SEP_SHARP_TOL:
                return "sharp mode"
            if SEP_BROAD_RANGE[0] <= t <= SEP_BROAD_RANGE[1]:
                return "broad basin"
            return "outside"
        return None
    df1["torsion_class"] = df1.apply(torsion_class, axis=1)
    df1["torsion_concordant"] = df1["torsion_class"].isin(
        ["dominant mode", "sharp mode", "broad basin"]
    )

    return df1


# ── Figure helpers ─────────────────────────────────────────────────────────

def savefig(fig, path: str, dpi: int = 200) -> None:
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {path}")


# ── Analysis 1: pLDDT stratification ──────────────────────────────────────

def fig_plddt_stratification(df: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("pLDDT stratification of PhosphoFill predictions", fontsize=13, fontweight="bold")

    for ax, (metric, label) in zip(axes, [
        ("confidence",  "PhosphoFill confidence score"),
        ("n_basic",     "Predicted basic contacts (4Å)"),
        ("salt_bridge", "Salt bridge formed (fraction)"),
    ]):
        data_high = df[df["plddt_class"] == "high"][metric].dropna()
        data_low  = df[df["plddt_class"] == "low"][metric].dropna()

        if metric == "salt_bridge":
            vals = [data_high.mean() * 100, data_low.mean() * 100]
            ax.bar(["pLDDT ≥70", "pLDDT <70"], vals,
                   color=[COLORS["high_plddt"], COLORS["low_plddt"]],
                   width=0.5, edgecolor="white")
            ax.set_ylabel("Sites with salt bridge (%)")
            for i, v in enumerate(vals):
                ax.text(i, v + 0.5, f"{v:.1f}%", ha="center", fontsize=9)
        else:
            ax.violinplot([data_high, data_low], positions=[0, 1],
                          showmedians=True, showextrema=False)
            for patch, col in zip(ax.collections[:2],
                                   [COLORS["high_plddt"], COLORS["low_plddt"]]):
                patch.set_facecolor(col)
                patch.set_alpha(0.6)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["pLDDT ≥70", "pLDDT <70"])
            ax.set_ylabel(label)

            # Mann-Whitney p-value
            if len(data_high) > 5 and len(data_low) > 5:
                _, p = stats.mannwhitneyu(data_high, data_low, alternative="two-sided")
                ax.set_title(f"p = {p:.2e}" if p < 0.05 else "n.s.", fontsize=9)

        n_h = (df["plddt_class"] == "high").sum()
        n_l = (df["plddt_class"] == "low").sum()
        ax.set_xlabel(f"pLDDT ≥70: n={n_h:,}   pLDDT <70: n={n_l:,}", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    savefig(fig, str(outdir / "fig_s_plddt_stratification.png"))


# ── Analysis 2: Residue type distribution ────────────────────────────────

def fig_residue_type(df: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Residue type distribution and torsion modes (proteome-wide)",
                 fontsize=13, fontweight="bold")

    # Panel 1: site counts
    counts = df["mod_type"].value_counts()
    colors = [COLORS.get(m, "#999") for m in counts.index]
    axes[0].bar(counts.index, counts.values, color=colors, edgecolor="white")
    axes[0].set_ylabel("Number of sites")
    axes[0].set_title("Site counts by residue type")
    for i, (k, v) in enumerate(counts.items()):
        axes[0].text(i, v + 200, f"{v:,}", ha="center", fontsize=8)
    axes[0].spines[["top", "right"]].set_visible(False)

    # Panel 2: confidence by type
    for i, mod in enumerate(["SEP", "TPO", "PTR"]):
        sub = df[df["mod_type"] == mod]["confidence"].dropna()
        axes[1].violinplot([sub], positions=[i], showmedians=True, showextrema=False)
    axes[1].collections[0].set_facecolor(COLORS["SEP"]); axes[1].collections[0].set_alpha(0.6)
    axes[1].collections[1].set_facecolor(COLORS["TPO"]); axes[1].collections[1].set_alpha(0.6)
    axes[1].collections[2].set_facecolor(COLORS["PTR"]); axes[1].collections[2].set_alpha(0.6)
    axes[1].set_xticks([0, 1, 2])
    axes[1].set_xticklabels(["SEP", "TPO", "PTR"])
    axes[1].set_ylabel("PhosphoFill confidence score")
    axes[1].set_title("Confidence by residue type")
    axes[1].spines[["top", "right"]].set_visible(False)

    # Panel 3: empirical-prior concordance. PTR is omitted because no PTR
    # torsion prior is applied. SEP's disjoint broad and sharp regions are
    # shown as a stacked bar so the two populations remain interpretable.
    tpo = df[(df["mod_type"] == "TPO") & df["torsion"].notna()]
    sep = df[(df["mod_type"] == "SEP") & df["torsion"].notna()]
    tpo_pct = 100 * (tpo["torsion_class"] == "dominant mode").mean()
    sep_broad_pct = 100 * (sep["torsion_class"] == "broad basin").mean()
    sep_sharp_pct = 100 * (sep["torsion_class"] == "sharp mode").mean()
    axes[2].bar(0, tpo_pct, color=COLORS["TPO"], edgecolor="white")
    axes[2].bar(1, sep_broad_pct, color="#78b87a", edgecolor="white", label="SEP broad basin")
    axes[2].bar(1, sep_sharp_pct, bottom=sep_broad_pct, color=COLORS["SEP"],
                edgecolor="white", label="SEP sharp mode")
    axes[2].text(0, tpo_pct + 1, f"{tpo_pct:.1f}%", ha="center", fontsize=8)
    axes[2].text(1, sep_broad_pct / 2, f"{sep_broad_pct:.1f}%", ha="center",
                 va="center", fontsize=8, color="white")
    axes[2].text(1, sep_broad_pct + sep_sharp_pct / 2, f"{sep_sharp_pct:.1f}%",
                 ha="center", va="center", fontsize=8, color="white")
    axes[2].set_xticks([0, 1])
    axes[2].set_xticklabels(["TPO dominant mode", "SEP supported regions"])
    axes[2].set_ylabel("Concordant sites (%)")
    axes[2].set_title("Concordance with empirical torsion criteria")
    axes[2].set_ylim(0, 110)
    axes[2].legend(frameon=False, fontsize=7)
    axes[2].spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    savefig(fig, str(outdir / "fig_s_residue_type_distribution.png"))


# ── Analysis 3: Secondary structure ──────────────────────────────────────

def fig_secondary_structure(df: pd.DataFrame, outdir: Path) -> None:
    # ss_class already set from DSSP ss3 in load_data — do not overwrite
    if "ss_class" not in df.columns or df["ss_class"].notna().sum() < 10:
        print("  [skip] No secondary structure data available")
        return

    ss_order  = [s for s in ["helix", "sheet", "loop"] if (df["ss_class"] == s).sum() > 5]
    mods      = ["SEP", "TPO", "PTR"]
    mod_colors = [COLORS["SEP"], COLORS["TPO"], COLORS["PTR"]]
    x = np.arange(len(ss_order))
    w = 0.25

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Secondary structure context of phosphosites — split by residue type",
                 fontsize=13, fontweight="bold")

    # Panel 1: stacked bar — site counts by SS and residue type
    bottom = np.zeros(len(ss_order))
    for mod, col in zip(mods, mod_colors):
        vals = np.array([len(df[(df["ss_class"]==s) & (df["mod_type"]==mod)])
                         for s in ss_order], dtype=float)
        axes[0].bar(ss_order, vals, bottom=bottom, color=col,
                    label=mod, edgecolor="white", alpha=0.85)
        bottom += vals
    axes[0].set_ylabel("Number of sites")
    axes[0].set_title("Site counts by SS and residue type")
    axes[0].legend(fontsize=8, frameon=False)
    axes[0].spines[["top","right"]].set_visible(False)

    # Panel 2: grouped bar — salt bridge % by SS and residue type
    for i, (mod, col) in enumerate(zip(mods, mod_colors)):
        vals = [df[(df["ss_class"]==s) & (df["mod_type"]==mod)]["salt_bridge"].mean()*100
                for s in ss_order]
        bars = axes[1].bar(x + i*w, vals, w, label=mod, color=col,
                           edgecolor="white", alpha=0.85)
        for bar, v in zip(bars, vals):
            if v > 1:
                axes[1].text(bar.get_x() + bar.get_width()/2, v + 0.3,
                             f"{v:.0f}%", ha="center", va="bottom", fontsize=6)
    axes[1].set_xticks(x + w); axes[1].set_xticklabels(ss_order)
    axes[1].set_ylabel("Sites with salt bridge (%)")
    axes[1].set_title("Salt bridge by SS and residue type")
    axes[1].legend(fontsize=8, frameon=False)
    axes[1].spines[["top","right"]].set_visible(False)

    # Panel 3: grouped bar — median confidence by SS and residue type
    for i, (mod, col) in enumerate(zip(mods, mod_colors)):
        vals = [df[(df["ss_class"]==s) & (df["mod_type"]==mod)]["confidence"].median()
                for s in ss_order]
        axes[2].bar(x + i*w, vals, w, label=mod, color=col,
                    edgecolor="white", alpha=0.85)
    axes[2].set_xticks(x + w); axes[2].set_xticklabels(ss_order)
    axes[2].set_ylabel("Median PhosphoFill confidence")
    axes[2].set_title("Confidence by SS and residue type")
    axes[2].legend(fontsize=8, frameon=False)
    axes[2].spines[["top","right"]].set_visible(False)

    fig.tight_layout()
    savefig(fig, str(outdir / "fig_s_secondary_structure.png"))


# ── Analysis 4: Salt bridge fraction ─────────────────────────────────────

def fig_salt_bridge(df: pd.DataFrame, outdir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Predicted salt bridge formation — proteome-wide",
                 fontsize=13, fontweight="bold")

    # Panel 1: by residue type + pLDDT class
    rows = []
    for mod in ["SEP", "TPO", "PTR"]:
        for pclass in ["high", "low"]:
            sub = df[(df["mod_type"] == mod) & (df["plddt_class"] == pclass)]
            rows.append({
                "mod_type": mod, "plddt_class": pclass,
                "sb_frac": sub["salt_bridge"].mean() * 100 if len(sub) > 0 else 0,
                "n": len(sub),
            })
    sb_df = pd.DataFrame(rows)

    x    = np.arange(3)
    w    = 0.35
    mods = ["SEP", "TPO", "PTR"]
    for i, pclass in enumerate(["high", "low"]):
        vals = [sb_df[(sb_df.mod_type == m) & (sb_df.plddt_class == pclass)]["sb_frac"].values[0]
                for m in mods]
        color = COLORS["high_plddt"] if pclass == "high" else COLORS["low_plddt"]
        axes[0].bar(x + i * w, vals, w, label=f"pLDDT {'≥70' if pclass=='high' else '<70'}",
                    color=color, alpha=0.8, edgecolor="white")
    axes[0].set_xticks(x + w / 2)
    axes[0].set_xticklabels(mods)
    axes[0].set_ylabel("Sites with ≥1 predicted salt bridge (%)")
    axes[0].set_title("Salt bridge fraction by residue type and pLDDT")
    axes[0].legend(fontsize=9, frameon=False)
    axes[0].spines[["top", "right"]].set_visible(False)

    # Panel 2: n_basic_contacts distribution
    for mod in ["SEP", "TPO", "PTR"]:
        sub = df[df["mod_type"] == mod]["n_basic"].dropna()
        counts = sub.value_counts().sort_index()
        axes[1].plot(counts.index, counts.values / counts.sum() * 100,
                     marker="o", label=mod, color=COLORS[mod], linewidth=1.5)
    axes[1].set_xlabel("Number of basic contacts (4Å)")
    axes[1].set_ylabel("Fraction of sites (%)")
    axes[1].set_title("Basic contact count distribution")
    axes[1].legend(fontsize=9, frameon=False)
    axes[1].spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    savefig(fig, str(outdir / "fig_s_salt_bridge_fraction.png"))


# ── Analysis 5: Multi-site crowding ──────────────────────────────────────


def fig_rsa(df: pd.DataFrame, outdir: Path) -> None:
    """RSA distribution and correlation with PhosphoFill metrics, split by residue type."""
    if "rsa" not in df.columns or df["rsa"].notna().sum() < 10:
        print("  [skip] No RSA data available")
        return

    mods = ["SEP", "TPO", "PTR"]
    mod_colors = [COLORS["SEP"], COLORS["TPO"], COLORS["PTR"]]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle("Solvent accessibility of phosphosites (proteome-wide)",
                 fontsize=13, fontweight="bold")

    # Panel 1: RSA violin by residue type
    data = [df[df["mod_type"]==m]["rsa"].dropna() for m in mods]
    vp = axes[0].violinplot(data, showmedians=True, showextrema=False)
    for patch, col in zip(vp["bodies"], mod_colors):
        patch.set_facecolor(col); patch.set_alpha(0.65)
    axes[0].set_xticks([1,2,3]); axes[0].set_xticklabels(mods)
    axes[0].set_ylabel("RSA")
    axes[0].set_title("RSA distribution by residue type")
    axes[0].axhline(0.25, color="#999", lw=1, ls="--", label="exposed (0.25)")
    axes[0].legend(fontsize=7, frameon=False)
    axes[0].spines[["top","right"]].set_visible(False)

    # Panel 2: salt bridge fraction by RSA bin and residue type
    bins   = [0, 0.1, 0.25, 0.5, 0.75, 1.01]
    blabels = ["0–0.1","0.1–0.25","0.25–0.5","0.5–0.75","0.75–1.0"]
    df = df.copy()
    df["rsa_bin"] = pd.cut(pd.to_numeric(df["rsa"], errors="coerce"),
                            bins=bins, labels=blabels, include_lowest=True)
    x = np.arange(len(blabels)); w = 0.25
    for i, (mod, col) in enumerate(zip(mods, mod_colors)):
        vals = [df[(df["rsa_bin"]==b) & (df["mod_type"]==mod)]["salt_bridge"].mean()*100
                for b in blabels]
        axes[1].bar(x + i*w, vals, w, label=mod, color=col,
                    edgecolor="white", alpha=0.85)
    axes[1].set_xticks(x + w); axes[1].set_xticklabels(blabels, rotation=25, ha="right", fontsize=8)
    axes[1].set_ylabel("Sites with salt bridge (%)")
    axes[1].set_title("Salt bridge fraction by RSA bin")
    axes[1].legend(fontsize=8, frameon=False)
    axes[1].spines[["top","right"]].set_visible(False)

    # Panel 3: binned mean confidence vs RSA per residue type
    for mod, col in zip(mods, mod_colors):
        sub = df[df["mod_type"]==mod][["rsa","confidence"]].dropna()
        if len(sub) < 10: continue
        sub = sub.copy()
        sub["rsa_bin20"] = pd.cut(pd.to_numeric(sub["rsa"], errors="coerce"), bins=20)
        binned = sub.groupby("rsa_bin20", observed=True)["confidence"].mean()
        mids = [(iv.left + iv.right)/2 for iv in binned.index]
        axes[2].plot(mids, binned.values, color=col, lw=2, label=mod)
    axes[2].set_xlabel("RSA"); axes[2].set_ylabel("Median PhosphoFill confidence")
    axes[2].set_title("Confidence vs solvent accessibility")
    axes[2].legend(fontsize=8, frameon=False)
    axes[2].spines[["top","right"]].set_visible(False)

    fig.tight_layout()
    savefig(fig, str(outdir / "fig_s_rsa.png"))

def fig_multisite(df: pd.DataFrame, outdir: Path) -> None:
    sites_per_protein = df.groupby("acc_id")["resseq"].nunique()
    df["n_sites_protein"] = df["acc_id"].map(sites_per_protein)
    df["site_class"] = pd.cut(df["n_sites_protein"],
                               bins=[0, 1, 3, 7, 15, 9999],
                               labels=["1", "2-3", "4-7", "8-15", ">15"])

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Multi-site proteins: crowding effect on PhosphoFill confidence",
                 fontsize=13, fontweight="bold")

    classes = ["1", "2-3", "4-7", "8-15", ">15"]
    medians = [df[df["site_class"] == c]["confidence"].median() for c in classes]
    counts  = [len(df[df["site_class"] == c]) for c in classes]

    axes[0].bar(classes, medians, color="#546e7a", edgecolor="white")
    axes[0].set_xlabel("Sites per protein")
    axes[0].set_ylabel("Median confidence score")
    axes[0].set_title("Confidence vs number of phosphosites per protein")
    for i, (m, n) in enumerate(zip(medians, counts)):
        axes[0].text(i, (m or 0) + 0.005, f"n={n:,}", ha="center", fontsize=7)
    axes[0].spines[["top", "right"]].set_visible(False)

    sb_frac = [df[df["site_class"] == c]["salt_bridge"].mean() * 100 for c in classes]
    axes[1].bar(classes, sb_frac, color="#37474f", edgecolor="white")
    axes[1].set_xlabel("Sites per protein")
    axes[1].set_ylabel("Sites with predicted salt bridge (%)")
    axes[1].set_title("Salt bridge formation vs protein phosphorylation density")
    axes[1].spines[["top", "right"]].set_visible(False)

    fig.tight_layout()
    savefig(fig, str(outdir / "fig_s_multisite_crowding.png"))


# ── Summary table ──────────────────────────────────────────────────────────

def write_summary(df: pd.DataFrame, outdir: Path) -> None:
    rows = []
    for mod in ["SEP", "TPO", "PTR", "ALL"]:
        sub = df if mod == "ALL" else df[df["mod_type"] == mod]
        rows.append({
            "mod_type":          mod,
            "n_sites":           len(sub),
            "n_proteins":        sub["acc_id"].nunique(),
            "median_confidence": sub["confidence"].median(),
            "sb_fraction_pct":   sub["salt_bridge"].mean() * 100,
            "plddt_high_pct":    (sub["plddt_class"] == "high").mean() * 100,
            "torsion_concordant_pct": (
                sub["torsion_concordant"].mean() * 100
                if mod in {"SEP", "TPO"} else np.nan
            ),
            "median_n_basic":    sub["n_basic"].median(),
        })
    pd.DataFrame(rows).round(2).to_csv(
        outdir / "proteome_summary_table.tsv", sep="\t", index=False)
    print(f"  Summary table: {outdir / 'proteome_summary_table.tsv'}")


# ── Main ───────────────────────────────────────────────────────────────────

def main(results_path: str, sites_path: str, outdir: str) -> None:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)

    print("Loading results...")
    df = load_data(results_path)
    print(f"  {len(df):,} site-pose rows loaded")
    print(f"  {df.acc_id.nunique():,} proteins")
    print(f"  mod_type counts: {df.mod_type.value_counts().to_dict()}")

    print("\nGenerating figures...")
    fig_plddt_stratification(df, out)
    fig_residue_type(df, out)
    fig_secondary_structure(df, out)
    fig_salt_bridge(df, out)
    fig_multisite(df, out)
    fig_rsa(df, out)
    write_summary(df, out)

    print(f"\nAll outputs written to: {out}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True)
    ap.add_argument("--sites",   required=False, default=None)
    ap.add_argument("--outdir",  default="proteome_run/analysis")
    args = ap.parse_args()
    main(args.results, args.sites, args.outdir)
