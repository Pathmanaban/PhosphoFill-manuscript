#!/usr/bin/env python3
"""Rebuild proteome context Supplementary Figures S1 and S2.

Figure S1 joins top-ranked PhosphoFill sites to the DSSP/RSA annotation table.
Figure S2 uses all successfully processed top-ranked sites and does not require
the DSSP/RSA join.  All plotted values are also written to long-form TSV files.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RESIDUES = ["SEP", "TPO", "PTR"]
RESIDUE_COLORS = {"SEP": "#2e7d32", "TPO": "#1565c0", "PTR": "#6a1b9a"}
SS_ORDER = ["Helix", "Sheet", "Loop"]
SS_COLORS = {"Helix": "#E76F51", "Sheet": "#F4A261", "Loop": "#74A9CF"}
RSA_LABELS = ["0-0.1", "0.1-0.25", "0.25-0.5", "0.5-0.75", "0.75-1.0"]
DENSITY_LABELS = ["1", "2-3", "4-7", "8-15", ">15"]
CONTACT_LABELS = ["0", "1", "2", "3+"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True, help="Complete three-rank proteome result TSV")
    parser.add_argument("--dssp", type=Path, required=True, help="DSSP/RSA annotation table (Structure_AF.txt)")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory")
    return parser.parse_args()


def configure_plotting() -> None:
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.dpi": 300,
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.14, 1.08, f"({label})", transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def polish(ax: plt.Axes, grid: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid:
        ax.grid(axis="y", color="#D9D9D9", linewidth=0.6, alpha=0.7)
        ax.set_axisbelow(True)


def load_top1(results_path: Path) -> tuple[pd.DataFrame, dict]:
    required = [
        "acc_id",
        "chain_id",
        "resseq",
        "pose_rank",
        "new_resname",
        "status",
        "plddt_site",
        "salt_bridge_count_after",
        "nearby_basic_count_after",
        "phosphofill_confidence",
    ]
    df = pd.read_csv(results_path, sep="\t", usecols=required, low_memory=False)
    df["pose_rank"] = pd.to_numeric(df["pose_rank"], errors="coerce")
    df["resseq"] = pd.to_numeric(df["resseq"], errors="coerce")
    df["site_key"] = (
        df["acc_id"].astype(str)
        + "|"
        + df["chain_id"].astype(str)
        + "|"
        + df["resseq"].astype("Int64").astype(str)
    )
    rank_sets = df.groupby("site_key", sort=False)["pose_rank"].agg(lambda x: tuple(sorted(set(x.dropna().astype(int)))))
    invalid_rank_sets = int((rank_sets != (1, 2, 3)).sum())
    if invalid_rank_sets:
        raise ValueError(f"Found {invalid_rank_sets} sites without exactly ranks 1, 2 and 3")

    top = df[(df["pose_rank"] == 1) & (df["status"] == "OK")].copy()
    if top["site_key"].duplicated().any():
        raise ValueError("Duplicate rank-1 site records found")
    for col in ["plddt_site", "salt_bridge_count_after", "nearby_basic_count_after", "phosphofill_confidence"]:
        top[col] = pd.to_numeric(top[col], errors="coerce")
    top["salt_bridge"] = top["salt_bridge_count_after"].fillna(0).gt(0)
    top["plddt_group"] = np.where(top["plddt_site"] >= 70, "pLDDT >=70", "pLDDT <70")
    manifest = {
        "input_rows": int(len(df)),
        "input_unique_sites": int(df["site_key"].nunique()),
        "rank1_successful_sites": int(len(top)),
        "rank1_proteins": int(top["acc_id"].nunique()),
        "rank_set_validation_failures": invalid_rank_sets,
    }
    return top, manifest


def load_dssp(dssp_path: Path) -> pd.DataFrame:
    dssp = pd.read_csv(dssp_path, sep="\t", low_memory=False)
    needed = {"ACC_ID", "UniProt_pos", "SS8", "RSA", "Position_Label"}
    missing = sorted(needed - set(dssp.columns))
    if missing:
        raise ValueError(f"DSSP table is missing columns: {', '.join(missing)}")
    # Position_Label can contain multiple comma-separated annotations, for
    # example "Phosphorylation,Mutation".  These are still phosphosites and
    # must not be discarded by an exact-string filter.
    dssp = dssp[
        dssp["Position_Label"].astype(str).str.contains("Phosphorylation", regex=False, na=False)
    ].copy()
    dssp["UniProt_pos"] = pd.to_numeric(dssp["UniProt_pos"], errors="coerce")
    dssp["RSA"] = pd.to_numeric(dssp["RSA"], errors="coerce")
    key_cols = ["ACC_ID", "UniProt_pos"]
    conflicting = (
        dssp.groupby(key_cols, dropna=False)[["SS8", "RSA"]]
        .nunique(dropna=False)
        .gt(1)
        .any(axis=1)
        .sum()
    )
    if conflicting:
        raise ValueError(f"DSSP table contains {int(conflicting)} conflicting duplicate site keys")
    return dssp.drop_duplicates(key_cols)[["ACC_ID", "UniProt_pos", "SS8", "RSA"]]


def annotate_context(top: pd.DataFrame, dssp: pd.DataFrame) -> pd.DataFrame:
    context = top.merge(
        dssp,
        left_on=["acc_id", "resseq"],
        right_on=["ACC_ID", "UniProt_pos"],
        how="left",
        validate="many_to_one",
    )
    ss_map = {
        "H": "Helix",
        "G": "Helix",
        "I": "Helix",
        "E": "Sheet",
        "B": "Sheet",
        "T": "Loop",
        "S": "Loop",
        "C": "Loop",
        # DSSP 4 uses P for polyproline-II; retain it in the broad loop/coil
        # class used by the original three-state manuscript analysis.
        "P": "Loop",
        "-": "Loop",
    }
    context["SS3"] = context["SS8"].astype(str).map(ss_map)
    context["RSA"] = pd.to_numeric(context["RSA"], errors="coerce")
    context["rsa_bin"] = pd.cut(
        context["RSA"],
        bins=[-np.inf, 0.1, 0.25, 0.5, 0.75, np.inf],
        labels=RSA_LABELS,
        right=False,
    )
    return context


def record(rows: list[dict], panel: str, metric: str, value: float, unit: str, **groups: object) -> None:
    item = {"figure_panel": panel, "metric": metric, **groups, "value": value, "unit": unit}
    rows.append(item)


def build_s1(context: pd.DataFrame, outdir: Path, manifest: dict) -> list[dict]:
    usable = context[context["SS3"].notna() & context["RSA"].notna()].copy()
    if usable.empty:
        raise ValueError("No top-ranked sites could be joined to DSSP/RSA annotations")
    rows: list[dict] = []
    fig, axes = plt.subplots(2, 3, figsize=(13.2, 8.0), constrained_layout=True)

    # a, left: secondary-structure composition.
    ax = axes[0, 0]
    counts = usable.groupby(["new_resname", "SS3"]).size().unstack(fill_value=0).reindex(index=RESIDUES, columns=SS_ORDER, fill_value=0)
    bottom = np.zeros(len(RESIDUES))
    for ss in SS_ORDER:
        vals = counts[ss].to_numpy()
        ax.bar(RESIDUES, vals, bottom=bottom, color=SS_COLORS[ss], width=0.68, label=ss)
        for residue, value in zip(RESIDUES, vals):
            record(rows, "a-left", "site_count", int(value), "sites", residue_type=residue, category=ss, denominator_n=int(counts.loc[residue].sum()))
        bottom += vals
    ax.set_ylabel("Annotated phosphosites")
    ax.set_title("Secondary-structure context")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    panel_label(ax, "a")
    polish(ax)

    # a, middle: salt bridges by secondary structure.
    ax = axes[0, 1]
    x = np.arange(len(SS_ORDER))
    width = 0.23
    for idx, residue in enumerate(RESIDUES):
        values = []
        for ss in SS_ORDER:
            sub = usable[(usable["new_resname"] == residue) & (usable["SS3"] == ss)]
            value = 100 * sub["salt_bridge"].mean() if len(sub) else np.nan
            values.append(value)
            record(rows, "a-middle", "salt_bridge_frequency", value, "percent", residue_type=residue, category=ss, denominator_n=int(len(sub)))
        ax.bar(x + (idx - 1) * width, values, width, color=RESIDUE_COLORS[residue], label=residue)
    ax.set_xticks(x, SS_ORDER)
    ax.set_ylabel("Sites with predicted salt bridge (%)")
    ax.set_title("Salt bridges by secondary structure")
    ax.legend(frameon=False)
    polish(ax)

    # a, right: confidence by secondary structure.
    ax = axes[0, 2]
    for idx, residue in enumerate(RESIDUES):
        values = []
        for ss in SS_ORDER:
            sub = usable[(usable["new_resname"] == residue) & (usable["SS3"] == ss)]["phosphofill_confidence"].dropna()
            value = float(sub.median()) if len(sub) else np.nan
            values.append(value)
            record(rows, "a-right", "median_confidence", value, "score", residue_type=residue, category=ss, denominator_n=int(len(sub)))
        ax.plot(x, values, marker="o", linewidth=2, markersize=5, color=RESIDUE_COLORS[residue], label=residue)
    ax.set_xticks(x, SS_ORDER)
    ax.set_ylim(50, 92)
    ax.set_ylabel("Median PhosphoFill confidence")
    ax.set_title("Confidence by secondary structure")
    ax.legend(frameon=False)
    polish(ax)

    # b, left: RSA distributions.
    ax = axes[1, 0]
    violin_data = [usable.loc[usable["new_resname"] == residue, "RSA"].dropna().clip(0, 1).to_numpy() for residue in RESIDUES]
    parts = ax.violinplot(violin_data, positions=np.arange(1, 4), widths=0.72, showmeans=False, showmedians=True, showextrema=False)
    for body, residue in zip(parts["bodies"], RESIDUES):
        body.set_facecolor(RESIDUE_COLORS[residue])
        body.set_edgecolor("white")
        body.set_alpha(0.82)
        values = usable.loc[usable["new_resname"] == residue, "RSA"].dropna()
        record(rows, "b-left", "median_rsa", float(values.median()), "RSA", residue_type=residue, category="all", denominator_n=int(len(values)))
    parts["cmedians"].set_color("#202020")
    parts["cmedians"].set_linewidth(2)
    ax.set_xticks([1, 2, 3], RESIDUES)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Relative solvent accessibility")
    ax.set_title("RSA distribution")
    panel_label(ax, "b")
    polish(ax)

    # b, middle: salt bridges by RSA bin.
    ax = axes[1, 1]
    x = np.arange(len(RSA_LABELS))
    for residue in RESIDUES:
        values = []
        for rsa_bin in RSA_LABELS:
            sub = usable[(usable["new_resname"] == residue) & (usable["rsa_bin"].astype(str) == rsa_bin)]
            value = 100 * sub["salt_bridge"].mean() if len(sub) else np.nan
            values.append(value)
            record(rows, "b-middle", "salt_bridge_frequency", value, "percent", residue_type=residue, category=rsa_bin, denominator_n=int(len(sub)))
        ax.plot(x, values, marker="o", linewidth=2, markersize=5, color=RESIDUE_COLORS[residue], label=residue)
    ax.set_xticks(x, RSA_LABELS, rotation=25, ha="right")
    ax.set_ylabel("Sites with predicted salt bridge (%)")
    ax.set_title("Salt bridges by RSA")
    ax.legend(frameon=False)
    polish(ax)

    # b, right: confidence by RSA bin (median, not mean).
    ax = axes[1, 2]
    for residue in RESIDUES:
        values = []
        for rsa_bin in RSA_LABELS:
            sub = usable[(usable["new_resname"] == residue) & (usable["rsa_bin"].astype(str) == rsa_bin)]["phosphofill_confidence"].dropna()
            value = float(sub.median()) if len(sub) else np.nan
            values.append(value)
            record(rows, "b-right", "median_confidence", value, "score", residue_type=residue, category=rsa_bin, denominator_n=int(len(sub)))
        ax.plot(x, values, marker="o", linewidth=2, markersize=5, color=RESIDUE_COLORS[residue], label=residue)
    ax.set_xticks(x, RSA_LABELS, rotation=25, ha="right")
    ax.set_ylim(55, 92)
    ax.set_ylabel("Median PhosphoFill confidence")
    ax.set_title("Confidence by RSA")
    ax.legend(frameon=False)
    polish(ax)

    for ext in ["pdf", "png"]:
        fig.savefig(outdir / f"fig_s1_proteome_ss_rsa.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)

    manifest["figure_s1"] = {
        "annotated_sites": int(len(usable)),
        "annotated_proteins": int(usable["acc_id"].nunique()),
        "coverage_percent": float(100 * len(usable) / len(context)),
        "overall_secondary_structure_percent": {
            key.lower(): float(100 * usable["SS3"].eq(key).mean()) for key in SS_ORDER
        },
    }
    return rows


def density_class(value: int) -> str:
    if value == 1:
        return "1"
    if value <= 3:
        return "2-3"
    if value <= 7:
        return "4-7"
    if value <= 15:
        return "8-15"
    return ">15"


def build_s2(top: pd.DataFrame, outdir: Path, manifest: dict) -> list[dict]:
    rows: list[dict] = []
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.0), constrained_layout=True)

    # a, left: pLDDT and salt bridges.
    ax = axes[0, 0]
    x = np.arange(len(RESIDUES))
    width = 0.34
    plddt_groups = [("pLDDT >=70", "#3A9D5D"), ("pLDDT <70", "#D95F59")]
    for idx, (group, color) in enumerate(plddt_groups):
        values = []
        for residue in RESIDUES:
            sub = top[(top["new_resname"] == residue) & (top["plddt_group"] == group)]
            value = 100 * sub["salt_bridge"].mean() if len(sub) else np.nan
            values.append(value)
            record(rows, "a-left", "salt_bridge_frequency", value, "percent", residue_type=residue, category=group, denominator_n=int(len(sub)))
        bars = ax.bar(x + (idx - 0.5) * width, values, width, color=color, label=group)
        ax.bar_label(bars, labels=[f"{v:.1f}%" for v in values], padding=2, fontsize=8)
    ax.set_xticks(x, RESIDUES)
    ax.set_ylim(0, 55)
    ax.set_ylabel("Sites with predicted salt bridge (%)")
    ax.set_title("Structural confidence and salt bridges")
    ax.legend(frameon=False)
    panel_label(ax, "a")
    polish(ax)

    # a, right: number of nearby basic residues.
    ax = axes[0, 1]
    x = np.arange(len(CONTACT_LABELS))
    width = 0.23
    contact_class = top["nearby_basic_count_after"].fillna(0).clip(lower=0).astype(int).map(lambda v: str(v) if v < 3 else "3+")
    for idx, residue in enumerate(RESIDUES):
        residue_mask = top["new_resname"].eq(residue)
        denominator = int(residue_mask.sum())
        values = []
        for category in CONTACT_LABELS:
            n = int((residue_mask & contact_class.eq(category)).sum())
            value = 100 * n / denominator if denominator else np.nan
            values.append(value)
            record(rows, "a-right", "nearby_basic_contact_distribution", value, "percent", residue_type=residue, category=category, numerator_n=n, denominator_n=denominator)
        ax.bar(x + (idx - 1) * width, values, width, color=RESIDUE_COLORS[residue], label=residue)
    ax.set_xticks(x, CONTACT_LABELS)
    ax.set_xlabel("Nearby Arg/Lys residues within 4 A")
    ax.set_ylabel("Phosphosites (%)")
    ax.set_title("Local basic-residue contacts")
    ax.legend(frameon=False)
    polish(ax)

    # b: analysed-site density per protein.
    site_counts = top.groupby("acc_id")["resseq"].nunique()
    plot_df = top.copy()
    plot_df["analysed_sites_per_protein"] = plot_df["acc_id"].map(site_counts)
    plot_df["density_class"] = plot_df["analysed_sites_per_protein"].astype(int).map(density_class)
    x = np.arange(len(DENSITY_LABELS))

    ax = axes[1, 0]
    medians = []
    ns = []
    for category in DENSITY_LABELS:
        sub = plot_df[plot_df["density_class"] == category]
        medians.append(float(sub["phosphofill_confidence"].median()))
        ns.append(int(len(sub)))
        record(rows, "b-left", "median_confidence", medians[-1], "score", residue_type="all", category=category, denominator_n=ns[-1])
    ax.plot(x, medians, color="#264653", marker="o", linewidth=2.3, markersize=6)
    for xx, yy, n in zip(x, medians, ns):
        ax.text(xx, yy + 1.1, f"n={n:,}", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(x, DENSITY_LABELS)
    ax.set_xlabel("Analysed phosphosites per protein")
    ax.set_ylabel("Median PhosphoFill confidence")
    ax.set_ylim(65, 84)
    ax.set_title("Confidence by phosphosite density")
    panel_label(ax, "b")
    polish(ax)

    ax = axes[1, 1]
    frequencies = []
    for category, n in zip(DENSITY_LABELS, ns):
        sub = plot_df[plot_df["density_class"] == category]
        value = float(100 * sub["salt_bridge"].mean())
        frequencies.append(value)
        record(rows, "b-right", "salt_bridge_frequency", value, "percent", residue_type="all", category=category, denominator_n=n)
    bars = ax.bar(x, frequencies, color="#5B8E7D", width=0.68)
    for bar, value, n in zip(bars, frequencies, ns):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 0.55, f"{value:.1f}%\nn={n:,}", ha="center", va="bottom", fontsize=7.5)
    ax.set_xticks(x, DENSITY_LABELS)
    ax.set_xlabel("Analysed phosphosites per protein")
    ax.set_ylabel("Sites with predicted salt bridge (%)")
    ax.set_ylim(0, max(frequencies) + 6)
    ax.set_title("Salt bridges by phosphosite density")
    polish(ax)

    for ext in ["pdf", "png"]:
        fig.savefig(outdir / f"fig_s2_proteome_saltbridge_multisite.{ext}", bbox_inches="tight", dpi=300)
    plt.close(fig)

    manifest["figure_s2"] = {
        "analysed_sites": int(len(top)),
        "analysed_proteins": int(top["acc_id"].nunique()),
        "density_site_counts": dict(zip(DENSITY_LABELS, ns)),
    }
    return rows


def write_captions(outdir: Path, manifest: dict) -> None:
    annotated = manifest["figure_s1"]["annotated_sites"]
    coverage = manifest["figure_s1"]["coverage_percent"]
    if abs(coverage - 100.0) < 1e-9:
        coverage_text = f"All {annotated:,} successfully processed sites had matching"
    else:
        coverage_text = f"Of 102,300 successfully processed sites, {annotated:,} ({coverage:.1f}\\%) could be matched to"
    s1 = rf"""\caption*{{\textbf{{Supplementary Figure S1.}} Secondary-structure and solvent-accessibility context of top-ranked proteome-wide PhosphoFill predictions. {coverage_text} DSSP-derived secondary-structure and relative-solvent-accessibility (RSA) annotations for the input AlphaFold models. \textbf{{(a, left)}} Secondary-structure composition by phosphoresidue type; approximately 75\% of annotated sites occur in loop regions, whereas PTR is relatively enriched in sheets (18\%). \textbf{{(a, middle)}} Fraction of sites with at least one predicted phosphate-mediated salt bridge by secondary-structure class and phosphoresidue type. Frequencies are highest in sheet contexts (PTR 48\%, SEP 36\%, and TPO 29\%). \textbf{{(a, right)}} Median PhosphoFill confidence score by secondary-structure class and phosphoresidue type. \textbf{{(b, left)}} RSA distributions by phosphoresidue type; horizontal bars mark medians. \textbf{{(b, middle)}} Fraction of sites with at least one predicted salt bridge by RSA interval. Predicted salt bridges are most frequent at low RSA ($<0.25$) and decline with increasing accessibility. \textbf{{(b, right)}} Median PhosphoFill confidence score by RSA interval. Confidence scores and contacts are model-derived PhosphoFill annotations rather than experimentally measured quantities.}}
"""
    s2 = r"""\caption*{\textbf{Supplementary Figure S2.} Association of structural confidence, local interaction environment, and analysed phosphosite density with top-ranked proteome-wide PhosphoFill predictions. All panels use 102,300 successfully processed sites across 15,049 canonical human proteins. \textbf{(a, left)} Fraction of sites with at least one predicted phosphate-mediated salt bridge, stratified by phosphoresidue type and AlphaFold pLDDT. Sites with pLDDT $\geq70$ show higher predicted salt-bridge frequencies than sites with pLDDT $<70$ for all three residue types. \textbf{(a, right)} Distribution of nearby Arg/Lys residues within 4~\AA\ of the reconstructed phosphate group. Most sites have no nearby basic residue, while a subset has one or more contacts. \textbf{(b, left)} Median PhosphoFill confidence score as a function of the number of successfully analysed phosphosites per protein. \textbf{(b, right)} Fraction of sites with at least one predicted salt bridge across the same density categories; the descriptive frequency decreases from 20.8\% in proteins with one analysed site to 13.0\% in proteins with more than 15 analysed sites. Values above points or bars give the number of sites in each category. These associations may reflect differences in residue composition, structural confidence, and accessibility and should not be interpreted as causal effects of phosphorylation density. Confidence scores and contacts are model-derived PhosphoFill annotations rather than experimentally measured quantities.}
"""
    (outdir / "fig_s1_caption.tex").write_text(s1, encoding="utf-8")
    (outdir / "fig_s2_caption.tex").write_text(s2, encoding="utf-8")


def main() -> None:
    args = parse_args()
    for path in [args.results, args.dssp]:
        if not path.exists():
            raise FileNotFoundError(path)
    args.outdir.mkdir(parents=True, exist_ok=True)
    configure_plotting()

    top, manifest = load_top1(args.results)
    dssp = load_dssp(args.dssp)
    context = annotate_context(top, dssp)
    s1_rows = build_s1(context, args.outdir, manifest)
    s2_rows = build_s2(top, args.outdir, manifest)

    pd.DataFrame(s1_rows).to_csv(args.outdir / "fig_s1_source_values.tsv", sep="\t", index=False)
    pd.DataFrame(s2_rows).to_csv(args.outdir / "fig_s2_source_values.tsv", sep="\t", index=False)
    write_captions(args.outdir, manifest)
    (args.outdir / "supplementary_context_statistics.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    print(args.outdir / "fig_s1_proteome_ss_rsa.pdf")
    print(args.outdir / "fig_s2_proteome_saltbridge_multisite.pdf")


if __name__ == "__main__":
    main()
