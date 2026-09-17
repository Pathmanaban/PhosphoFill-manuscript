#!/usr/bin/env python3
"""Freeze Figures 5/6 values, paired statistics, Table 4, and captions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


KEYS = ["acc_id", "pdb_key", "restype_3", "chain_key", "position_key"]
METHODS = [
    ("PyTMs default", "pytms_default_rmsd"),
    ("PTM-Psi", "ptmpsi_rmsd"),
    ("PyTMs optimised", "pytms_optimised_rmsd"),
    ("PhosphoFill Stage 0 seed", "pf_stage0_seed_rmsd"),
    ("PhosphoFill top-1", "pf_top1_rmsd"),
    ("PhosphoFill top-3", "pf_top3_rmsd"),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--phosphofill-tsv", nargs=3, required=True)
    p.add_argument("--pytms-default", required=True)
    p.add_argument("--pytms-optimised", required=True)
    p.add_argument("--ptmpsi", required=True)
    p.add_argument("--audit-tsv", nargs=3, required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def norm_restype(value: object) -> str:
    return {"T": "TPO", "S": "SEP", "Y": "PTR"}.get(
        str(value).strip().upper(), str(value).strip().upper()
    )


def add_pf_keys(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["pdb_key"] = out["pdb_id"].astype(str).str.lower()
    out["chain_key"] = out["chain"].astype(str)
    out["position_key"] = pd.to_numeric(out["position"], errors="raise").astype(int)
    return out


def add_comparator_keys(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["pdb_key"] = out["reference_path"].astype(str).map(lambda value: Path(value).stem.lower())
    out["chain_key"] = out["chain"].astype(str)
    out["position_key"] = pd.to_numeric(out["position"], errors="raise").astype(int)
    return out


def load_pf(paths: list[str]) -> pd.DataFrame:
    frames = []
    for raw in paths:
        path = Path(raw)
        primary = next(rt for rt in ("TPO", "SEP", "PTR") if rt in path.stem.upper())
        df = pd.read_csv(path, sep="\t")
        df = df[(df["status"] == "OK") & (df["context_n_sites"] == 1)].copy()
        df["restype_3"] = df["restype"].map(norm_restype)
        df = df[df["restype_3"] == primary]
        frames.append(df)
    return add_pf_keys(pd.concat(frames, ignore_index=True))


def load_common(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    full_pf = load_pf(args.phosphofill_tsv)
    default = add_comparator_keys(pd.read_csv(args.pytms_default, sep="\t"))
    optimised = add_comparator_keys(pd.read_csv(args.pytms_optimised, sep="\t"))
    ptmpsi = add_comparator_keys(pd.read_csv(args.ptmpsi, sep="\t"))

    common = full_pf[KEYS].drop_duplicates()
    for frame in (default, optimised, ptmpsi):
        common = common.merge(frame.loc[frame["status"] == "OK", KEYS].drop_duplicates(), on=KEYS)

    wide = full_pf.merge(common, on=KEYS)
    wide = wide.merge(
        default.loc[default["status"] == "OK", KEYS + ["pytms_phosphate_sym_rmsd"]]
        .rename(columns={"pytms_phosphate_sym_rmsd": "pytms_default_rmsd"}),
        on=KEYS,
    )
    wide = wide.merge(
        optimised.loc[optimised["status"] == "OK", KEYS + ["pytms_phosphate_sym_rmsd"]]
        .rename(columns={"pytms_phosphate_sym_rmsd": "pytms_optimised_rmsd"}),
        on=KEYS,
    )
    wide = wide.merge(
        ptmpsi.loc[ptmpsi["status"] == "OK", KEYS + ["ptmpsi_phosphate_sym_rmsd"]]
        .rename(columns={"ptmpsi_phosphate_sym_rmsd": "ptmpsi_rmsd"}),
        on=KEYS,
    )
    wide["pf_stage0_seed_rmsd"] = pd.to_numeric(wide["stage0_seed_rmsd"], errors="coerce")
    wide["pf_top1_rmsd"] = pd.to_numeric(wide["rank1_postmin_rmsd"], errors="coerce")
    wide["pf_top3_rmsd"] = pd.to_numeric(wide["top3_best_postmin_rmsd"], errors="coerce")
    for _, column in METHODS[:3]:
        wide[column] = pd.to_numeric(wide[column], errors="coerce")
    return full_pf, wide


def rank_biserial_phosphofill_advantage(pf: np.ndarray, comparator: np.ndarray) -> float:
    difference = comparator - pf
    nonzero = difference != 0
    difference = difference[nonzero]
    ranks = stats.rankdata(np.abs(difference), method="average")
    positive = ranks[difference > 0].sum()
    negative = ranks[difference < 0].sum()
    return float((positive - negative) / (positive + negative))


def paired_statistics(data: pd.DataFrame, comparator_column: str, label: str) -> dict:
    pf = data["pf_top1_rmsd"].to_numpy(float)
    comparator = data[comparator_column].to_numpy(float)
    wilcoxon = stats.wilcoxon(pf, comparator, alternative="two-sided", method="auto")
    pf_good = pf <= 1.0
    comparator_good = comparator <= 1.0
    both = int(np.sum(pf_good & comparator_good))
    pf_only = int(np.sum(pf_good & ~comparator_good))
    comparator_only = int(np.sum(~pf_good & comparator_good))
    neither = int(np.sum(~pf_good & ~comparator_good))
    discordant = pf_only + comparator_only
    chi2 = (abs(pf_only - comparator_only) - 1) ** 2 / discordant
    return {
        "comparison": label,
        "n": int(len(data)),
        "phosphofill_lower_n": int(np.sum(pf < comparator)),
        "comparator_lower_n": int(np.sum(comparator < pf)),
        "ties_n": int(np.sum(comparator == pf)),
        "phosphofill_median_A": float(np.median(pf)),
        "comparator_median_A": float(np.median(comparator)),
        "wilcoxon_statistic": float(wilcoxon.statistic),
        "wilcoxon_p": float(wilcoxon.pvalue),
        "matched_rank_biserial_r": rank_biserial_phosphofill_advantage(pf, comparator),
        "mcnemar_threshold_A": 1.0,
        "both_recovered_n": both,
        "phosphofill_only_n": pf_only,
        "comparator_only_n": comparator_only,
        "neither_recovered_n": neither,
        "mcnemar_continuity_corrected_chi2": float(chi2),
        "mcnemar_continuity_corrected_p": float(stats.chi2.sf(chi2, 1)),
    }


def bootstrap_summary(data: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.RandomState(42)
    rows = []
    for restype in ("TPO", "SEP", "PTR", "Combined"):
        subset = data if restype == "Combined" else data[data["restype_3"] == restype]
        for method, column in METHODS:
            values = subset[column].dropna().to_numpy(float)
            boot = np.empty(10000, dtype=float)
            for index in range(10000):
                boot[index] = np.median(rng.choice(values, size=len(values), replace=True))
            low, high = np.percentile(boot, [2.5, 97.5])
            rows.append(
                {
                    "restype": restype,
                    "method": method,
                    "n": len(values),
                    "median_RMSD_A": np.median(values),
                    "bootstrap_95CI_low_A": low,
                    "bootstrap_95CI_high_A": high,
                }
            )
    return pd.DataFrame(rows)


def method_summary(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for restype in ("TPO", "SEP", "PTR", "Combined"):
        subset = data if restype == "Combined" else data[data["restype_3"] == restype]
        for method, column in METHODS:
            values = subset[column].dropna().to_numpy(float)
            rows.append(
                {
                    "restype": restype,
                    "method": method,
                    "n": len(values),
                    "median_RMSD_A": float(np.median(values)),
                    "recovery_le_1A_pct": float(np.mean(values <= 1.0) * 100),
                    "recovery_le_1_5A_pct": float(np.mean(values <= 1.5) * 100),
                }
            )
    return pd.DataFrame(rows)


def environment_matches(data: pd.DataFrame, paths: list[str]) -> pd.DataFrame:
    frames = []
    for raw in paths:
        path = Path(raw)
        restype = next(rt for rt in ("TPO", "SEP", "PTR") if rt in path.stem.upper())
        frame = pd.read_csv(path, sep="\t")
        if "status" in frame:
            frame = frame[frame["status"] == "ok"]
        frame["restype_3"] = restype
        frame["pdb_key"] = frame["pdb_id"].astype(str).str.lower()
        frame["chain_key"] = frame["chain_id"].astype(str)
        frame["position_key"] = pd.to_numeric(frame["target_position"], errors="raise").astype(int)
        frames.append(frame)
    audit = pd.concat(frames, ignore_index=True)
    return data.merge(audit[KEYS + ["interaction_class"]], on=KEYS, how="inner")


def table4_tex(summary: pd.DataFrame) -> str:
    lines = [
        r"\begin{table}[h!]",
        r"\caption{Comparison of phosphate-placement methods on the common complete-case cohort of 785 single-site phosphosites. PyTMs default uses its unoptimised placement. PyTMs optimised uses the program's built-in van der Waals optimisation with 30$^{\circ}$ angular sampling. PTM-Psi uses NERF internal-coordinate construction. PhosphoFill Stage~0 is the residue-specific internal-coordinate seed before conformational scanning and minimisation. PhosphoFill top-3 reports the minimum RMSD across three independently minimised ranked poses. RMSD is permutation-invariant phosphate RMSD against experimental coordinates.}",
        r"\label{tab:toolcomp}",
        r"\centering",
        r"\small",
        r"\begin{tabular}{llrrrr}",
        r"\toprule",
        r"\textbf{Type} & \textbf{Method} & \textbf{$n$} & \textbf{Median (\AA)} & \textbf{$\leq$1.0 (\%)} & \textbf{$\leq$1.5 (\%)} \\",
        r"\midrule",
    ]
    for restype in ("TPO", "SEP", "PTR", "Combined"):
        subset = summary[summary["restype"] == restype]
        for index, row in enumerate(subset.itertuples(index=False)):
            type_label = restype if index == 0 else ""
            method = str(row.method).replace("Naive", r'Na\"ive')
            if method.startswith("PhosphoFill"):
                method = r"\textbf{" + method + "}"
            lines.append(
                f"{type_label} & {method} & {row.n} & {row.median_RMSD_A:.2f} & "
                rf"{row.recovery_le_1A_pct:.1f}\% & {row.recovery_le_1_5A_pct:.1f}\% \\" 
            )
        if restype != "Combined":
            lines.append(r"\addlinespace")
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}"])
    return "\n".join(lines) + "\n"


def write_captions(outdir: Path, stats_payload: dict) -> None:
    rho = stats_payload["spearman_scoring_gap_common_785"]["rho"]
    rho_p = stats_payload["spearman_scoring_gap_common_785"]["p"]
    coefficient, exponent = f"{rho_p:.1e}".split("e")
    rho_p_tex = rf"{coefficient}\times10^{{{int(exponent)}}}"
    figure5 = r"""\caption{\textbf{Comparison with existing phosphorylation grafting tools.} \textbf{(a)} Median phosphate RMSD by residue type for six method configurations: PyTMs default (unoptimised placement), PTM-Psi (NERF), PyTMs optimised (built-in van der Waals optimisation with 30$^{\circ}$ angular sampling), the PhosphoFill Stage~0 residue-specific internal-coordinate seed, and PhosphoFill top-1 and top-3 outputs. All methods were evaluated on the same 785 successfully processed single-site phosphosites (328 TPO, 241 SEP, and 216 PTR). \textbf{(b)} Cumulative recovery curves showing the fraction of the common benchmark sites placed within each RMSD threshold. At 1.0~\AA, PhosphoFill top-3 recovers 84.5\% of sites, compared with 66.8\% for PhosphoFill top-1, 49.4\% for the PhosphoFill Stage~0 seed, 30.2\% for PyTMs optimised, 6.6\% for PyTMs default, and 5.5\% for PTM-Psi.}
"""
    figure6 = rf"""\caption{{\textbf{{Statistical comparison of PhosphoFill with existing phosphorylation grafting methods.}} All comparative analyses use the same 785 successfully processed single-site phosphosites. \textbf{{(a,b)}} Paired site-level comparison of PhosphoFill top-1 RMSD with PyTMs optimised and PTM-Psi, respectively. PhosphoFill produces the lower-RMSD placement for 619/785 sites (79\%) relative to PyTMs optimised and 670/785 sites (85\%) relative to PTM-Psi. \textbf{{(c,d)}} McNemar contingency analyses at a 1.0~\AA\ recovery threshold. PhosphoFill uniquely recovers 336 sites compared with 49 uniquely recovered by PyTMs optimised ($p=4.0\times10^{{-48}}$), and uniquely recovers 511 sites compared with 30 uniquely recovered by PTM-Psi ($p=1.3\times10^{{-94}}$). \textbf{{(e)}} Median RMSD and bootstrap 95\% confidence intervals by residue type and method, calculated using 10,000 bootstrap resamples. \textbf{{(f)}} Association between the PhosphoFill scoring gap and final top-1 RMSD (Spearman $\rho={rho:.2f}$, $p={rho_p_tex}$). \textbf{{(g)}} Final top-1 RMSD stratified by the structure-specific experimental phosphoresidue environment; 773 paired sites were matched to the environment audit, with salt-bridge-like and no-contact classes shown. \textbf{{(h)}} Overall median RMSD across the paired benchmark cohort for all six method configurations.}}
"""
    (outdir / "figure5_caption.tex").write_text(figure5, encoding="utf-8")
    (outdir / "figure6_caption.tex").write_text(figure6, encoding="utf-8")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    full_pf, common = load_common(args)
    summary = method_summary(common)
    bootstrap = bootstrap_summary(common)
    env = environment_matches(common, args.audit_tsv)

    kruskal = stats.kruskal(
        *[full_pf.loc[full_pf["restype_3"] == rt, "rank1_postmin_rmsd"].astype(float) for rt in ("TPO", "SEP", "PTR")]
    )
    rho_all = stats.spearmanr(full_pf["stage2_scoring_gap"].astype(float), full_pf["rank1_postmin_rmsd"].astype(float))
    rho_common = stats.spearmanr(common["stage2_scoring_gap"].astype(float), common["pf_top1_rmsd"].astype(float))
    paired = [
        paired_statistics(common, "pytms_optimised_rmsd", "PhosphoFill top-1 vs PyTMs optimised"),
        paired_statistics(common, "ptmpsi_rmsd", "PhosphoFill top-1 vs PTM-Psi"),
    ]
    payload = {
        "phosphofill_validation_n": int(len(full_pf)),
        "complete_case_n": int(len(common)),
        "complete_case_residue_counts": common["restype_3"].value_counts().sort_index().to_dict(),
        "environment_matches_n": int(len(env)),
        "kruskal_wallis_top1_by_residue_all_790": {"H": float(kruskal.statistic), "p": float(kruskal.pvalue)},
        "spearman_scoring_gap_all_790": {"rho": float(rho_all.statistic), "p": float(rho_all.pvalue)},
        "spearman_scoring_gap_common_785": {"rho": float(rho_common.statistic), "p": float(rho_common.pvalue)},
        "paired_comparisons": paired,
        "method_definition_checkpoint": {
            "PhosphoFill Stage 0 seed": "Residue-specific internal-coordinate phosphate seed before conformational scanning and minimisation; the legacy TSV column name is stage0_kabsch_rmsd",
            "PyTMs default": "PyTMs 1.2 phosphorylate(selection=site), optimize=0 default",
            "PyTMs optimised": "PyTMs 1.2 phosphorylate(selection=site, optimize=1, interval=30, remove_radius=0); the built-in loop may stop early at zero VDW strain",
            "PTM-Psi": "PTM-Psi 0.1 Protein.modify(site, 'phosphorylation'); NERF internal-coordinate construction with a 180-degree attachment torsion",
        },
    }
    summary.to_csv(outdir / "table4_tool_comparison_values.tsv", sep="\t", index=False)
    bootstrap.to_csv(outdir / "figure6_bootstrap_confidence_intervals.tsv", sep="\t", index=False)
    common.to_csv(outdir / "common_785_site_metrics.tsv", sep="\t", index=False)
    (outdir / "comparison_statistics.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (outdir / "table4_tool_comparison.tex").write_text(table4_tex(summary), encoding="utf-8")
    write_captions(outdir, payload)
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
