#!/usr/bin/env python3
"""Figure 5: Tool comparison. Nature-style."""
import argparse, os
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESTYPE_MAP = {"T": "TPO", "S": "SEP", "Y": "PTR"}
TYPE_COLORS = {"TPO": "#1565c0", "SEP": "#2e7d32", "PTR": "#6a1b9a"}
METHOD_COLORS = {
    "PyTMs default": "#E8703A", "PTM-Psi": "#C2185B", "PyTMs optimised": "#F5A627",
    "PF Stage 0 seed": "#9E9E9E", "PF top-1": "#1565c0", "PF top-3": "#2e7d32",
}
RC = {"axes.titlesize":13,"axes.titleweight":"bold","axes.labelsize":14,"axes.labelweight":"bold",
      "xtick.labelsize":12,"ytick.labelsize":12,"legend.fontsize":12,"legend.frameon":False,
      "figure.facecolor":"white","axes.facecolor":"white","savefig.dpi":300}

def norm_rt(x): return RESTYPE_MAP.get(str(x).strip().upper(), str(x).strip().upper())

KEYS = ["acc_id", "pdb_key", "restype_3", "chain", "position_key"]

def add_keys(df, reference_col=None):
    df = df.copy()
    if reference_col:
        df["pdb_key"] = df[reference_col].astype(str).map(lambda x: Path(x).stem.lower())
    else:
        df["pdb_key"] = df["pdb_id"].astype(str).str.lower()
    df["position_key"] = pd.to_numeric(df["position"], errors="raise").astype(int)
    return df

def load_single(paths):
    frames = []
    for p in paths:
        df = pd.read_csv(p, sep="\t")
        ok = df[df["status"]=="OK"].copy()
        ok["restype_3"] = ok["restype"].apply(norm_rt)
        primary = "TPO" if "TPO" in Path(p).stem.upper() else ("SEP" if "SEP" in Path(p).stem.upper() else "PTR")
        ok = ok[(ok["restype_3"]==primary) & (ok["context_n_sites"]==1)]
        frames.append(ok)
    return add_keys(pd.concat(frames, ignore_index=True))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phosphofill-tsv", nargs="+", required=True)
    ap.add_argument("--pytms-default-detail", required=True)
    ap.add_argument("--pytms-detail", required=True)
    ap.add_argument("--ptmpsi-detail", required=True)
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    plt.rcParams.update(RC)

    single = load_single(args.phosphofill_tsv)
    pytms_default_ok = pd.read_csv(args.pytms_default_detail, sep="\t"); pytms_default_ok = add_keys(pytms_default_ok[pytms_default_ok["status"]=="OK"], "reference_path")
    pytms_ok = pd.read_csv(args.pytms_detail, sep="\t"); pytms_ok = add_keys(pytms_ok[pytms_ok["status"]=="OK"], "reference_path")
    ptmpsi_ok = pd.read_csv(args.ptmpsi_detail, sep="\t"); ptmpsi_ok = add_keys(ptmpsi_ok[ptmpsi_ok["status"]=="OK"], "reference_path")

    # Fair comparison: retain only sites successfully evaluated by every method.
    common = single[KEYS].drop_duplicates()
    for df in [pytms_default_ok, pytms_ok, ptmpsi_ok]:
        common = common.merge(df[KEYS].drop_duplicates(), on=KEYS, how="inner")
    single = single.merge(common, on=KEYS, how="inner")
    pytms_default_ok = pytms_default_ok.merge(common, on=KEYS, how="inner")
    pytms_ok = pytms_ok.merge(common, on=KEYS, how="inner")
    ptmpsi_ok = ptmpsi_ok.merge(common, on=KEYS, how="inner")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(16, 6))

    # (a) Grouped bar
    methods_data = {}
    for res in ["TPO","SEP","PTR"]:
        sub_py = pytms_ok[pytms_ok["restype_3"]==res]
        sub_pd = pytms_default_ok[pytms_default_ok["restype_3"]==res]
        sub_pp = ptmpsi_ok[ptmpsi_ok["restype_3"]==res]
        sub_pf = single[single["restype_3"]==res]
        methods_data[res] = {
            "PyTMs default": pd.to_numeric(sub_pd["pytms_phosphate_sym_rmsd"], errors="coerce").median(),
            "PTM-Psi": pd.to_numeric(sub_pp["ptmpsi_phosphate_sym_rmsd"], errors="coerce").median(),
            "PyTMs optimised": pd.to_numeric(sub_py["pytms_phosphate_sym_rmsd"], errors="coerce").median(),
            "PF Stage 0 seed": pd.to_numeric(sub_pf["stage0_seed_rmsd"], errors="coerce").median(),
            "PF top-1": pd.to_numeric(sub_pf["stage3_post_minimization_rmsd"], errors="coerce").median(),
            "PF top-3": pd.to_numeric(sub_pf["top3_best_rmsd"], errors="coerce").median(),
        }

    method_names = list(list(methods_data.values())[0].keys())
    x = np.arange(3)
    width = 0.13
    for j, meth in enumerate(method_names):
        vals = [methods_data[res][meth] for res in ["TPO","SEP","PTR"]]
        offset = (j - len(method_names)/2 + 0.5) * width
        ax_a.bar(x+offset, vals, width, label=meth, color=METHOD_COLORS[meth], edgecolor="white")
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(["TPO","SEP","PTR"], fontsize=13, fontweight="bold")
    ax_a.set_ylabel("Median phosphate RMSD (Å)")
    ax_a.legend(fontsize=9, loc="upper left", ncol=2)
    ax_a.axhline(1.0, color="red", alpha=0.3, ls="--")
    ax_a.set_ylim(0, 3.2)
    ax_a.grid(axis="y", alpha=0.15, ls="--")
    ax_a.text(-0.08, 1.05, "a", transform=ax_a.transAxes, fontsize=20, fontweight="bold", color="black")

    # (b) Cumulative recovery
    thresholds = np.arange(0, 4.01, 0.05)
    datasets = [
        ("PyTMs default", pytms_default_ok["pytms_phosphate_sym_rmsd"], METHOD_COLORS["PyTMs default"], (0, (3, 1, 1, 1))),
        ("PTM-Psi", ptmpsi_ok["ptmpsi_phosphate_sym_rmsd"], METHOD_COLORS["PTM-Psi"], "--"),
        ("PyTMs optimised", pytms_ok["pytms_phosphate_sym_rmsd"], METHOD_COLORS["PyTMs optimised"], "-."),
        ("PF Stage 0 seed", single["stage0_seed_rmsd"], METHOD_COLORS["PF Stage 0 seed"], ":"),
        ("PF top-1", single["stage3_post_minimization_rmsd"], METHOD_COLORS["PF top-1"], "-"),
        ("PF top-3", single["top3_best_rmsd"], METHOD_COLORS["PF top-3"], "-"),
    ]
    for label, col, color, ls in datasets:
        vals = pd.to_numeric(col, errors="coerce").dropna().values
        recovery = [(vals <= t).sum()/len(vals)*100 for t in thresholds]
        ax_b.plot(thresholds, recovery, color=color, ls=ls, lw=2.5, label=label)
    ax_b.axvline(1.0, color="red", alpha=0.3, ls="--")
    ax_b.axvline(1.5, color="orange", alpha=0.3, ls="--")
    ax_b.set_xlabel("RMSD threshold (Å)")
    ax_b.set_ylabel("Sites recovered (%)")
    ax_b.legend(fontsize=10, loc="lower right")
    ax_b.set_xlim(0, 4); ax_b.set_ylim(0, 100)
    ax_b.grid(alpha=0.15, ls="--")
    ax_b.text(-0.08, 1.05, "b", transform=ax_b.transAxes, fontsize=20, fontweight="bold", color="black")

    fig.tight_layout()
    for ext in ["png","pdf"]:
        fig.savefig(f"{args.outdir}/fig5_tool_comparison.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure 5 saved from {len(common)} common sites.")

if __name__ == "__main__":
    main()
