#!/usr/bin/env python3
"""Figure 4: residue-filtered phosphosite environment comparisons."""
import argparse, os
import pandas as pd, numpy as np
from pathlib import Path
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

TYPE_COLORS = {"TPO": "#1565c0", "SEP": "#2e7d32", "PTR": "#6a1b9a"}
UNMOD_COLOR = "#9E9E9E"
RC = {"axes.titlesize":13,"axes.titleweight":"bold","axes.labelsize":14,"axes.labelweight":"bold",
      "xtick.labelsize":12,"ytick.labelsize":12,"legend.fontsize":12,"legend.frameon":False,
      "figure.facecolor":"white","axes.facecolor":"white","savefig.dpi":300}

def load_audit(paths, type_map):
    frames = []
    for af in paths:
        adf = pd.read_csv(af, sep="\t")
        matches = [res for res, key in type_map.items() if key in Path(af).stem.upper()]
        if len(matches) != 1:
            raise ValueError(f"Cannot identify one expected residue type from {af}")
        expected = matches[0]
        if "status" not in adf or "actual_resname" not in adf:
            raise ValueError(f"Missing status or actual_resname column in {af}")
        adf = adf.loc[
            adf["status"].eq("ok")
            & adf["actual_resname"].astype(str).str.upper().eq(expected)
        ].copy()
        adf["restype_3"] = expected
        frames.append(adf)
    return pd.concat(frames, ignore_index=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mod-audit", nargs="+", required=True)
    ap.add_argument("--unmod-audit", nargs="+", required=True)
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    plt.rcParams.update(RC)

    mod = load_audit(args.mod_audit, {"TPO":"TPO","SEP":"SEP","PTR":"PTR"})
    unmod = load_audit(args.unmod_audit, {"THR":"THR","SER":"SER","TYR":"TYR"})

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))

    for i, (unmod_res, mod_res) in enumerate([("THR","TPO"),("SER","SEP"),("TYR","PTR")]):
        ax = axes[i]
        u = pd.to_numeric(unmod[unmod["restype_3"]==unmod_res]["nearest_basic_dist"], errors="coerce").dropna()
        m = pd.to_numeric(mod[mod["restype_3"]==mod_res]["nearest_basic_dist"], errors="coerce").dropna()
        _, pval = stats.mannwhitneyu(u, m, alternative="two-sided", method="asymptotic")

        clip_max = 18.0
        u_c = u.clip(upper=clip_max)
        m_c = m.clip(upper=clip_max)

        # Violin
        vp = ax.violinplot([u_c.values, m_c.values], positions=[1, 2],
                           showmeans=False, showmedians=False, showextrema=False, widths=0.7)
        vp["bodies"][0].set_facecolor(UNMOD_COLOR); vp["bodies"][0].set_alpha(0.35)
        vp["bodies"][1].set_facecolor(TYPE_COLORS[mod_res]); vp["bodies"][1].set_alpha(0.35)

        # Thin box
        bp = ax.boxplot([u_c.values, m_c.values], positions=[1, 2], widths=0.15,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color="black", lw=2),
                        whiskerprops=dict(lw=1), capprops=dict(lw=1))
        bp["boxes"][0].set_facecolor(UNMOD_COLOR); bp["boxes"][0].set_alpha(0.8)
        bp["boxes"][1].set_facecolor(TYPE_COLORS[mod_res]); bp["boxes"][1].set_alpha(0.8)

        # Median annotations — above whisker tops
        u_q75 = u_c.quantile(0.75)
        u_iqr = u_q75 - u_c.quantile(0.25)
        u_wtop = min(u_q75 + 1.5 * u_iqr, u_c.max())
        m_q75 = m_c.quantile(0.75)
        m_iqr = m_q75 - m_c.quantile(0.25)
        m_wtop = min(m_q75 + 1.5 * m_iqr, m_c.max())
        ax.text(1, u_wtop + 0.3, f"{u.median():.1f} Å", ha="center", fontsize=10, fontweight="bold", color="#555555")
        ax.text(2, m_wtop + 0.3, f"{m.median():.1f} Å", ha="center", fontsize=10, fontweight="bold", color="black")

        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"{unmod_res}\n(n={len(u)})", f"{mod_res}\n(n={len(m)})"], fontsize=11, fontweight="bold")
        ax.set_ylabel("Distance to nearest basic-residue atom (Å)" if i == 0 else "")

        pexp = int(np.floor(np.log10(max(pval, 1e-300))))
        pcoef = pval / (10 ** pexp)
        ax.set_title(f"{unmod_res} → {mod_res}   ($p$ = {pcoef:.1f} × 10$^{{{pexp}}}$)")
        ax.set_ylim(0, min(clip_max + 1, u.quantile(0.95) + 3))
        ax.grid(axis="y", alpha=0.15, ls="--")
        print(f"{unmod_res}/{mod_res}: n={len(u)}/{len(m)}, "
              f"medians={u.median():.4f}/{m.median():.4f} Å, p={pval:.4e}")

        # Panel label
        ax.text(-0.12, 1.05, chr(97+i), transform=ax.transAxes, fontsize=20, fontweight="bold", color="black")

    fig.tight_layout()
    fig.savefig(f"{args.outdir}/fig4_environment_corrected.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Corrected Figure 4 saved.")

if __name__ == "__main__":
    main()
