#!/usr/bin/env python3
"""Figure 6: Statistical comparison. Nature-style. 2x4 layout."""
import argparse, os
import numpy as np, pandas as pd
from pathlib import Path
from scipy import stats
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

RESTYPE_MAP = {"T": "TPO", "S": "SEP", "Y": "PTR"}
TYPE_COLORS = {"TPO": "#1565c0", "SEP": "#2e7d32", "PTR": "#6a1b9a"}
METHOD_COLORS = {
    "PyTMs default": "#E8703A", "PTM-Psi": "#C2185B", "PyTMs opt.": "#F5A627",
    "PF Stage 0 seed": "#9E9E9E", "PF top-1": "#1565c0", "PF top-3": "#2e7d32",
}
RC = {"axes.titlesize":13,"axes.titleweight":"bold","axes.labelsize":14,"axes.labelweight":"bold",
      "xtick.labelsize":12,"ytick.labelsize":12,"legend.fontsize":12,"legend.frameon":False,
      "figure.facecolor":"white","axes.facecolor":"white","savefig.dpi":300}

def norm_rt(x): return RESTYPE_MAP.get(str(x).strip().upper(), str(x).strip().upper())

KEYS = ["acc_id", "pdb_key", "chain_key", "pos_key", "restype_3"]

def add_keys(df, reference_col=None):
    df = df.copy()
    if reference_col:
        df["pdb_key"] = df[reference_col].astype(str).map(lambda x: Path(x).stem.lower())
    else:
        df["pdb_key"] = df["pdb_id"].astype(str).str.lower()
    df["chain_key"] = df["chain"].astype(str)
    df["pos_key"] = pd.to_numeric(df["position"], errors="raise").astype(int)
    return df

def p_math(pval):
    exponent = int(np.floor(np.log10(max(pval, 1e-300))))
    coefficient = pval / 10**exponent
    return f"$p$ = {coefficient:.1f} × 10$^{{{exponent}}}$"

def draw_group_bracket(ax, y0, y1, label, color, x=-0.58, arm=0.035):
    """Draw a two-level y-axis group label outside the plotting area."""
    transform = ax.get_yaxis_transform()
    line = dict(transform=transform, color="#555555", lw=1.0, clip_on=False)
    ax.plot([x, x], [y0, y1], **line)
    ax.plot([x, x + arm], [y0, y0], **line)
    ax.plot([x, x + arm], [y1, y1], **line)
    ax.text(
        x - 0.035,
        (y0 + y1) / 2,
        label,
        transform=transform,
        ha="right",
        va="center",
        fontsize=9,
        fontweight="bold",
        color=color,
        clip_on=False,
    )

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

def draw_scatter(ax, pf_col, comp_col, data, comp_name, panel_label):
    for rt in ["TPO", "SEP", "PTR"]:
        sub = data[data["restype_3"]==rt]
        pf = pd.to_numeric(sub[pf_col], errors="coerce")
        comp = pd.to_numeric(sub[comp_col], errors="coerce")
        v = pf.notna() & comp.notna()
        ax.scatter(comp[v], pf[v], c=TYPE_COLORS[rt], marker="o",
                   label=rt, alpha=0.5, s=18, edgecolors="white", linewidth=0.3)
    ax.plot([0,5],[0,5],"k--",alpha=0.4,lw=1)
    ax.set_xlim(0,5); ax.set_ylim(0,5); ax.set_aspect("equal")
    ax.set_xlabel(f"{comp_name} RMSD (Å)")
    ax.set_ylabel("PF top-1 RMSD (Å)")
    ax.legend(fontsize=9, loc="upper left")
    # Stats top-right
    pf_a = pd.to_numeric(data[pf_col], errors="coerce")
    c_a = pd.to_numeric(data[comp_col], errors="coerce")
    v = pf_a.notna() & c_a.notna()
    pf_w = int((pf_a[v]<c_a[v]).sum()); c_w = int((c_a[v]<pf_a[v]).sum()); n = int(v.sum())
    ax.text(0.95, 0.95,
            f"PF: {pf_a[v].median():.2f} Å ({pf_w}/{n}, {pf_w/n*100:.0f}%)\n"
            f"{comp_name}: {c_a[v].median():.2f} Å ({c_w}/{n}, {c_w/n*100:.0f}%)",
            transform=ax.transAxes, fontsize=8, ha="right", va="top",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))
    ax.text(-0.12, 1.05, panel_label, transform=ax.transAxes, fontsize=20, fontweight="bold", color="black")
    ax.set_title(f"PF vs {comp_name}", pad=10)

def draw_mcnemar(ax, pf_col, comp_col, data, comp_name, panel_label, threshold=1.0):
    pf = pd.to_numeric(data[pf_col], errors="coerce")
    comp = pd.to_numeric(data[comp_col], errors="coerce")
    v = pf.notna() & comp.notna()
    pf_g = pf[v] <= threshold; c_g = comp[v] <= threshold
    both = int((pf_g&c_g).sum()); pf_only = int((pf_g&~c_g).sum())
    c_only = int((~pf_g&c_g).sum()); neither = int((~pf_g&~c_g).sum())
    # Green = PF wins, Red = PF loses, Grey = both/neither
    colors = [["#E0E0E0","#A5D6A7"],["#EF9A9A","#E0E0E0"]]
    mat = np.array([[both,pf_only],[c_only,neither]])
    ax.set_xlim(-0.5,1.5); ax.set_ylim(-0.5,1.5); ax.invert_yaxis()
    for i in range(2):
        for j in range(2):
            rect = plt.Rectangle((j-0.5,i-0.5),1,1,facecolor=colors[i][j],edgecolor="white",lw=2)
            ax.add_patch(rect)
            ax.text(j,i,f"{mat[i,j]}",ha="center",va="center",fontsize=18,fontweight="bold",color="#333333")
    ax.set_xticks([0,1]); ax.set_xticklabels([f"{comp_name} ≤1Å",f"{comp_name} >1Å"],fontsize=10)
    ax.set_yticks([0,1]); ax.set_yticklabels(["PF ≤1Å","PF >1Å"],fontsize=10)
    ratio = pf_only/max(c_only,1)
    disc = pf_only+c_only
    if disc > 0:
        chi2 = (abs(pf_only-c_only)-1)**2/disc
        pval = stats.chi2.sf(chi2,1)
        ax.set_title(f"PF advantage {pf_only}:{c_only} ({ratio:.1f}:1)\n{p_math(pval)}")
    ax.set_aspect("equal")
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.tick_params(length=0)
    ax.text(-0.12, 1.05, panel_label, transform=ax.transAxes, fontsize=20, fontweight="bold", color="black")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phosphofill-tsv", nargs="+", required=True)
    ap.add_argument("--pytms-default-detail", required=True)
    ap.add_argument("--pytms-detail", required=True)
    ap.add_argument("--ptmpsi-detail", required=True)
    ap.add_argument("--audit-tsv", nargs="*", default=[])
    ap.add_argument("--outdir", default="figures")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    plt.rcParams.update(RC)

    single = load_single(args.phosphofill_tsv)
    pytms_default_ok = pd.read_csv(args.pytms_default_detail, sep="\t"); pytms_default_ok = add_keys(pytms_default_ok[pytms_default_ok["status"]=="OK"], "reference_path")
    pytms_ok = pd.read_csv(args.pytms_detail, sep="\t"); pytms_ok = add_keys(pytms_ok[pytms_ok["status"]=="OK"], "reference_path")
    ptmpsi_ok = pd.read_csv(args.ptmpsi_detail, sep="\t"); ptmpsi_ok = add_keys(ptmpsi_ok[ptmpsi_ok["status"]=="OK"], "reference_path")

    # Use the identical successfully evaluated cohort for every panel.
    common = single[KEYS].drop_duplicates()
    for df in [pytms_default_ok, pytms_ok, ptmpsi_ok]:
        common = common.merge(df[KEYS].drop_duplicates(), on=KEYS, how="inner")
    single = single.merge(common, on=KEYS, how="inner")
    pytms_default_ok = pytms_default_ok.merge(common, on=KEYS, how="inner")
    pytms_ok = pytms_ok.merge(common, on=KEYS, how="inner")
    ptmpsi_ok = ptmpsi_ok.merge(common, on=KEYS, how="inner")

    merged_audit = pd.DataFrame()
    if args.audit_tsv:
        afs = []
        for af in args.audit_tsv:
            adf = pd.read_csv(af, sep="\t")
            adf = adf[adf["status"]=="ok"] if "status" in adf.columns else adf
            for res in ["TPO","SEP","PTR"]:
                if res in Path(af).stem.upper(): adf["restype_3"] = res
            adf["pdb_key"] = adf["pdb_id"].astype(str).str.lower()
            adf["chain_key"] = adf["chain_id"].astype(str)
            adf["pos_key"] = pd.to_numeric(adf["target_position"], errors="raise").astype(int)
            afs.append(adf)
        audit = pd.concat(afs, ignore_index=True)
        merged_audit = single.merge(audit[KEYS + ["interaction_class"]], on=KEYS, how="inner")

    fig = plt.figure(figsize=(24, 12))
    gs = gridspec.GridSpec(2, 4, hspace=0.4, wspace=0.78)

    # (a) Scatter PF vs PyTMs
    draw_scatter(fig.add_subplot(gs[0,0]), "stage3_post_minimization_rmsd", "pytms_phosphate_sym_rmsd", pytms_ok, "PyTMs opt.", "a")

    # (b) Scatter PF vs PTM-Psi
    ax_b = fig.add_subplot(gs[0,1])
    draw_scatter(ax_b, "stage3_post_minimization_rmsd", "ptmpsi_phosphate_sym_rmsd", ptmpsi_ok, "PTM-Psi", "b")

    # (c) McNemar vs PyTMs
    draw_mcnemar(fig.add_subplot(gs[0,2]), "stage3_post_minimization_rmsd", "pytms_phosphate_sym_rmsd", pytms_ok, "PyTMs", "c")

    # (d) McNemar vs PTM-Psi
    ax_d = fig.add_subplot(gs[0,3])
    draw_mcnemar(ax_d, "stage3_post_minimization_rmsd", "ptmpsi_phosphate_sym_rmsd", ptmpsi_ok, "PTM-Psi", "d")

    # (e) Bootstrap CI
    ax_e = fig.add_subplot(gs[1,0])
    rng = np.random.RandomState(42)
    y_pos=0.0; y_ticks=[]; y_labels=[]; e_groups=[]
    ci_list = [
        ("PyTMs default", pytms_default_ok, "pytms_phosphate_sym_rmsd", METHOD_COLORS["PyTMs default"]),
        ("PyTMs opt.", pytms_ok, "pytms_phosphate_sym_rmsd", METHOD_COLORS["PyTMs opt."]),
        ("PTM-Psi", ptmpsi_ok, "ptmpsi_phosphate_sym_rmsd", METHOD_COLORS["PTM-Psi"]),
        ("PF Stage 0 seed", single, "stage0_seed_rmsd", METHOD_COLORS["PF Stage 0 seed"]),
        ("PF top-1", single, "stage3_post_minimization_rmsd", METHOD_COLORS["PF top-1"]),
        ("PF top-3", single, "top3_best_rmsd", METHOD_COLORS["PF top-3"]),
    ]
    for rt in ["TPO","SEP","PTR"]:
        group_start = y_pos
        for mn, data, col, color in ci_list:
            sub = data[data["restype_3"]==rt]
            vals = pd.to_numeric(sub[col], errors="coerce").dropna().values
            if len(vals)<5: continue
            boot = [np.median(rng.choice(vals,len(vals),True)) for _ in range(10000)]
            lo,hi = np.percentile(boot,[2.5,97.5])
            med = np.median(vals)
            ax_e.errorbar(med,y_pos,xerr=[[med-lo],[hi-med]],fmt="o",color=color,markersize=5,capsize=3,lw=1.5)
            y_ticks.append(y_pos); y_labels.append(f"{rt} {mn}"); y_pos+=1
        group_end = y_pos - 1
        e_groups.append((rt, group_start - 0.35, group_end + 0.35))
        y_pos+=0.8
    method_labels = [label.split(" ", 1)[1] for label in y_labels]
    ax_e.set_yticks(y_ticks); ax_e.set_yticklabels(method_labels, fontsize=8)
    ax_e.set_xlabel("Median RMSD (Å)"); ax_e.axvline(1.0,color="red",alpha=0.3,ls="--")
    ax_e.set_ylim(-0.7, y_pos - 0.45); ax_e.invert_yaxis(); ax_e.grid(axis="x",alpha=0.15,ls="--")
    ax_e.tick_params(axis="y", length=0, pad=4)
    for rt, y0, y1 in e_groups:
        draw_group_bracket(ax_e, y0, y1, rt, TYPE_COLORS[rt], x=-0.58)
    ax_e.text(-0.12, 1.05, "e", transform=ax_e.transAxes, fontsize=20, fontweight="bold", color="black")

    # (f) Scoring gap
    ax_f = fig.add_subplot(gs[1,1])
    gap_all = pd.to_numeric(single["stage2_scoring_gap"],errors="coerce")
    rmsd_all = pd.to_numeric(single["stage3_post_minimization_rmsd"],errors="coerce")
    valid = gap_all.notna() & rmsd_all.notna()
    ax_f.hexbin(gap_all[valid],rmsd_all[valid],gridsize=25,cmap="Greys",alpha=0.3,mincnt=1)
    for rt in ["TPO","SEP","PTR"]:
        sub = single[single["restype_3"]==rt]
        g = pd.to_numeric(sub["stage2_scoring_gap"],errors="coerce")
        r = pd.to_numeric(sub["stage3_post_minimization_rmsd"],errors="coerce")
        v = g.notna()&r.notna()
        ax_f.scatter(g[v],r[v],c=TYPE_COLORS[rt],marker="o",label=rt,alpha=0.5,s=15,edgecolors="white",linewidth=0.3)
    rho, rho_p = stats.spearmanr(gap_all[valid],rmsd_all[valid])
    ax_f.set_xlabel("Scoring gap (Å)"); ax_f.set_ylabel("Final RMSD (Å)")
    ax_f.legend(fontsize=9, loc="upper left")
    ax_f.text(0.95,0.95,f"ρ = {rho:.2f}\n{p_math(rho_p)}",transform=ax_f.transAxes,fontsize=10,ha="right",va="top",
              bbox=dict(boxstyle="round",facecolor="lightyellow",alpha=0.8))
    ax_f.set_xlim(-0.1,4); ax_f.set_ylim(-0.1,4); ax_f.plot([0,4],[0,4],"k--",alpha=0.3)
    ax_f.grid(alpha=0.15, ls="--")
    ax_f.text(-0.12, 1.05, "f", transform=ax_f.transAxes, fontsize=20, fontweight="bold", color="black")

    # (g) Environment
    ax_g = fig.add_subplot(gs[1,2])
    if len(merged_audit)>0:
        dp=[]; le=[]; ce=[]; positions=[]; g_groups=[]
        y_pos = 1.0
        for rt in ["TPO","SEP","PTR"]:
            group_start = y_pos
            sub = merged_audit[merged_audit["restype_3"]==rt]
            for ic,lab in [("salt_bridge_like","Salt bridge"),("none_obvious","No contact")]:
                g = sub[sub["interaction_class"]==ic]
                vals = pd.to_numeric(g["stage3_post_minimization_rmsd"],errors="coerce").dropna()
                if len(vals)>=5:
                    dp.append(vals.values); le.append(f"{lab} (n={len(vals)})"); ce.append(TYPE_COLORS[rt]); positions.append(y_pos); y_pos += 1.0
            group_end = y_pos - 1.0
            g_groups.append((rt, group_start - 0.35, group_end + 0.35))
            y_pos += 0.8
        bp = ax_g.boxplot(dp, patch_artist=True, showfliers=False, widths=0.6,
                          positions=positions, orientation="horizontal")
        for patch, color in zip(bp["boxes"], ce):
            patch.set_facecolor(color); patch.set_alpha(0.5)
        ax_g.set_yticks(positions); ax_g.set_yticklabels(le, fontsize=8)
        ax_g.set_xlabel("RMSD (Å)"); ax_g.axvline(1.0,color="red",alpha=0.3,ls="--")
        ax_g.set_ylim(0.3, y_pos - 0.45); ax_g.invert_yaxis(); ax_g.grid(axis="x",alpha=0.15,ls="--")
        ax_g.tick_params(axis="y", length=0, pad=4)
        for rt, y0, y1 in g_groups:
            draw_group_bracket(ax_g, y0, y1, rt, TYPE_COLORS[rt], x=-0.58)
    ax_g.text(-0.12, 1.05, "g", transform=ax_g.transAxes, fontsize=20, fontweight="bold", color="black")

    # (h) All methods bar
    ax_h = fig.add_subplot(gs[1,3])
    mo = ["PyTMs default","PTM-Psi (NERF)","PyTMs optimised","PF Stage 0 seed","PF top-1","PF top-3"]
    meds = [pd.to_numeric(pytms_default_ok["pytms_phosphate_sym_rmsd"],errors="coerce").median(),
            pd.to_numeric(ptmpsi_ok["ptmpsi_phosphate_sym_rmsd"],errors="coerce").median(),
            pd.to_numeric(pytms_ok["pytms_phosphate_sym_rmsd"],errors="coerce").median(),
            pd.to_numeric(single["stage0_seed_rmsd"],errors="coerce").median(),
            pd.to_numeric(single["stage3_post_minimization_rmsd"],errors="coerce").median(),
            pd.to_numeric(single["top3_best_rmsd"],errors="coerce").median()]
    cf = [METHOD_COLORS[k] for k in ["PyTMs default","PTM-Psi","PyTMs opt.","PF Stage 0 seed","PF top-1","PF top-3"]]
    bars = ax_h.barh(range(len(mo)),meds,color=cf,edgecolor="white")
    ax_h.bar_label(bars,fmt="%.2f",fontsize=9,padding=3)
    ax_h.set_yticks(range(len(mo))); ax_h.set_yticklabels(mo,fontsize=9)
    ax_h.set_xlabel("Median RMSD (Å)"); ax_h.axvline(1.0,color="red",alpha=0.3,ls="--")
    ax_h.set_xlim(0,2.6); ax_h.invert_yaxis(); ax_h.grid(axis="x",alpha=0.15,ls="--")
    ax_h.text(-0.12, 1.05, "h", transform=ax_h.transAxes, fontsize=20, fontweight="bold", color="black")

    for ext in ["png","pdf"]:
        fig.savefig(f"{args.outdir}/fig6_statistical_panels.{ext}",dpi=300,bbox_inches="tight")
    plt.close(fig)
    print(f"Figure 6 saved from {len(common)} common sites; environment matches: {len(merged_audit)}.")

if __name__ == "__main__":
    main()
