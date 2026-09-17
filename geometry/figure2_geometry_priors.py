from pathlib import Path
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


ROOT = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description="Rebuild manuscript Figure 2 from the empirical geometry tables.")
parser.add_argument(
    "--data-dir",
    type=Path,
    default=ROOT / "mod_geometry_prior_all",
    help="Directory containing TPO/SEP/PTR_geometry_prior.tsv.",
)
parser.add_argument(
    "--outdir",
    type=Path,
    default=ROOT / "output",
    help="Output directory for PDF and PNG files.",
)
parser.add_argument("--single", action="store_true", help="Plot only residue-type single-site rows")
args = parser.parse_args()
DATA = args.data_dir.resolve()
OUT = args.outdir.resolve()
OUT.mkdir(parents=True, exist_ok=True)

SPECS = {
    "TPO": dict(title="pThr (TPO)", color="#1565c0", n=508,
                torsion="torsion_CG2_CB_OG1_P", torsion_label="CG2-CB-OG1-P (°)",
                production=(1.613, 118.8, -58.9), empirical=(1.611, 118.0, -56.1)),
    "SEP": dict(title="pSer (SEP)", color="#2e7d32", n=388,
                torsion="torsion_CA_CB_OG_P", torsion_label="CA-CB-OG-P (°)",
                production=(1.614, 115.4, 66.3), empirical=(1.611, 114.1, None)),
    "PTR": dict(title="pTyr (PTR)", color="#6a1b9a", n=501,
                torsion="torsion_CE1_CZ_OH_P", torsion_label="CE1-CZ-OH-P (°)",
                production=(1.609, 125.1, None), empirical=(1.607, 124.9, None)),
}


def kde(ax, values, color, xmin, xmax):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    x = np.linspace(xmin, xmax, 700)
    y = gaussian_kde(values)(x)
    ax.fill_between(x, y, color=color, alpha=.25)
    ax.plot(x, y, color=color, lw=2.2)


def circular_kde(values, bandwidth=10.0):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    x = np.linspace(-180, 180, 1441)
    delta = (x[:, None] - values[None, :] + 180.0) % 360.0 - 180.0
    y = np.exp(-0.5 * (delta / bandwidth) ** 2).mean(axis=1)
    y /= bandwidth * np.sqrt(2.0 * np.pi)
    return x, y


fig, axes = plt.subplots(3, 3, figsize=(15.8, 12.6))
for col, (code, spec) in enumerate(SPECS.items()):
    df = pd.read_csv(DATA / f"{code}_geometry_prior.tsv", sep="\t")
    df = df[df["status"].eq("ok") & df["resname"].eq(code)]
    if args.single:
        df = df[df["subset_class"].eq(f"{code}_single")]
    color = spec["color"]
    axes[0, col].set_title(f'{spec["title"]}\n$n$ = {len(df)}', color=color,
                           fontsize=20, fontweight="bold", pad=13)

    length = df["d_anchor_p"].to_numpy(float)
    angle = df["angle_base_anchor_p"].to_numpy(float)
    torsion = df[spec["torsion"]].to_numpy(float)

    kde(axes[0, col], length, color, 1.38, 1.86)
    length_median = float(np.median(length))
    angle_median = float(np.median(angle))
    subset_label = "single-site" if args.single else "all-site"
    axes[0, col].axvline(length_median, color="#c0392b", ls="--", lw=1.8)
    axes[0, col].text(length_median + .008, axes[0, col].get_ylim()[1] * .88,
                      f'{length_median:.3f} Å {subset_label} median', color="#c0392b", fontweight="bold")
    axes[0, col].set_xlabel("Anchor-P (Å)", fontweight="bold")

    kde(axes[1, col], angle, color, 82, 166)
    axes[1, col].axvline(angle_median, color="#c0392b", ls="--", lw=1.8)
    axes[1, col].text(angle_median + 1.5, axes[1, col].get_ylim()[1] * .88,
                      f'{angle_median:.1f}° {subset_label} median', color="#c0392b", fontweight="bold")
    axes[1, col].set_xlabel("Base-anchor-P (°)", fontweight="bold")

    ax = axes[2, col]
    bins = np.arange(-180, 181, 10)
    ax.hist(torsion, bins=bins, density=True, color=color, alpha=.13, edgecolor=color, linewidth=.8)
    x, y = circular_kde(torsion)
    ax.plot(x, y, color=color, lw=2.2)
    ax.set_xlim(-180, 180)
    ax.set_xticks([-180, -90, 0, 90, 180])
    ax.set_xlabel(spec["torsion_label"], fontweight="bold")

    if code == "TPO":
        peak = float(x[np.argmax(y)])
        ax.axvline(peak, color="#c0392b", ls="--", lw=1.8)
        ax.text(.36, .92, "−56°", transform=ax.transAxes,
                color="#c0392b", fontweight="bold")
    elif code == "SEP":
        ax.axvspan(-60, 60, color=color, alpha=.08)
        ax.axvline(66.0, color="#c0392b", ls="--", lw=1.8)
        ax.text(.68, .92, "+66° sharp mode", transform=ax.transAxes,
                color="#c0392b", fontweight="bold")
        ax.text(.08, .80, "Broad basin: −60° to +60°", transform=ax.transAxes,
                color=color, fontweight="bold")
    else:
        ax.text(0, ax.get_ylim()[1] * .90, "No preferred orientation", ha="center",
                color=color, fontweight="bold", style="italic")

for row in range(3):
    axes[row, 0].set_ylabel("Density", fontweight="bold")
for ax in axes.flat:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=11)
    ax.grid(False)

fig.text(.018, .77, "Bond length", rotation=90, va="center", fontsize=16, fontweight="bold")
fig.text(.018, .48, "Valence angle", rotation=90, va="center", fontsize=16, fontweight="bold")
fig.text(.018, .19, "Key dihedral", rotation=90, va="center", fontsize=16, fontweight="bold")
fig.tight_layout(rect=(.045, .025, 1, .995), h_pad=2.0, w_pad=2.4)

suffix = "single_site" if args.single else "corrected"
pdf = OUT / f"fig2_geometry_priors_{suffix}.pdf"
png = OUT / f"fig2_geometry_priors_{suffix}.png"
fig.savefig(pdf, bbox_inches="tight")
fig.savefig(png, dpi=180, bbox_inches="tight")
print(pdf)
