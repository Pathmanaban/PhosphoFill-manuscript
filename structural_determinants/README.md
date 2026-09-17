# Structural determinants: reproducible correction audit

This directory audits the manuscript's structural-determinants section without
editing the original source files or `phosphoFill_final(2).tex`. The starting
point was `C:/Users/drpat/Downloads/stat_test_4_environment.py`; its outputs
were **not** treated as validated results. The corrected analysis and all of its
source-file SHA-256 hashes are in `outputs/audit_summary.json`.

## Sources and scope

- `../geometry/mod_context_audit_all/{TPO,SEP,PTR}_context_audit.tsv`
  and `../geometry/unmod_context_audit_all/{THR,SER,TYR}_context_audit.tsv`:
  1,397 actual modified and 4,513 actual unmodified structure records after
  requiring `status=ok` and `actual_resname` to match the intended residue.
- Corresponding `*_context_audit_single` TSVs: single-site orientation summaries.
- Corresponding `*_geometry_prior_{single,all}/*.tsv` and `*.summary.json`:
  side-chain torsions. The raw TSV counts are checked against JSON histograms.
- `../benchmark_validation/outputs/raw/{TPO,SEP,PTR}_postmin.tsv`:
  790 successful primary-type single-site strip-and-regraft results.
- `comparator_validation/outputs/analysis/common_785_site_metrics.tsv`:
  the common 785-site comparator cohort used in Figure 6.
- The external general PDB cache: 30 direct-coordinate torsion checks via
  Bio.PDB. These checks show that the extractor/production dihedral angle
  equals `wrap(180 degrees - Bio.PDB dihedral)`. The context-audit script's
  dihedral output uses the opposite sign from the extractor. This does not
  invalidate the production seeds, which use the extractor convention, but
  the old chi1 **rotamer names** were incorrect.

All default source paths are derived from this folder's location. They can be
overridden with `--phospho-root`, `--benchmark-dir`, and `--common-785`.

## Reproduce

From this directory in the user's Ubuntu Python environment with `numpy`,
`pandas`, `scipy`, `scikit-learn`, `matplotlib`, and `biopython` installed:

```bash
python scripts/audit_structural_determinants.py
python scripts/verify_torsion_convention.py
python scripts/analyze_contact_rmsd_and_s4.py
```

To rebuild Figure 6 with explicit contact-class labels, run from
`../comparator_validation`:

```bash
python ../structural_determinants_audit/scripts/fig6_statistical_panels_environment_corrected.py \
  --phosphofill-tsv ../benchmark_validation/outputs/raw/TPO_postmin.tsv ../benchmark_validation/outputs/raw/SEP_postmin.tsv ../benchmark_validation/outputs/raw/PTR_postmin.tsv \
  --pytms-default-detail inputs/comparators/pytms_detail.tsv \
  --pytms-detail inputs/comparators/pytms_detail_optimized.tsv \
  --ptmpsi-detail inputs/comparators/ptmpsi_detail.tsv \
  --audit-tsv inputs/environment/TPO_context_audit.tsv inputs/environment/SEP_context_audit.tsv inputs/environment/PTR_context_audit.tsv \
  --outdir ../structural_determinants_audit/figures
```

The existing residue-filtered Figure 4 is
`../figures/main/figure4_environment.pdf`. Its numerical medians and
two-sided Mann--Whitney p-values agree with this audit, so it was not redrawn.
The earlier label-only Figure 6 is retained for provenance at
`figures/fig6_statistical_panels_environment_corrected.pdf`. Panel g retains
the same 773 exact matches and boxplots; only the contact-class labels were
made operationally precise. Other Figure 6 panels were not changed.

The final contact-analysis Figure 6 is retained with the script-generated
outputs in this folder's `figures/` directory. Its panel g uses all
775 exact context matches from the 790-site PhosphoFill benchmark and compares
`salt_bridge_like` with the other three audit classes combined. Panels a--f
and h still use the common 785-site tool-comparison cohort. Panel f's x-axis
now identifies its actual quantity, the retrospective S2--S1 RMSD difference;
the plotted 785-site Spearman result is rho=0.85322318,
p=1.4580068e-223. Rebuild this PDF using the Figure 6 command above with
`--environment-contrast salt_vs_all_other --output-basename
fig6_statistical_panels_contact_contrast`.

`outputs/contact_rmsd_and_s4_new_analysis.json` records the new retrospective
association analysis. The primary contrast is salt-bridge-like versus all
other annotated classes, separately for TPO, SEP and PTR. The narrow
salt-versus-none and salt-versus-none-or-water comparisons are sensitivity
checks. Median RMSD differences and 1-A recovery differences have 95%
percentile intervals from 5,000 bootstrap resamples of protein/accession
clusters (seed 20260914); every structure from a sampled protein is retained
together. No causal or prospective-prediction interpretation is assigned to
these reference-derived contact classes. Top-1 and retrospective best-of-three
RMSD are both recorded; Figure 6(g) displays top-1.

| Type | Salt-like / other structures | Top-1 medians, salt / other (A) | Median difference, salt minus other (A) | 95% protein-cluster bootstrap interval (A) |
| --- | ---: | ---: | ---: | ---: |
| TPO | 264 / 55 | 0.410 / 0.425 | -0.015 | [-0.429, 0.415] |
| SEP | 153 / 86 | 0.558 / 0.517 | +0.041 | [-0.878, 0.236] |
| PTR | 180 / 37 | 1.002 / 1.527 | -0.526 | [-1.624, -0.053] |

Only PTR has an interval excluding zero in the primary top-1 contrast. Its
interval includes zero when the other group is restricted to no-obvious-contact
plus water-proximal records, and its top-3-best interval also includes zero.
These are exploratory associations conditional on the experimental contact
labels, not evidence of a universal salt-bridge effect.

The same JSON audits Supplementary Table S4 nesting. The single-site and
all-site selections are not strictly nested: TPO 275/276 and SEP 227/228
single-site exact keys occur in their all-site audits, whereas PTR is 208/208.
In both mismatches the all-site audit selected a different PDB chain. Shared
records have identical measured orientation fields. All-site-only records
also include some records with `n_target_sites=1`, so subtracting single-site
percentages from all-site percentages does **not** isolate a pure multi-site
effect.

## What the original statistics script got wrong

1. It accepted `status=ok` rows without checking `actual_resname`. That adds
   three PTR residues from the TPO audit, 14 TYR residues from the THR audit,
   and two THR residues from the TYR audit to the wrong cohorts.
2. It joined the benchmark to the modified audit on accession, position, and
   type but not PDB and chain. An audit class from one structure could be
   broadcast to another structure at the same protein position. The exact
   structure-level key yields 775/790 matches. Of 15 unmatched benchmark
   records, 12 have no audit record and three are residue mismatches; all are
   listed in `outputs/benchmark_unmatched_environment.tsv`.
3. Its logistic-regression code listed four features but merged only two,
   silently omitting `n_basic_within6` and `nonwater_heavy_within6`. The
   manuscript's `+1.28` coefficient cannot be reproduced from that script.
4. The old chi1 labels used nonstandard bins, calling the strong gauche-minus
   population “trans.” The conventional chi1 analysis is in
   `outputs/table_s3_rotamers.tsv`; the old-bin result is retained only as
   `outputs/legacy_rotamer_classification.tsv` for traceability.

## Definitions and interpretation limits

The context audit searches selected Arg, Lys, and His side-chain atoms within
8 A. Its basic-atom set contains Arg NE/NH1/NH2/CZ, Lys NZ/CE, and His
ND1/NE2/CE1/CD2; some selected atoms are carbon. In modified residues,
nearest-basic distance is the minimum distance to any of P/O1P/O2P/O3P; in
unmodified residues it is the distance to the parent hydroxyl oxygen. Thus
the unmodified-versus-modified violin comparison has **different query
geometry** and is descriptive, not evidence that phosphorylation recruits
basic residues. Only records with a selected basic atom within 8 A appear in
Figure 4's distance distributions. Pooled PDB entries can repeat the same
protein position; a unique protein-position sensitivity analysis is supplied.

`salt_bridge_like` is a distance class: the nearest selected basic atom is
within 3.5 A of a phosphate atom. It is not an energetic or chemically
validated salt bridge. `hbond_like` is selected polar-sidechain proximity;
backbone donor geometry was not tested in this audit. `water_mediated_like`
means nearby water, not demonstrated water mediation. The four classes are
assigned hierarchically and are mutually exclusive. The audit does not
establish metal dependence or a fundamental accuracy limit in solvent-free
models.

“Toward basic” means the angle between anchor-to-P and anchor-to-nearest-basic
atom vectors is at most 60 degrees, conditional on a measurable nearest basic
atom. It is not a phosphate-oxygen donor interaction. All-site values are
71.0% TPO, 53.8% SEP, and 41.5% PTR; dedicated single-site values are 75.6%,
66.4%, and 50.5%, respectively.

The corrected contact/RMSD analysis is retrospective: classes are measured
from the **experimental phospho structure**, then joined to top-1 RMSD against
that same reference. The “no obvious contact” PTR group has only seven
benchmark structures. Salt-bridge-like versus no-obvious-contact medians are
0.41/0.32 A (TPO), 0.56/1.26 A (SEP), and 1.00/2.38 A (PTR). When water-like
records are added to the comparison group, the SEP medians are 0.56/0.47 A
with `p=0.159`. No uniform accuracy advantage should be claimed across types.
The exact four-feature logistic fit has 744 complete cases and a +1.49
standardized coefficient for `n_basic_within4`, but feature collinearity
(`r=0.836` between 4- and 6-A basic counts), reference-derived predictors,
and no external validation make “strongest predictor” inappropriate for the
manuscript.

## Output map

- `table_s2_environment.tsv`: all-site counts, distance medians, contact
  classes, and exposure classes.
- `nearest_basic_mannwhitney.tsv`,
  `unique_site_nearest_basic_sensitivity.tsv`,
  `shared_protein_position_sensitivity.tsv`, `exposure_chi2.tsv`: comparison
  statistics and sensitivity checks.
- `table_s3_rotamers.tsv`, `torsion_convention_examples.tsv`: conventional
  chi1 classification and direct-coordinate checks.
- `table_s4_orientation.tsv`: measured-direction denominators and fractions.
- `benchmark_contact_rmsd.tsv`: all three explicit comparison definitions for
  the 790-site and common 785-site cohorts.
- `benchmark_joined_environment.tsv`, `benchmark_unmatched_environment.tsv`:
  exact-match provenance.
- `logistic_coefficients.tsv`: exploratory fit, not a validated prediction.

See `manuscript_replacements.md` for paste-ready section, caption, and table
revisions. The source manuscript and original scripts were not overwritten.

## Two additional sensitivity analyses (same-site anchors and multi-site RMSD)

Run `scripts/analyze_paired_and_multisite.py`, then
`scripts/plot_new_additions.py`. The analysis reads the original context audit,
cached experimental mmCIF structures, and frozen post-minimisation benchmark
TSVs. Its new outputs are in `outputs/new_additions/`, including structure-level
join files, coordinate-read errors, all source SHA-256 hashes, and a JSON
summary. The review figure is retained in this folder's `figures/`
directory.
Neither the existing Figure 4 nor Figure 6 PDF was overwritten.

The first check matched accession and protein position between modified and
unmodified cohorts, then recalculated nearest selected basic-atom distance
from the **same bridging oxygen** (OG1, OG, or OH) in both states. Distances
were searched across the same chain without an 8-A cap. Each protein-position
contributed one median per state, regardless of the number of PDB structures.
All 159 shared positions had complete distances and zero coordinate-read
errors; 2,658 distinct cached PDB files were parsed. Six independent
unmodified-record spot checks reproduced the original audit's nearest-basic
distance exactly. The groups were not paired before/after crystal structures:
they were separate PDB observations of the same protein positions.

| Pair | Shared positions | Median modified-minus-unmodified anchor distance (A) | 95% protein-cluster bootstrap CI (A) | Both-state <=8 A sensitivity: n, median, CI (A) |
| --- | ---: | ---: | ---: | --- |
| TPO/THR | 47 | -0.542 | [-1.370, -0.090] | 32, -0.192, [-1.018, 0.100] |
| SEP/SER | 56 | -0.474 | [-0.790, -0.040] | 38, -0.458, [-0.788, -0.040] |
| PTR/TYR | 56 | -0.141 | [-0.652, 0.084] | 44, -0.029, [-0.303, 0.162] |

The common-marker differences are considerably smaller than Figure 4's
pooled phosphate-versus-hydroxyl differences. SEP has the clearest consistent
same-site signal; TPO loses a zero-excluding interval in the stricter
within-8-A subset, and PTR has no clear paired difference. This is a
sensitivity analysis, not proof that phosphorylation recruits a basic residue.

The second check extended the exact accession/PDB/chain/position/type join to
the 494 successful primary-type multi-site benchmark records. It matched 490:
160/162 TPO, 57/59 SEP, and 273/273 PTR. The four unmatched records are absent
from the corresponding audit, not replaced using another structure. All 494
source best-of-three RMSDs were independently checked to equal the minimum of
their three post-minimisation rank values.

| Type | Salt-like / all other | Top-1 medians (A) | Salt-minus-other median CI (A) | Best-of-3 medians (A) |
| --- | ---: | ---: | ---: | ---: |
| TPO | 123 / 37 | 0.933 / 1.283 | [-0.951, -0.009] | 0.576 / 0.788 |
| SEP | 30 / 27 | 0.833 / 1.456 | [-2.132, -0.049] | 0.584 / 0.533 |
| PTR | 213 / 60 | 0.964 / 1.559 | [-1.579, 0.050] | 0.600 / 0.717 |

Protein/accession-cluster bootstrap intervals used 5,000 resamples. The
top-1 intervals exclude zero for TPO and SEP, but not PTR; no best-of-three
interval excludes zero. Contacts are assigned from the *experimental*
phosphate structures and the four audit classes are distance-only proxies.
Consequently this is retrospective association, not a deployable contact
predictor or a causal contribution of the scoring function. See the optional
paste-ready paragraphs and supplementary caption at the end of
`manuscript_replacements.md`.
