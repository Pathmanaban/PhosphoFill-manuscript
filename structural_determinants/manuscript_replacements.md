# Paste-ready manuscript changes (source `.tex` left untouched)

These changes were checked against `phosphoFill_final(2).tex` in Downloads.
The numerical source is the TSV/JSON output in this directory. Only the
blocks below need replacing for the structural-determinants correction.

## 1. Results subsection: replace from `\subsection{Structural determinants...}` through the end of the chi1 paragraph

Keep the next subsection, `\subsection{Comparison with existing phosphorylation grafting tools}`, in place.

```tex
\subsection{Structural determinants of phosphate orientation}

To investigate phosphate orientation, we analysed 1,397 experimental phosphosites and compared them with 4,513 unmodified residues of the corresponding parent types (1,185 THR, 1,204 SER, and 2,124 TYR) from structures in the same PDB resolution range (Supplementary Tables~S2--S4). The comparisons use pooled structure records rather than matched phospho--unmodified structures; some protein positions are represented in both cohorts.

The phosphate group in modified structures lay closer to a selected Arg, Lys, or His side-chain atom than did the hydroxyl oxygen in the unmodified cohorts (Fig.~\ref{fig:environment}). Median nearest-basic distances were 2.77~\AA\ for TPO versus 6.00~\AA\ for THR, 3.08~\AA\ for SEP versus 5.31~\AA\ for SER, and 2.91~\AA\ for PTR versus 4.69~\AA\ for TYR (two-sided Mann--Whitney $U$ tests: $p=8.8\times10^{-132}$, $2.1\times10^{-56}$, and $1.3\times10^{-142}$, respectively). These distances are measured from different query atoms---the phosphate atoms for modified residues and the parent hydroxyl oxygen for unmodified residues---and include repeated PDB structures of some protein positions. They describe the observed local environments but do not establish that phosphorylation itself recruits basic residues.

Using the audit's distance-based categories, 402/508 TPO (79.1\%), 210/388 SEP (54.1\%), and 404/501 PTR (80.6\%) phosphosites were salt-bridge-like (Supplementary Table~S2). Smaller groups were classified as polar-contact-like or water-proximal, while 11.2\%, 23.7\%, and 8.8\%, respectively, had no obvious contact under these criteria. These categories report proximity, not interaction energy or a demonstrated dependence on water or metal ions.

Among phosphosites with a measurable nearby basic atom, the anchor-to-phosphorus vector pointed toward it in 71.0\% of all-site TPO, 53.8\% of SEP, and 41.5\% of PTR structures. The corresponding single-site fractions were 75.6\%, 66.4\%, and 50.5\% (Supplementary Table~S4). Basic-atom proximity therefore provides a possible directional cue for many sites, but is not a universal constraint on phosphate orientation.

We matched the experimental environment audit to the single-site strip-and-regraft benchmark by accession, PDB identifier, chain, residue position, and phosphoresidue type; 775 of 790 benchmark records had an exact match. In a retrospective comparison of salt-bridge-like sites with all other annotated classes, median top-1 RMSDs were 0.41 versus 0.43~\AA\ for TPO (264 versus 55 structures), 0.56 versus 0.52~\AA\ for SEP (153 versus 86), and 1.00 versus 1.53~\AA\ for PTR (180 versus 37; Fig.~\ref{fig:stats}g). Protein-cluster bootstrap 95\% intervals for the salt-like-minus-other median difference were [$-0.43$, $+0.42$]~\AA\ for TPO, [$-0.88$, $+0.24$]~\AA\ for SEP, and [$-1.62$, $-0.05$]~\AA\ for PTR. The PTR association was not uniform across the narrower alternative contact-group definitions, and the experimental phosphate coordinates define these classes; this analysis does not establish a causal effect or a predictor available for new structures.

The sidechain $\chi_1$ distributions also differed between phosphorylated and unmodified cohorts (Supplementary Table~S3). After assigning conventional $\chi_1$ rotamer names, the gauche$-$ fraction was 55.9\% for single-site THR and 87.3\% for TPO. Single-site SER was predominantly gauche$+$ (61.8\%), whereas SEP had its largest fraction in gauche$-$ (53.9\%). Four-class tests on the all-site cohorts gave $\chi^2=178.0$, $p=2.4\times10^{-38}$, Cram\'{e}r's $V=0.32$ for THR/TPO and $\chi^2=87.5$, $p=7.6\times10^{-19}$, $V=0.23$ for SER/SEP. These are differences between structural cohorts, not observed before-and-after rotamer switches at matched sites.
```

## 2. Figure 4 caption: replace the complete caption

The existing `Benchmark_combined_final/fig4_environment.pdf` already has the
correct residue-filtered data and the p-values given above. Only its caption
needs revision.

```tex
\caption{\textbf{Basic-atom proximity in phosphorylated and unmodified PDB cohorts.} Distributions of the nearest distance to a selected Arg, Lys, or His side-chain atom. For modified residues, the distance is the minimum from P, O1P, O2P, or O3P; for unmodified residues, it is measured from the parent hydroxyl oxygen. Only records with a selected basic atom within 8~\AA\ are shown. Medians and two-sided Mann--Whitney $U$ test $p$-values compare pooled structure records; the different query-atom sets and repeated PDB entries preclude interpreting these comparisons as a direct effect of phosphorylation.}
```

## 3. Supplementary Table S2: replace the complete table

Percentages use all actual-type, status-ok records in each column. The nearest
distance median uses only records with a basic atom within the 8-A search
radius. “Buried” and “exposed” are the audit's RSA/fallback classifications.

```tex
\begin{table}[h!]
\caption*{\textbf{Supplementary Table S2.} Distance-based structural context of modified phosphosites and unmodified parent residues in the all-site PDB cohorts. Interaction classes apply only to phosphate-containing residues and are mutually exclusive, hierarchical proximity labels rather than validated salt bridges or hydrogen bonds. Nearest-basic medians use only records with a selected Arg, Lys, or His side-chain atom within 8~\AA; class percentages use all records of the indicated residue type.}
\centering
\small
\begin{tabularx}{\textwidth}{l c c c c c c}
\toprule
\textbf{Metric} & \textbf{THR} & \textbf{TPO} & \textbf{SER} & \textbf{SEP} & \textbf{TYR} & \textbf{PTR} \\
\midrule
$n$ & 1,185 & 508 & 1,204 & 388 & 2,124 & 501 \\
Nearest basic measured ($n$) & 975 & 486 & 904 & 338 & 1,572 & 492 \\
Nearest basic median (\AA) & 6.00 & 2.77 & 5.31 & 3.08 & 4.69 & 2.91 \\
Salt-bridge-like (\%) & -- & 79.1 & -- & 54.1 & -- & 80.6 \\
Polar-contact-like (\%) & -- & 5.5 & -- & 6.2 & -- & 2.6 \\
Water-proximal (\%) & -- & 4.1 & -- & 16.0 & -- & 8.0 \\
None obvious (\%) & -- & 11.2 & -- & 23.7 & -- & 8.8 \\
Buried (\%) & 37.5 & 12.0 & 38.9 & 9.5 & 38.5 & 2.8 \\
Exposed (\%) & 16.7 & 35.2 & 18.8 & 33.2 & 9.0 & 31.3 \\
\bottomrule
\end{tabularx}
\end{table}
```

## 4. Supplementary Table S3: replace the complete table

This is the substantive rotamer correction. The old `trans 56% -> 87%`
statement is really gauche-minus enrichment. Raw-coordinate checks on 30
experimental torsions established the extractor-to-standard transformation.
Values below are single-site / all-site.

```tex
\begin{table}[h!]
\caption*{\textbf{Supplementary Table S3.} Sidechain $\chi_1$ rotamers in unmodified and phosphorylated PDB cohorts. Extractor torsions were converted to the conventional dihedral convention before assigning gauche$-$ ($-90^{\circ}$ to $-30^{\circ}$), gauche$+$ ($+30^{\circ}$ to $+90^{\circ}$), and trans ($|\chi_1|\geq150^{\circ}$); other values fall outside these intervals. Values are single-site / all-site. Four-class chi-squared tests on the all-site cohorts give THR/TPO $\chi^2=178.0$, $p=2.4\times10^{-38}$, $V=0.32$ and SER/SEP $\chi^2=87.5$, $p=7.6\times10^{-19}$, $V=0.23$. These are unpaired structural cohorts.}
\centering
\small
\begin{tabularx}{\textwidth}{l c c c c c}
\toprule
\textbf{Residue} & \textbf{$n$} & \textbf{Gauche$-$ (\%)} & \textbf{Gauche$+$ (\%)} & \textbf{Trans (\%)} & \textbf{Other (\%)} \\
\midrule
THR & 490 / 1,185 & 55.9 / 52.0 & 32.4 / 38.3 & 9.8 / 6.5 & 1.8 / 3.2 \\
TPO & 276 / 508 & 87.3 / 67.7 & 3.3 / 8.3 & 6.2 / 17.3 & 3.3 / 6.7 \\
\addlinespace
SER & 714 / 1,204 & 20.4 / 29.0 & 61.8 / 55.6 & 14.7 / 12.1 & 3.1 / 3.2 \\
SEP & 228 / 388 & 53.9 / 39.9 & 21.1 / 34.0 & 15.4 / 13.1 & 9.6 / 12.9 \\
\bottomrule
\end{tabularx}
\end{table}
```

## 5. Supplementary Table S4: replace the complete table

The percentages are conditional on a measurable nearest basic atom, not all
phosphosites. “Same side” is a separate side-of-plane metric; it does not sum
with toward and lateral. Values are single-site / all-site. The cohorts were
selected independently and should not be subtracted to infer a multi-site
effect; one TPO and one SEP single-site chain is absent from the corresponding
all-site selection.

```tex
\begin{table}[h!]
\caption*{\textbf{Supplementary Table S4.} Phosphate orientation relative to the nearest selected basic-residue atom. “Toward” means an angle $\leq60^{\circ}$ between the anchor-to-P and anchor-to-basic-atom vectors; “lateral” means $60^{\circ}$--$120^{\circ}$. Opposite orientations account for the remaining evaluable sites. “Same side” is an independent side-of-plane measure. Fractions use the measurable-orientation denominator and are reported as single-site / all-site. The two cohorts were selected separately and are not a paired comparison.}
\centering
\small
\begin{tabularx}{\textwidth}{l c c c c c}
\toprule
\textbf{Residue} & \textbf{All $n$} & \textbf{Evaluable $n$} & \textbf{Toward (\%)} & \textbf{Lateral (\%)} & \textbf{Same side (\%)} \\
\midrule
TPO & 276 / 508 & 270 / 486 & 75.6 / 71.0 & 15.2 / 23.5 & 84.4 / 82.7 \\
SEP & 228 / 388 & 211 / 338 & 66.4 / 53.8 & 32.2 / 44.4 & 78.2 / 73.1 \\
PTR & 208 / 501 & 204 / 492 & 50.5 / 41.5 & 49.5 / 57.5 & 77.5 / 73.6 \\
\bottomrule
\end{tabularx}
\end{table}
```

## 6. Figure 6 and nearby Results text

The final PDF is `Benchmarkv2/output/pdf/fig6_statistical_panels_contact_contrast.pdf`
(also reproducible under this audit folder's `figures/` directory).
When applying these changes, use it in place of the old
`media/fig6_statistical_panels.pdf`. Panels a--e and h retain the common 785
comparative cohort. Panel f also retains those 785 sites, but its plotted
Spearman statistic is $\rho=0.8532$, $p=1.4580\times10^{-223}$, not the
caption's older $0.73$ and $8.0\times10^{-134}$; the x-axis is specifically the
retrospective S2--S1 RMSD difference. Panel g now shows the new primary contact
contrast on 775 exact-matched sites from the full 790-site benchmark.

Replace the sentence in the comparison subsection beginning `To ensure a fair
direct comparison` with:

```tex
To ensure a fair direct comparison, Figures~\ref{fig:toolcomp} and~\ref{fig:stats}a--f,h use the 785 sites successfully processed by every method (328 TPO, 241 SEP, and 216 PTR), matched by accession, PDB identifier, chain, residue position, and phosphoresidue type. Figure~\ref{fig:stats}g instead uses 775 of the 790 single-site PhosphoFill benchmark records that matched the experimental environment audit at the same structure-level key.
```

Replace the entire Figure 6 caption with:

```tex
\caption{\textbf{Statistical comparison of PhosphoFill with existing phosphorylation grafting methods.} Panels a--f and h use the same 785 single-site phosphosites successfully processed by every method; panel g uses 775 exact structure-level matches from the full 790-site PhosphoFill benchmark. \textbf{(a,b)} Paired site-level comparison of PhosphoFill top-1 RMSD with PyTMs optimised and PTM-Psi, respectively. PhosphoFill produces the lower-RMSD placement for 619/785 sites (79\%) relative to PyTMs optimised and 670/785 sites (85\%) relative to PTM-Psi. \textbf{(c,d)} Continuity-corrected McNemar analyses at a 1.0~\AA\ recovery threshold. PhosphoFill uniquely recovers 336 sites compared with 49 uniquely recovered by PyTMs optimised ($p=4.0\times10^{-48}$) and 511 sites compared with 30 uniquely recovered by PTM-Psi ($p=1.3\times10^{-94}$). \textbf{(e)} Median RMSD and bootstrap 95\% confidence intervals by residue type and method, calculated using 10,000 bootstrap resamples. \textbf{(f)} Association between the retrospective S2--S1 RMSD difference and final top-1 RMSD (Spearman $\rho=0.85$, $p=1.5\times10^{-223}$); S1 uses experimental phosphate coordinates and is not available for new predictions. \textbf{(g)} Top-1 RMSD stratified by the experimental structure's distance-based contact class: salt-bridge-like versus the combined polar-contact-like, water-proximal, and no-obvious-contact classes. These classes are descriptive and do not establish interaction energy or causation. \textbf{(h)} Overall median RMSD for PyTMs default, PTM-Psi, PyTMs optimised, and PhosphoFill Stage~0 seed, top-1, and top-3 outputs.}
```

In the comparison-with-tools Results paragraph beginning `With VDW
optimisation enabled`, replace the claim that SEP and PTR orientations
*depend* on electrostatic interactions with:

```tex
However, VDW clash avoidance alone is insufficient for SEP and PTR, where local structural context can influence orientation beyond steric constraints (Table~\ref{tab:toolcomp}).
```

At the end of that paragraph, replace the two sentences beginning `The
contribution of local interaction environments...` with:

```tex
Figure~\ref{fig:stats}g stratifies top-1 RMSD by contact classes measured in the corresponding experimental phospho structures. Salt-bridge-like sites had a lower median than the other classes for PTR, but not consistently for TPO or SEP. These reference-derived classes cannot be known before phosphate grafting and do not establish a causal contribution to placement accuracy.
```

The visual check of Figure 6(h) also exposed a linked, independent Table 4
rounding/label inconsistency in the current TeX. The common-cohort top-3
median is **0.4873~\AA**, which rounds to **0.49~\AA**, as Figure 6(h) shows.
Change the `Combined` PhosphoFill top-3 median in Table~\ref{tab:toolcomp}
from `0.48` to `0.49`, and use `0.49~[0.46, 0.51]~\AA` in the adjacent
bootstrap-CI sentence. The 2.7-fold comparison is unchanged. The current
TeX also calls the fourth method `na\"ive Kabsch`; the source metric is
PhosphoFill's residue-specific **Stage~0 internal-coordinate seed**, not an
independent Kabsch-only comparator. Use `PhosphoFill Stage~0 seed` in the
Figure 5 caption, Table~\ref{tab:toolcomp}, and the two adjacent comparator
paragraphs to match Figures 5 and 6. The 0.48~\AA\ in Table~\ref{tab:multisite}
for SEP with 3+ sites is a different subgroup value and should not change.

## 7. Methods: replace `Phosphosite environment annotation` paragraph

The current paragraph describes a different orientation vector and unsupported
water/metal dependence classes. Replace it with the implementation used for
Supplementary Tables S2--S4 and Figure 4:

```tex
\subsection{Phosphosite environment annotation}

Experimental PDB records were retained when the context audit completed successfully and the observed residue type matched the intended modified or unmodified type. For modified residues, the audit measured nearest selected Arg, Lys, and His side-chain atoms from the P and three phosphate oxygens; for unmodified parent residues it used the hydroxyl oxygen. The basic-atom set comprised Arg NE, NH1, NH2, and CZ; Lys NZ and CE; and His ND1, NE2, CE1, and CD2. Atoms beyond 8~\AA\ were not counted as measurable nearest-basic contacts. Modified sites were assigned hierarchically to salt-bridge-like (nearest selected basic atom $\leq3.5$~\AA), polar-contact-like (nearest selected polar sidechain atom $\leq3.5$~\AA), water-proximal (nearest water $\leq3.5$~\AA), or no-obvious-contact classes. These proximity classes do not establish interaction energies or hydrogen-bond geometries. When DSSP-derived relative solvent accessibility (RSA) was available, sites with RSA $\geq0.25$ were classified as exposed and the others as buried; otherwise the audit used local water and non-water heavy-atom counts to assign exposed, buried, or intermediate classes. The orientation angle was calculated between the anchor-to-P vector and the anchor-to-nearest-basic-atom vector; angles $\leq60^{\circ}$ were called toward, $60^{\circ}<\theta<120^{\circ}$ lateral, and $\geq120^{\circ}$ opposite. A separate same-side measure compared the side of the local sidechain or aromatic-ring plane occupied by P and by the nearest basic atom. Benchmark RMSD was joined to the audit using accession, PDB identifier, chain, residue position, and phosphoresidue type, excluding unmatched records without substitution from another structure.
```

## 8. Methods: torsion convention and statistical analysis

Append this sentence to `Internal-coordinate geometry derivation` to make
Figure 2's and the production seed angles reproducible:

```tex
The phosphate torsions in Figure~\ref{fig:geometry} and the PhosphoFill priors use the extractor's signed convention; relative to Bio.PDB's conventional dihedral, $\phi_{\mathrm{Bio.PDB}}=\mathrm{wrap}(180^{\circ}-\phi_{\mathrm{PhosphoFill}})$, where wrap maps angles to [$-180^{\circ}$, $+180^{\circ}$). Conventional $\chi_1$ rotamer labels in Supplementary Table~S3 were assigned after this conversion.
```

In the `Statistical analysis` Methods paragraph, replace the sentences
beginning `Spearman's rank correlation...`, `Scoring function
discrimination...`, `Mann--Whitney U-test compared...`, `Chi-squared tests...`,
and `Logistic regression...` with:

```tex
Spearman's rank correlation quantified the association between the retrospective S2--S1 RMSD difference and final top-1 placement RMSD. S1 is the lowest experimental RMSD among scanned orientations and therefore requires experimental phosphate coordinates. Two-sided Mann--Whitney $U$ tests compared measurable nearest-basic distances between pooled modified and unmodified PDB cohorts. Chi-squared tests with Cram\'{e}r's $V$ assessed four-category conventional $\chi_1$ rotamer distributions. For contact-stratified benchmark accuracy, the primary retrospective contrast compared salt-bridge-like sites with all other audit classes after exact structure-level matching. Median RMSD differences were assigned percentile 95\% confidence intervals from 5,000 bootstrap resamples of proteins, retaining all PDB records from each resampled accession together. Comparisons with no-obvious-contact alone and with no-obvious-contact plus water-proximal sites were sensitivity checks. Experimental phosphate coordinates define these contact classes, so they were not used as prospective predictors of placement success.
```

Delete the old AUC sentence rather than retaining it elsewhere: a diagnostic
derived from S1 cannot predict accuracy before the experimental phosphate
position is known, and no AUC result is reported in this section.

The current opening sentence of that Methods paragraph says `SciPy~1.11 and
scikit-learn~1.3`; that does not describe this corrected run. Replace it with:

```tex
Statistical analyses were performed in Python; exact package versions for each reproducibility run are recorded with the associated scripts and analysis outputs.
```

For this specific audit, the recorded versions are SciPy 1.18.1 and
scikit-learn 1.9.1 (`outputs/audit_summary.json`). The comparator-analysis
manifest also records SciPy 1.18.1. This wording avoids assigning one
unverified version pair to every figure in the study.

The Figure 2 caption can retain its empirical numerical modes, but append:

```tex
Angles are reported in the signed torsion convention used by the geometry extractor and PhosphoFill core.
```

## 9. Discussion: replace the paragraph beginning `The structural environment surrounding phosphorylation sites`

```tex
The experimental context audit shows that phosphate atoms often lie near selected basic-residue atoms and that phosphate direction relative to those atoms varies by residue type. These observations motivate the local-environment terms in PhosphoFill's scoring function. The unmodified-versus-phosphorylated distance comparison, however, uses different query-atom sets and pooled structural cohorts, while contact-stratified benchmark accuracy does not show a uniform salt-bridge-like advantage. The evidence therefore supports local contacts as one source of orientation information, rather than establishing that phosphorylation recruits basic residues or that any single contact type determines the experimental phosphate pose.
```

## 10. Claims removed rather than silently revised

- “Salt bridge formation is the dominant interaction determining phosphate
  orientation”: the audit shows frequent proximity, not energetic dominance.
- Non-basic **backbone-donor** hydrogen bonds: not measured by this audit.
- Water/metal **dependence** or a fundamental solvent-free accuracy ceiling:
  nearby waters are counted, but necessity and metal effects were not tested.
- Old SEP `0.46 vs 0.94 A, p<10^-6` contact/RMSD result: fails exact-key
  reproduction; the exact-key analysis is group-definition-sensitive.
- `+1.28` as the strongest validated logistic predictor: source script omitted
  two declared features and used a non-structure-specific join. The corrected
  exploratory four-feature coefficient is `+1.49` on 744 complete cases, but
  its collinearity and reference-derived features make it unsuitable as the
  manuscript's predictive claim.
- A phosphorylation-induced trans rotamer switch: the old label corresponds
  to conventional gauche-minus, and these are unpaired cohorts.

## 11. Optional additions: same-site anchor comparison and multi-site context

These analyses are new sensitivity checks, not corrections to the existing
Figure 4 or Figure 6. The review PDF is
`Benchmarkv2/output/pdf/structural_determinants_additions_review.pdf`. It can
become one new supplementary figure after the final S-number is assigned.

After the paragraph beginning `The phosphate group in modified structures
lay closer`, optionally add:

```tex
We also compared experimental structures of the same protein positions using the same bridging oxygen (OG1, OG, or OH) as the distance marker in both modified and unmodified residues. After summarising each position by the median across available structures in each state, the modified-minus-unmodified nearest-basic distance was $-0.54$~\AA\ for TPO/THR ($n=47$ positions), $-0.47$~\AA\ for SEP/SER ($n=56$), and $-0.14$~\AA\ for PTR/TYR ($n=56$; Supplementary Fig.~\ref{fig:environment_sensitivity}a). Protein-cluster bootstrap 95\% intervals excluded zero for TPO and SEP, but not PTR. The TPO interval included zero when analysis was restricted to positions with a selected basic atom within 8~\AA\ in both states. These same-marker differences were substantially smaller than the pooled Figure~\ref{fig:environment} contrasts; separate experimental structures at each protein position still cannot establish that phosphorylation itself recruited a basic residue.
```

After the paragraph beginning `We matched the experimental environment
audit`, optionally add:

```tex
Extending the exact structure-level match to the multi-site benchmark yielded 490 of 494 successful records (160 TPO, 57 SEP, and 273 PTR; Supplementary Fig.~\ref{fig:environment_sensitivity}b). For top-1 poses, salt-bridge-like sites had median RMSDs of 0.93, 0.83, and 0.96~\AA\ for TPO, SEP, and PTR, compared with 1.28, 1.46, and 1.56~\AA\ in the other audit classes. Protein-cluster bootstrap intervals for the salt-like-minus-other median difference excluded zero for TPO and SEP, but not PTR. No best-of-three difference had an interval excluding zero. As in the single-site analysis, the contact classes were assigned from experimental phosphate coordinates and are retrospective associations rather than prospective scoring predictors.
```

Optional new supplementary caption (assign its final S-number after deciding
the order of the other new supplementary figures):

```tex
\caption*{\textbf{Supplementary Figure Sx.} Same-site basic-atom proximity and multi-site contact sensitivity analyses. \textbf{(a)} Nearest selected Arg, Lys, or His side-chain atom distance recalculated from the same bridging oxygen in modified and unmodified experimental PDB structures. Each point is one protein position, represented by the median across available structures in each state. Negative modified-minus-unmodified values indicate a shorter distance in modified structures. Black diamonds denote the median paired difference and horizontal bars the 95\% protein-cluster bootstrap interval. The cohorts contain 47 TPO/THR, 56 SEP/SER, and 56 PTR/TYR shared positions; they are not before-and-after measurements in the same crystal structure. \textbf{(b)} Top-1 PhosphoFill phosphate RMSD for 490 multi-site benchmark records matched exactly to their experimental environment audit, stratified by residue type and distance-defined salt-bridge-like versus other contact classes. Sample sizes appear below each box. These reference-derived contact labels describe association, not interaction energy or causation.}
\label{fig:environment_sensitivity}
```

Append to `Phosphosite environment annotation` if the new supplement is used:

```tex
For the same-site sensitivity analysis, modified and unmodified records were linked by accession and protein position. Nearest selected basic-atom distance was recalculated from the bridging oxygen in both states using the original mmCIF coordinates and all selected basic atoms in the same chain, without an 8-\AA\ search cap. Each state contributed one median across its available PDB structures per protein position; differences were calculated within positions. A restricted analysis retained positions whose median distances were within 8~\AA\ in both states. The multi-site benchmark was linked to the experimental context audit using the same exact structure-level key as the single-site analysis, with four unmatched records excluded.
```

Append to `Statistical analysis` if the new supplement is used:

```tex
Median paired bridging-oxygen distance differences and multi-site salt-like-versus-other top-1 and best-of-three RMSD differences were assigned percentile 95\% confidence intervals from 5,000 bootstrap resamples of protein accessions, retaining all positions and structures from each sampled accession together. The paired-site and multi-site comparisons were exploratory sensitivity analyses.
```

If included, add this sentence after the Discussion paragraph in section 9:

```tex
Using the same bridging-oxygen distance marker at protein positions represented in both structural cohorts attenuated the pooled proximity contrast, and the multi-site contact analysis remained sensitive to residue type and whether top-1 or best-of-three placement was assessed.
```
