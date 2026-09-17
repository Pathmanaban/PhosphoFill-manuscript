# Comparator validation and corrected benchmark figures

This directory completes benchmark checkpoints 5 and 6. It combines the
corrected PhosphoFill benchmark, in which ranks 1--3 were independently
minimised, with the archived PyTMs and PTM-Psi strip-and-regraft results.

## Rerun decision

The PyTMs and PTM-Psi calculations do **not** require rerunning. Their requested
790 site keys match the corrected PhosphoFill benchmark, their grafts are
independent of the PhosphoFill minimisation change, and all 2,361 archived
successful comparator outputs reproduce the recorded phosphate RMSDs to within
0.000052 A.

The common complete-case cohort contains 785 sites: 328 TPO, 241 SEP, and 216
PTR. PyTMs default succeeded for 789/790 requests. PyTMs optimised and PTM-Psi
each succeeded for 786/790 requests. The recorded failures are retained in
`outputs/comparator_audit.json`; they are tool/input mapping failures, not
silently discarded measurements.

## Comparator definitions

- PhosphoFill Stage 0 seed: residue-specific internal-coordinate construction
  using the empirical anchor--P distance, base--anchor--P angle, and initial
  torsion before conformational scanning or minimisation. Kabsch alignment is
  used only to transfer the terminal phosphate oxygens onto the constructed
  anchor--oxygen--phosphorus triad. The historical TSV field
  `stage0_kabsch_rmsd` is retained for compatibility but is not a naive Kabsch
  baseline.
- PyTMs default: PyTMs 1.2 `phosphorylate(selection=site)` with `optimize=0`.
  Coordinate audit shows its fixed output torsions at approximately -57.7
  degrees for TPO, -180.0 degrees for SEP, and +90.0 degrees for PTR.
- PyTMs optimised: PyTMs 1.2
  `phosphorylate(selection=site, optimize=1, interval=30, remove_radius=0)`.
  This is the package's built-in van der Waals optimisation with 30-degree
  angular sampling; the loop may stop early when its VDW strain reaches zero,
  so it should not be described as an unconditional 12-pose scan.
- PTM-Psi: PTM-Psi 0.1.0, commit
  `eb5422b6ff1ae00e8e214a5afe91147794865a8e`, using
  `Protein.modify(site, "phosphorylation")`. The package constructs the group
  by NERF internal coordinates with a 180-degree attachment torsion.

Every comparator starts from the corresponding experimental phosphorylated
structure, records the experimental phosphate coordinates, strips the
phosphate and restores the parent SER/THR/TYR residue, grafts with the selected
tool, and calculates permutation-invariant phosphate RMSD against the saved
experimental coordinates.

Software versions are frozen in `outputs/analysis/software_versions.json`.
The three comparator input TSVs and source scripts used here are copied into
`inputs/comparators/` and `scripts/frozen/`.

## Principal outputs

- `outputs/figures/fig5_tool_comparison.pdf`
- `outputs/figures/fig6_statistical_panels_contact_contrast.pdf`
- `outputs/case_study/figures/fig8_case_study_corrected.pdf`
- `outputs/analysis/table4_tool_comparison.tex`
- `outputs/analysis/figure5_caption.tex`
- `outputs/analysis/figure6_caption.tex`
- `outputs/analysis/figure8_caption.tex`
- `outputs/analysis/comparison_statistics.json`
- `outputs/comparator_audit.json`
- `outputs/comparator_coordinate_audit.tsv`

Figure 8 panels A--C were recalculated from the archived, independently
minimised PhosphoFill ranks. Panels D--G reuse the existing manually prepared
PyMOL images. The largest substantive correction is ERK2 pTyr187 in 4IZA:
the old best-of-three value of 0.26 A came from pre-correction rank handling;
the corrected independently minimised best-of-three RMSD is 0.94 A.

## Key frozen statistics

- At 1.0 A on the 785-site common cohort: PhosphoFill top-3 84.5%, top-1
  66.8%, PhosphoFill Stage 0 seed 49.4%, PyTMs optimised 30.2%, PyTMs default 6.6%, and
  PTM-Psi 5.5%.
- PhosphoFill top-1 versus PyTMs optimised: 619/785 lower RMSD; Wilcoxon
  p = 1.37e-56; McNemar cells both 188, PhosphoFill only 336, PyTMs only 49,
  neither 212; p = 4.00e-48.
- PhosphoFill top-1 versus PTM-Psi: 670/785 lower RMSD; Wilcoxon
  p = 6.88e-91; McNemar cells both 13, PhosphoFill only 511, PTM-Psi only 30,
  neither 231; p = 1.28e-94.
- Scoring gap versus final top-1 RMSD: Spearman rho = 0.853 on the common 785
  sites (p = 1.46e-223), or rho = 0.854 on all 790 PhosphoFill sites
  (p = 1.18e-225).

## Reproduction commands

Run from `PhosphoFill-manuscript/comparator_validation` after installing
the packages listed in `outputs/analysis/software_versions.json`.

```bash
python scripts/audit_comparator_outputs.py \
  --phosphofill-tsv ../benchmark_validation/outputs/raw/TPO_postmin.tsv ../benchmark_validation/outputs/raw/SEP_postmin.tsv ../benchmark_validation/outputs/raw/PTR_postmin.tsv \
  --pytms-default inputs/comparators/pytms_detail.tsv \
  --pytms-optimised inputs/comparators/pytms_detail_optimized.tsv \
  --ptmpsi inputs/comparators/ptmpsi_detail.tsv \
  --comparator-root /path/to/full/comparator_structure_archive \
  --output-json outputs/comparator_audit.json \
  --output-tsv outputs/comparator_coordinate_audit.tsv
```

```bash
python scripts/compute_comparison_statistics.py \
  --phosphofill-tsv ../benchmark_validation/outputs/raw/TPO_postmin.tsv ../benchmark_validation/outputs/raw/SEP_postmin.tsv ../benchmark_validation/outputs/raw/PTR_postmin.tsv \
  --pytms-default inputs/comparators/pytms_detail.tsv \
  --pytms-optimised inputs/comparators/pytms_detail_optimized.tsv \
  --ptmpsi inputs/comparators/ptmpsi_detail.tsv \
  --audit-tsv inputs/environment/TPO_context_audit.tsv inputs/environment/SEP_context_audit.tsv inputs/environment/PTR_context_audit.tsv \
  --outdir outputs/analysis
```

```bash
python scripts/frozen/fig5_tool_comparison_fair.py \
  --phosphofill-tsv ../benchmark_validation/outputs/raw/TPO_postmin.tsv ../benchmark_validation/outputs/raw/SEP_postmin.tsv ../benchmark_validation/outputs/raw/PTR_postmin.tsv \
  --pytms-default-detail inputs/comparators/pytms_detail.tsv \
  --pytms-detail inputs/comparators/pytms_detail_optimized.tsv \
  --ptmpsi-detail inputs/comparators/ptmpsi_detail.tsv \
  --outdir outputs/figures
```

```bash
python ../structural_determinants/scripts/fig6_statistical_panels_environment_corrected.py \
  --phosphofill-tsv ../benchmark_validation/outputs/raw/TPO_postmin.tsv ../benchmark_validation/outputs/raw/SEP_postmin.tsv ../benchmark_validation/outputs/raw/PTR_postmin.tsv \
  --pytms-default-detail inputs/comparators/pytms_detail.tsv \
  --pytms-detail inputs/comparators/pytms_detail_optimized.tsv \
  --ptmpsi-detail inputs/comparators/ptmpsi_detail.tsv \
  --audit-tsv inputs/environment/TPO_context_audit.tsv inputs/environment/SEP_context_audit.tsv inputs/environment/PTR_context_audit.tsv \
  --outdir outputs/figures \
  --environment-contrast salt_vs_all_other \
  --output-basename fig6_statistical_panels_contact_contrast
```

This final Figure 6 command uses 785 common comparator sites in panels a--f
and h, and 775 exact structure-level environment matches in panel g. The
frozen `fig6_statistical_panels_fair.py` and its 773-match historical output
are retained for provenance but are not the final manuscript Figure 6.

```bash
python scripts/rebuild_case_study_figure.py \
  --benchmark-root case_study_assets \
  --subset-root ../validation_structure_subset \
  --output-dir outputs/reproduced_case_study
```

The Figure 8 rebuild uses the six focused ranks in the sibling
`validation_structure_subset/` and the compact comparator/reference structures
in `case_study_assets/`. The first coordinate-audit command above is optional
and requires the full 785-site raw comparator structure archive; it is not
required to rebuild Figures 5, 6, or 8 from the included audited tables.

Rebuild Supplementary Figure S3 from the corrected case-study recovery tables:

```bash
python scripts/build_contact_recovery_supplement.py
```

The script reads the corrected environment- and basic-contact-recovery TSVs in
`outputs/case_study/` by default. It writes the PDF, PNG and TIFF heatmap,
long-form source values, a detailed contact table and a manuscript-ready
caption. Display labels use `PhosphoFill top-1` and `PyTMs optimised`;
best-of-3 refers to the independently minimised rank with the lowest phosphate
RMSD to the corresponding experimental coordinates.
