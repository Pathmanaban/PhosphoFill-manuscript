# Figure reproduction map

Run commands from the repository root unless a section says otherwise.

## Figure 2 — empirical phosphate geometry

Inputs:

- `geometry/mod_geometry_prior_all/*_geometry_prior.tsv`
- single-site and unmodified counterparts are retained beside the all-site data
  for Supplementary Table S1 and sensitivity checks.

Command:

```bash
python geometry/figure2_geometry_priors.py \
  --data-dir geometry/mod_geometry_prior_all \
  --outdir figures/reproduced
```

Add `--single` to generate the single-site-only diagnostic plot.

## Figure 3 — strip-and-regraft stages

Inputs, scripts, outputs, and verification are under
`benchmark_validation/`. The frozen analysis script is
`benchmark_validation/scripts/frozen/analyse_benchmark_postmin.py`; the
reference output is
`benchmark_validation/plots/fig3_validation_stages_postmin.pdf`.

To verify the frozen result first:

```bash
cd benchmark_validation
python scripts/verify_frozen_benchmark.py
```

Rebuild Figure 3 and the source tables:

```bash
python scripts/frozen/analyse_benchmark_postmin.py \
  --phosphofill-tsv \
    outputs/raw/TPO_postmin.tsv \
    outputs/raw/SEP_postmin.tsv \
    outputs/raw/PTR_postmin.tsv \
  --outdir outputs/reproduced
```

The same bundle contains the source values and LaTeX for manuscript Tables 2
and 3.

## Figure 4 — experimental structural environment

Command:

```bash
python structural_determinants/figure4_environment.py \
  --mod-audit \
    geometry/mod_context_audit_all/TPO_context_audit.tsv \
    geometry/mod_context_audit_all/SEP_context_audit.tsv \
    geometry/mod_context_audit_all/PTR_context_audit.tsv \
  --unmod-audit \
    geometry/unmod_context_audit_all/THR_context_audit.tsv \
    geometry/unmod_context_audit_all/SER_context_audit.tsv \
    geometry/unmod_context_audit_all/TYR_context_audit.tsv \
  --outdir figures/reproduced
```

The corrected audit, sensitivity analyses, Supplementary Tables S2--S4 source
values, and interpretation limits are documented in
`structural_determinants/README.md`.

## Figures 5 and 6 — comparator benchmark

The common 785-site data, PyTMs and PTM-Psi input tables, statistical outputs,
and plotting scripts are under `comparator_validation/`. Use the Figure 5 and
final corrected Figure 6 commands in
`comparator_validation/README.md` from that directory. Final Figure 6
uses the exact structure-level 775-site environment match in panel g.

The comparator calculations themselves do not need to be rerun to rebuild the
figures: the archived per-site RMSDs were coordinate-audited, and the complete
input tables are included.

## Figure 7 — proteome-scale analysis

Script:

```bash
python proteome_analysis/06_finalize_proteome_outputs.py \
  --results proteome_analysis/data/proteome_results_complete_canonical.tsv.gz \
  --outdir figures/reproduced_proteome
```

The included compressed table contains all 306,900 rows for 102,300 sites and
three independently minimised ranks. The manuscript reference PDF is retained
at `figures/main/figure7_proteome.pdf`.

## Figure 8 — CDK2/ERK2 case study

The numerical panel data, focused reference/comparator structures, manual PyMOL
panels, and final outputs are under `comparator_validation/`. The six corrected
PhosphoFill rank archives are retained under
`validation_structure_subset/outputs/archives/`.

```bash
cd comparator_validation
python scripts/rebuild_case_study_figure.py \
  --benchmark-root case_study_assets \
  --subset-root ../validation_structure_subset \
  --output-dir outputs/reproduced_case_study
python scripts/build_contact_recovery_supplement.py
```

## Supplementary figures

- S1 and S2: `proteome_analysis/08_build_supplementary_context_figures.py`;
  the included three-rank table and external DSSP/RSA annotation table are
  required.
- S3: `comparator_validation/scripts/build_contact_recovery_supplement.py`.
- S4: `multisite_validation/scripts/plot_multisite_sensitivity.py`
  rebuilds the figure directly from the included source-value TSV.
- S5: generated with `proteome_analysis/06_finalize_proteome_outputs.py`.
- S6: `proteome_analysis/07_compare_proteome_ranks.py`; use the included
  compressed three-rank table.
