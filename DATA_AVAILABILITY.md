# Data availability and exclusions

## Included

- Experimental geometry and environment TSV/JSON files for modified and
  unmodified, all-site and single-site cohorts.
- Frozen TPO, SEP, and PTR strip-and-regraft definitions and result tables.
- Per-site PyTMs default, PyTMs optimised, and PTM-Psi comparator tables.
- Processed data, statistics, captions, table source values, and reference PDFs.
- The complete canonical three-rank proteome result table, compressed as
  `proteome_analysis/data/proteome_results_complete_canonical.tsv.gz`.
- Compact outputs for all four multi-site controls.
- The focused CDK2/ERK2 structures required for Figure 8.
- SHA-256 manifests for excluded general structure caches.

## Not included

1. **Per-protein AlphaFold/PhosphoFill production artifacts.** Coordinate
   structures, rich JSON/TSV reports, and HTML visualisations are distributed
   through Scop3P. The consolidated three-rank table required by the analysis
   scripts is included here.
2. **General benchmark PDB cache.** Missing structures are downloaded by the
   benchmark runner. The original cache inventory and checksums are retained in
   `benchmark_validation/inputs/manifests/`.
3. **Full structural archives for all multi-site controls.** Compact per-site
   results, summaries, verification, and analysis outputs are included. Fresh
   coordinate archives can be generated with
   `multisite_validation/scripts/run_multisite_controls.py`.

These exclusions affect storage, not the frozen numerical results included in
the repository.

## Licensing

Original PhosphoFill data, figures, tables, and documentation are provided
under CC BY 4.0. Analysis and reproduction code is provided under Apache 2.0.
Third-party structures and comparator materials retain their source terms and
are not relicensed. See `LICENSE_SCOPE.md`.
