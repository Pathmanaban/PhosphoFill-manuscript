# PhosphoFill manuscript data and figure reproduction

This repository contains the scripts, input tables, processed outputs, audit
records, and reference PDFs used for the PhosphoFill manuscript figures from
Figure 2 onward. It is a curated copy: the original working directories were
not modified.

## Scope

- Figures 2--6 and 8 are accompanied by their local source data, analysis
  scripts, processed outputs, and reference PDFs.
- The strip-and-regraft benchmark includes the frozen TPO, SEP, and PTR results,
  all three independently minimised ranks, comparator results, and the
  multi-site sensitivity analysis.
- Figure 8 includes the six focused CDK2/ERK2 structure sets needed for the
  case-study rebuild, but not the unrelated archived validation structures.
- Figure 7 and the proteome supplementary figures include their analysis
  scripts, final PDFs, and the compressed canonical three-rank result table.
  Per-protein AlphaFold/PhosphoFill structures and HTML/rich reports are
  intentionally excluded because they are distributed through Scop3P.
- The 2.3 GB multi-site coordinate archive and the 0.78 GB general PDB cache are
  intentionally excluded. Their compact result tables, checksums, manifests,
  verification reports, and regeneration scripts are included.

See [FIGURE_MAP.md](FIGURE_MAP.md) for the exact source and command for each
figure, and [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) for included and
external assets.

## Repository layout

```text
geometry/                       Figure 2 data and plotting script
benchmark_validation/           Figure 3 and Tables 2--3
structural_determinants/        Figure 4 and corrected environment audits
comparator_validation/          Figures 5, 6, 8 and Supplementary Figure S3
multisite_validation/           multi-site control analysis and Supplementary Figure S4
validation_structure_subset/    compact case-study validation subset
proteome_analysis/              Figure 7 and proteome supplementary scripts
figures/main/                   manuscript reference PDFs
figures/supplementary/          supplementary reference PDFs
```

## Environment

Python 3.10 or newer is recommended. For the plotting and tabular analyses:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

The structural benchmark additionally requires OpenMM. Exact versions used by
the frozen benchmark are in
`benchmark_validation/environment/requirements.txt`; comparator versions are
recorded in
`comparator_validation/outputs/analysis/software_versions.json`.

## Quick verification

The following commands verify frozen benchmark outputs without rerunning
OpenMM:

```bash
cd benchmark_validation
python scripts/verify_frozen_benchmark.py
cd ../multisite_validation
python scripts/verify_multisite_run.py
python scripts/plot_multisite_sensitivity.py
```

The plot-only multi-site command reads the included frozen source-value TSV.
Recalculating its physical-distance columns from structures requires the
regenerable input structure cache. Full benchmark and comparator commands are
documented in the README within each analysis directory.

## Proteome data boundary

The complete canonical three-rank PhosphoFill table is included as
`proteome_analysis/data/proteome_results_complete_canonical.tsv.gz`.
It contains 306,900 rows representing 102,300 sites and can be passed directly
to the analysis scripts as `--results`. The corresponding per-protein
coordinate models, rich reports, and HTML files are not duplicated here; those
storage-heavy artifacts will be available through Scop3P.

## Integrity

The frozen benchmark and structural-audit subdirectories retain their original
SHA-256 manifests and verification reports. Newly curated repository-level
files are documented through the source map above.
