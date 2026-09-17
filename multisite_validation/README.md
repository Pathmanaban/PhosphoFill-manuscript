# Multi-site order and coupling validation

This directory contains a self-contained sensitivity analysis for the PhosphoFill multi-site benchmark. It compares scan order, minimization order, independent-site placement, and true joint OpenMM minimization without changing the completed reference benchmark.

## Important implementation finding

The frozen benchmark implementation used for the current manuscript results is:

1. sites sorted N-to-C;
2. context-aware phosphate scan on the accumulating grafted model;
3. all score-ranked candidates assembled for each output rank;
4. modified residues minimized one at a time in N-to-C order, retaining coordinate updates.

Therefore, its multi-site minimization is sequential, not one joint OpenMM call containing all phosphates. The `joint_min` control below implements that distinct joint call.

## Four controls

| Control | Grafting/prescan context | Minimization |
|---|---|---|
| `n_to_c` | accumulating model, N→C | sequential N→C; current benchmark reference |
| `c_to_n` | accumulating model, C→N | sequential C→N |
| `independent` | each site on the original fully stripped structure | each site minimized separately |
| `joint_min` | accumulating model, N→C | all modified sites in one OpenMM call per pose rank |

Every control produces three score-ranked poses. The comparison reports top-1 RMSD and the minimum experimental RMSD among the three independently minimized ranks.

## Dataset

The input selection rule is `n_target_sites > 1` applied to the frozen TPO, SEP, and PTR benchmark definitions.

| Primary benchmark group | Rows | Requested target sites |
|---|---:|---:|
| TPO | 80 | 169 |
| SEP | 32 | 76 |
| PTR | 136 | 274 |
| **Total** | **248** | **519** |

Mixed-context benchmark labels are retained exactly as in the source TSVs. Analysis uses the explicit `primary_type` field rather than inferring the benchmark group from the mixed `MODRES_1` string. The grafting code determines the actual residue type from the stripped structure itself.

## Recommended Ubuntu run

Activate the same Python/OpenMM environment used for the completed benchmark, then enter this directory:

```bash
cd /path/to/PhosphoFill-manuscript/multisite_validation
```

First confirm the generated inputs and commands without running OpenMM:

```bash
python scripts/run_multisite_controls.py --prepare-only
```

Run all four controls and all three residue groups:

```bash
python scripts/run_multisite_controls.py
```

If interrupted, retain successful checkpoints and retry `ERROR` cases:

```bash
python scripts/run_multisite_controls.py --resume
```

The runner executes the 12 mode/type jobs sequentially. A single control can be run separately:

```bash
python scripts/run_multisite_controls.py --modes n_to_c
python scripts/run_multisite_controls.py --modes c_to_n
python scripts/run_multisite_controls.py --modes independent
python scripts/run_multisite_controls.py --modes joint_min
```

Likewise, residue groups can be restricted using, for example, `--types TPO SEP`.

## Structural archive policy

The published compact repository retains the completed per-site result TSVs,
summaries, analyses, and verification reports, but omits the multi-gigabyte
coordinate archives and input structure cache. A fresh full run downloads or
rebuilds the missing cache and can archive reports, coordinates, prescan frames,
and metadata while omitting duplicated stage CIF files. Saving every structural
stage across 248 rows and four controls would require several gigabytes.

To archive all reference, stripped, S0, S1, S2, S3, and top-1/2/3 structures, start a fresh run with:

```bash
python scripts/run_multisite_controls.py --save-stage-structures
```

Do not add `--resume` when changing from metrics-only to full-stage archiving because completed checkpoint rows would be skipped. The preferred workflow is to finish the metrics-first comparison, identify sensitive cases, and rerun only those cases with structures if PyMOL inspection is needed.

## Verification and analysis

After all controls finish:

```bash
python scripts/verify_multisite_run.py
python scripts/analyze_multisite_controls.py
```

For the compact published repository, rebuild the manuscript figure directly
from the included frozen source-value table without the omitted structure
cache:

```bash
python scripts/plot_multisite_sensitivity.py
```

The verifier checks:

- all expected mode/type result files;
- 519 requested sites per control, allowing explicit case-level errors;
- recorded scan order, graft mode, and minimization coupling;
- independent rank-1, rank-2, and rank-3 RMSDs;
- exact agreement between `top3_best_postmin_rmsd` and the minimum of the three ranks.

The analysis generates:

- `outputs/analysis/multisite_per_site_comparison.tsv`;
- `outputs/analysis/multisite_mode_summary.tsv`;
- `outputs/analysis/multisite_paired_summary.tsv`;
- `outputs/analysis/multisite_failures.tsv`;
- `outputs/analysis/multisite_order_sensitive_sites.tsv`;
- `outputs/analysis/multisite_distance_sensitivity.tsv`;
- `outputs/analysis/multisite_within_case_environment.tsv`;
- `outputs/analysis/multisite_sensitivity_summary.json`;
- `outputs/plots/multi_site_sensitivity.pdf` and `.png`.

The manuscript figure shows the complete top-1 and best-of-three RMSD distributions across all four controls, paired site-level RMSD changes relative to the N-to-C reference, and order sensitivity stratified by nearest experimental phosphorus--phosphorus distance. Recovery percentages are retained in the TSV/JSON summaries rather than plotted because N-to-C, C-to-N, and independent-site recovery are identical and joint minimisation differs by only one site at the 1.0-A threshold. The phosphorus--phosphorus distance is used for the physical-proximity panel because sequence or Cα separation alone can substantially overstate direct phosphate coupling.

## Reported comparisons

For every site and control, the analysis includes:

- top-1 and best-of-3 RMSD;
- change relative to `n_to_c`;
- top-1 and best-of-3 recovery at 1.0 and 1.5 Å;
- failure or missing status;
- score-selected prescan angle;
- experimentally best rank among the three outputs;
- nearest sequence separation, Cα distance, and experimental phosphate distance to another target site when available.

An order-sensitive site is flagged if C→N versus N→C produces any of:

- an absolute top-1 or best-of-3 RMSD change of at least 0.5 Å;
- a recovery change at 1.0 or 1.5 Å;
- a change in the experimentally best output rank;
- a score-selected prescan-angle change of at least 30°.

A close cluster is defined as a nearest target-site sequence separation of at most five residues or a nearest target-site Cα distance of at most 10 Å. Both raw distances remain in the per-site table so these thresholds can be changed later without rerunning OpenMM.

## Frozen placement priors and benchmark settings

The controls retain the same benchmark parameters:

- TPO: 1.613 Å, 118.8°, −58.9° seed;
- SEP: 1.614 Å, 115.4°, +66.3° primary seed;
- PTR: 1.609 Å, 125.1°, neutral 0° start with no torsion penalty;
- 12 scan orientations at 30° increments;
- prescan contact weight 0.0 and packing weight 0.2;
- three independently minimized score-ranked outputs.

## Directory layout

```text
multi-site-validation/
├── environment/requirements.txt
├── inputs/
│   ├── source_benchmarks/
│   └── cache/
├── scripts/
│   ├── frozen/
│   ├── run_multisite_controls.py
│   ├── verify_multisite_run.py
│   └── analyze_multisite_controls.py
├── outputs/
│   ├── generated_inputs/
│   ├── runs/<control>/
│   ├── analysis/
│   └── plots/
└── verification/
```
