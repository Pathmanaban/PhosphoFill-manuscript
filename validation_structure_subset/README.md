# Archived structural validation subset

This directory is the structural spot-check for the frozen PhosphoFill
benchmark. It preserves all 14 source benchmark rows, exact crystal inputs and
CCD templates, frozen benchmark scripts, run configuration, numerical results,
verification, and a PyMOL loader. To keep the repository compact, archived
stage/rank structures are included only for the six CDK2/ERK2 manuscript
case-study targets. The remaining eight archives can be regenerated with the
included runner.

The subset is intentionally small enough to inspect structure by structure while covering the important benchmark behaviors:

- all six manuscript case-study targets;
- TPO, SEP, and PTR;
- single-site and joint multi-site grafting;
- a SEP site where the empirical seed is already close to the crystal pose;
- a SEP site where the 12-angle scan rescues the initial orientation;
- records containing valid sites alongside sites without experimental phosphate coordinates.

## Case selection

| Group | PDB / target | Role |
|---|---|---|
| TPO | 2CCH A160 | CDK2 pThr160 case study |
| TPO | 4EOJ A160 | CDK2 pThr160 case study and PyMOL-panel source |
| TPO | 4IZA A185 | ERK2 pThr185 case study |
| TPO | 5V62 A185 | ERK2 pThr185 case study |
| PTR | 4IZA A187 | ERK2 pTyr187 case study |
| PTR | 5V62 A187 | ERK2 pTyr187 case study |
| TPO | 4BN1 A287,A288 | clean two-site joint-grafting control |
| SEP | 8R27 A261 | empirical-seed match control |
| SEP | 1UU3 A241 | 12-angle scan-rescue control |
| SEP | 5N6N C60,C83 | clean two-site joint-grafting control |
| PTR | 2J0L A576,A577 | clean two-site joint-grafting control |
| SEP | 1R0Z C422,C659,C660,C670 | archive regression; C659 and C670 lack reference phosphate coordinates |
| SEP | 1R0Z D659,D660,D670 | archive regression; D659 lacks reference phosphate coordinates |
| PTR | 2ZM3 A1161,A1165,A1166 | archive regression; A1161 lacks reference phosphate coordinates |

The 14 rows contain 24 requested target sites. The expected outcome is 20 `OK` sites and four explicitly retained `NO_COORDS` sites.

### Case-study residue-label correction

The general PTR benchmark rows for ERK2 4IZA/5V62 carry the mixed source label `T,Y` even though the requested target is residue 187. The subset explicitly assigns those target-only rows to `Y`, matching pTyr187, the experimental residue, and the case-study dataset. The corresponding TPO rows are explicitly assigned to `T` at residue 185. This prevents the primary target from inheriting the first label in an unrelated mixed-residue list.

## Archived stages

For every row the run saves:

1. `input/reference_phospho.cif` — the experimental phosphoprotein structure;
2. `input/stripped_input.cif` — the same structure with phosphate atoms removed at the requested sites;
3. `stage0/structure_s0.cif` — internal-coordinate seed with no torsion scan and no minimization;
4. `stage1/structure_s1.cif` — geometric oracle assembled from the best of the 12 scanned orientations against the experimental coordinates;
5. `stage2/structure_s2.cif` — environment-score-selected prescan orientation;
6. `stage3/phosphoFill1.cif`, `phosphoFill2.cif`, and `phosphoFill3.cif` — the three score-ranked candidates, each independently minimized with OpenMM;
7. rank-specific OpenMM reports, prescan frames, coordinate snapshots, and per-site metrics.

Stage 1 is an evaluation-only geometric ceiling and is not used to select production poses.

## Frozen empirical priors

| Residue | Anchor--P | Valence angle | Initial torsion seed |
|---|---:|---:|---:|
| TPO | OG1--P = 1.613 Å | CB--OG1--P = 118.8° | CG2--CB--OG1--P = -58.9° |
| SEP | OG--P = 1.614 Å | CB--OG--P = 115.4° | CA--CB--OG--P = +66.3° |
| PTR | OH--P = 1.609 Å | CZ--OH--P = 125.1° | neutral 0° starting orientation; no torsion-prior penalty |

The scan evaluates 12 orientations in 30° increments. The benchmark configuration uses prescan contact weight `0.0`, packing weight `0.2`, joint grafting within each benchmark row, and `openmm_local` minimization.

## Run commands

Run inside the same Python/OpenMM environment used for the complete benchmark. From a terminal in this directory:

```bash
python scripts/run_validation_subset.py
python scripts/verify_validation_subset.py
```

If a run is interrupted, resume from checkpointed cases:

```bash
python scripts/run_validation_subset.py --resume
```

Resume also retries any checkpointed benchmark row whose previous status was
`ERROR`; successful rows are retained and skipped.

To generate and inspect the selected input TSVs without running OpenMM:

```bash
python scripts/run_validation_subset.py --prepare-only
```

## Verification performed

`scripts/verify_validation_subset.py` checks that:

- all 14 selected rows and all six case-study cases are present;
- all 24 requested site rows are retained, including the four expected `NO_COORDS` records;
- phosphate atoms are absent from stripped inputs;
- Stage 0 reproduces the configured bond length, valence angle, and TPO/SEP torsion seed;
- each site has exactly 12 prescan frames at 30° increments;
- reference, stripped, Stage 0, Stage 1, Stage 2, and Stage 3 structures exist;
- ranks 1, 2, and 3 each have their own structure and OpenMM report with `relax_status=OK`;
- rank metadata explicitly marks all three structures as independently minimized;
- phosphate RMSD recomputed from each archived rank structure matches the result TSV;
- `top3_best_postmin_rmsd` and its rank are the minimum of the three recomputed values.

The detailed machine-readable result is written to `verification/spotcheck_report.json`. Exact commands, software versions, and SHA-256 hashes of the frozen inputs/scripts are written to the other files in `verification/`.

## Visual inspection

After verification, load `visual_checks/load_archived_subset.pml` in PyMOL. The loader creates objects for the reference, stripped input, Stage 0, Stage 1, Stage 2, and all three minimized ranks for every case, then disables them. Enable one case at a time, for example:

```pml
@visual_checks/load_archived_subset.pml
enable case_cdk2_4eoj_t160_*
zoom case_cdk2_4eoj_t160_*
```

Object colors are: reference gray, Stage 0 cyan, Stage 1 marine, Stage 2 orange, top-1 blue, top-2 violet, and top-3 salmon.

## Directory layout

```text
validation_structure_subset/
├── environment/requirements.txt
├── inputs/
│   ├── subset_cases.json
│   ├── source_benchmarks/
│   └── cache/
├── scripts/
│   ├── frozen/
│   ├── run_validation_subset.py
│   └── verify_validation_subset.py
├── outputs/
│   ├── generated_inputs/
│   ├── results/
│   ├── archives/
│   └── work/
├── verification/
└── visual_checks/
```
