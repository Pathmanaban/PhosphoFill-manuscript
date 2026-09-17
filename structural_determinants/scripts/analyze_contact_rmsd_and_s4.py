#!/usr/bin/env python3
"""Audit retrospective contact/RMSD contrasts and single/all-site S4 nesting.

Uses only the Python standard library. The main contrast is fixed as
salt_bridge_like versus every other context-audit class. Narrower contrasts
are sensitivity checks, not independent confirmatory tests. Protein-cluster
bootstrap intervals respect repeated structures from the same accession.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import statistics
from collections import defaultdict
from pathlib import Path


RESIDUES = ("TPO", "SEP", "PTR")
CLASSES = ("salt_bridge_like", "hbond_like", "water_mediated_like", "none_obvious")
CONTRASTS = {
    "primary_salt_vs_all_other": CLASSES,
    "sensitivity_salt_vs_none": ("salt_bridge_like", "none_obvious"),
    "sensitivity_salt_vs_none_or_water":
        ("salt_bridge_like", "none_obvious", "water_mediated_like"),
}
METRICS = ("rank1_postmin_rmsd", "top3_best_postmin_rmsd")
ORIENTATION = ("toward", "lateral", "opposite")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def percentile(values: list[float], q: float) -> float:
    ordered = sorted(values)
    if not ordered:
        raise ValueError("Cannot calculate a percentile of no values")
    index = (len(ordered) - 1) * q
    lower = math.floor(index)
    upper = math.ceil(index)
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (index - lower)


def group_rows(rows: list[dict], contrast: str) -> tuple[list[dict], list[dict]]:
    included = set(CONTRASTS[contrast])
    selected = [row for row in rows if row["interaction_class"] in included]
    salt = [row for row in selected if row["interaction_class"] == "salt_bridge_like"]
    other = [row for row in selected if row["interaction_class"] != "salt_bridge_like"]
    return salt, other


def summaries(salt: list[dict], other: list[dict], metric: str) -> dict:
    a = [row[metric] for row in salt]
    b = [row[metric] for row in other]
    if not a or not b:
        raise ValueError("Both comparison groups must be nonempty")
    a_median = statistics.median(a)
    b_median = statistics.median(b)
    return {
        "salt_n": len(a),
        "other_n": len(b),
        "salt_proteins_n": len({row["acc_id"] for row in salt}),
        "other_proteins_n": len({row["acc_id"] for row in other}),
        "salt_protein_positions_n":
            len({(row["acc_id"], row["position_key"]) for row in salt}),
        "other_protein_positions_n":
            len({(row["acc_id"], row["position_key"]) for row in other}),
        "salt_median_angstrom": a_median,
        "other_median_angstrom": b_median,
        "median_difference_salt_minus_other_angstrom": a_median - b_median,
        "salt_recovery_le_1A_pct": 100 * sum(x <= 1.0 for x in a) / len(a),
        "other_recovery_le_1A_pct": 100 * sum(x <= 1.0 for x in b) / len(b),
        "recovery_difference_salt_minus_other_pp":
            100 * (sum(x <= 1.0 for x in a) / len(a)
                   - sum(x <= 1.0 for x in b) / len(b)),
    }


def cluster_bootstrap(rows: list[dict], replicates: int, seed: int) -> dict:
    clusters = defaultdict(list)
    for row in rows:
        clusters[row["acc_id"]].append(row)
    proteins = sorted(clusters)
    rng = random.Random(seed)
    distributions = {contrast: {metric: {"median_difference": [],
                                          "recovery_difference": []}
                                 for metric in METRICS}
                     for contrast in CONTRASTS}
    skipped = {contrast: 0 for contrast in CONTRASTS}
    for _ in range(replicates):
        sample = []
        for _ in proteins:
            sample.extend(clusters[rng.choice(proteins)])
        for contrast in CONTRASTS:
            salt, other = group_rows(sample, contrast)
            if not salt or not other:
                skipped[contrast] += 1
                continue
            for metric in METRICS:
                a = [row[metric] for row in salt]
                b = [row[metric] for row in other]
                distributions[contrast][metric]["median_difference"].append(
                    statistics.median(a) - statistics.median(b))
                distributions[contrast][metric]["recovery_difference"].append(
                    100 * (sum(x <= 1.0 for x in a) / len(a)
                           - sum(x <= 1.0 for x in b) / len(b)))
    intervals = {}
    for contrast in CONTRASTS:
        intervals[contrast] = {"valid_replicates": replicates - skipped[contrast],
                               "skipped_replicates": skipped[contrast], "metrics": {}}
        for metric in METRICS:
            sample = distributions[contrast][metric]
            intervals[contrast]["metrics"][metric] = {
                "median_difference_95pct_cluster_bootstrap_ci_angstrom":
                    [percentile(sample["median_difference"], 0.025),
                     percentile(sample["median_difference"], 0.975)],
                "recovery_difference_95pct_cluster_bootstrap_ci_pp":
                    [percentile(sample["recovery_difference"], 0.025),
                     percentile(sample["recovery_difference"], 0.975)],
            }
    return intervals


def context_key(row: dict[str, str]) -> tuple[str, str, str, int, str]:
    return (row["acc_id"].strip().upper(), row["pdb_id"].strip().lower(),
            row["chain_id"].strip(), int(float(row["target_position"])),
            row["actual_resname"].strip().upper())


def usable_context(rows: list[dict[str, str]], residue: str) -> dict[tuple, dict]:
    result = {}
    for row in rows:
        if row["status"].strip().lower() != "ok":
            continue
        if row["actual_resname"].strip().upper() != residue:
            continue
        key = context_key(row)
        if key in result:
            raise ValueError(f"Duplicate context key: {key}")
        result[key] = row
    return result


def audit_s4(context_root: Path) -> tuple[dict, list[Path]]:
    output = {}
    sources = []
    for residue in RESIDUES:
        single_path = context_root / "mod_context_audit_single" / f"{residue}_context_audit.tsv"
        all_path = context_root / "mod_context_audit_all" / f"{residue}_context_audit.tsv"
        sources.extend((single_path, all_path))
        single = usable_context(read_tsv(single_path), residue)
        all_sites = usable_context(read_tsv(all_path), residue)
        missing = sorted(set(single) - set(all_sites))
        overlap = set(single) & set(all_sites)
        fields = ("orientation_toward_basic", "phosphate_to_basic_angle_deg",
                  "nearest_basic_dist", "same_side_as_nearest_basic")
        changed = {field: sum(single[key][field] != all_sites[key][field] for key in overlap)
                   for field in fields}
        additional = [all_sites[key] for key in set(all_sites) - set(single)]
        valid = [row for row in additional if row["orientation_toward_basic"] in ORIENTATION]
        target_counts = defaultdict(int)
        for row in additional:
            target_counts[row["n_target_sites"]] += 1
        output[residue] = {
            "single_n": len(single), "all_n": len(all_sites),
            "single_missing_from_all_n": len(missing),
            "single_missing_from_all_examples": missing[:5],
            "overlap_field_mismatch_n": changed,
            "additional_n": len(additional),
            "additional_n_target_sites_distribution": dict(sorted(target_counts.items())),
            "additional_measurable_orientation_n": len(valid),
            "additional_toward_n":
                sum(row["orientation_toward_basic"] == "toward" for row in valid),
            "additional_toward_pct":
                100 * sum(row["orientation_toward_basic"] == "toward" for row in valid)
                / len(valid) if valid else None,
        }
    return output, sources


def main() -> None:
    audit_dir = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--joined", type=Path,
                        default=audit_dir / "outputs" / "benchmark_joined_environment.tsv")
    parser.add_argument("--context-root", type=Path,
                        default=audit_dir.parent / "geometry")
    parser.add_argument("--out", type=Path,
                        default=audit_dir / "outputs" / "contact_rmsd_and_s4_new_analysis.json")
    parser.add_argument("--bootstrap-replicates", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260914)
    args = parser.parse_args()
    if args.bootstrap_replicates < 100:
        raise ValueError("Use at least 100 bootstrap replicates")

    raw = read_tsv(args.joined)
    if len(raw) != 775:
        raise ValueError(f"Expected 775 exact-joined benchmark rows; found {len(raw)}")
    rows = []
    for row in raw:
        if row["restype_3"] not in RESIDUES or row["interaction_class"] not in CLASSES:
            raise ValueError(f"Unexpected residue/contact class: {row}")
        parsed = dict(row)
        for metric in METRICS:
            parsed[metric] = float(row[metric])
            if not math.isfinite(parsed[metric]):
                raise ValueError(f"Non-finite RMSD for {metric}")
        rows.append(parsed)

    results = {}
    for residue_index, residue in enumerate(RESIDUES):
        subset = [row for row in rows if row["restype_3"] == residue]
        if not subset:
            raise ValueError(f"No benchmark records for {residue}")
        boot = cluster_bootstrap(subset, args.bootstrap_replicates,
                                 args.seed + residue_index)
        contrasts = {}
        for contrast in CONTRASTS:
            salt, other = group_rows(subset, contrast)
            contrasts[contrast] = {
                "class_definition": list(CONTRASTS[contrast]),
                "metrics": {metric: summaries(salt, other, metric) for metric in METRICS},
                "cluster_bootstrap": boot[contrast],
            }
        results[residue] = {
            "exact_joined_structure_records_n": len(subset),
            "unique_proteins_n": len({row["acc_id"] for row in subset}),
            "interaction_class_counts": {
                label: sum(row["interaction_class"] == label for row in subset)
                for label in CLASSES},
            "contrasts": contrasts,
        }

    s4, context_sources = audit_s4(args.context_root)
    report = {
        "question": "Retrospective association of experimental contact class with reconstructed phosphate RMSD",
        "unit": "PDB structure record; bootstrap resampling cluster is accession/protein",
        "primary_contrast": "salt_bridge_like versus hbond_like + water_mediated_like + none_obvious",
        "sensitivity_contrasts": [key for key in CONTRASTS if key != "primary_salt_vs_all_other"],
        "interpretation_limit": "Experimental phosphate coordinates define the contact class; this is not a prospective predictor or causal test.",
        "bootstrap": {"replicates": args.bootstrap_replicates, "seed": args.seed,
                      "interval": "2.5th and 97.5th percentiles of protein-cluster bootstrap contrasts"},
        "rmsd": results,
        "s4_single_all_nesting": s4,
        "inputs": [{"path": str(path.resolve()), "sha256": sha256(path)}
                   for path in [args.joined, *context_sources]],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Wrote {args.out}")
    for residue in RESIDUES:
        primary = results[residue]["contrasts"]["primary_salt_vs_all_other"]
        top1 = primary["metrics"]["rank1_postmin_rmsd"]
        interval = primary["cluster_bootstrap"]["metrics"]["rank1_postmin_rmsd"][
            "median_difference_95pct_cluster_bootstrap_ci_angstrom"]
        print(residue, "salt/other n", top1["salt_n"], top1["other_n"],
              "median difference", round(top1["median_difference_salt_minus_other_angstrom"], 3),
              "95% cluster CI", [round(x, 3) for x in interval],
              "S4 nested", s4[residue]["single_missing_from_all_n"] == 0)


if __name__ == "__main__":
    main()
