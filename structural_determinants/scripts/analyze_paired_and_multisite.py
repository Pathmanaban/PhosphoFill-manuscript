#!/usr/bin/env python3
"""Two exploratory structural-environment controls from unmodified source data.

1. Match modified and unmodified records at the same accession/position and
   recalculate nearest selected basic-atom distance from the same anchor O.
2. Join successful multi-site benchmark rows to exact experimental structures
   and compare RMSD across audited contact classes.

Source files are read only. Output files are written below --outdir.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import statistics
from collections import OrderedDict, defaultdict
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser
from scipy import stats


PAIRS = (("TPO", "THR", "OG1"), ("SEP", "SER", "OG"), ("PTR", "TYR", "OH"))
BASIC_ATOMS = {
    "ARG": {"NE", "NH1", "NH2", "CZ"},
    "LYS": {"NZ", "CE"},
    "HIS": {"ND1", "NE2", "CE1", "CD2"},
}
CLASSES = {"salt_bridge_like", "hbond_like", "water_mediated_like", "none_obvious"}
METRICS = ("rank1_postmin_rmsd", "top3_best_postmin_rmsd")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def percentile(data: list[float], q: float) -> float:
    values = sorted(data)
    index = (len(values) - 1) * q
    lower, upper = math.floor(index), math.ceil(index)
    return values[lower] + (values[upper] - values[lower]) * (index - lower)


def key(row: dict, *, benchmark: bool = False) -> tuple[str, str, str, int, str]:
    residue = row["restype_3"] if benchmark else row["actual_resname"]
    position = row["position"] if benchmark else row["target_position"]
    chain = row["chain"] if benchmark else row["chain_id"]
    return (row["acc_id"].strip().upper(), row["pdb_id"].strip().lower(),
            chain.strip(), int(float(position)), residue.strip().upper())


def audit_rows(path: Path, residue: str) -> list[dict[str, str]]:
    return [row for row in read_tsv(path)
            if row["status"].strip().lower() == "ok"
            and row["actual_resname"].strip().upper() == residue]


def anchor_nearest_basic(row: dict, anchor: str, cache_dir: Path,
                         parser: MMCIFParser, parsed: OrderedDict[str, object],
                         seen_pdb: set[str]) -> float | None:
    pdb_id = row["pdb_id"].strip().lower()
    if pdb_id not in parsed:
        cif = cache_dir / f"{pdb_id}.cif"
        if not cif.exists():
            raise FileNotFoundError(cif)
        parsed[pdb_id] = next(parser.get_structure(pdb_id, str(cif)).get_models())
        seen_pdb.add(pdb_id)
        if len(parsed) > 3:
            parsed.popitem(last=False)
    else:
        parsed.move_to_end(pdb_id)
    model = parsed[pdb_id]
    chain_id = row["chain_id"].strip()
    if chain_id not in model:
        raise ValueError(f"Chain absent: {pdb_id} {chain_id}")
    chain = model[chain_id]
    position = int(float(row["target_position"]))
    target = next((res for res in chain.get_residues()
                   if res.id[1] == position and res.get_resname().upper()
                   == row["actual_resname"].strip().upper()), None)
    if target is None or anchor not in target:
        raise ValueError(f"Anchor residue absent: {pdb_id} {chain_id} {position}")
    query = np.asarray(target[anchor].get_coord(), dtype=float)
    best = math.inf
    for residue in chain.get_residues():
        if residue is target:
            continue
        names = BASIC_ATOMS.get(residue.get_resname().upper())
        if not names:
            continue
        for atom in residue.get_atoms():
            if atom.get_name().strip().upper() in names:
                distance = float(np.linalg.norm(query - atom.get_coord()))
                best = min(best, distance)
    return best if math.isfinite(best) else None


def paired_bootstrap(rows: list[dict], replicates: int, seed: int) -> list[float]:
    clusters = defaultdict(list)
    for row in rows:
        clusters[row["acc_id"]].append(row["difference_A"])
    proteins = sorted(clusters)
    rng = random.Random(seed)
    distribution = []
    for _ in range(replicates):
        sample = []
        for _ in proteins:
            sample.extend(clusters[rng.choice(proteins)])
        distribution.append(statistics.median(sample))
    return [percentile(distribution, 0.025), percentile(distribution, 0.975)]


def paired_analysis(context_root: Path, output: Path, replicates: int,
                    source_paths: set[Path]) -> dict:
    parser = MMCIFParser(QUIET=True)
    parsed = OrderedDict()
    seen_pdb = set()
    cache_dir = context_root / "pdb_cache"
    output_rows = []
    errors = []
    results = {}
    for modified, unmodified, anchor in PAIRS:
        path_m = context_root / "mod_context_audit_all" / f"{modified}_context_audit.tsv"
        path_u = context_root / "unmod_context_audit_all" / f"{unmodified}_context_audit.tsv"
        source_paths.update((path_m, path_u))
        mod_rows = audit_rows(path_m, modified)
        unmod_rows = audit_rows(path_u, unmodified)
        m_sites = defaultdict(list)
        u_sites = defaultdict(list)
        for row in mod_rows:
            m_sites[(row["acc_id"].strip().upper(), int(float(row["target_position"])))].append(row)
        for row in unmod_rows:
            u_sites[(row["acc_id"].strip().upper(), int(float(row["target_position"])))].append(row)
        shared = sorted(set(m_sites) & set(u_sites))
        n_complete = 0
        for acc_id, position in shared:
            distances = {}
            for state, records in (("modified", m_sites[(acc_id, position)]),
                                   ("unmodified", u_sites[(acc_id, position)])):
                measured = []
                for row in records:
                    try:
                        distance = anchor_nearest_basic(row, anchor, cache_dir, parser,
                                                        parsed, seen_pdb)
                        if distance is not None:
                            measured.append(distance)
                        source_paths.add(cache_dir / f"{row['pdb_id'].strip().lower()}.cif")
                    except (OSError, ValueError, KeyError) as exc:
                        errors.append({"comparison": f"{modified}/{unmodified}",
                                       "state": state, "acc_id": acc_id,
                                       "pdb_id": row["pdb_id"], "chain_id": row["chain_id"],
                                       "position": position, "error": str(exc)})
                distances[state] = measured
            if not distances["modified"] or not distances["unmodified"]:
                continue
            n_complete += 1
            m_med = statistics.median(distances["modified"])
            u_med = statistics.median(distances["unmodified"])
            output_rows.append({"comparison": f"{modified}/{unmodified}",
                                "acc_id": acc_id, "position": position,
                                "modified_structures_n": len(distances["modified"]),
                                "unmodified_structures_n": len(distances["unmodified"]),
                                "modified_anchor_basic_median_A": m_med,
                                "unmodified_anchor_basic_median_A": u_med,
                                "difference_A": m_med - u_med,
                                "both_state_medians_within8": m_med <= 8 and u_med <= 8})
        subset = [row for row in output_rows if row["comparison"] == f"{modified}/{unmodified}"]
        values = [row["difference_A"] for row in subset]
        if len(values) != n_complete or not values:
            raise ValueError(f"No complete paired observations for {modified}")
        wilcoxon = stats.wilcoxon(values, alternative="two-sided", zero_method="wilcox")
        restricted_rows = [row for row in subset if row["both_state_medians_within8"]]
        restricted = [row["difference_A"] for row in restricted_rows]
        results[modified] = {
            "unmodified": unmodified, "shared_protein_positions_n": len(shared),
            "complete_anchor_distance_pairs_n": len(values),
            "protein_clusters_n": len({row["acc_id"] for row in subset}),
            "modified_anchor_distance_median_A": statistics.median(
                row["modified_anchor_basic_median_A"] for row in subset),
            "unmodified_anchor_distance_median_A": statistics.median(
                row["unmodified_anchor_basic_median_A"] for row in subset),
            "paired_difference_median_A": statistics.median(values),
            "paired_difference_95pct_protein_cluster_bootstrap_ci_A": paired_bootstrap(
                subset, replicates, 20260914 + len(shared)),
            "fraction_modified_anchor_closer_pct": 100 * sum(x < 0 for x in values) / len(values),
            "wilcoxon_two_sided_p_exploratory": float(wilcoxon.pvalue),
            "both_state_medians_within8_n": len(restricted),
            "within8_sensitivity_difference_median_A": statistics.median(restricted) if restricted else None,
            "within8_sensitivity_95pct_protein_cluster_bootstrap_ci_A":
                paired_bootstrap(restricted_rows, replicates, 20260914 + len(restricted))
                if restricted else None,
        }
    write_tsv(output / "paired_anchor_basic_sites.tsv", output_rows,
              ["comparison", "acc_id", "position", "modified_structures_n",
               "unmodified_structures_n", "modified_anchor_basic_median_A",
               "unmodified_anchor_basic_median_A", "difference_A",
               "both_state_medians_within8"])
    write_tsv(output / "paired_anchor_coordinate_errors.tsv", errors,
              ["comparison", "state", "acc_id", "pdb_id", "chain_id", "position", "error"])
    return {"by_residue": results, "coordinate_errors_n": len(errors),
            "parsed_pdb_files_n": len(seen_pdb), "site_rows_written_n": len(output_rows)}


def bootstrap_contact(rows: list[dict], metric: str, replicates: int,
                      seed: int) -> list[float] | None:
    clusters = defaultdict(list)
    for row in rows:
        clusters[row["acc_id"]].append(row)
    proteins = sorted(clusters)
    rng = random.Random(seed)
    differences = []
    for _ in range(replicates):
        sample = []
        for _ in proteins:
            sample.extend(clusters[rng.choice(proteins)])
        salt = [row[metric] for row in sample if row["contact_group"] == "salt_like"]
        other = [row[metric] for row in sample if row["contact_group"] == "other"]
        if salt and other:
            differences.append(statistics.median(salt) - statistics.median(other))
    if len(differences) < replicates * 0.9:
        return None
    return [percentile(differences, 0.025), percentile(differences, 0.975)]


def multisite_analysis(context_root: Path, benchmark_root: Path,
                       output: Path, replicates: int, source_paths: set[Path]) -> dict:
    joined = []
    unmatched = []
    totals = {}
    for modified, _, _ in PAIRS:
        audit_path = context_root / "mod_context_audit_all" / f"{modified}_context_audit.tsv"
        benchmark_path = benchmark_root / f"{modified}_postmin.tsv"
        source_paths.update((audit_path, benchmark_path))
        audit_all = read_tsv(audit_path)
        audit_filtered = [row for row in audit_all
                          if row["status"].strip().lower() == "ok"
                          and row["actual_resname"].strip().upper() == modified]
        contexts = {key(row): row for row in audit_filtered}
        other_audit = defaultdict(list)
        for audit_row in audit_all:
            other_audit[(audit_row["acc_id"].strip().upper(),
                         audit_row["pdb_id"].strip().lower(),
                         audit_row["chain_id"].strip(),
                         int(float(audit_row["target_position"])))].append(audit_row)
        if len(contexts) != len(audit_filtered):
            raise ValueError(f"Duplicate context keys for {modified}")
        raw = read_tsv(benchmark_path)
        mapping = {"TPO": "T", "SEP": "S", "PTR": "Y"}
        selected = [row for row in raw if row["status"] == "OK"
                    and int(float(row["context_n_sites"])) >= 2
                    and row["restype"].strip().upper() in (modified, mapping[modified])]
        if len({key({**row, "restype_3": modified}, benchmark=True) for row in selected}) != len(selected):
            raise ValueError(f"Duplicate multi-site benchmark keys for {modified}")
        matched = 0
        for row in selected:
            row_key = key({**row, "restype_3": modified}, benchmark=True)
            context = contexts.get(row_key)
            if context is None:
                alternatives = other_audit.get(row_key[:4], [])
                unmatched.append({"restype_3": modified, "acc_id": row["acc_id"],
                                  "pdb_id": row["pdb_id"], "chain": row["chain"],
                                  "position": row["position"],
                                  "context_n_sites": row["context_n_sites"],
                                  "reason": "audit_wrong_residue_or_status" if alternatives
                                            else "audit_record_absent",
                                  "audit_status": ",".join(sorted({a["status"] for a in alternatives})),
                                  "audit_actual_resname": ",".join(sorted(
                                      {a["actual_resname"] for a in alternatives}))})
                continue
            cls = context["interaction_class"]
            if cls not in CLASSES:
                raise ValueError(f"Unknown interaction class {cls}")
            matched += 1
            item = {"restype_3": modified, "acc_id": row["acc_id"].strip().upper(),
                    "pdb_id": row["pdb_id"].strip().lower(),
                    "chain": row["chain"].strip(),
                    "position": int(float(row["position"])),
                    "context_n_sites": int(float(row["context_n_sites"])),
                    "interaction_class": cls,
                    "contact_group": "salt_like" if cls == "salt_bridge_like" else "other"}
            for metric in METRICS:
                item[metric] = float(row[metric])
            joined.append(item)
        totals[modified] = {"successful_primary_type_multi_n": len(selected),
                            "exact_context_matches_n": matched,
                            "unmatched_n": len(selected) - matched}
    write_tsv(output / "multisite_joined_environment.tsv", joined,
              ["restype_3", "acc_id", "pdb_id", "chain", "position", "context_n_sites",
               "interaction_class", "contact_group", *METRICS])
    write_tsv(output / "multisite_unmatched_environment.tsv", unmatched,
              ["restype_3", "acc_id", "pdb_id", "chain", "position",
               "context_n_sites", "reason", "audit_status", "audit_actual_resname"])
    summary = {}
    for modified, _, _ in PAIRS:
        subset = [row for row in joined if row["restype_3"] == modified]
        salt = [row for row in subset if row["contact_group"] == "salt_like"]
        other = [row for row in subset if row["contact_group"] == "other"]
        item = {**totals[modified], "salt_n": len(salt), "other_n": len(other),
                "proteins_n": len({row["acc_id"] for row in subset}),
                "structures_n": len({(row["pdb_id"], row["chain"]) for row in subset}),
                "interaction_class_counts": {
                    cls: sum(row["interaction_class"] == cls for row in subset)
                    for cls in sorted(CLASSES)}}
        for metric in METRICS:
            if salt and other:
                sa = [row[metric] for row in salt]
                ot = [row[metric] for row in other]
                item[metric] = {"salt_median_A": statistics.median(sa),
                                "other_median_A": statistics.median(ot),
                                "median_difference_salt_minus_other_A":
                                    statistics.median(sa) - statistics.median(ot),
                                "difference_95pct_protein_cluster_bootstrap_ci_A":
                                    bootstrap_contact(subset, metric, replicates,
                                                      20260914 + len(subset) + len(metric)),
                                "salt_recovery_le1A_pct": 100 * sum(x <= 1 for x in sa) / len(sa),
                                "other_recovery_le1A_pct": 100 * sum(x <= 1 for x in ot) / len(ot)}
        summary[modified] = item
    return {"by_residue": summary, "all_successful_multi_n": sum(
        item["successful_primary_type_multi_n"] for item in totals.values()),
            "all_exact_matches_n": len(joined), "unmatched_n": len(unmatched)}


def main() -> None:
    repository = Path(__file__).resolve().parents[2]
    default_benchmark = repository / "benchmark_validation" / "outputs" / "raw"
    default_context = repository / "geometry"
    default_output = Path(__file__).resolve().parents[1] / "outputs" / "new_additions"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--context-root", type=Path, default=default_context)
    parser.add_argument("--benchmark-root", type=Path, default=default_benchmark)
    parser.add_argument("--outdir", type=Path, default=default_output)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--analysis", choices=("both", "paired", "multisite"),
                        default="both")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    inputs = set()
    paired = (paired_analysis(args.context_root, args.outdir, args.bootstrap, inputs)
              if args.analysis in ("both", "paired") else None)
    multi = (multisite_analysis(args.context_root, args.benchmark_root,
                                args.outdir, args.bootstrap, inputs)
             if args.analysis in ("both", "multisite") else None)
    payload = {"definitions": {
        "paired_key": "accession and residue position; structure medians within state",
        "paired_query": "same OG1/OG/OH bridging oxygen in both states, nearest selected Arg/Lys/His atom in the same chain; no distance cap",
        "multi_join": "accession, PDB ID, chain, residue position, modified residue type",
        "multi_contrast": "experimental salt_bridge_like audit class versus all other classes",
        "bootstrap": f"{args.bootstrap} accession-cluster resamples, percentile 95% intervals"},
        "paired": paired, "multisite": multi,
        "source_sha256": {str(path): sha256(path) for path in sorted(inputs) if path.exists()}}
    summary_path = args.outdir / ("paired_and_multisite_summary.json" if args.analysis == "both"
                                  else f"{args.analysis}_summary.json")
    summary_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"paired": paired["by_residue"] if paired else None,
                      "multisite": multi["by_residue"] if multi else None,
                      "summary_file": str(summary_path)}, indent=2))


if __name__ == "__main__":
    main()
