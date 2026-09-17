#!/usr/bin/env python3
"""Reproduce the phosphosite-environment statistics with structure-level joins.

Source TSVs and JSON files are read without modification. Outputs are placed in
the directory given by --outdir. This script intentionally reports alternative
contact-group definitions instead of choosing a favourable comparison silently.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

try:
    from Bio import __version__ as BIOPYTHON_VERSION
except ImportError:
    BIOPYTHON_VERSION = None


MOD_TYPES = ("TPO", "SEP", "PTR")
UNMOD_TYPES = ("THR", "SER", "TYR")
PAIRS = (("THR", "TPO"), ("SER", "SEP"), ("TYR", "PTR"))
KEY = ["acc_key", "pdb_key", "chain_key", "position_key", "restype_3"]
ROTAMER_PAIRS = (("THR", "TPO", "torsion_N_CA_CB_OG1"),
                 ("SER", "SEP", "torsion_N_CA_CB_OG"))
INTERACTION_CLASSES = ("salt_bridge_like", "hbond_like",
                       "water_mediated_like", "none_obvious")
FEATURES = ("nearest_basic_dist", "n_basic_within4", "n_basic_within6",
            "nonwater_heavy_within6")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add_keys(frame: pd.DataFrame, *, audit: bool) -> pd.DataFrame:
    out = frame.copy()
    out["acc_key"] = out["acc_id"].astype(str).str.strip().str.upper()
    out["pdb_key"] = out["pdb_id"].astype(str).str.strip().str.lower()
    out["chain_key"] = out["chain_id" if audit else "chain"].astype(str).str.strip()
    pos_col = "target_position" if audit else "position"
    out["position_key"] = pd.to_numeric(out[pos_col], errors="raise").astype(int)
    return out


def assert_unique(frame: pd.DataFrame, label: str) -> None:
    duplicated = frame.duplicated(KEY, keep=False)
    if duplicated.any():
        sample = frame.loc[duplicated, KEY].head(10).to_dict("records")
        raise ValueError(f"{label} has {int(duplicated.sum())} duplicate exact keys: {sample}")


def load_audits(folder: Path, residues: tuple[str, ...], input_paths: list[Path]
                ) -> tuple[pd.DataFrame, pd.DataFrame]:
    kept = []
    qc = []
    for residue in residues:
        path = folder / f"{residue}_context_audit.tsv"
        input_paths.append(path)
        raw = pd.read_csv(path, sep="\t", low_memory=False)
        if not {"status", "actual_resname"}.issubset(raw.columns):
            raise ValueError(f"Audit lacks status/actual_resname: {path}")
        status_ok = raw["status"].astype(str).str.lower().eq("ok")
        actual_ok = raw["actual_resname"].astype(str).str.upper().eq(residue)
        qc.append({"subset": folder.name, "expected_residue": residue,
                   "raw_rows": len(raw), "status_ok": int(status_ok.sum()),
                   "status_ok_actual_match": int((status_ok & actual_ok).sum()),
                   "excluded_wrong_actual": int((status_ok & ~actual_ok).sum()),
                   "excluded_bad_status": int((~status_ok).sum())})
        frame = raw.loc[status_ok & actual_ok].copy()
        frame["restype_3"] = residue
        kept.append(add_keys(frame, audit=True))
    result = pd.concat(kept, ignore_index=True)
    assert_unique(result, folder.name)
    return result, pd.DataFrame(qc)


def load_benchmark(benchmark_dir: Path, input_paths: list[Path]) -> pd.DataFrame:
    frames = []
    rt_map = {"T": "TPO", "S": "SEP", "Y": "PTR"}
    for residue in MOD_TYPES:
        path = benchmark_dir / f"{residue}_postmin.tsv"
        input_paths.append(path)
        raw = pd.read_csv(path, sep="\t", low_memory=False)
        actual = raw["restype"].astype(str).str.upper().map(lambda x: rt_map.get(x, x))
        subset = raw.loc[raw["status"].eq("OK")
                         & pd.to_numeric(raw["context_n_sites"], errors="coerce").eq(1)
                         & actual.eq(residue)].copy()
        subset["restype_3"] = residue
        frames.append(add_keys(subset, audit=False))
    result = pd.concat(frames, ignore_index=True)
    assert_unique(result, "single-site benchmark")
    if len(result) != 790:
        raise ValueError(f"Expected 790 primary-type single-site results; found {len(result)}")
    return result


def rotamer_counts(path: Path, residue: str, torsion: str) -> tuple[int, list[int]]:
    with path.open(encoding="utf-8") as stream:
        payload = json.load(stream)
    detail = payload["by_resname"][residue]["torsions"][torsion]
    hist = detail["hist_30deg"]
    n = int(detail["n"])
    if sum(int(v) for v in hist.values()) != n:
        raise ValueError(f"Histogram sum != n in {path}")
    trans = sum(int(hist.get(k, 0)) for k in ("[-180,-150)", "[-150,-120)", "[-120,-90)"))
    gp = sum(int(hist.get(k, 0)) for k in ("[90,120)", "[120,150)", "[150,180)"))
    gm = sum(int(hist.get(k, 0)) for k in ("[-30,0)", "[0,30)"))
    return n, [trans, gp, gm, n - trans - gp - gm]


def check_raw_rotamers(path: Path, residue: str, torsion: str,
                       expected_n: int, expected_counts: list[int]) -> None:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    subset = raw.loc[raw["status"].eq("ok") & raw["resname"].eq(residue)]
    values = pd.to_numeric(subset[torsion], errors="coerce").dropna().to_numpy()
    trans = int(((values >= -180) & (values < -90)).sum())
    gp = int(((values >= 90) & (values <= 180)).sum())
    gm = int(((values >= -30) & (values < 30)).sum())
    observed = [trans, gp, gm, len(values) - trans - gp - gm]
    if len(values) != expected_n or observed != expected_counts:
        raise ValueError(f"Raw geometry disagrees with summary JSON: {path}: "
                         f"raw n={len(values)}, bins={observed}; "
                         f"JSON n={expected_n}, bins={expected_counts}")


def calc_rotamers(geometry_root: Path, input_paths: list[Path]) -> pd.DataFrame:
    rows = []
    for subset in ("single", "all"):
        for unmod, mod, torsion in ROTAMER_PAIRS:
            upath = geometry_root / f"unmod_geometry_prior_{subset}" / f"{unmod}_geometry_prior.summary.json"
            mpath = geometry_root / f"mod_geometry_prior_{subset}" / f"{mod}_geometry_prior.summary.json"
            input_paths.extend((upath, mpath))
            un, uc = rotamer_counts(upath, unmod, torsion)
            mn, mc = rotamer_counts(mpath, mod, torsion)
            uraw = upath.with_name(f"{unmod}_geometry_prior.tsv")
            mraw = mpath.with_name(f"{mod}_geometry_prior.tsv")
            input_paths.extend((uraw, mraw))
            check_raw_rotamers(uraw, unmod, torsion, un, uc)
            check_raw_rotamers(mraw, mod, torsion, mn, mc)
            test = stats.chi2_contingency(np.array([uc, mc]), correction=False)
            cramer_v = float(np.sqrt(test.statistic / (un + mn)))
            row = {"subset": subset, "comparison": f"{unmod}->{mod}",
                   "unmodified_n": un, "modified_n": mn,
                   "chi2": float(test.statistic), "df": int(test.dof),
                   "p_value": float(test.pvalue), "cramers_v": cramer_v}
            for label, count in zip(("trans", "gauche_plus", "gauche_minus", "other"), uc):
                row[f"unmodified_{label}_n"] = count
                row[f"unmodified_{label}_pct"] = 100 * count / un
            for label, count in zip(("trans", "gauche_plus", "gauche_minus", "other"), mc):
                row[f"modified_{label}_n"] = count
                row[f"modified_{label}_pct"] = 100 * count / mn
            rows.append(row)
    return pd.DataFrame(rows)


def standard_chi1_counts(path: Path, residue: str, torsion: str) -> tuple[int, list[int]]:
    """Convert the extractor convention to standard Bio.PDB chi1 angles.

    The extractor's dihedral is wrap(180 - Bio.PDB dihedral), verified against
    the N-CA-CB-OG1 coordinates of TPO 1E9H chain A residue 160.
    """
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    values = pd.to_numeric(raw.loc[raw["status"].eq("ok") & raw["resname"].eq(residue), torsion],
                           errors="coerce").dropna().to_numpy()
    standard = ((180.0 - values + 180.0) % 360.0) - 180.0
    trans = int((np.abs(standard) >= 150).sum())
    gp = int(((standard >= 30) & (standard < 90)).sum())
    gm = int(((standard >= -90) & (standard < -30)).sum())
    return len(values), [trans, gp, gm, len(values) - trans - gp - gm]


def calc_standard_rotamers(geometry_root: Path) -> pd.DataFrame:
    rows = []
    for subset in ("single", "all"):
        for unmod, mod, torsion in ROTAMER_PAIRS:
            upath = (geometry_root / f"unmod_geometry_prior_{subset}"
                     / f"{unmod}_geometry_prior.tsv")
            mpath = (geometry_root / f"mod_geometry_prior_{subset}"
                     / f"{mod}_geometry_prior.tsv")
            un, uc = standard_chi1_counts(upath, unmod, torsion)
            mn, mc = standard_chi1_counts(mpath, mod, torsion)
            test = stats.chi2_contingency(np.array([uc, mc]), correction=False)
            row = {"subset": subset, "comparison": f"{unmod}->{mod}",
                   "unmodified_n": un, "modified_n": mn,
                   "chi2": float(test.statistic), "df": int(test.dof),
                   "p_value": float(test.pvalue),
                   "cramers_v": float(np.sqrt(test.statistic / (un + mn)))}
            for label, count in zip(("trans", "gauche_plus", "gauche_minus", "other"), uc):
                row[f"unmodified_{label}_n"] = count
                row[f"unmodified_{label}_pct"] = 100 * count / un
            for label, count in zip(("trans", "gauche_plus", "gauche_minus", "other"), mc):
                row[f"modified_{label}_n"] = count
                row[f"modified_{label}_pct"] = 100 * count / mn
            rows.append(row)
    return pd.DataFrame(rows)


def mw_summary(left: pd.Series, right: pd.Series) -> dict:
    x = pd.to_numeric(left, errors="coerce").dropna()
    y = pd.to_numeric(right, errors="coerce").dropna()
    result = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
    return {"left_n": len(x), "right_n": len(y),
            "left_median": float(x.median()), "right_median": float(y.median()),
            "u": float(result.statistic), "p_value": float(result.pvalue)}


def calc_environment(mod_all: pd.DataFrame, unmod_all: pd.DataFrame
                     ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    table = []
    distances = []
    for residue in (*UNMOD_TYPES, *MOD_TYPES):
        source = unmod_all if residue in UNMOD_TYPES else mod_all
        sub = source.loc[source["restype_3"].eq(residue)]
        dists = pd.to_numeric(sub["nearest_basic_dist"], errors="coerce").dropna()
        row = {"residue": residue, "n": len(sub), "nearest_basic_measured_n": len(dists),
               "nearest_basic_median_angstrom": float(dists.median()),
               "nearest_basic_mean_angstrom": float(dists.mean())}
        for label in ("buried", "intermediate", "exposed"):
            count = int(sub["exposure_class"].eq(label).sum())
            row[f"{label}_n"] = count
            row[f"{label}_pct"] = 100 * count / len(sub)
        if residue in MOD_TYPES:
            for label in INTERACTION_CLASSES:
                count = int(sub["interaction_class"].eq(label).sum())
                row[f"{label}_n"] = count
                row[f"{label}_pct"] = 100 * count / len(sub)
            salt = sub.loc[sub["interaction_class"].eq("salt_bridge_like")]
            row["salt_like_nearest_carbon_n"] = int(salt["nearest_basic_atom"].astype(str).str.startswith("C").sum())
        table.append(row)
    for unmod, mod in PAIRS:
        u = unmod_all.loc[unmod_all["restype_3"].eq(unmod)]
        m = mod_all.loc[mod_all["restype_3"].eq(mod)]
        row = mw_summary(u["nearest_basic_dist"], m["nearest_basic_dist"])
        row["comparison"] = f"{unmod}->{mod}"
        distances.append(row)
    # Same protein-position is a sensitivity check, not a structure-matched pair.
    pairs = []
    for unmod, mod in PAIRS:
        u = unmod_all.loc[unmod_all["restype_3"].eq(unmod)]
        m = mod_all.loc[mod_all["restype_3"].eq(mod)]
        group = ["acc_key", "position_key"]
        umed = u.groupby(group)["nearest_basic_dist"].median().rename("unmodified_median")
        mmed = m.groupby(group)["nearest_basic_dist"].median().rename("modified_median")
        shared = pd.concat([umed, mmed], axis=1, join="inner").dropna()
        diff = shared["modified_median"] - shared["unmodified_median"]
        test = stats.wilcoxon(diff) if len(diff) and not np.allclose(diff, 0) else None
        pairs.append({"comparison": f"{unmod}->{mod}", "shared_protein_positions_n": len(shared),
                      "median_modified_minus_unmodified_angstrom": float(diff.median()) if len(diff) else np.nan,
                      "fraction_modified_closer": float((diff < 0).mean()) if len(diff) else np.nan,
                      "wilcoxon_p_value": float(test.pvalue) if test else np.nan})
    return pd.DataFrame(table), pd.DataFrame(distances), pd.DataFrame(pairs)


def calc_exposure_chi2(environment: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for unmod, mod in PAIRS:
        u = environment.loc[environment["residue"].eq(unmod)].iloc[0]
        m = environment.loc[environment["residue"].eq(mod)].iloc[0]
        table = np.array([[u[f"{label}_n"] for label in ("buried", "intermediate", "exposed")],
                          [m[f"{label}_n"] for label in ("buried", "intermediate", "exposed")]],
                         dtype=int)
        test = stats.chi2_contingency(table, correction=False)
        rows.append({"comparison": f"{unmod}->{mod}", "chi2": float(test.statistic),
                     "df": int(test.dof), "p_value": float(test.pvalue),
                     "cramers_v": float(np.sqrt(test.statistic / table.sum()))})
    return pd.DataFrame(rows)


def calc_unique_site_distances(mod_all: pd.DataFrame, unmod_all: pd.DataFrame) -> pd.DataFrame:
    """Sensitivity analysis after collapsing repeated PDB entries by protein-position."""
    rows = []
    for unmod, mod in PAIRS:
        u = unmod_all.loc[unmod_all["restype_3"].eq(unmod)]
        m = mod_all.loc[mod_all["restype_3"].eq(mod)]
        key = ["acc_key", "position_key"]
        ud = u.groupby(key)["nearest_basic_dist"].median().dropna()
        md = m.groupby(key)["nearest_basic_dist"].median().dropna()
        row = mw_summary(ud, md)
        row["comparison"] = f"{unmod}->{mod}"
        rows.append(row)
    return pd.DataFrame(rows)


def calc_orientation(mod_single: pd.DataFrame, mod_all: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for subset, source in (("single", mod_single), ("all", mod_all)):
        for residue in MOD_TYPES:
            sub = source.loc[source["restype_3"].eq(residue)]
            orientation = sub["orientation_toward_basic"].astype("string")
            valid = orientation.isin(["toward", "lateral", "opposite"])
            valid_n = int(valid.sum())
            row = {"subset": subset, "residue": residue, "all_sites_n": len(sub),
                   "measurable_basic_orientation_n": valid_n}
            for label in ("toward", "lateral", "opposite"):
                count = int(orientation.eq(label).sum())
                row[f"{label}_n"] = count
                row[f"{label}_pct"] = 100 * count / valid_n
            side = sub["same_side_as_nearest_basic"].astype("string")
            side_valid = side.isin(["same", "opposite"])
            row["same_side_valid_n"] = int(side_valid.sum())
            row["same_side_n"] = int(side.eq("same").sum())
            row["same_side_pct"] = 100 * row["same_side_n"] / row["same_side_valid_n"]
            rows.append(row)
    return pd.DataFrame(rows)


def exact_join(benchmark: pd.DataFrame, mod_all: pd.DataFrame
               ) -> tuple[pd.DataFrame, pd.DataFrame]:
    fields = KEY + ["actual_resname", "interaction_class", *FEATURES,
                    "nearest_basic_atom", "orientation_toward_basic"]
    merged = benchmark.merge(mod_all[fields], on=KEY, how="left", validate="one_to_one",
                             indicator=True)
    missing = merged.loc[merged["_merge"].eq("left_only"),
                         KEY + ["benchmark_case_key"]].copy()
    joined = merged.loc[merged["_merge"].eq("both")].drop(columns="_merge")
    return joined, missing


def explain_unmatched(missing: pd.DataFrame, audit_dir: Path) -> pd.DataFrame:
    raw_frames = []
    for expected in MOD_TYPES:
        raw = pd.read_csv(audit_dir / f"{expected}_context_audit.tsv", sep="\t", low_memory=False)
        raw["restype_3"] = expected
        raw_frames.append(add_keys(raw, audit=True))
    raw_audit = pd.concat(raw_frames, ignore_index=True)
    if raw_audit.duplicated(KEY).any():
        raise ValueError("Raw modified audit has duplicate exact keys")
    explained = missing.merge(raw_audit[KEY + ["status", "actual_resname"]],
                              on=KEY, how="left", validate="one_to_one")
    explained["unmatched_reason"] = np.where(
        explained["status"].isna(), "no_audit_record",
        np.where(explained["status"].ne("ok"), "audit_status_not_ok",
                 "actual_residue_does_not_match_expected"))
    return explained


def calc_contact_rmsd(joined: pd.DataFrame, cohort: str) -> pd.DataFrame:
    rows = []
    groups = {
        "none_obvious_only": ("none_obvious",),
        "none_or_water": ("none_obvious", "water_mediated_like"),
        "all_non_salt": ("none_obvious", "water_mediated_like", "hbond_like"),
    }
    for residue in MOD_TYPES:
        sub = joined.loc[joined["restype_3"].eq(residue)]
        salt = sub.loc[sub["interaction_class"].eq("salt_bridge_like"), "rank1_postmin_rmsd"]
        for name, labels in groups.items():
            other = sub.loc[sub["interaction_class"].isin(labels), "rank1_postmin_rmsd"]
            row = mw_summary(salt, other)
            row.update({"cohort": cohort, "residue": residue,
                        "left_group": "salt_bridge_like", "right_group": name,
                        "left_recovery_1A_pct": 100 * pd.to_numeric(salt).le(1.0).mean(),
                        "right_recovery_1A_pct": 100 * pd.to_numeric(other).le(1.0).mean()})
            rows.append(row)
    return pd.DataFrame(rows)


def calc_logistic(joined: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    cols = list(FEATURES)
    frame = joined[cols + ["rank1_postmin_rmsd"]].copy()
    for col in (*cols, "rank1_postmin_rmsd"):
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    frame = frame.dropna()
    x = StandardScaler().fit_transform(frame[cols].to_numpy())
    y = frame["rank1_postmin_rmsd"].le(1.0).astype(int).to_numpy()
    model = LogisticRegression(max_iter=1000, random_state=0)
    model.fit(x, y)
    coefficients = pd.DataFrame({"feature": cols, "standardized_log_odds_coefficient": model.coef_[0],
                                 "odds_ratio_per_1sd": np.exp(model.coef_[0])})
    correlation = frame[cols].corr().round(4).to_dict()
    info = {"complete_case_n": len(frame), "success_n": int(y.sum()),
            "success_threshold_angstrom": 1.0, "regularization": "sklearn LogisticRegression default L2, C=1",
            "feature_correlation": correlation}
    return coefficients, info


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    repository = Path(__file__).resolve().parents[2]
    parser.add_argument("--context-root", type=Path, default=repository / "geometry")
    parser.add_argument("--benchmark-dir", type=Path,
                        default=repository / "benchmark_validation" / "outputs" / "raw")
    parser.add_argument("--common-785", type=Path,
                        default=repository / "comparator_validation"
                        / "outputs" / "analysis" / "common_785_site_metrics.tsv")
    parser.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parents[1] / "outputs")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    input_paths: list[Path] = []
    source = args.context_root
    mod_all, qc_ma = load_audits(source / "mod_context_audit_all", MOD_TYPES, input_paths)
    mod_single, qc_ms = load_audits(source / "mod_context_audit_single", MOD_TYPES, input_paths)
    unmod_all, qc_ua = load_audits(source / "unmod_context_audit_all", UNMOD_TYPES, input_paths)
    unmod_single, qc_us = load_audits(source / "unmod_context_audit_single", UNMOD_TYPES, input_paths)
    audit_qc = pd.concat([qc_ma, qc_ms, qc_ua, qc_us], ignore_index=True)
    benchmark = load_benchmark(args.benchmark_dir, input_paths)
    legacy_rotamers = calc_rotamers(source, input_paths)
    rotamers = calc_standard_rotamers(source)
    environment, distances, shared_sites = calc_environment(mod_all, unmod_all)
    unique_site_distances = calc_unique_site_distances(mod_all, unmod_all)
    exposure_tests = calc_exposure_chi2(environment)
    orientation = calc_orientation(mod_single, mod_all)
    joined, missing = exact_join(benchmark, mod_all)
    missing = explain_unmatched(missing, source / "mod_context_audit_all")
    contact = calc_contact_rmsd(joined, "all_790")
    coefficients, logistic_info = calc_logistic(joined)

    input_paths.append(args.common_785)
    common = pd.read_csv(args.common_785, sep="\t", low_memory=False)
    if len(common) != 785:
        raise ValueError(f"Common comparator file has {len(common)} rows, not 785")
    common = add_keys(common, audit=False)
    assert_unique(common, "common comparator cohort")
    common_joined, common_missing = exact_join(common, mod_all)
    common_missing = explain_unmatched(common_missing, source / "mod_context_audit_all")
    contact = pd.concat([contact, calc_contact_rmsd(common_joined, "common_785")], ignore_index=True)

    outputs = {"audit_filter_qc.tsv": audit_qc,
               "table_s2_environment.tsv": environment,
               "nearest_basic_mannwhitney.tsv": distances,
               "unique_site_nearest_basic_sensitivity.tsv": unique_site_distances,
               "exposure_chi2.tsv": exposure_tests,
               "shared_protein_position_sensitivity.tsv": shared_sites,
               "table_s3_rotamers.tsv": rotamers,
               "legacy_rotamer_classification.tsv": legacy_rotamers,
               "table_s4_orientation.tsv": orientation,
               "benchmark_contact_rmsd.tsv": contact,
               "benchmark_joined_environment.tsv": joined,
               "benchmark_unmatched_environment.tsv": missing,
               "common_785_unmatched_environment.tsv": common_missing,
               "logistic_coefficients.tsv": coefficients}
    for name, frame in outputs.items():
        frame.to_csv(args.outdir / name, sep="\t", index=False)
    summary = {"modified_all_n": len(mod_all), "unmodified_all_n": len(unmod_all),
               "benchmark_single_site_n": len(benchmark), "benchmark_exact_join_n": len(joined),
               "benchmark_unmatched_n": len(missing), "common_785_exact_join_n": len(common_joined),
               "common_785_unmatched_n": len(common_missing),
               "logistic": logistic_info,
               "software": {"python": sys.version.split()[0],
                            "python_dependencies": {"numpy": np.__version__, "pandas": pd.__version__,
                                                    "scipy": scipy.__version__,
                                                    "scikit_learn": sklearn.__version__,
                                                    "matplotlib": metadata.version("matplotlib"),
                                                    "biopython": BIOPYTHON_VERSION}},
               "inputs": [{"path": str(path.resolve()), "sha256": sha256(path)}
                          for path in sorted(set(input_paths))]}
    (args.outdir / "audit_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("modified_all_n", "unmodified_all_n",
                                          "benchmark_single_site_n", "benchmark_exact_join_n",
                                          "benchmark_unmatched_n", "common_785_exact_join_n",
                                          "common_785_unmatched_n")}, indent=2))
    print(f"Outputs: {args.outdir}")


if __name__ == "__main__":
    main()
