from __future__ import annotations

import argparse
import copy
import csv
import math
import os
import shutil
import tempfile
import urllib.request
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import numpy as np
from Bio.PDB import MMCIF2Dict, MMCIFParser, NeighborSearch, PDBIO, PDBParser, Select
from Bio.PDB.Atom import Atom
from Bio.PDB.mmcifio import MMCIFIO


PHOSPHO_MAP = {
    "SER": "SEP",
    "THR": "TPO",
    "TYR": "PTR",
}

# Residue-specific Kabsch anchor sets used only to align the parent residue frame.
# These must include the true bridging oxygen so the phosphate is seeded from the
# correct attachment atom rather than from an arbitrary all-heavy-atom fit.
KABSCH_ALIGN_ATOMS = {
    "SER": ["N", "CA", "CB", "OG"],
    "THR": ["N", "CA", "CB", "OG1"],
    "TYR": ["N", "CA", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
}

BRIDGING_OXYGEN = {
    "SER": "OG",
    "THR": "OG1",
    "TYR": "OH",
}

# === PDB-wide geometry priors (v5) ===
# Derived from experimental phospho structures at ≤2.5 Å resolution.
# TPO: n=276 single-site, SEP: n=228 single-site, PTR: n=208 single-site.

# TPO: unimodal torsion, pre-organized by methyl constraint
TPO_PRIOR_OG1_P = 1.613          # median from n=276 single-site
TPO_PRIOR_CB_OG1_P_DEG = 118.8   # median from n=276 single-site
TPO_PRIOR_CG2_CB_OG1_P_DEG = -58.9  # confirmed: 90% within ±30° at single-site

# SEP: bimodal torsion — sharp spike at +66.3° plus broad context-dependent basin
SEP_PRIOR_OG_P = 1.614            # median from n=228 single-site
SEP_PRIOR_CB_OG_P_DEG = 115.4     # median from n=228 single-site
SEP_PRIOR_SHARP_MODE_DEG = 66.3   # locked mode, σ≈3°, ~31% of sites (GMM fit)
SEP_PRIOR_SHARP_MODE_SIGMA = 5.0  # tight sigma for the sharp spike
SEP_PRIOR_BROAD_RANGE = (-60.0, 50.0)  # broad basin, ~69% of sites, no torsion penalty
# Legacy compatibility — used only by seeding logic in prescan_phospho_rotamer
# to generate starting angles. NOT used for torsion penalty scoring.
_SEP_LEGACY_SEED_ANGLES = [66.3, 12.6, -45.0]  # reordered: dominant first
SEP_PRIOR_CA_CB_OG_P_MODES_DEG = _SEP_LEGACY_SEED_ANGLES
SEP_PRIOR_CA_CB_OG_P_DEG = 66.3   # primary seed angle

# PTR: NO torsion prior (uniform distribution confirmed at n=501 and n=208)
# Only distance and angle priors are used.
PTR_PRIOR_OH_P = 1.609            # median from n=208 single-site
PTR_PRIOR_CZ_OH_P_DEG = 125.1    # median from n=208 single-site
PTR_PRIOR_CE1_CZ_OH_P_MODES_DEG = []   # empty: torsion is flat, no modes
PTR_PRIOR_CE2_CZ_OH_P_MODES_DEG = []   # empty: torsion is flat, no modes
PTR_PRIOR_DEFAULT_SEED_MODE_DEG = 0.0   # neutral seed, environment-driven

EXPECTED_PHOSPHO_BOND_MAX = 2.0  # Å, safety cutoff for seeded P-anchor distance

VDW_RADII = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
}

# Residue types whose sidechains are allowed to flex during minimization
# when they are near a modified residue.  Default: all standard amino acids.
# Any neighbor packed against the phosphosite may need to shift to
# accommodate the larger PO3 group, regardless of its charge or polarity.
FLEXIBLE_NEIGHBOR_RESTYPES_ALL = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
    "TYR", "VAL",
    # Also include the phospho-residues in case two sites are nearby
    "SEP", "TPO", "PTR",
}
BACKBONE_NAMES_RIGID = {"N", "CA", "C", "O", "OXT"}

# Optional CPU-friendly fallback.  Rotates entire sidechain (chi1/chi2).
RELAX_BONDS = {
    "SEP": [("CA", "CB", ["OG", "P", "O1P", "O2P", "O3P", "H1P", "H2P", "H3P"])],
    "TPO": [("CA", "CB", ["OG1", "CG2", "P", "O1P", "O2P", "O3P", "H1P", "H2P", "H3P"])],
    "PTR": [
        ("CA", "CB", ["CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH", "P", "O1P", "O2P", "O3P", "H1P", "H2P", "H3P"]),
        ("CB", "CG", ["CD1", "CD2", "CE1", "CE2", "CZ", "OH", "P", "O1P", "O2P", "O3P", "H1P", "H2P", "H3P"]),
    ],
}

# Prescan bonds — ONLY rotates the PO3 group around the bridging bond.
# The sidechain orientation (chi1/chi2 from AlphaFold) is preserved;
# only the phosphate orientation around the O–P axis is optimized.
#   SEP: P rotates around CB–OG    (serine hydroxyl direction preserved)
#   TPO: P rotates around CB–OG1   (threonine hydroxyl direction preserved)
#   PTR: P rotates around CZ–OH    (tyrosine hydroxyl direction preserved)
PRESCAN_BONDS = {
    "SEP": [("CB", "OG", ["P", "O1P", "O2P", "O3P", "H1P", "H2P", "H3P"])],
    "TPO": [("CB", "OG1", ["P", "O1P", "O2P", "O3P", "H1P", "H2P", "H3P"])],
    "PTR": [("CZ", "OH", ["P", "O1P", "O2P", "O3P", "H1P", "H2P", "H3P"])],
}

# Minimal bond graphs for the modified residue itself.
PHOSPHO_BONDS = {
    "SEP": [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"), ("CB", "OG"), ("OG", "P"), ("P", "O1P"), ("P", "O2P"), ("P", "O3P")],
    "TPO": [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"), ("CB", "OG1"), ("CB", "CG2"), ("OG1", "P"), ("P", "O1P"), ("P", "O2P"), ("P", "O3P")],
    "PTR": [
        ("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"), ("CB", "CG"),
        ("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"),
        ("CZ", "OH"), ("OH", "P"), ("P", "O1P"), ("P", "O2P"), ("P", "O3P"),
    ],
}

TEMPLATE_URL = "https://files.rcsb.org/ligands/download/{code}.cif"
PLDDT_BINS = [
    (90.0, "very_high"),
    (70.0, "confident"),
    (50.0, "low"),
    (-float("inf"), "very_low"),
]


@dataclass(frozen=True)
class SiteSpec:
    chain_id: Optional[str]
    resseq: int
    icode: str = " "


@dataclass
class TemplateAtom:
    name: str
    element: str
    coord: np.ndarray


@dataclass
class TemplateGeometry:
    """CCD-derived ideal bond lengths (Å) and angles (radians) for a
    modified residue.  Used as equilibrium values in OpenMM harmonic
    forces instead of measuring from the (potentially distorted)
    post-graft coordinates."""
    bond_lengths: Dict[Tuple[str, str], float]   # (atomA, atomB) -> Å
    bond_angles: Dict[Tuple[str, str, str], float]  # (a, center, c) -> radians


class GraftingError(Exception):
    pass


class AltlocPruneSelect(Select):
    def __init__(self, residue_altloc_map: Dict[Tuple, Optional[str]]):
        super().__init__()
        self.residue_altloc_map = residue_altloc_map

    def accept_atom(self, atom):
        if not atom.is_disordered():
            return True
        residue = atom.get_parent()
        preferred = self.residue_altloc_map.get(residue.full_id, None)
        altloc = (atom.get_altloc() or " ").strip()
        if not altloc:
            return True
        if preferred is None:
            return altloc == "A"
        return altloc == preferred


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Graft phosphate groups onto SER/THR/TYR residues in a PDB or mmCIF structure, "
            "with optional OpenMM local minimization, protonation hints, pLDDT QC, altloc handling, "
            "and multi-model support."
        )
    )
    parser.add_argument("input_structure", help="Input PDB or mmCIF file")
    parser.add_argument("output_structure", help="Output PDB or mmCIF file")
    parser.add_argument("--sites", nargs="+", required=True, help="Sites to phosphorylate. Examples: 206 350 A:450 A:451B")
    parser.add_argument("--report", default=None, help="Optional TSV report file. Default: <output>.report.tsv")
    parser.add_argument("--cache-dir", default=None, help="Optional directory for cached CCD template CIF files")
    parser.add_argument("--clash-threshold", type=float, default=2.0, help="Distance threshold (Å) used for simple close-contact counts")
    parser.add_argument("--model-indices", default="0", help="Comma-separated model indices to modify, or 'all'. Default: 0")
    parser.add_argument("--altloc-mode", choices=["highest_occupancy", "A", "first"], default="highest_occupancy")
    parser.add_argument("--prune-altlocs", action="store_true", help="When saving, keep only the selected altloc conformer")
    parser.add_argument("--plddt-source", choices=["bfactor", "none"], default="bfactor")
    parser.add_argument("--plddt-warn-below", type=float, default=70.0, help="Flag sites below this pLDDT mean as low-confidence")
    parser.add_argument("--protonation", choices=["0", "1", "2"], default="0", help="How many phosphate protons to add as explicit hydrogen atoms")
    parser.add_argument("--relax-mode", choices=["none", "torsion_scan", "openmm_local"], default="openmm_local")
    parser.add_argument(
        "--site-order",
        choices=["n_to_c", "c_to_n", "input"],
        default="n_to_c",
        help="Order used for context-aware grafting/prescan and sequential minimization.",
    )
    parser.add_argument(
        "--minimization-coupling",
        choices=["sequential", "joint"],
        default="sequential",
        help="Minimize modified sites one at a time or together in one OpenMM call.",
    )
    parser.add_argument("--torsion-step-deg", type=float, default=30.0)
    parser.add_argument("--prescan-rotamer", action="store_true", default=True, help="Pre-scan phospho rotamer before OpenMM minimization (default: on)")
    parser.add_argument("--no-prescan-rotamer", action="store_false", dest="prescan_rotamer", help="Disable phospho rotamer pre-scan")
    parser.add_argument("--prescan-step-deg", type=float, default=30.0, help="Step size (degrees) for rotamer pre-scan")
    parser.add_argument("--prescan-contact-weight", type=float, default=0.5, help="Weight for generic contact-preservation bonus in prescan scoring (0=pure clash, higher=prefer contacts)")
    parser.add_argument("--prescan-basic-weight", type=float, default=2.0, help="Extra weight for phosphate contacts to basic residues (Lys/Arg/His)")
    parser.add_argument("--prescan-hbond-weight", type=float, default=0.5, help="Extra weight for phosphate H-bond-like contacts to donor atoms")
    parser.add_argument("--prescan-polar-clash-scale", type=float, default=0.5, help="Scale factor for clash penalties involving phosphate oxygens and polar/basic donor atoms")
    parser.add_argument("--prescan-pack-weight", type=float, default=0.0, help="Weight for local protein packing bonus around phosphate oxygens")
    parser.add_argument("--prescan-geom-prior-weight", type=float, default=2.0, help="Weight for local empirical geometry prior during prescan scoring (currently strongest for TPO)")
    parser.add_argument("--prescan-tpo-torsion-sigma", type=float, default=20.0, help="Sigma in degrees for the TPO branch-relative torsion prior during prescan scoring")
    parser.add_argument("--prescan-keep-kabsch-if-good", action="store_true", help="Preserve the initial 0-degree graft if it is already clash-free and chemically supported")
    parser.add_argument("--flex-neighbor-radius", type=float, default=6.0, help="Radius (Å) to find flexible neighbor sidechains. 0 = disabled.")
    parser.add_argument("--restraint-k-neighbor", type=float, default=100.0, help="Positional restraint for flexible neighbor sidechain atoms")
    parser.add_argument("--openmm-platform", default="auto", help="OpenMM platform: auto, CPU, CUDA, OpenCL, HIP, Reference")
    parser.add_argument("--minimize-max-iterations", type=int, default=200)
    parser.add_argument("--minimize-tolerance", type=float, default=10.0, help="OpenMM minimization tolerance in kJ/mol/nm")
    parser.add_argument("--repulsion-scale", type=float, default=0.60, help="Scale factor applied to sum of VDW radii in the local overlap penalty")
    parser.add_argument("--repulsion-k", type=float, default=5000.0, help="OpenMM local overlap penalty constant in kJ/mol/nm^2")
    parser.add_argument("--restraint-k-rigid", type=float, default=5000.0, help="Positional restraint for atoms outside modified residues")
    parser.add_argument("--restraint-k-backbone", type=float, default=1000.0, help="Positional restraint for modified residue backbone atoms")
    parser.add_argument("--restraint-k-sidechain", type=float, default=50.0, help="Positional restraint for modified residue sidechain/template atoms")
    parser.add_argument("--bond-k", type=float, default=20000.0, help="Harmonic bond constant for modified residue geometry")
    parser.add_argument("--angle-k", type=float, default=200.0, help="Harmonic angle constant for modified residue geometry")
    parser.add_argument("--n-poses", type=int, default=1, choices=[1, 2, 3],
                        help="Number of independently minimized ranked poses to write (default: 1)")
    return parser.parse_args()


def parse_structure(path: str):
    suffix = Path(path).suffix.lower()
    parser = MMCIFParser(QUIET=True) if suffix in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    return parser.get_structure("input", path)


def save_structure(structure, output_path: str, residue_altloc_map: Optional[Dict[Tuple, Optional[str]]] = None, prune_altlocs: bool = False) -> None:
    suffix = Path(output_path).suffix.lower()
    io = MMCIFIO() if suffix in {".cif", ".mmcif"} else PDBIO()
    io.set_structure(structure)
    if prune_altlocs and residue_altloc_map is not None:
        io.save(output_path, select=AltlocPruneSelect(residue_altloc_map))
    else:
        io.save(output_path)


def _pose_output_path(base_path: str, rank: int) -> str:
    """Return the output path for a ranked pose."""
    path = Path(base_path)
    return str(path.parent / f"{path.stem}_phosphoFill_{rank}{path.suffix}")


def parse_site_token(token: str) -> SiteSpec:
    token = token.strip()
    chain_id: Optional[str] = None
    site_part = token
    if ":" in token:
        chain_id, site_part = token.split(":", 1)
        chain_id = chain_id.strip() or None
    if not site_part:
        raise ValueError(f"Invalid site token: {token!r}")
    i = len(site_part)
    while i > 0 and site_part[i - 1].isalpha():
        i -= 1
    num_part = site_part[:i]
    icode = site_part[i:] if i < len(site_part) else " "
    if not num_part or not num_part.lstrip("-").isdigit():
        raise ValueError(f"Invalid residue number in site token: {token!r}")
    return SiteSpec(chain_id=chain_id, resseq=int(num_part), icode=icode or " ")


def _ensure_list(value) -> List[str]:
    if isinstance(value, list):
        return value
    return [value]


def fetch_template_cif(code: str, cache_dir: Optional[str] = None) -> Tuple[str, Optional[str]]:
    code = code.upper()
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
        local_path = os.path.join(cache_dir, f"{code}.cif")
        if not os.path.exists(local_path):
            urllib.request.urlretrieve(TEMPLATE_URL.format(code=code), local_path)
        return local_path, None
    tmp_dir = tempfile.mkdtemp(prefix=f"{code}_ccd_")
    local_path = os.path.join(tmp_dir, f"{code}.cif")
    urllib.request.urlretrieve(TEMPLATE_URL.format(code=code), local_path)
    return local_path, tmp_dir


def load_template_atoms(code: str, cache_dir: Optional[str] = None, include_hydrogens: bool = False) -> Dict[str, TemplateAtom]:
    cif_path, tmp_dir = fetch_template_cif(code, cache_dir=cache_dir)
    try:
        d = MMCIF2Dict.MMCIF2Dict(cif_path)
        atom_ids = _ensure_list(d["_chem_comp_atom.atom_id"])
        elements = _ensure_list(d["_chem_comp_atom.type_symbol"])
        x_key_candidates = ["_chem_comp_atom.pdbx_model_Cartn_x_ideal", "_chem_comp_atom.model_Cartn_x", "_chem_comp_atom.model_cartn_x_ideal", "_chem_comp_atom.model_cartn_x"]
        y_key_candidates = ["_chem_comp_atom.pdbx_model_Cartn_y_ideal", "_chem_comp_atom.model_Cartn_y", "_chem_comp_atom.model_cartn_y_ideal", "_chem_comp_atom.model_cartn_y"]
        z_key_candidates = ["_chem_comp_atom.pdbx_model_Cartn_z_ideal", "_chem_comp_atom.model_Cartn_z", "_chem_comp_atom.model_cartn_z_ideal", "_chem_comp_atom.model_cartn_z"]

        def pick_key(candidates: Sequence[str]) -> str:
            for key in candidates:
                if key in d:
                    return key
            raise KeyError(f"No coordinate key found in template {code}; checked {candidates}")

        x_vals = _ensure_list(d[pick_key(x_key_candidates)])
        y_vals = _ensure_list(d[pick_key(y_key_candidates)])
        z_vals = _ensure_list(d[pick_key(z_key_candidates)])
        atoms: Dict[str, TemplateAtom] = {}
        leaving_flags = _ensure_list(d.get("_chem_comp_atom.pdbx_leaving_atom_flag", ["N"] * len(atom_ids)))

        for name, element, x, y, z, leaving in zip(atom_ids, elements, x_vals, y_vals, z_vals, leaving_flags):
            if not include_hydrogens and (str(element).upper() == "H" or str(name).upper().startswith("H")):
                continue
            if str(leaving).upper() == "Y":
                continue
            atoms[str(name).strip()] = TemplateAtom(str(name).strip(), str(element).strip().upper(), np.array([float(x), float(y), float(z)], dtype=float))
        return atoms
    finally:
        if tmp_dir is not None:
            shutil.rmtree(tmp_dir, ignore_errors=True)


def load_template_geometry(
    code: str,
    template_atoms: Dict[str, TemplateAtom],
    cache_dir: Optional[str] = None,
) -> TemplateGeometry:
    """Extract ideal bond lengths from the CCD ``_chem_comp_bond`` table
    and compute ideal angles from the ideal coordinates in *template_atoms*.

    Returns a :class:`TemplateGeometry` with reference values that should
    be used as equilibrium parameters in OpenMM harmonic forces, instead
    of measuring from the (potentially distorted) post-graft coordinates.
    """
    bond_lengths: Dict[Tuple[str, str], float] = {}
    bond_angles: Dict[Tuple[str, str, str], float] = {}

    # --- Bond lengths from CCD _chem_comp_bond table ---
    cif_path, tmp_dir = fetch_template_cif(code, cache_dir=cache_dir)
    try:
        d = MMCIF2Dict.MMCIF2Dict(cif_path)
        if "_chem_comp_bond.atom_id_1" in d and "_chem_comp_bond.atom_id_2" in d:
            ids1 = _ensure_list(d["_chem_comp_bond.atom_id_1"])
            ids2 = _ensure_list(d["_chem_comp_bond.atom_id_2"])
            # Try value_dist first (ideal), fall back to value_dist_nucleus
            dist_key = None
            for candidate in ["_chem_comp_bond.value_dist", "_chem_comp_bond.value_dist_nucleus"]:
                if candidate in d:
                    dist_key = candidate
                    break
            if dist_key is not None:
                dists = _ensure_list(d[dist_key])
                for a1, a2, dist_str in zip(ids1, ids2, dists):
                    a1 = str(a1).strip()
                    a2 = str(a2).strip()
                    try:
                        dist_val = float(dist_str)
                    except (ValueError, TypeError):
                        continue
                    # Skip hydrogen bonds
                    if a1.startswith("H") or a2.startswith("H"):
                        continue
                    # Store both orderings for easy lookup
                    bond_lengths[(a1, a2)] = dist_val
                    bond_lengths[(a2, a1)] = dist_val
    finally:
        if tmp_dir is not None:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    # --- Ideal angles from ideal coordinates ---
    # Build bond graph from the PHOSPHO_BONDS table (which lists the
    # bonds we actually restrain in OpenMM) and compute angles from the
    # CCD ideal coordinates.
    phospho_bonds = PHOSPHO_BONDS.get(code.upper(), [])
    nbrs: Dict[str, Set[str]] = defaultdict(set)
    for a, b in phospho_bonds:
        nbrs[a].add(b)
        nbrs[b].add(a)

    for center, neighbors in nbrs.items():
        if center not in template_atoms:
            continue
        coord_center = template_atoms[center].coord
        neigh_list = sorted(neighbors)
        for i in range(len(neigh_list)):
            for j in range(i + 1, len(neigh_list)):
                a_name, c_name = neigh_list[i], neigh_list[j]
                if a_name not in template_atoms or c_name not in template_atoms:
                    continue
                va = template_atoms[a_name].coord - coord_center
                vc = template_atoms[c_name].coord - coord_center
                norm_product = np.linalg.norm(va) * np.linalg.norm(vc)
                if norm_product < 1e-8:
                    continue
                cos_angle = float(np.clip(np.dot(va, vc) / norm_product, -1.0, 1.0))
                theta = math.acos(cos_angle)
                bond_angles[(a_name, center, c_name)] = theta
                bond_angles[(c_name, center, a_name)] = theta  # symmetric

    return TemplateGeometry(bond_lengths=bond_lengths, bond_angles=bond_angles)


def parse_model_indices(structure, text: str) -> List[int]:
    n_models = len(list(structure.get_models()))
    if text.strip().lower() == "all":
        return list(range(n_models))
    out = []
    for chunk in text.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        idx = int(chunk)
        if idx < 0 or idx >= n_models:
            raise ValueError(f"Model index {idx} out of range; structure has {n_models} model(s)")
        out.append(idx)
    if not out:
        raise ValueError("No valid model indices parsed")
    return sorted(set(out))


def unique_chain_ids(model) -> List[str]:
    return [chain.id for chain in model]


def find_residue(model, site: SiteSpec):
    if site.chain_id is not None:
        candidate_chains = [site.chain_id]
    else:
        chains = unique_chain_ids(model)
        if len(chains) != 1:
            raise GraftingError(f"Site {site.resseq}{site.icode.strip() or ''} has no chain specified, but structure has multiple chains: {chains}")
        candidate_chains = chains
    for chain_id in candidate_chains:
        if chain_id not in model:
            continue
        chain = model[chain_id]
        residue_id = (" ", site.resseq, site.icode)
        if residue_id in chain:
            return chain[residue_id]
    raise GraftingError(f"Residue not found for site {site.chain_id or '?'}:{site.resseq}{site.icode.strip() or ''}")


def get_residue_altloc_label(residue, altloc_mode: str) -> Optional[str]:
    altloc_scores: Dict[str, float] = {}
    seen: List[str] = []
    for atom in residue.get_unpacked_list():
        altloc = (atom.get_altloc() or " ").strip()
        if not altloc:
            continue
        if altloc not in altloc_scores:
            altloc_scores[altloc] = 0.0
            seen.append(altloc)
        altloc_scores[altloc] += atom.get_occupancy() or 0.0
    if not altloc_scores:
        return None
    if altloc_mode == "A" and "A" in altloc_scores:
        return "A"
    if altloc_mode == "highest_occupancy":
        return max(altloc_scores.items(), key=lambda kv: (kv[1], kv[0] == "A", kv[0]))[0]
    return seen[0]


def residue_atom_dict(residue, altloc_mode: str) -> Dict[str, Atom]:
    selected_altloc = get_residue_altloc_label(residue, altloc_mode)
    atoms: Dict[str, Atom] = {}
    best_rank: Dict[str, Tuple] = {}
    for atom in residue.get_unpacked_list():
        element = (atom.element or "").strip().upper()
        name = atom.get_name().strip()
        if element == "H" or name.startswith("H"):
            continue
        altloc = (atom.get_altloc() or " ").strip()
        occupancy = atom.get_occupancy() or 0.0
        if selected_altloc is not None and altloc not in {"", selected_altloc}:
            continue
        if altloc_mode == "highest_occupancy":
            rank = (altloc == selected_altloc, altloc == "", occupancy)
        elif altloc_mode == "A":
            rank = (altloc == "A", altloc == "", occupancy)
        else:
            rank = (True, altloc == "", -len(best_rank))
        if name not in atoms or rank > best_rank[name]:
            atoms[name] = atom
            best_rank[name] = rank
    return atoms


def kabsch_fit(mobile_xyz: np.ndarray, target_xyz: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if mobile_xyz.shape != target_xyz.shape:
        raise ValueError("mobile_xyz and target_xyz must have the same shape")
    if mobile_xyz.shape[0] < 3:
        raise ValueError("Need at least 3 points for a stable rigid-body fit")
    mobile_centroid = mobile_xyz.mean(axis=0)
    target_centroid = target_xyz.mean(axis=0)
    mobile_centered = mobile_xyz - mobile_centroid
    target_centered = target_xyz - target_centroid
    h = mobile_centered.T @ target_centered
    u, _, vt = np.linalg.svd(h)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T
    t = target_centroid - (r @ mobile_centroid)
    return r, t


def transform_points(points: np.ndarray, rotation: np.ndarray, translation: np.ndarray) -> np.ndarray:
    return (rotation @ points.T).T + translation


def angle_from_points(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> Optional[float]:
    ba = a - b
    bc = c - b
    nba = np.linalg.norm(ba)
    nbc = np.linalg.norm(bc)
    if nba < 1e-8 or nbc < 1e-8:
        return None
    cosang = float(np.dot(ba, bc) / (nba * nbc))
    cosang = max(-1.0, min(1.0, cosang))
    return float(np.arccos(cosang))


def dihedral_deg(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> Optional[float]:
    b0 = np.asarray(a, dtype=float) - np.asarray(b, dtype=float)
    b1 = np.asarray(c, dtype=float) - np.asarray(b, dtype=float)
    b2 = np.asarray(d, dtype=float) - np.asarray(c, dtype=float)
    n1 = np.cross(b0, b1)
    n2 = np.cross(b1, b2)
    n1n = np.linalg.norm(n1)
    n2n = np.linalg.norm(n2)
    b1n = np.linalg.norm(b1)
    if n1n < 1e-8 or n2n < 1e-8 or b1n < 1e-8:
        return None
    n1 /= n1n
    n2 /= n2n
    b1u = b1 / b1n
    m1 = np.cross(n1, b1u)
    x = float(np.dot(n1, n2))
    y = float(np.dot(m1, n2))
    return float(np.degrees(np.arctan2(y, x)))


def place_atom_from_internal(a: np.ndarray, b: np.ndarray, c: np.ndarray, length: float, angle_deg: float, dihedral_deg_value: float) -> np.ndarray:
    """Place atom d given atoms a-b-c, with distance(c,d)=length, angle(b,c,d)=angle_deg, and dihedral(a,b,c,d)=dihedral_deg_value."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)
    bc = b - c
    bc_n = np.linalg.norm(bc)
    if bc_n < 1e-8:
        raise GraftingError('Cannot place atom from internal coordinates: B and C are coincident')
    e1 = bc / bc_n
    ab = a - b
    n = np.cross(ab, e1)
    nn = np.linalg.norm(n)
    if nn < 1e-8:
        raise GraftingError('Cannot place atom from internal coordinates: anchor atoms are collinear')
    e3 = n / nn
    e2 = np.cross(e3, e1)
    theta = math.radians(angle_deg)
    phi = math.radians(dihedral_deg_value)
    direction = math.cos(theta) * e1 + math.sin(theta) * (math.cos(phi) * e2 + math.sin(phi) * e3)
    return c + length * direction


def build_sep_seed_coords(template_atoms: Dict[str, TemplateAtom], residue_atoms: Dict[str, Atom], template_only_names: List[str], mode_deg: Optional[float] = None) -> Tuple[Dict[str, np.ndarray], Dict[str, float]]:
    required_template = ["CB", "OG", "P"]
    required_target = ["CB", "OG", "CA"]
    for name in required_template:
        if name not in template_atoms:
            raise GraftingError(f'SEP template missing required atom {name} for local-frame seed')
    for name in required_target:
        if name not in residue_atoms:
            raise GraftingError(f'Target SER missing required atom {name} for SEP local-frame seed')

    cb_t = np.asarray(template_atoms["CB"].coord, dtype=float)
    og_t = np.asarray(template_atoms["OG"].coord, dtype=float)
    p_t = np.asarray(template_atoms["P"].coord, dtype=float)

    cb = np.asarray(residue_atoms["CB"].coord, dtype=float)
    og = np.asarray(residue_atoms["OG"].coord, dtype=float)
    ca = np.asarray(residue_atoms["CA"].coord, dtype=float)

    chosen_mode_deg = SEP_PRIOR_CA_CB_OG_P_DEG if mode_deg is None else float(mode_deg)
    placed_p = place_atom_from_internal(
        a=ca, b=cb, c=og,
        length=SEP_PRIOR_OG_P,
        angle_deg=SEP_PRIOR_CB_OG_P_DEG,
        dihedral_deg_value=chosen_mode_deg,
    )

    template_triad = np.array([cb_t, og_t, p_t], dtype=float)
    target_triad = np.array([cb, og, placed_p], dtype=float)
    rotation, translation = kabsch_fit(template_triad, target_triad)

    transformed = {}
    for name in template_only_names:
        if name not in template_atoms:
            continue
        coord = np.asarray(template_atoms[name].coord, dtype=float)
        transformed[name] = transform_points(np.array([coord], dtype=float), rotation, translation)[0]
    transformed["P"] = placed_p

    metrics = {
        "seed_mode": "sep_local_internal_prior",
        "seed_sep_mode_deg": chosen_mode_deg,
        "seed_target_OG_P": float(np.linalg.norm(placed_p - og)),
        "seed_target_CB_OG_P_deg": float(np.degrees(angle_from_points(cb, og, placed_p))) if angle_from_points(cb, og, placed_p) is not None else float('nan'),
    }
    dih = dihedral_deg(ca, cb, og, placed_p)
    if dih is not None:
        metrics["seed_target_CA_CB_OG_P_deg"] = float(dih)
    return transformed, metrics


def build_tpo_seed_coords(template_atoms: Dict[str, TemplateAtom], residue_atoms: Dict[str, Atom], template_only_names: List[str]) -> Tuple[Dict[str, np.ndarray], Dict[str, float]]:
    required_template = ["CB", "OG1", "P"]
    required_target = ["CB", "OG1", "CG2", "CA"]
    for name in required_template:
        if name not in template_atoms:
            raise GraftingError(f'TPO template missing required atom {name} for local-frame seed')
    for name in required_target:
        if name not in residue_atoms:
            raise GraftingError(f'Target THR missing required atom {name} for TPO local-frame seed')

    cb_t = np.asarray(template_atoms["CB"].coord, dtype=float)
    og1_t = np.asarray(template_atoms["OG1"].coord, dtype=float)
    p_t = np.asarray(template_atoms["P"].coord, dtype=float)

    cb = np.asarray(residue_atoms["CB"].coord, dtype=float)
    og1 = np.asarray(residue_atoms["OG1"].coord, dtype=float)
    cg2 = np.asarray(residue_atoms["CG2"].coord, dtype=float)
    ca = np.asarray(residue_atoms["CA"].coord, dtype=float)

    placed_p = place_atom_from_internal(
        a=cg2, b=cb, c=og1,
        length=TPO_PRIOR_OG1_P,
        angle_deg=TPO_PRIOR_CB_OG1_P_DEG,
        dihedral_deg_value=TPO_PRIOR_CG2_CB_OG1_P_DEG,
    )

    template_triad = np.array([cb_t, og1_t, p_t], dtype=float)
    target_triad = np.array([cb, og1, placed_p], dtype=float)
    rotation, translation = kabsch_fit(template_triad, target_triad)

    transformed = {}
    for name in template_only_names:
        if name not in template_atoms:
            continue
        coord = np.asarray(template_atoms[name].coord, dtype=float)
        transformed[name] = transform_points(np.array([coord], dtype=float), rotation, translation)[0]
    # enforce exact P from internal coordinates
    transformed["P"] = placed_p

    metrics = {
        "seed_mode": "tpo_local_internal_prior",
        "seed_target_OG1_P": float(np.linalg.norm(placed_p - og1)),
        "seed_target_CB_OG1_P_deg": float(np.degrees(angle_from_points(cb, og1, placed_p))) if angle_from_points(cb, og1, placed_p) is not None else float('nan'),
    }
    dih1 = dihedral_deg(ca, cb, og1, placed_p)
    dih2 = dihedral_deg(cg2, cb, og1, placed_p)
    if dih1 is not None:
        metrics["seed_target_CA_CB_OG1_P_deg"] = float(dih1)
    if dih2 is not None:
        metrics["seed_target_CG2_CB_OG1_P_deg"] = float(dih2)
    return transformed, metrics


def build_ptr_seed_coords(template_atoms: Dict[str, TemplateAtom], residue_atoms: Dict[str, Atom], template_only_names: List[str], mode_deg: Optional[float] = None) -> Tuple[Dict[str, np.ndarray], Dict[str, float]]:
    required_template = ["CZ", "OH", "P"]
    required_target = ["CE1", "CZ", "OH"]
    for name in required_template:
        if name not in template_atoms:
            raise GraftingError(f'PTR template missing required atom {name} for local-frame seed')
    for name in required_target:
        if name not in residue_atoms:
            raise GraftingError(f'Target TYR missing required atom {name} for PTR local-frame seed')

    cz_t = np.asarray(template_atoms["CZ"].coord, dtype=float)
    oh_t = np.asarray(template_atoms["OH"].coord, dtype=float)
    p_t = np.asarray(template_atoms["P"].coord, dtype=float)

    ce1 = np.asarray(residue_atoms["CE1"].coord, dtype=float)
    cz = np.asarray(residue_atoms["CZ"].coord, dtype=float)
    oh = np.asarray(residue_atoms["OH"].coord, dtype=float)

    chosen_mode_deg = PTR_PRIOR_DEFAULT_SEED_MODE_DEG if mode_deg is None else float(mode_deg)
    placed_p = place_atom_from_internal(
        a=ce1, b=cz, c=oh,
        length=PTR_PRIOR_OH_P,
        angle_deg=PTR_PRIOR_CZ_OH_P_DEG,
        dihedral_deg_value=chosen_mode_deg,
    )

    template_triad = np.array([cz_t, oh_t, p_t], dtype=float)
    target_triad = np.array([cz, oh, placed_p], dtype=float)
    rotation, translation = kabsch_fit(template_triad, target_triad)

    transformed = {}
    for name in template_only_names:
        if name not in template_atoms:
            continue
        coord = np.asarray(template_atoms[name].coord, dtype=float)
        transformed[name] = transform_points(np.array([coord], dtype=float), rotation, translation)[0]
    transformed["P"] = placed_p

    metrics = {
        "seed_mode": "ptr_local_internal_prior",
        "seed_ptr_mode_deg": chosen_mode_deg,
        "seed_target_OH_P": float(np.linalg.norm(placed_p - oh)),
        "seed_target_CZ_OH_P_deg": float(np.degrees(angle_from_points(cz, oh, placed_p))) if angle_from_points(cz, oh, placed_p) is not None else float('nan'),
    }
    dih = dihedral_deg(ce1, cz, oh, placed_p)
    if dih is not None:
        metrics["seed_target_CE1_CZ_OH_P_deg"] = float(dih)
    return transformed, metrics


def pick_shared_atom_names(original_resname: str, residue_atoms: Dict[str, Atom], template_atoms: Dict[str, TemplateAtom]) -> List[str]:
    original_resname = original_resname.strip().upper()
    if original_resname not in KABSCH_ALIGN_ATOMS:
        raise GraftingError(f"Unsupported parent residue for phosphograft alignment: {original_resname}")
    align_names = KABSCH_ALIGN_ATOMS[original_resname]
    missing_residue = [name for name in align_names if name not in residue_atoms]
    missing_template = [name for name in align_names if name not in template_atoms]
    if missing_residue:
        raise GraftingError(
            f"Target residue {original_resname} missing required Kabsch anchor atoms: {','.join(missing_residue)}"
        )
    if missing_template:
        raise GraftingError(
            f"Template for {original_resname}->{PHOSPHO_MAP.get(original_resname, '?')} missing required Kabsch anchor atoms: {','.join(missing_template)}"
        )
    if len(align_names) < 3:
        raise GraftingError(f"Too few Kabsch anchor atoms for alignment: found {len(align_names)}")
    return align_names


def validate_seeded_phosphate_placement(residue, original_resname: str) -> Dict[str, float]:
    original_resname = original_resname.strip().upper()
    anchor_name = BRIDGING_OXYGEN.get(original_resname)
    if anchor_name is None:
        return {}
    atoms = residue_atom_by_name(residue)
    if anchor_name not in atoms:
        raise GraftingError(f"Target residue missing expected bridging oxygen {anchor_name}")
    if "P" not in atoms:
        raise GraftingError("Seeded phosphate is missing P atom")
    anchor_coord = np.asarray(atoms[anchor_name].coord, dtype=float)
    p_coord = np.asarray(atoms["P"].coord, dtype=float)
    metrics = {
        "anchor_atom": anchor_name,
        "anchor_p_distance": float(np.linalg.norm(anchor_coord - p_coord)),
    }
    if original_resname == "SER":
        if "CB" in atoms and "OG" in atoms:
            angle = angle_from_points(np.asarray(atoms["CB"].coord, dtype=float), anchor_coord, p_coord)
            if angle is not None:
                metrics["anchor_angle_deg"] = float(np.degrees(angle))
    elif original_resname == "THR":
        if "CB" in atoms and "OG1" in atoms:
            angle = angle_from_points(np.asarray(atoms["CB"].coord, dtype=float), anchor_coord, p_coord)
            if angle is not None:
                metrics["anchor_angle_deg"] = float(np.degrees(angle))
        if "CG2" in atoms:
            metrics["cg2_p_distance"] = float(np.linalg.norm(np.asarray(atoms["CG2"].coord, dtype=float) - p_coord))
    elif original_resname == "TYR":
        if "CZ" in atoms and "OH" in atoms:
            angle = angle_from_points(np.asarray(atoms["CZ"].coord, dtype=float), anchor_coord, p_coord)
            if angle is not None:
                metrics["anchor_angle_deg"] = float(np.degrees(angle))

    if metrics["anchor_p_distance"] > EXPECTED_PHOSPHO_BOND_MAX:
        raise GraftingError(
            f"Seeded phosphate failed validation for {original_resname}: P is {metrics['anchor_p_distance']:.3f} Å from {anchor_name} (expected ~1.6 Å)"
        )
    if original_resname == "THR" and "cg2_p_distance" in metrics and metrics["cg2_p_distance"] <= metrics["anchor_p_distance"]:
        raise GraftingError(
            f"Seeded TPO phosphate failed validation: CG2-P distance {metrics['cg2_p_distance']:.3f} Å is <= OG1-P distance {metrics['anchor_p_distance']:.3f} Å"
        )
    return metrics


def next_serial(structure) -> int:
    serials = []
    for atom in structure.get_atoms():
        try:
            serials.append(int(atom.get_serial_number()))
        except Exception:
            continue
    return max(serials, default=0) + 1


def atom_radius(atom) -> float:
    return VDW_RADII.get((atom.element or "").strip().upper(), 1.70)


def clash_penalty(dist: float, r1: float, r2: float) -> float:
    cutoff = 0.75 * (r1 + r2)
    if dist >= cutoff:
        return 0.0
    overlap = cutoff - dist
    return overlap * overlap


def count_new_atom_clashes(structure, residue, new_atom_names: Iterable[str], threshold: float) -> int:
    new_atoms = [residue[name] for name in new_atom_names if name in residue]
    if not new_atoms:
        return 0
    residue_atoms = set(residue.get_atoms())
    search_atoms = [atom for atom in structure.get_atoms() if atom not in residue_atoms]
    if not search_atoms:
        return 0
    ns = NeighborSearch(search_atoms)
    clashing_neighbors: Set[int] = set()
    for new_atom in new_atoms:
        for atom in ns.search(new_atom.coord, threshold, level="A"):
            clashing_neighbors.add(id(atom))
    return len(clashing_neighbors)


def plddt_category(value: Optional[float]) -> str:
    if value is None:
        return "not_available"
    for threshold, label in PLDDT_BINS:
        if value >= threshold:
            return label
    return "not_available"


def get_plddt_for_residue(residue, altloc_mode: str, source: str) -> Tuple[Optional[float], Optional[float], str]:
    if source == "none":
        return None, None, "not_available"
    atoms = list(residue_atom_dict(residue, altloc_mode).values())
    if not atoms:
        return None, None, "not_available"
    vals = [float(atom.get_bfactor()) for atom in atoms]
    mean_val = float(np.mean(vals))
    min_val = float(np.min(vals))
    return mean_val, min_val, plddt_category(mean_val)


def plddt_qc_flag(mean_val: Optional[float], warn_below: float) -> str:
    if mean_val is None:
        return "not_available"
    if mean_val < 50.0:
        return "very_low_confidence"
    if mean_val < warn_below:
        return "warn_low_confidence"
    return "ok"


def rotate_points_about_axis(points: np.ndarray, axis_point1: np.ndarray, axis_point2: np.ndarray, angle_deg: float) -> np.ndarray:
    if points.size == 0:
        return points
    axis = axis_point2 - axis_point1
    norm = np.linalg.norm(axis)
    if norm < 1e-8:
        return points.copy()
    axis = axis / norm
    theta = math.radians(angle_deg)
    c = math.cos(theta)
    s = math.sin(theta)
    ux, uy, uz = axis
    rot = np.array([
        [c + ux * ux * (1 - c), ux * uy * (1 - c) - uz * s, ux * uz * (1 - c) + uy * s],
        [uy * ux * (1 - c) + uz * s, c + uy * uy * (1 - c), uy * uz * (1 - c) - ux * s],
        [uz * ux * (1 - c) - uy * s, uz * uy * (1 - c) + ux * s, c + uz * uz * (1 - c)],
    ], dtype=float)
    shifted = points - axis_point1
    return (rot @ shifted.T).T + axis_point1




def apply_recorded_bond_rotations(residue, rotation_text: str) -> List[str]:
    """Apply semicolon-separated bond rotations like ``CB-OG:60.0`` to a
    residue in-place. Used to transfer the best independent prescan
    phosphate orientation onto the real multi-site model before joint
    minimization. Only bonds declared in PRESCAN_BONDS / RELAX_BONDS are
    considered valid. Returns the applied rotation specs.
    """
    rotation_text = (rotation_text or "").strip()
    if not rotation_text:
        return []

    resname = residue.get_resname().strip().upper()
    allowed = list(PRESCAN_BONDS.get(resname, [])) + list(RELAX_BONDS.get(resname, []))
    if not allowed:
        return []

    all_atoms = {atom.get_name().strip(): atom for atom in residue.get_atoms()}
    applied: List[str] = []

    for chunk in rotation_text.split(';'):
        chunk = chunk.strip()
        if not chunk or ':' not in chunk or '-' not in chunk:
            continue
        axis_text, angle_text = chunk.split(':', 1)
        axis_atom1, axis_atom2 = [part.strip() for part in axis_text.split('-', 1)]
        try:
            angle_deg = float(angle_text)
        except ValueError:
            continue

        moving_names = None
        for a1, a2, names in allowed:
            if a1 == axis_atom1 and a2 == axis_atom2:
                moving_names = names
                break
        if moving_names is None or axis_atom1 not in all_atoms or axis_atom2 not in all_atoms:
            continue

        moving_atoms = [all_atoms[name] for name in moving_names if name in all_atoms]
        if not moving_atoms:
            continue
        axis_p1 = all_atoms[axis_atom1].coord.copy()
        axis_p2 = all_atoms[axis_atom2].coord.copy()
        current_xyz = np.array([atom.coord.copy() for atom in moving_atoms], dtype=float)
        rotated_xyz = rotate_points_about_axis(current_xyz, axis_p1, axis_p2, angle_deg)
        for atom, coord in zip(moving_atoms, rotated_xyz):
            atom.coord = coord
        applied.append(f"{axis_atom1}-{axis_atom2}:{angle_deg:.1f}")

    return applied


# Legacy: not used in production sequential mode.  Kept for backward
# compatibility with older benchmark scripts that call it directly.
def independent_prescan_site(
    input_structure: str,
    model_index: int,
    site: SiteSpec,
    template_atoms: Dict[str, TemplateAtom],
    phospho_code: str,
    clash_threshold: float,
    altloc_mode: str,
    protonation: str,
    plddt_mean: Optional[float],
    prescan_step_deg: float,
    prescan_contact_weight: float,
    prescan_basic_weight: float,
    prescan_hbond_weight: float,
    prescan_polar_clash_scale: float,
    prescan_pack_weight: float,
    prescan_geom_prior_weight: float,
    prescan_tpo_torsion_sigma: float,
    prescan_keep_kabsch_if_good: bool,
) -> Dict[str, object]:
    """Run the phospho prescan for one site on a fresh copy of the original
    unmodified input structure so site scores are not affected by grafting
    order in multi-site jobs.
    """
    work_structure = parse_structure(input_structure)
    work_model = list(work_structure.get_models())[model_index]
    work_residue = find_residue(work_model, site)
    serial_counter = [next_serial(work_structure)]
    graft_one_site(
        structure=work_structure,
        residue=work_residue,
        template_atoms=template_atoms,
        phospho_code=phospho_code,
        clash_threshold=clash_threshold,
        altloc_mode=altloc_mode,
        protonation=protonation,
        plddt_mean=plddt_mean,
        serial_counter=serial_counter,
    )
    prescan_info = prescan_phospho_rotamer(
        work_structure,
        work_residue,
        phospho_code,
        step_deg=prescan_step_deg,
        contact_weight=prescan_contact_weight,
        basic_weight=prescan_basic_weight,
        hbond_weight=prescan_hbond_weight,
        polar_clash_scale=prescan_polar_clash_scale,
        pack_weight=prescan_pack_weight,
        geom_prior_weight=prescan_geom_prior_weight,
        tpo_torsion_sigma=prescan_tpo_torsion_sigma,
        keep_kabsch_if_good=prescan_keep_kabsch_if_good,
    )
    prescan_info['prescan_context'] = 'independent_unmodified_model'
    return prescan_info

def score_moved_atoms(structure, residue, moved_atoms: Sequence[Atom]) -> float:
    residue_atoms = set(residue.get_atoms())
    score = 0.0
    for atom in structure.get_atoms():
        if atom in residue_atoms:
            continue
        r2 = atom_radius(atom)
        coord2 = atom.coord
        for moved in moved_atoms:
            d = float(np.linalg.norm(moved.coord - coord2))
            score += clash_penalty(d, atom_radius(moved), r2)
    return score


def relax_residue_by_torsion_scan(structure, residue, phospho_code: str, step_deg: float) -> Dict[str, object]:
    if phospho_code not in RELAX_BONDS:
        return {"relax_status": "SKIPPED", "relax_score_before": "", "relax_score_after": "", "relax_rotations": ""}
    all_atoms = {atom.get_name().strip(): atom for atom in residue.get_atoms()}
    score_before = score_moved_atoms(structure, residue, list(residue.get_atoms()))
    applied_rotations: List[str] = []
    for axis_atom1, axis_atom2, moving_names in RELAX_BONDS[phospho_code]:
        if axis_atom1 not in all_atoms or axis_atom2 not in all_atoms:
            continue
        moving_atoms = [all_atoms[name] for name in moving_names if name in all_atoms]
        if not moving_atoms:
            continue
        axis_p1 = all_atoms[axis_atom1].coord.copy()
        axis_p2 = all_atoms[axis_atom2].coord.copy()
        current_xyz = np.array([atom.coord.copy() for atom in moving_atoms], dtype=float)
        best_xyz = current_xyz.copy()
        best_angle = 0.0
        best_score = score_moved_atoms(structure, residue, moving_atoms)
        for angle in np.arange(0.0, 360.0, max(1.0, step_deg)):
            trial_xyz = rotate_points_about_axis(current_xyz, axis_p1, axis_p2, float(angle))
            for atom, coord in zip(moving_atoms, trial_xyz):
                atom.coord = coord
            trial_score = score_moved_atoms(structure, residue, moving_atoms)
            if trial_score < best_score:
                best_score = trial_score
                best_xyz = trial_xyz.copy()
                best_angle = float(angle)
        for atom, coord in zip(moving_atoms, best_xyz):
            atom.coord = coord
        if abs(best_angle) > 1e-6:
            applied_rotations.append(f"{axis_atom1}-{axis_atom2}:{best_angle:.1f}")
    score_after = score_moved_atoms(structure, residue, list(residue.get_atoms()))
    return {
        "relax_status": "OK",
        "relax_score_before": f"{score_before:.4f}",
        "relax_score_after": f"{score_after:.4f}",
        "relax_rotations": ";".join(applied_rotations),
    }



_BASIC_ATOMS = {
    ("LYS", "NZ"),
    ("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"),
    ("HIS", "ND1"), ("HIS", "NE2"), ("HID", "ND1"), ("HIE", "NE2"),
}
_HBOND_DONOR_ATOMS = _BASIC_ATOMS | {
    ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"),
    ("ASN", "ND2"), ("GLN", "NE2"), ("TRP", "NE1"),
}

# Same-residue atoms that should still contribute steric penalties during
# prescan scoring. The previous logic skipped the whole modified residue,
# which allowed unrealistic self-contacts to win.
_SAME_RESIDUE_CLASH_ATOMS = {
    "SER": {"N", "CA", "C", "O", "CB"},
    "SEP": {"N", "CA", "C", "O", "CB"},
    "THR": {"N", "CA", "C", "O", "CB", "CG2"},
    "TPO": {"N", "CA", "C", "O", "CB", "CG2"},
    "TYR": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PTR": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
}

def _atom_resname_atomname(atom: Atom) -> Tuple[str, str]:
    parent = atom.get_parent()
    resname = parent.get_resname().strip().upper() if parent is not None else ""
    return resname, atom.get_name().strip().upper()

def _is_basic_donor_atom(atom: Atom) -> bool:
    resname, aname = _atom_resname_atomname(atom)
    return (resname, aname) in _BASIC_ATOMS

def _is_hbond_donor_atom(atom: Atom) -> bool:
    resname, aname = _atom_resname_atomname(atom)
    if (resname, aname) in _HBOND_DONOR_ATOMS:
        return True
    # Backbone N donors on standard amino acids
    return aname == "N" and resname in THREE_TO_ONE

def _is_packable_protein_atom(atom: Atom) -> bool:
    resname, _ = _atom_resname_atomname(atom)
    elem = (atom.element or "").strip().upper()
    return elem != "H" and (resname in THREE_TO_ONE or resname in {"SEP", "TPO", "PTR"})

def _is_soft_polar_pair(moved: Atom, neighbor: Atom) -> bool:
    moved_elem = (moved.element or "").strip().upper()
    neigh_elem = (neighbor.element or "").strip().upper()
    return moved_elem == "O" and neigh_elem in {"N", "O", "S"} and (_is_basic_donor_atom(neighbor) or _is_hbond_donor_atom(neighbor) or neigh_elem in {"O", "S"})

def _distance_bonus(d: float, lo: float, peak: float, hi: float) -> float:
    if d < lo or d > hi:
        return 0.0
    if d <= peak:
        # gentle increase toward the optimum
        return 0.5 + 0.5 * ((d - lo) / max(peak - lo, 1e-6))
    return max(0.0, 1.0 - (d - peak) / max(hi - peak, 1e-6))


# --- Angular H-bond geometry (used for PTR scoring only) ---
_DONOR_PARENT_ATOM = {
    ("LYS", "NZ"):  "CE",
    ("ARG", "NE"):  "CD",
    ("ARG", "NH1"): "CZ",
    ("ARG", "NH2"): "CZ",
    ("HIS", "ND1"): "CG",
    ("HIS", "NE2"): "CD2",
    ("HID", "ND1"): "CG",
    ("HIE", "NE2"): "CD2",
    ("SER", "OG"):  "CB",
    ("THR", "OG1"): "CB",
    ("TYR", "OH"):  "CZ",
    ("ASN", "ND2"): "CG",
    ("GLN", "NE2"): "CD",
    ("TRP", "NE1"): "CE2",
}


def _get_donor_parent_coord(atom: Atom) -> Optional[np.ndarray]:
    """Return the coordinate of the parent heavy atom for an H-bond donor."""
    resname, aname = _atom_resname_atomname(atom)
    parent_name = _DONOR_PARENT_ATOM.get((resname, aname))
    if parent_name is None:
        if aname == "N":
            parent_name = "CA"
        else:
            return None
    parent_res = atom.get_parent()
    if parent_res is None:
        return None
    for a in parent_res.get_atoms():
        if a.get_name().strip().upper() == parent_name:
            return np.asarray(a.coord, dtype=float)
    return None


def _angular_hbond_factor(donor_parent_coord: np.ndarray, donor_coord: np.ndarray,
                          acceptor_coord: np.ndarray) -> float:
    """Angular scaling for H-bond quality: 1.0 at ≥140°, 0.0 at ≤80°."""
    v1 = donor_parent_coord - donor_coord
    v2 = acceptor_coord - donor_coord
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < 1e-8 or n2 < 1e-8:
        return 0.5
    cos_angle = float(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))
    angle_deg = float(np.degrees(np.arccos(cos_angle)))
    if angle_deg >= 140.0:
        return 1.0
    elif angle_deg <= 80.0:
        return 0.0
    else:
        return (angle_deg - 80.0) / 60.0


def _same_residue_scoring_atoms(residue, moved_atoms: Sequence[Atom]) -> List[Atom]:
    allowed = _SAME_RESIDUE_CLASH_ATOMS.get(residue.get_resname().strip().upper(), set())
    moved_set = set(moved_atoms)
    out = []
    for atom in residue.get_atoms():
        if atom in moved_set:
            continue
        if atom.get_name().strip().upper() in allowed:
            out.append(atom)
    return out


def _wrap_deg(x: float) -> float:
    return ((float(x) + 180.0) % 360.0) - 180.0


def _zero_geom_detail() -> Dict[str, float]:
    return {
        "geom_prior_penalty": 0.0,
        "geom_prior_dist_penalty": 0.0,
        "geom_prior_angle_penalty": 0.0,
        "geom_prior_torsion_penalty": 0.0,
        "geom_prior_mode": "none",
    }


def _tpo_geometry_prior_detail(residue, moved_atoms: Sequence[Atom], torsion_sigma: float = 20.0) -> Dict[str, float]:
    atoms = {a.get_name().strip().upper(): a for a in residue.get_atoms()}
    if residue.get_resname().strip().upper() not in {"THR", "TPO"}:
        return _zero_geom_detail()
    required = ["CB", "OG1", "CG2", "P"]
    if any(name not in atoms for name in required):
        return _zero_geom_detail()
    cb = np.asarray(atoms["CB"].coord, dtype=float)
    og1 = np.asarray(atoms["OG1"].coord, dtype=float)
    cg2 = np.asarray(atoms["CG2"].coord, dtype=float)
    p = np.asarray(atoms["P"].coord, dtype=float)
    d = float(np.linalg.norm(p - og1))
    ang = angle_from_points(cb, og1, p)
    dih = dihedral_deg(cg2, cb, og1, p)
    dist_sigma = 0.08
    angle_sigma = 8.0
    torsion_sigma = max(1.0, float(torsion_sigma))
    dist_pen = ((d - TPO_PRIOR_OG1_P) / dist_sigma) ** 2
    if ang is None:
        angle_pen = 0.0
    else:
        angle_deg = float(np.degrees(ang))
        angle_pen = ((angle_deg - TPO_PRIOR_CB_OG1_P_DEG) / angle_sigma) ** 2
    if dih is None:
        torsion_pen = 0.0
    else:
        torsion_pen = (_wrap_deg(dih - TPO_PRIOR_CG2_CB_OG1_P_DEG) / torsion_sigma) ** 2
    total = 0.5 * dist_pen + 0.5 * angle_pen + 1.5 * torsion_pen
    return {
        "geom_prior_penalty": float(total),
        "geom_prior_dist_penalty": float(dist_pen),
        "geom_prior_angle_penalty": float(angle_pen),
        "geom_prior_torsion_penalty": float(torsion_pen),
        "geom_prior_mode": "TPO:-58.9",
    }


def _sep_geometry_prior_detail(residue, moved_atoms: Sequence[Atom], torsion_sigma: float = 20.0) -> Dict[str, float]:
    """SEP geometry prior: bimodal torsion model.

    Sharp spike at +66.3° (σ≈3°, ~31% of experimental sites) — a locked
    salt-bridge/kinase-bound conformation.  Broad basin from -60° to +50°
    (~69% of sites) — context-dependent, no torsion penalty applied.

    Strategy: if the current torsion is near the sharp mode, reward it
    with a tight prior.  If it's in the broad basin, apply zero torsion
    penalty and let the scoring function (basic_weight, clash, hbond) decide.
    Outside both regions, apply a soft penalty to guide toward the basin.
    """
    atoms = {a.get_name().strip().upper(): a for a in residue.get_atoms()}
    if residue.get_resname().strip().upper() not in {"SER", "SEP"}:
        return _zero_geom_detail()
    required = ["CA", "CB", "OG", "P"]
    if any(name not in atoms for name in required):
        return _zero_geom_detail()
    ca = np.asarray(atoms["CA"].coord, dtype=float)
    cb = np.asarray(atoms["CB"].coord, dtype=float)
    og = np.asarray(atoms["OG"].coord, dtype=float)
    p = np.asarray(atoms["P"].coord, dtype=float)
    d = float(np.linalg.norm(p - og))
    ang = angle_from_points(cb, og, p)
    dih = dihedral_deg(ca, cb, og, p)
    dist_sigma = 0.08
    angle_sigma = 8.0
    dist_pen = ((d - SEP_PRIOR_OG_P) / dist_sigma) ** 2
    if ang is None:
        angle_pen = 0.0
    else:
        angle_deg = float(np.degrees(ang))
        angle_pen = ((angle_deg - SEP_PRIOR_CB_OG_P_DEG) / angle_sigma) ** 2
    if dih is None:
        torsion_pen = 0.0
        best_mode = "SEP:none"
    else:
        # Check sharp spike at +66.3°
        spike_diff = abs(_wrap_deg(dih - SEP_PRIOR_SHARP_MODE_DEG))
        spike_sigma = max(1.0, SEP_PRIOR_SHARP_MODE_SIGMA)
        # Check broad basin [-60°, +50°]
        broad_lo, broad_hi = SEP_PRIOR_BROAD_RANGE
        in_broad = broad_lo <= _wrap_deg(dih) <= broad_hi

        if spike_diff <= 3.0 * spike_sigma:
            # Near sharp spike: apply tight torsion prior (reward locked state)
            torsion_pen = (spike_diff / spike_sigma) ** 2
            best_mode = f"SEP:spike:{SEP_PRIOR_SHARP_MODE_DEG:+.1f}"
        elif in_broad:
            # In broad basin: no torsion penalty, let scoring decide
            torsion_pen = 0.0
            best_mode = f"SEP:broad:{dih:+.1f}"
        else:
            # Outside both regions: soft penalty toward nearest boundary
            dist_to_broad_lo = abs(_wrap_deg(dih - broad_lo))
            dist_to_broad_hi = abs(_wrap_deg(dih - broad_hi))
            dist_to_spike = spike_diff
            nearest_dist = min(dist_to_broad_lo, dist_to_broad_hi, dist_to_spike)
            outer_sigma = max(1.0, float(torsion_sigma))
            torsion_pen = (nearest_dist / outer_sigma) ** 2
            best_mode = f"SEP:outer:{dih:+.1f}"
    total = 0.5 * dist_pen + 0.5 * angle_pen + 1.5 * torsion_pen
    return {
        "geom_prior_penalty": float(total),
        "geom_prior_dist_penalty": float(dist_pen),
        "geom_prior_angle_penalty": float(angle_pen),
        "geom_prior_torsion_penalty": float(torsion_pen),
        "geom_prior_mode": best_mode,
    }




def _ptr_geometry_prior_detail(residue, moved_atoms: Sequence[Atom], torsion_sigma: float = 20.0) -> Dict[str, float]:
    """PTR geometry prior: distance and angle only, NO torsion.

    PDB-wide analysis (n=501 all, n=208 single-site) confirmed the
    CE1-CZ-OH-P torsion is uniformly distributed — the tyrosyl hydroxyl
    rotates freely.  Phosphate orientation is environment-driven (83%
    salt-bridge, but 50/50 toward vs lateral).  Only the P-OH bond length
    and CZ-OH-P angle carry geometric information.
    """
    atoms = {a.get_name().strip().upper(): a for a in residue.get_atoms()}
    if residue.get_resname().strip().upper() not in {"TYR", "PTR"}:
        return _zero_geom_detail()
    required = ["CZ", "OH", "P"]
    if any(name not in atoms for name in required):
        return _zero_geom_detail()
    cz = np.asarray(atoms["CZ"].coord, dtype=float)
    oh = np.asarray(atoms["OH"].coord, dtype=float)
    p = np.asarray(atoms["P"].coord, dtype=float)
    d = float(np.linalg.norm(p - oh))
    ang = angle_from_points(cz, oh, p)
    dist_sigma = 0.08
    angle_sigma = 8.0
    dist_pen = ((d - PTR_PRIOR_OH_P) / dist_sigma) ** 2
    if ang is None:
        angle_pen = 0.0
    else:
        angle_deg = float(np.degrees(ang))
        angle_pen = ((angle_deg - PTR_PRIOR_CZ_OH_P_DEG) / angle_sigma) ** 2

    # No torsion penalty — orientation is environment-driven
    torsion_pen = 0.0
    best_mode = "PTR:env"

    total = 0.5 * dist_pen + 0.5 * angle_pen
    return {
        "geom_prior_penalty": float(total),
        "geom_prior_dist_penalty": float(dist_pen),
        "geom_prior_angle_penalty": float(angle_pen),
        "geom_prior_torsion_penalty": float(torsion_pen),
        "geom_prior_mode": best_mode,
    }


def _geometry_prior_detail(residue, moved_atoms: Sequence[Atom], torsion_sigma: float = 20.0) -> Dict[str, float]:
    resname = residue.get_resname().strip().upper()
    if resname in {"THR", "TPO"}:
        return _tpo_geometry_prior_detail(residue, moved_atoms, torsion_sigma=torsion_sigma)
    if resname in {"SER", "SEP"}:
        return _sep_geometry_prior_detail(residue, moved_atoms, torsion_sigma=torsion_sigma)
    if resname in {"TYR", "PTR"}:
        return _ptr_geometry_prior_detail(residue, moved_atoms, torsion_sigma=torsion_sigma)
    return _zero_geom_detail()


def score_moved_atoms_with_contacts(
    structure, residue, moved_atoms: Sequence[Atom],
    contact_weight: float = 0.5,
    basic_weight: float = 2.0,
    hbond_weight: float = 0.5,
    polar_clash_scale: float = 0.5,
    pack_weight: float = 0.0,
    geom_prior_weight: float = 0.0,
    tpo_torsion_sigma: float = 20.0,
) -> float:
    detail = score_moved_atoms_detailed(
        structure, residue, moved_atoms,
        contact_weight=contact_weight,
        clash_threshold=2.0,
        basic_weight=basic_weight,
        hbond_weight=hbond_weight,
        polar_clash_scale=polar_clash_scale,
        pack_weight=pack_weight,
    )
    return detail["total_score"]


def score_moved_atoms_detailed(
    structure, residue, moved_atoms: Sequence[Atom],
    contact_weight: float = 0.5,
    clash_threshold: float = 2.0,
    basic_weight: float = 2.0,
    hbond_weight: float = 0.5,
    polar_clash_scale: float = 0.5,
    pack_weight: float = 0.0,
    geom_prior_weight: float = 0.0,
    tpo_torsion_sigma: float = 20.0,
) -> Dict[str, float]:
    """Like score_moved_atoms_with_contacts but returns breakdown."""
    POLAR_ELEMS = {"N", "O", "S"}
    resname_score = residue.get_resname().strip().upper()
    residue_atoms = set(residue.get_atoms())
    same_residue_atoms = set(_same_residue_scoring_atoms(residue, moved_atoms))
    ptr_pack_best_by_residue = {}   # residue-id -> best packing contribution for PTR only
    ptr_basic_best_by_residue = {}  # residue-id -> best basic contribution for PTR only
    clash_score = 0.0
    self_clash_score = 0.0
    contact_bonus = 0.0
    basic_bonus = 0.0
    hbond_bonus = 0.0
    pack_bonus = 0.0
    geom_detail = _geometry_prior_detail(residue, moved_atoms, torsion_sigma=tpo_torsion_sigma)
    geom_prior_penalty = geom_detail["geom_prior_penalty"]
    n_clashing = 0
    n_self_clashing = 0

    polar_moved = [a for a in moved_atoms if (a.element or "").strip().upper() in POLAR_ELEMS]
    phosphate_oxy_moved = [a for a in moved_atoms if (a.element or "").strip().upper() == "O"]

    for atom in structure.get_atoms():
        in_same_residue = atom in residue_atoms
        if in_same_residue and atom not in same_residue_atoms:
            continue

        r2 = atom_radius(atom)
        coord2 = atom.coord
        elem2 = (atom.element or "").strip().upper()
        is_polar2 = elem2 in POLAR_ELEMS
        atom_clashes = False

        for moved in moved_atoms:
            d = float(np.linalg.norm(moved.coord - coord2))
            cp = clash_penalty(d, atom_radius(moved), r2)
            if not in_same_residue and cp > 0.0 and _is_soft_polar_pair(moved, atom):
                cp *= polar_clash_scale
            clash_score += cp
            if in_same_residue:
                self_clash_score += cp
            if d < clash_threshold and not atom_clashes:
                n_clashing += 1
                if in_same_residue:
                    n_self_clashing += 1
                atom_clashes = True

        # Same-residue atoms only contribute steric penalties, not favorable
        # chemistry bonuses, to avoid selecting self-collapsed phosphate poses.
        if in_same_residue:
            continue

        if is_polar2 and polar_moved:
            for pm in polar_moved:
                d = float(np.linalg.norm(pm.coord - coord2))
                if 2.3 <= d <= 4.5:
                    strength = max(0.0, 1.0 - (d - 2.3) / 2.2)
                    contact_bonus += strength

        if phosphate_oxy_moved:
            for pm in phosphate_oxy_moved:
                d = float(np.linalg.norm(pm.coord - coord2))
                if _is_basic_donor_atom(atom):
                    bc = _distance_bonus(d, 2.4, 3.0, 4.2)
                    if resname_score in {"TYR", "PTR"} and bc > 0.0:
                        # PTR: angular H-bond factor scales distance bonus by
                        # donor geometry quality.  Suppresses contributions from
                        # basic atoms that are nearby but poorly oriented.
                        parent_coord = _get_donor_parent_coord(atom)
                        if parent_coord is not None:
                            bc *= _angular_hbond_factor(parent_coord, coord2, pm.coord)
                        # PTR: per-residue cap for basic bonus
                        rp = atom.get_parent()
                        rk = (rp.get_parent().id, rp.id[1], rp.id[2].strip(), rp.get_resname().strip())
                        prev = ptr_basic_best_by_residue.get(rk, 0.0)
                        if bc > prev:
                            ptr_basic_best_by_residue[rk] = bc
                    else:
                        # SEP/TPO: v5.0 per-atom basic bonus (proven)
                        basic_bonus += bc
                elif _is_hbond_donor_atom(atom):
                    hb = _distance_bonus(d, 2.4, 2.9, 3.6)
                    if resname_score in {"TYR", "PTR"} and hb > 0.0:
                        parent_coord = _get_donor_parent_coord(atom)
                        if parent_coord is not None:
                            hb *= _angular_hbond_factor(parent_coord, coord2, pm.coord)
                    hbond_bonus += hb
                if _is_packable_protein_atom(atom):
                    pack_contrib = _distance_bonus(d, 2.8, 3.6, 5.0)
                    if resname_score in {"TYR", "PTR"}:
                        rp = atom.get_parent()
                        rk = (rp.get_parent().id, rp.id[1], rp.id[2].strip(), rp.get_resname().strip())
                        prev = ptr_pack_best_by_residue.get(rk, 0.0)
                        if pack_contrib > prev:
                            ptr_pack_best_by_residue[rk] = pack_contrib
                    else:
                        pack_bonus += pack_contrib

    if resname_score in {"TYR", "PTR"}:
        if ptr_pack_best_by_residue:
            pack_bonus += sum(ptr_pack_best_by_residue.values())
        if ptr_basic_best_by_residue:
            basic_bonus += sum(ptr_basic_best_by_residue.values())

    eff_contact_weight = contact_weight
    eff_basic_weight = basic_weight
    eff_hbond_weight = hbond_weight
    eff_pack_weight = pack_weight
    clash_count_penalty = 0.0
    if resname_score in {"SER", "SEP"}:
        # SEP scoring (v5, PDB-wide audit: n=228 single-site)
        # 64% salt-bridge, 66% toward-basic orientation.
        # Packing drops 24→19 upon phosphorylation → reduce pack_bonus.
        # Contact support still important for the 36% without basic contacts.
        eff_contact_weight = max(contact_weight, 1.0)
        eff_basic_weight = basic_weight * 0.8     # moderate-high (was 1.0)
        eff_pack_weight = pack_weight * 0.4       # reduced (was 0.5): packing drops upon phospho
        clash_count_penalty = 2.0 * float(n_clashing)
    elif resname_score in {"TYR", "PTR"}:
        # PTR scoring (v5, PDB-wide audit: n=208 single-site)
        # 83% salt-bridge but 50/50 toward vs lateral orientation.
        # Packing unchanged (15→15) → keep pack low but not penalizing.
        # basic_weight effective because distance matters, not directionality.
        eff_contact_weight = contact_weight
        eff_basic_weight = basic_weight * 0.8     # moderate-high (was 0.85)
        eff_hbond_weight = hbond_weight * 0.85
        eff_pack_weight = pack_weight * 0.3       # low (was 0.6): environment-driven, not packing
        clash_count_penalty = 2.5 * float(n_clashing)
    # TPO: no rebalancing needed — all weights at face value
    # (81% salt-bridge, 76% toward-basic, packing preserved 21→23)

    total = (
        clash_score
        + clash_count_penalty
        - eff_contact_weight * contact_bonus
        - eff_basic_weight * basic_bonus
        - eff_hbond_weight * hbond_bonus
        - eff_pack_weight * pack_bonus
        + geom_prior_weight * geom_prior_penalty
    )
    return {
        "geom_prior_penalty": geom_prior_penalty,
        "geom_prior_dist_penalty": geom_detail["geom_prior_dist_penalty"],
        "geom_prior_angle_penalty": geom_detail["geom_prior_angle_penalty"],
        "geom_prior_torsion_penalty": geom_detail["geom_prior_torsion_penalty"],
        "geom_prior_mode": geom_detail["geom_prior_mode"],
        "clash_score": clash_score,
        "self_clash_score": self_clash_score,
        "contact_bonus": contact_bonus,
        "basic_bonus": basic_bonus,
        "hbond_bonus": hbond_bonus,
        "pack_bonus": pack_bonus,
        "effective_contact_weight": eff_contact_weight,
        "effective_pack_weight": eff_pack_weight,
        "clash_count_penalty": clash_count_penalty,
        "total_score": total,
        "n_clashing": n_clashing,
        "n_self_clashing": n_self_clashing,
    }


def _neighborhood_smoothed_scores(candidates: List[Dict], step_deg: float = 30.0) -> List[float]:
    """Compute neighborhood-smoothed scores for prescan candidates.

    For each candidate at angle θ, the smoothed score is:
      0.5 * score(θ) + 0.25 * score(θ-step) + 0.25 * score(θ+step)

    Rewards candidates in consistently good angular regions.  Isolated
    low-score frames (likely scoring noise) get penalized.
    """
    if len(candidates) <= 2:
        return [float(c["detail"]["total_score"]) for c in candidates]
    angle_score = {}
    for c in candidates:
        angle = round(float(c["angle"]) % 360.0, 1)
        angle_score[angle] = float(c["detail"]["total_score"])
    smoothed = []
    for c in candidates:
        angle = round(float(c["angle"]) % 360.0, 1)
        score_center = angle_score.get(angle, 0.0)
        angle_prev = round((angle - step_deg) % 360.0, 1)
        angle_next = round((angle + step_deg) % 360.0, 1)
        score_prev = angle_score.get(angle_prev, score_center)
        score_next = angle_score.get(angle_next, score_center)
        smoothed.append(0.5 * score_center + 0.25 * score_prev + 0.25 * score_next)
    return smoothed


def _choose_best_prescan_candidate(phospho_code: str, candidates: List[Dict]) -> Optional[Dict]:
    """Select the best prescan candidate.

    - TPO/SEP: lowest total_score (v5.0 behavior, validated).
    - PTR: consensus-aware tiered selection with neighborhood-smoothed
      scores and angular H-bond geometry (v5.3 improvement).
    """
    if not candidates:
        return None

    if phospho_code != "PTR":
        # TPO/SEP: simple lowest-score selection (v5.0, proven)
        return min(candidates, key=lambda c: float(c["detail"]["total_score"]))

    # PTR: consensus-aware tiered selection with smoothed scores
    smoothed = _neighborhood_smoothed_scores(candidates)

    def tier_members(max_clashes: int, max_geom: float):
        return [
            (i, c) for i, c in enumerate(candidates)
            if int(c["detail"].get("n_clashing", 999)) <= max_clashes
            and float(c["detail"].get("geom_prior_penalty", 999.0)) <= max_geom
        ]

    tiers = [
        tier_members(1, 2.0),
        tier_members(0, 4.0),
        tier_members(1, 4.0),
        [(i, c) for i, c in enumerate(candidates) if int(c["detail"].get("n_clashing", 999)) == 0],
        [(i, c) for i, c in enumerate(candidates) if int(c["detail"].get("n_clashing", 999)) <= 1],
        list(enumerate(candidates)),
    ]

    for tier in tiers:
        if tier:
            best_ic = min(tier, key=lambda ic: smoothed[ic[0]])
            return best_ic[1]
    return min(candidates, key=lambda c: float(c["detail"]["total_score"]))


def _rank_prescan_candidates(phospho_code: str, candidates: List[Dict], k: int = 3) -> List[Dict]:
    """Rank prescan candidates while guaranteeing that rank 1 is the
    candidate selected by the production scoring rule.

    TPO and SEP use the raw composite score. PTR retains its consensus-aware
    rank-1 selection and orders the remaining candidates by the smoothed
    neighbourhood score.
    """
    if not candidates or k <= 0:
        return []

    selected = _choose_best_prescan_candidate(phospho_code, candidates)
    if phospho_code == "PTR":
        smoothed = _neighborhood_smoothed_scores(candidates)
        ordered = [candidate for _, candidate in sorted(
            enumerate(candidates), key=lambda item: smoothed[item[0]])
        ]
    else:
        ordered = sorted(candidates, key=lambda c: float(c["detail"]["total_score"]))

    ranked: List[Dict] = []
    if selected is not None:
        ranked.append(selected)
    ranked.extend(candidate for candidate in ordered if candidate is not selected)
    return ranked[:k]


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "SEP": "pS", "TPO": "pT", "PTR": "pY",
}


def _build_2d_basis(axis_p1: np.ndarray, axis_p2: np.ndarray):
    """Return (origin, axis_norm, e1, e2) for projecting 3D points onto
    the plane perpendicular to the rotation axis at axis_p2."""
    axis = axis_p2 - axis_p1
    norm = np.linalg.norm(axis)
    if norm < 1e-8:
        return axis_p2, np.array([0, 0, 1.0]), np.array([1, 0, 0.0]), np.array([0, 1, 0.0])
    axis_norm = axis / norm
    # Pick a vector not parallel to axis
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(axis_norm, ref)) > 0.9:
        ref = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(axis_norm, ref)
    e1 = e1 / np.linalg.norm(e1)
    e2 = np.cross(axis_norm, e1)
    e2 = e2 / np.linalg.norm(e2)
    return axis_p2, axis_norm, e1, e2


def _project_2d(point: np.ndarray, origin: np.ndarray, e1: np.ndarray, e2: np.ndarray) -> Tuple[float, float]:
    """Project a 3D point onto the 2D plane defined by origin, e1, e2."""
    v = point - origin
    return float(np.dot(v, e1)), float(np.dot(v, e2))


def prescan_phospho_rotamer(
    structure, residue, phospho_code: str, step_deg: float = 30.0,
    contact_weight: float = 0.5,
    clash_threshold: float = 2.0,
    basic_weight: float = 2.0,
    hbond_weight: float = 0.5,
    polar_clash_scale: float = 0.5,
    pack_weight: float = 0.0,
    geom_prior_weight: float = 2.0,
    tpo_torsion_sigma: float = 20.0,
    keep_kabsch_if_good: bool = False,
) -> Dict[str, object]:
    """Rotate only the PO3 group around the bridging O–P bond to find the
    best phosphate orientation before OpenMM minimization.

    Records the full score distribution, per-frame contacts with neighbor
    residues, and 2D projected positions for visualization.
    """
    import json as _json

    empty_result = {
        "prescan_status": "SKIPPED",
        "prescan_score_before": "",
        "prescan_score_after": "",
        "prescan_rotations": "",
        "prescan_score_min": "",
        "prescan_score_max": "",
        "prescan_score_mean": "",
        "prescan_n_zero_clash": "",
        "prescan_zero_clash_arc": "",
        "prescan_frames_json": "",
        "prescan_top_k_candidates": [],
    }
    if phospho_code not in PRESCAN_BONDS:
        return empty_result

    all_atoms = {atom.get_name().strip(): atom for atom in residue.get_atoms()}
    all_atom_list = list(residue.get_atoms())
    score_before = score_moved_atoms_with_contacts(
        structure, residue, all_atom_list, contact_weight,
        basic_weight=basic_weight, hbond_weight=hbond_weight, polar_clash_scale=polar_clash_scale, pack_weight=pack_weight,
        geom_prior_weight=geom_prior_weight, tpo_torsion_sigma=tpo_torsion_sigma,
    )
    applied_rotations: List[str] = []
    frames: List[Dict] = []
    top_k_candidates: List[Dict] = []

    # Build NeighborSearch for fast contact queries (non-self, non-H atoms).
    residue_atoms_set = set(residue.get_atoms())
    search_atoms = [a for a in structure.get_atoms()
                    if a not in residue_atoms_set
                    and (a.element or "").strip().upper() != "H"]
    ns = NeighborSearch(search_atoms) if search_atoms else None
    CONTACT_RADIUS = 5.0

    # Collect all neighbor metadata across all frames.
    all_neighbor_keys: Dict[str, Dict] = {}  # "RES chain:resseq" -> info

    for axis_atom1, axis_atom2, moving_names in PRESCAN_BONDS[phospho_code]:
        if axis_atom1 not in all_atoms or axis_atom2 not in all_atoms:
            continue
        moving_atoms = [all_atoms[name] for name in moving_names if name in all_atoms]
        if not moving_atoms:
            continue
        moving_atom_names = [a.get_name().strip() for a in moving_atoms]
        axis_p1 = all_atoms[axis_atom1].coord.copy()
        axis_p2 = all_atoms[axis_atom2].coord.copy()
        current_xyz = np.array([atom.coord.copy() for atom in moving_atoms], dtype=float)
        best_xyz = current_xyz.copy()
        best_angle = 0.0
        best_score_detail = score_moved_atoms_detailed(
            structure, residue, moving_atoms, contact_weight, clash_threshold,
            basic_weight=basic_weight, hbond_weight=hbond_weight, polar_clash_scale=polar_clash_scale, pack_weight=pack_weight,
            geom_prior_weight=geom_prior_weight, tpo_torsion_sigma=tpo_torsion_sigma,
        )
        best_score = best_score_detail["total_score"]
        prescan_candidates: List[Dict] = []

        if keep_kabsch_if_good and best_score_detail["n_clashing"] == 0 and (best_score_detail["basic_bonus"] > 0.0 or best_score_detail["hbond_bonus"] > 0.0):
            angles = [0.0]
        else:
            angles = list(np.arange(0.0, 360.0, max(1.0, step_deg)))

        # 2D projection basis from rotation axis.
        origin_2d, _, e1, e2 = _build_2d_basis(axis_p1, axis_p2)

        for angle in angles:
            trial_xyz = rotate_points_about_axis(current_xyz, axis_p1, axis_p2, float(angle))
            for atom, coord in zip(moving_atoms, trial_xyz):
                atom.coord = coord
            detail = score_moved_atoms_detailed(
                structure, residue, moving_atoms, contact_weight, clash_threshold,
                basic_weight=basic_weight, hbond_weight=hbond_weight, polar_clash_scale=polar_clash_scale, pack_weight=pack_weight,
                geom_prior_weight=geom_prior_weight, tpo_torsion_sigma=tpo_torsion_sigma,
            )

            # Collect contacts for this frame.
            frame_contacts: List[Dict] = []
            if ns is not None:
                seen_pairs: Set[str] = set()
                for po3_name, po3_coord in zip(moving_atom_names, trial_xyz):
                    for neighbor in ns.search(po3_coord, CONTACT_RADIUS, level="A"):
                        parent = neighbor.get_parent()
                        parent_chain = parent.get_parent().id if parent.get_parent() else "?"
                        resname = parent.get_resname().strip()
                        resseq = parent.id[1]
                        res_key = f"{resname} {parent_chain}:{resseq}"
                        pair_key = f"{res_key}:{neighbor.get_name().strip()}:{po3_name}"
                        if pair_key in seen_pairs:
                            continue
                        seen_pairs.add(pair_key)
                        d = float(np.linalg.norm(po3_coord - neighbor.coord))
                        is_clash = d < clash_threshold
                        frame_contacts.append({
                            "res": res_key,
                            "n_atom": neighbor.get_name().strip(),
                            "po3": po3_name,
                            "dist": round(d, 2),
                            "clash": is_clash,
                        })
                        # Record neighbor metadata for 2D positioning.
                        if res_key not in all_neighbor_keys:
                            one_letter = THREE_TO_ONE.get(resname, resname[:3])
                            all_neighbor_keys[res_key] = {
                                "label": res_key,
                                "short": f"{one_letter}{resseq}",
                                "coord_3d": neighbor.coord.copy(),
                                "min_dist": d,
                            }
                        elif d < all_neighbor_keys[res_key]["min_dist"]:
                            all_neighbor_keys[res_key]["coord_3d"] = neighbor.coord.copy()
                            all_neighbor_keys[res_key]["min_dist"] = d
            # Sort contacts by distance.
            frame_contacts.sort(key=lambda c: c["dist"])

            # Project PO3 atoms to 2D.
            po3_2d = {}
            for name, coord in zip(moving_atom_names, trial_xyz):
                x2, y2 = _project_2d(coord, origin_2d, e1, e2)
                po3_2d[name] = [round(x2, 3), round(y2, 3)]

            frame = {
                "angle": round(float(angle), 1),
                "clash_score": round(detail["clash_score"], 4),
                "contact_bonus": round(detail["contact_bonus"], 4),
                "basic_bonus": round(detail["basic_bonus"], 4),
                "hbond_bonus": round(detail["hbond_bonus"], 4),
                "pack_bonus": round(detail["pack_bonus"], 4),
                "geom_prior_penalty": round(detail["geom_prior_penalty"], 4),
                "geom_prior_torsion_penalty": round(detail["geom_prior_torsion_penalty"], 4),
                "geom_prior_mode": detail.get("geom_prior_mode", "none"),
                "total_score": round(detail["total_score"], 4),
                "n_clashing": detail["n_clashing"],
                "coords": {
                    name: [round(float(c), 3) for c in coord]
                    for name, coord in zip(moving_atom_names, trial_xyz)
                },
                "po3_2d": po3_2d,
                "contacts": frame_contacts[:15],  # cap to keep JSON small
            }
            frames.append(frame)
            prescan_candidates.append({
                "angle": float(angle),
                "xyz": trial_xyz.copy(),
                "detail": detail,
            })

        chosen_candidate = _choose_best_prescan_candidate(phospho_code, prescan_candidates)
        if chosen_candidate is not None:
            best_xyz = chosen_candidate["xyz"].copy()
            best_angle = float(chosen_candidate["angle"])
            best_score_detail = chosen_candidate["detail"]
            best_score = best_score_detail["total_score"]

        top_k_candidates = []
        for candidate in _rank_prescan_candidates(phospho_code, prescan_candidates, k=3):
            top_k_candidates.append({
                "angle": float(candidate["angle"]),
                "total_score": float(candidate["detail"]["total_score"]),
                "coords": {
                    name: [float(value) for value in coord]
                    for name, coord in zip(moving_atom_names, candidate["xyz"])
                },
            })
        rank_by_angle = {
            round(float(candidate["angle"]), 1): rank
            for rank, candidate in enumerate(top_k_candidates, start=1)
        }
        for frame in frames:
            frame_angle = round(float(frame.get("angle", 0.0)), 1)
            if frame_angle in rank_by_angle:
                frame["selection_rank"] = rank_by_angle[frame_angle]

        for atom, coord in zip(moving_atoms, best_xyz):
            atom.coord = coord
        if abs(best_angle) > 1e-6:
            applied_rotations.append(f"{axis_atom1}-{axis_atom2}:{best_angle:.1f}")

    score_after = score_moved_atoms_with_contacts(
        structure, residue, all_atom_list, contact_weight,
        basic_weight=basic_weight, hbond_weight=hbond_weight, polar_clash_scale=polar_clash_scale, pack_weight=pack_weight,
        geom_prior_weight=geom_prior_weight, tpo_torsion_sigma=tpo_torsion_sigma,
    )

    # Compute 2D positions for all neighbor residues.
    neighbors_2d: List[Dict] = []
    if all_neighbor_keys:
        for key, info in all_neighbor_keys.items():
            x2, y2 = _project_2d(info["coord_3d"], origin_2d, e1, e2)
            neighbors_2d.append({
                "label": info["label"],
                "short": info["short"],
                "x": round(x2, 3),
                "y": round(y2, 3),
            })

    # Compute distribution stats from frames.
    if frames:
        all_scores = [f["total_score"] for f in frames]
        all_clashes = [f["n_clashing"] for f in frames]
        n_zero = sum(1 for c in all_clashes if c == 0)
        zero_flags = [1 if c == 0 else 0 for c in all_clashes]
        doubled = zero_flags + zero_flags
        max_run = 0
        run = 0
        for flag in doubled:
            if flag:
                run += 1
                max_run = max(max_run, run)
            else:
                run = 0
        max_run = min(max_run, len(zero_flags))
        arc_deg = round(max_run * max(1.0, step_deg), 1)
    else:
        all_scores = []
        n_zero = 0
        arc_deg = 0.0

    # Bundle frames with neighbor data for the JSON sidecar.
    frames_bundle = {
        "neighbors": neighbors_2d,
        "residue": phospho_code,
        "frames": frames,
    }

    result = {
        "prescan_status": "OK",
        "prescan_score_before": f"{score_before:.4f}",
        "prescan_score_after": f"{score_after:.4f}",
        "prescan_rotations": ";".join(applied_rotations),
        "prescan_score_min": f"{min(all_scores):.4f}" if all_scores else "",
        "prescan_score_max": f"{max(all_scores):.4f}" if all_scores else "",
        "prescan_score_mean": f"{float(np.mean(all_scores)):.4f}" if all_scores else "",
        "prescan_n_zero_clash": n_zero,
        "prescan_n_total_steps": len(frames),
        "prescan_zero_clash_arc": f"{arc_deg:.1f}",
        "prescan_frames_json": _json.dumps(frames_bundle),
        # In-memory absolute coordinates used to construct each ranked pose.
        # DictWriter ignores this field, while the frame sidecar retains all
        # candidates for auditability.
        "prescan_top_k_candidates": top_k_candidates,
    }
    return result


def apply_prescan_candidate_coords(residue, candidate: Dict[str, object]) -> None:
    """Apply one prescan candidate using its absolute saved coordinates.

    Using absolute coordinates avoids compounding a rank-2/rank-3 scan angle
    on top of the already-applied rank-1 orientation.
    """
    atom_map = {atom.get_name().strip(): atom for atom in residue.get_atoms()}
    coords = candidate.get("coords", {}) if isinstance(candidate, dict) else {}
    if not isinstance(coords, dict):
        return
    for atom_name, xyz in coords.items():
        atom = atom_map.get(str(atom_name))
        if atom is not None:
            atom.coord = np.asarray(xyz, dtype=float)


def find_flexible_neighbors(
    model,
    modified_residues: Sequence,
    radius: float = 6.0,
    residue_types: Optional[Set[str]] = None,
) -> List:
    """Return non-modified residues near any modified site whose sidechains
    should receive soft positional restraints during minimization.

    Any residue packed against the phosphosite may need to shift to
    accommodate the larger PO3 group — not just charged/polar ones.
    A MET sulfur at 3.4 Å, a LEU blocking the phosphate, or a SER
    hydroxyl in the way all need room to move.

    Parameters
    ----------
    model : Bio.PDB Model
        The model being processed.
    modified_residues : sequence of Bio.PDB Residue
        The phosphorylated residues.
    radius : float
        Distance cutoff in Angstroms.  Any residue with at least one
        heavy atom within *radius* of a modified residue atom is included.
    residue_types : set of str or None
        Three-letter residue names to consider.  Defaults to all standard
        amino acids (including SEP/TPO/PTR for adjacent phosphosites).

    Returns
    -------
    list of Bio.PDB Residue
        Unique neighbor residues, in no particular order.
    """
    if residue_types is None:
        residue_types = FLEXIBLE_NEIGHBOR_RESTYPES_ALL
    modified_ids = {r.full_id for r in modified_residues}

    # Collect all non-modified atoms for the neighbor search.
    non_mod_atoms = [
        a for a in model.get_atoms()
        if a.get_parent().full_id not in modified_ids
    ]
    if not non_mod_atoms:
        return []

    ns = NeighborSearch(non_mod_atoms)
    found_ids: Set[Tuple] = set()
    neighbors: List = []

    for residue in modified_residues:
        for atom in residue.get_atoms():
            nearby_residues = ns.search(atom.coord, radius, level="R")
            for r in nearby_residues:
                rid = r.full_id
                if rid in modified_ids or rid in found_ids:
                    continue
                resname = r.get_resname().strip().upper()
                if resname in residue_types:
                    found_ids.add(rid)
                    neighbors.append(r)
    return neighbors


def choose_protonated_oxygens(residue, n_protons: int) -> List[str]:
    if n_protons <= 0 or "P" not in residue:
        return []
    oxygen_names = [name for name in ["O1P", "O2P", "O3P"] if name in residue]
    if not oxygen_names:
        return []
    p = residue["P"].coord
    candidates = []
    for name in oxygen_names:
        o_atom = residue[name]
        direction = o_atom.coord - p
        norm = np.linalg.norm(direction)
        if norm < 1e-8:
            continue
        direction = direction / norm
        trial_h = o_atom.coord + 0.98 * direction
        penalty = 0.0
        structure = residue.get_parent().get_parent().get_parent()
        for atom in structure.get_atoms():
            if atom is o_atom or atom.get_parent() is residue:
                continue
            penalty += clash_penalty(float(np.linalg.norm(trial_h - atom.coord)), VDW_RADII["H"], atom_radius(atom))
        candidates.append((penalty, name))
    candidates.sort()
    return [name for _, name in candidates[:n_protons]]


def add_phosphate_hydrogens(residue, protonation: str, bfactor: float, serial_counter: List[int]) -> List[str]:
    n_protons = int(protonation)
    if n_protons <= 0 or "P" not in residue:
        return []
    protonated_oxygens = choose_protonated_oxygens(residue, n_protons)
    p = residue["P"].coord
    added = []
    for idx, oxy_name in enumerate(protonated_oxygens, start=1):
        if oxy_name not in residue:
            continue
        o = residue[oxy_name]
        direction = o.coord - p
        norm = np.linalg.norm(direction)
        if norm < 1e-8:
            continue
        direction = direction / norm
        h_coord = o.coord + 0.98 * direction
        h_name = f"H{idx}P"
        if h_name in residue:
            continue
        atom = Atom(name=h_name, coord=h_coord, bfactor=bfactor, occupancy=1.0, altloc=" ", fullname=f"{h_name:>4}", serial_number=serial_counter[0], element="H")
        residue.add(atom)
        serial_counter[0] += 1
        added.append(h_name)
    return added


def graft_one_site(
    structure,
    residue,
    template_atoms: Dict[str, TemplateAtom],
    phospho_code: str,
    clash_threshold: float,
    altloc_mode: str,
    protonation: str,
    plddt_mean: Optional[float],
    serial_counter: List[int],
    sep_seed_mode_deg: Optional[float] = None,
    ptr_seed_mode_deg: Optional[float] = None,
) -> Dict[str, object]:
    original_resname = residue.get_resname().strip().upper()
    residue_atoms = residue_atom_dict(residue, altloc_mode)
    shared_names = pick_shared_atom_names(original_resname, residue_atoms, template_atoms)
    template_xyz = np.array([template_atoms[name].coord for name in shared_names], dtype=float)
    residue_xyz = np.array([residue_atoms[name].coord for name in shared_names], dtype=float)
    rotation, translation = kabsch_fit(template_xyz, residue_xyz)
    template_only_names = [name for name in template_atoms if name not in residue_atoms]
    added_names: List[str] = []
    bfactor = float(plddt_mean) if plddt_mean is not None else 20.0
    seed_debug: Dict[str, object] = {}
    if template_only_names:
        if original_resname == "SER" and phospho_code == "SEP":
            placed_coords, seed_debug = build_sep_seed_coords(template_atoms, residue_atoms, template_only_names, mode_deg=sep_seed_mode_deg)
            for name in template_only_names:
                if name not in placed_coords:
                    continue
                coord = placed_coords[name]
                t_atom = template_atoms[name]
                atom = Atom(name=t_atom.name, coord=coord, bfactor=bfactor, occupancy=1.0, altloc=" ", fullname=f"{t_atom.name:>4}", serial_number=serial_counter[0], element=t_atom.element)
                residue.add(atom)
                added_names.append(t_atom.name)
                serial_counter[0] += 1
        elif original_resname == "THR" and phospho_code == "TPO":
            placed_coords, seed_debug = build_tpo_seed_coords(template_atoms, residue_atoms, template_only_names)
            for name in template_only_names:
                if name not in placed_coords:
                    continue
                coord = placed_coords[name]
                t_atom = template_atoms[name]
                atom = Atom(name=t_atom.name, coord=coord, bfactor=bfactor, occupancy=1.0, altloc=" ", fullname=f"{t_atom.name:>4}", serial_number=serial_counter[0], element=t_atom.element)
                residue.add(atom)
                added_names.append(t_atom.name)
                serial_counter[0] += 1
        elif original_resname == "TYR" and phospho_code == "PTR":
            placed_coords, seed_debug = build_ptr_seed_coords(template_atoms, residue_atoms, template_only_names, mode_deg=ptr_seed_mode_deg)
            for name in template_only_names:
                if name not in placed_coords:
                    continue
                coord = placed_coords[name]
                t_atom = template_atoms[name]
                atom = Atom(name=t_atom.name, coord=coord, bfactor=bfactor, occupancy=1.0, altloc=" ", fullname=f"{t_atom.name:>4}", serial_number=serial_counter[0], element=t_atom.element)
                residue.add(atom)
                added_names.append(t_atom.name)
                serial_counter[0] += 1
        else:
            transformed_xyz = transform_points(np.array([template_atoms[name].coord for name in template_only_names], dtype=float), rotation, translation)
            for name, coord in zip(template_only_names, transformed_xyz):
                t_atom = template_atoms[name]
                atom = Atom(name=t_atom.name, coord=coord, bfactor=bfactor, occupancy=1.0, altloc=" ", fullname=f"{t_atom.name:>4}", serial_number=serial_counter[0], element=t_atom.element)
                residue.add(atom)
                added_names.append(t_atom.name)
                serial_counter[0] += 1
    residue.resname = phospho_code
    seed_metrics = validate_seeded_phosphate_placement(residue, original_resname)
    proton_h_names = add_phosphate_hydrogens(residue, protonation, bfactor=bfactor, serial_counter=serial_counter)
    if proton_h_names:
        added_names.extend(proton_h_names)
    clashes = count_new_atom_clashes(structure, residue, added_names, threshold=clash_threshold)
    result = {
        "shared_atom_names": ",".join(shared_names),
        "added_atom_names": ",".join(added_names),
        "n_shared_atoms": len(shared_names),
        "n_added_atoms": len(added_names),
        "clashes_below_threshold": clashes,
    }
    result.update(seed_debug)
    result.update(seed_metrics)
    return result


def residue_atom_by_name(residue) -> Dict[str, Atom]:
    return {atom.get_name().strip(): atom for atom in residue.get_atoms()}


def residue_bonds(residue) -> List[Tuple[str, str]]:
    bonds = PHOSPHO_BONDS.get(residue.get_resname().strip().upper(), [])
    atoms = residue_atom_by_name(residue)
    result = [(a, b) for a, b in bonds if a in atoms and b in atoms]
    # C-terminal OXT: must be bonded to C to survive minimization
    if "OXT" in atoms and "C" in atoms and ("C", "OXT") not in result:
        result.append(("C", "OXT"))
    return result


def angle_triples_from_bonds(bonds: List[Tuple[str, str]]) -> List[Tuple[str, str, str]]:
    nbrs: Dict[str, Set[str]] = defaultdict(set)
    for a, b in bonds:
        nbrs[a].add(b)
        nbrs[b].add(a)
    triples: Set[Tuple[str, str, str]] = set()
    for center, neigh in nbrs.items():
        neigh = sorted(neigh)
        for i in range(len(neigh)):
            for j in range(i + 1, len(neigh)):
                triples.add((neigh[i], center, neigh[j]))
    return sorted(triples)


def _import_openmm():
    try:
        import openmm  # type: ignore
        from openmm import unit  # type: ignore
        return openmm, unit
    except Exception:
        from simtk import openmm, unit  # type: ignore
        return openmm, unit


def _mass_from_element(openmm, symbol: str) -> float:
    symbol = (symbol or "C").strip().upper()
    masses = {
        "H": 1.008,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "P": 30.974,
        "S": 32.06,
    }
    return masses.get(symbol, 12.011)


def _choose_platform(openmm, platform_name: str):
    if platform_name.lower() == "auto":
        return None
    try:
        return openmm.Platform.getPlatformByName(platform_name)
    except Exception as exc:
        raise RuntimeError(f"Requested OpenMM platform '{platform_name}' is not available: {exc}")


def openmm_local_minimize_model(
    model,
    modified_residues: Sequence,
    residue_plddt: Dict[Tuple, Optional[float]],
    reference_geometry: Optional[Dict[str, TemplateGeometry]] = None,
    flexible_neighbors: Optional[Sequence] = None,
    platform_name: str = "auto",
    max_iterations: int = 200,
    tolerance_kj_per_mol_nm: float = 10.0,
    repulsion_scale: float = 0.60,
    repulsion_k: float = 5000.0,
    restraint_k_rigid: float = 5000.0,
    restraint_k_backbone: float = 1000.0,
    restraint_k_sidechain: float = 50.0,
    restraint_k_neighbor: float = 100.0,
    bond_k: float = 20000.0,
    angle_k: float = 200.0,
) -> Dict[str, object]:
    openmm, unit = _import_openmm()
    model_atoms = list(model.get_atoms())
    atom_to_index = {atom: i for i, atom in enumerate(model_atoms)}
    system = openmm.System()
    for atom in model_atoms:
        system.addParticle(_mass_from_element(openmm, atom.element))

    # Coordinates.
    positions = [openmm.Vec3(*(atom.coord / 10.0)) for atom in model_atoms] * unit.nanometer

    modified_set = set(modified_residues)
    modified_atom_indices: Set[int] = set()
    backbone_indices: Set[int] = set()

    for residue in modified_residues:
        plddt = residue_plddt.get(residue.full_id, None)
        pscale = 1.0 if plddt is None else max(0.25, min(1.0, float(plddt) / 100.0))
        atom_map = residue_atom_by_name(residue)
        backbone_names = {"N", "CA", "C", "O", "OXT"}
        for atom in residue.get_atoms():
            idx = atom_to_index[atom]
            modified_atom_indices.add(idx)
            if atom.get_name().strip() in backbone_names:
                backbone_indices.add(idx)

    # Flexible neighbor sidechain atom indices.
    # These get soft positional restraints so nearby Arg/Lys/etc. can
    # rearrange toward or away from the phosphate during minimization.
    flex_neighbor_sc_indices: Set[int] = set()
    flex_neighbor_bb_indices: Set[int] = set()
    if flexible_neighbors:
        for residue in flexible_neighbors:
            for atom in residue.get_atoms():
                idx = atom_to_index.get(atom)
                if idx is None:
                    continue
                name = atom.get_name().strip()
                if name in BACKBONE_NAMES_RIGID:
                    flex_neighbor_bb_indices.add(idx)
                else:
                    flex_neighbor_sc_indices.add(idx)

    # Positional restraints.
    restraint = openmm.CustomExternalForce("0.5*k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")
    restraint.addPerParticleParameter("k")
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    for atom in model_atoms:
        idx = atom_to_index[atom]
        residue = atom.get_parent()
        plddt = residue_plddt.get(residue.full_id, None)
        pscale = 1.0 if plddt is None else max(0.25, min(1.0, float(plddt) / 100.0))
        if idx in backbone_indices:
            k = restraint_k_backbone * pscale
        elif idx in modified_atom_indices:
            k = restraint_k_sidechain * pscale
        elif idx in flex_neighbor_sc_indices:
            k = restraint_k_neighbor * pscale
        elif idx in flex_neighbor_bb_indices:
            k = restraint_k_backbone  # keep neighbor backbone firm
        else:
            k = restraint_k_rigid
        x0, y0, z0 = atom.coord / 10.0
        restraint.addParticle(idx, [k, float(x0), float(y0), float(z0)])
    system.addForce(restraint)

    # Harmonic bonds and angles for modified residues only.
    bond_force = openmm.HarmonicBondForce()
    angle_force = openmm.HarmonicAngleForce()
    residue_internal_pairs: Set[Tuple[int, int]] = set()

    for residue in modified_residues:
        atoms = residue_atom_by_name(residue)
        bonds = residue_bonds(residue)
        resname = residue.get_resname().strip().upper()
        geom = reference_geometry.get(resname) if reference_geometry else None

        for a, b in bonds:
            ia = atom_to_index[atoms[a]]
            ib = atom_to_index[atoms[b]]
            pair = tuple(sorted((ia, ib)))
            residue_internal_pairs.add(pair)
            # Prefer CCD ideal bond length; fall back to measured.
            ideal_length = None
            if geom is not None:
                ideal_length = geom.bond_lengths.get((a, b))
            if ideal_length is not None:
                length_nm = ideal_length / 10.0
            else:
                length_nm = float(np.linalg.norm(atoms[a].coord - atoms[b].coord) / 10.0)
            bond_force.addBond(ia, ib, length_nm * unit.nanometer, bond_k * unit.kilojoule_per_mole / unit.nanometer**2)

        for a, b, c in angle_triples_from_bonds(bonds):
            # Prefer CCD ideal angle; fall back to measured.
            ideal_angle = None
            if geom is not None:
                ideal_angle = geom.bond_angles.get((a, b, c))
            if ideal_angle is not None:
                theta = ideal_angle
            else:
                va = atoms[a].coord - atoms[b].coord
                vc = atoms[c].coord - atoms[b].coord
                cosa = np.dot(va, vc) / (np.linalg.norm(va) * np.linalg.norm(vc))
                cosa = float(np.clip(cosa, -1.0, 1.0))
                theta = math.acos(cosa)
            angle_force.addAngle(atom_to_index[atoms[a]], atom_to_index[atoms[b]], atom_to_index[atoms[c]], theta * unit.radian, angle_k * unit.kilojoule_per_mole / unit.radian**2)

        # Exclude all internal residue pairs from clash force.
        residue_atoms = [atom_to_index[a] for a in residue.get_atoms()]
        for i in range(len(residue_atoms)):
            for j in range(i + 1, len(residue_atoms)):
                residue_internal_pairs.add(tuple(sorted((residue_atoms[i], residue_atoms[j]))))

    system.addForce(bond_force)
    system.addForce(angle_force)

    # Also exclude intra-residue atom pairs for flexible neighbor residues
    # so the repulsive force doesn't act within a single residue.
    if flexible_neighbors:
        for residue in flexible_neighbors:
            res_atoms = [atom_to_index[a] for a in residue.get_atoms() if a in atom_to_index]
            for i in range(len(res_atoms)):
                for j in range(i + 1, len(res_atoms)):
                    residue_internal_pairs.add(tuple(sorted((res_atoms[i], res_atoms[j]))))

    # Local overlap penalty.  The interaction group includes both the
    # modified-residue atoms AND flexible-neighbor sidechain atoms, so
    # that the repulsive force properly captures interactions when those
    # neighbor sidechains move during minimization.
    repulse = openmm.CustomNonbondedForce("0.5*k_rep*step(r0-r)*(r-r0)^2; r0=scale*(radius1+radius2)")
    repulse.addGlobalParameter("k_rep", repulsion_k)
    repulse.addGlobalParameter("scale", repulsion_scale)
    repulse.addPerParticleParameter("radius")
    for atom in model_atoms:
        repulse.addParticle([atom_radius(atom) / 10.0])
    movable_list = sorted(modified_atom_indices | flex_neighbor_sc_indices)
    all_list = list(range(len(model_atoms)))
    repulse.addInteractionGroup(movable_list, all_list)
    for i, j in sorted(residue_internal_pairs):
        repulse.addExclusion(i, j)
    system.addForce(repulse)

    integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
    platform = _choose_platform(openmm, platform_name)
    context = openmm.Context(system, integrator, platform) if platform is not None else openmm.Context(system, integrator)
    context.setPositions(positions)
    state_before = context.getState(getEnergy=True, getPositions=True)
    energy_before = state_before.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    openmm.LocalEnergyMinimizer.minimize(context, tolerance_kj_per_mol_nm * unit.kilojoule_per_mole / unit.nanometer, max_iterations)
    state_after = context.getState(getEnergy=True, getPositions=True)
    energy_after = state_after.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    out_positions = state_after.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    for atom, pos in zip(model_atoms, out_positions):
        atom.coord = np.array(pos, dtype=float) * 10.0

    # --- Per-residue repulsive-energy decomposition ---
    # Build a lightweight evaluation system with one force per modified residue,
    # each in its own force group.  Evaluate with pre- and post-minimization
    # positions to get a per-site delta_E.
    per_residue_energies: Dict[Tuple, Dict[str, object]] = {}
    n_modified = len(modified_residues)
    if n_modified > 0 and n_modified <= 32:          # OpenMM force-group limit
        positions_before_mm = state_before.getPositions()
        positions_after_mm = state_after.getPositions()

        eval_system = openmm.System()
        for atom in model_atoms:
            eval_system.addParticle(_mass_from_element(openmm, atom.element))

        # OpenMM requires all CustomNonbondedForce objects in one System to
        # share an identical exclusion list.  Collect every intra-residue pair
        # across ALL modified residues first, then apply the full set to each.
        all_exclusion_pairs: list = []
        for residue in modified_residues:
            res_indices = sorted(atom_to_index[a] for a in residue.get_atoms())
            for i in range(len(res_indices)):
                for j in range(i + 1, len(res_indices)):
                    all_exclusion_pairs.append((res_indices[i], res_indices[j]))

        fg_map: list = []                             # [(residue, force_group_idx)]
        for fg_idx, residue in enumerate(modified_residues):
            res_force = openmm.CustomNonbondedForce(
                "0.5*k_rep*step(r0-r)*(r-r0)^2; r0=scale*(radius1+radius2)"
            )
            res_force.addGlobalParameter("k_rep", repulsion_k)
            res_force.addGlobalParameter("scale", repulsion_scale)
            res_force.addPerParticleParameter("radius")
            for atom in model_atoms:
                res_force.addParticle([atom_radius(atom) / 10.0])

            res_indices = sorted(atom_to_index[a] for a in residue.get_atoms())
            all_indices = list(range(len(model_atoms)))
            res_force.addInteractionGroup(res_indices, all_indices)

            # Apply the unified exclusion set so all forces are identical.
            for ai, aj in all_exclusion_pairs:
                res_force.addExclusion(ai, aj)

            res_force.setForceGroup(fg_idx)
            eval_system.addForce(res_force)
            fg_map.append((residue, fg_idx))

        eval_integrator = openmm.VerletIntegrator(0.001 * unit.picoseconds)
        if platform is not None:
            eval_context = openmm.Context(eval_system, eval_integrator, platform)
        else:
            eval_context = openmm.Context(eval_system, eval_integrator)

        # Evaluate with pre-minimization positions.
        eval_context.setPositions(positions_before_mm)
        eb_map: Dict[Tuple, float] = {}
        for residue, fg in fg_map:
            e = eval_context.getState(getEnergy=True, groups={fg}).getPotentialEnergy()
            eb_map[residue.full_id] = e.value_in_unit(unit.kilojoule_per_mole)

        # Evaluate with post-minimization positions.
        eval_context.setPositions(positions_after_mm)
        for residue, fg in fg_map:
            ea = eval_context.getState(getEnergy=True, groups={fg}).getPotentialEnergy()
            ea_val = ea.value_in_unit(unit.kilojoule_per_mole)
            eb_val = eb_map[residue.full_id]
            n_atoms = len(list(residue.get_atoms()))
            delta = ea_val - eb_val
            per_residue_energies[residue.full_id] = {
                "site_repulsion_energy_before": round(eb_val, 4),
                "site_repulsion_energy_after": round(ea_val, 4),
                "site_repulsion_energy_delta": round(delta, 4),
                "site_n_atoms": n_atoms,
                "site_repulsion_energy_delta_per_atom": round(delta / n_atoms, 4) if n_atoms > 0 else None,
            }

        del eval_context

    model_info = {
        "relax_status": "OK",
        "relax_backend": "openmm_local",
        "relax_energy_before": f"{energy_before:.4f}",
        "relax_energy_after": f"{energy_after:.4f}",
        "relax_rotations": "",
    }
    return model_info, per_residue_energies


def run_grafting(
    input_structure: str,
    output_structure: str,
    site_tokens: Sequence[str],
    report_path: Optional[str] = None,
    cache_dir: Optional[str] = None,
    clash_threshold: float = 2.0,
    model_indices_text: str = "0",
    altloc_mode: str = "highest_occupancy",
    prune_altlocs: bool = False,
    plddt_source: str = "bfactor",
    plddt_warn_below: float = 70.0,
    protonation: str = "0",
    relax_mode: str = "openmm_local",
    torsion_step_deg: float = 30.0,
    prescan_rotamer: bool = True,
    prescan_step_deg: float = 30.0,
    prescan_contact_weight: float = 0.5,
    prescan_basic_weight: float = 2.0,
    prescan_hbond_weight: float = 0.5,
    prescan_polar_clash_scale: float = 0.5,
    prescan_pack_weight: float = 0.0,
    prescan_geom_prior_weight: float = 2.0,
    prescan_tpo_torsion_sigma: float = 20.0,
    prescan_keep_kabsch_if_good: bool = False,
    flex_neighbor_radius: float = 6.0,
    restraint_k_neighbor: float = 100.0,
    openmm_platform: str = "auto",
    minimize_max_iterations: int = 200,
    minimize_tolerance: float = 10.0,
    repulsion_scale: float = 0.60,
    repulsion_k: float = 5000.0,
    restraint_k_rigid: float = 5000.0,
    restraint_k_backbone: float = 1000.0,
    restraint_k_sidechain: float = 50.0,
    bond_k: float = 20000.0,
    angle_k: float = 200.0,
    n_poses: int = 1,
    site_order: str = "n_to_c",
    minimization_coupling: str = "sequential",
) -> None:
    structure = parse_structure(input_structure)
    models = list(structure.get_models())
    selected_model_indices = parse_model_indices(structure, model_indices_text)
    sites = [parse_site_token(tok) for tok in site_tokens]
    template_cache: Dict[str, Dict[str, TemplateAtom]] = {}
    geometry_cache: Dict[str, TemplateGeometry] = {}
    report_rows: List[Dict[str, object]] = []
    residue_altloc_map: Dict[Tuple, Optional[str]] = {}
    serial_counter: List[int] = [next_serial(structure)]
    pose_output_paths: List[str] = []
    pose_report_rows_all: List[List[Dict[str, object]]] = []
    actual_n_poses = max(1, min(int(n_poses), 3))

    for model_index in selected_model_indices:
        model = models[model_index]
        modified_residues_for_model = []
        residue_plddt: Dict[Tuple, Optional[float]] = {}
        per_residue_rows: Dict[Tuple, Dict[str, object]] = {}

        for chain in model:
            for residue in chain:
                residue_altloc_map[residue.full_id] = get_residue_altloc_label(residue, altloc_mode)

        # Explicit experimental control over scan/minimization order.
        if site_order == "n_to_c":
            sites_sorted = sorted(sites, key=lambda s: (s.chain_id or "", s.resseq))
        elif site_order == "c_to_n":
            sites_sorted = sorted(sites, key=lambda s: (s.chain_id or "", s.resseq), reverse=True)
        else:
            sites_sorted = list(sites)

        # Phase 1: graft and prescan every site on the accumulating model.
        # Later sites therefore see all earlier rank-1 phosphates, but no site
        # is minimized until the full scan phase is complete.
        site_prescan_cache: Dict[Tuple[str, int, str], Dict[str, object]] = {}
        site_residue_cache: Dict[Tuple[str, int, str], object] = {}

        for site in sites_sorted:
            site_key = (site.chain_id or "", site.resseq, site.icode.strip())
            row: Dict[str, object] = {
                "model_index": model_index,
                "model_id": getattr(model, "id", model_index),
                "chain_id": site.chain_id or "",
                "resseq": site.resseq,
                "icode": site.icode.strip(),
                "status": "OK",
                "message": "",
                "altloc_mode": altloc_mode,
                "selected_altloc": "",
                "plddt_mean": "",
                "plddt_min": "",
                "plddt_category": "not_available",
                "plddt_qc": "not_available",
                "protonation": protonation,
                "relax_backend": relax_mode,
                "relax_status": "SKIPPED",
                "relax_energy_before": "",
                "relax_energy_after": "",
                "relax_rotations": "",
                "prescan_context": "",
                "site_order": site_order,
                "minimization_coupling": minimization_coupling,
                "applied_initial_rotations": "",
            }
            report_rows.append(row)
            try:
                residue = find_residue(model, site)
                row["selected_altloc"] = get_residue_altloc_label(residue, altloc_mode) or ""
                original_resname = residue.get_resname().strip().upper()
                row["original_resname"] = original_resname
                plddt_mean, plddt_min, plddt_cat = get_plddt_for_residue(residue, altloc_mode, plddt_source)
                row["plddt_mean"] = "" if plddt_mean is None else f"{plddt_mean:.2f}"
                row["plddt_min"] = "" if plddt_min is None else f"{plddt_min:.2f}"
                row["plddt_category"] = plddt_cat
                row["plddt_qc"] = plddt_qc_flag(plddt_mean, plddt_warn_below)
                if original_resname not in PHOSPHO_MAP:
                    raise GraftingError(f"Residue is {original_resname}, not one of SER/THR/TYR")
                phospho_code = PHOSPHO_MAP[original_resname]
                row["new_resname"] = phospho_code
                if phospho_code not in template_cache:
                    template_cache[phospho_code] = load_template_atoms(phospho_code, cache_dir=cache_dir, include_hydrogens=False)
                    geometry_cache[phospho_code] = load_template_geometry(phospho_code, template_cache[phospho_code], cache_dir=cache_dir)

                graft_info = graft_one_site(
                    structure=structure,
                    residue=residue,
                    template_atoms=template_cache[phospho_code],
                    phospho_code=phospho_code,
                    clash_threshold=clash_threshold,
                    altloc_mode=altloc_mode,
                    protonation=protonation,
                    plddt_mean=plddt_mean,
                    serial_counter=serial_counter,
                )
                row.update(graft_info)

                prescan_info: Dict[str, object] = {}
                if prescan_rotamer and relax_mode == "openmm_local":
                    prescan_info = prescan_phospho_rotamer(
                        structure,
                        residue,
                        phospho_code,
                        step_deg=prescan_step_deg,
                        contact_weight=prescan_contact_weight,
                        basic_weight=prescan_basic_weight,
                        hbond_weight=prescan_hbond_weight,
                        polar_clash_scale=prescan_polar_clash_scale,
                        pack_weight=prescan_pack_weight,
                        geom_prior_weight=prescan_geom_prior_weight,
                        tpo_torsion_sigma=prescan_tpo_torsion_sigma,
                        keep_kabsch_if_good=prescan_keep_kabsch_if_good,
                    )
                    prescan_info["prescan_context"] = f"sequential_grafted_model_premin_{site_order}"
                    row.update(prescan_info)

                site_prescan_cache[site_key] = prescan_info
                site_residue_cache[site_key] = residue
                modified_residues_for_model.append(residue)
                residue_plddt[residue.full_id] = plddt_mean
                per_residue_rows[residue.full_id] = row

            except Exception as exc:
                row["status"] = "FAILED"
                row["message"] = str(exc)

        # Phase 2: capture the all-grafted, rank-1, pre-minimization model.
        baseline_coords = {id(atom): atom.coord.copy() for atom in model.get_atoms()}

        # Phase 3: reconstruct each rank directly from its absolute prescan
        # coordinates, then minimize sites sequentially with updates retained.
        for pose_rank in range(actual_n_poses):
            for atom in model.get_atoms():
                saved = baseline_coords.get(id(atom))
                if saved is not None:
                    atom.coord = saved.copy()

            for site in sites_sorted:
                site_key = (site.chain_id or "", site.resseq, site.icode.strip())
                residue = site_residue_cache.get(site_key)
                prescan_info = site_prescan_cache.get(site_key, {})
                candidates = prescan_info.get("prescan_top_k_candidates", [])
                if residue is not None and isinstance(candidates, list) and pose_rank < len(candidates):
                    apply_prescan_candidate_coords(residue, candidates[pose_rank])

            pose_rows = copy.deepcopy(report_rows)
            pose_per_residue: Dict[Tuple, Dict[str, object]] = {}
            for pose_row in pose_rows:
                pose_row["pose_rank"] = pose_rank + 1
                row_key = (pose_row.get("chain_id", ""), pose_row.get("resseq", ""), str(pose_row.get("icode", "")).strip())
                for residue_id, original_row in per_residue_rows.items():
                    original_key = (original_row.get("chain_id", ""), original_row.get("resseq", ""), str(original_row.get("icode", "")).strip())
                    if row_key == original_key:
                        pose_per_residue[residue_id] = pose_row
                        break

            # True joint-minimization control: all successfully grafted
            # phosphoresidues and their shared flexible-neighbour set enter one
            # OpenMM minimization call for this pose rank.
            if relax_mode == "openmm_local" and minimization_coupling == "joint":
                joint_residues = []
                for residue in modified_residues_for_model:
                    joint_row = pose_per_residue.get(residue.full_id)
                    if joint_row is not None and joint_row.get("status") != "FAILED":
                        joint_residues.append(residue)
                if joint_residues:
                    joint_flex_neighbors: List = []
                    if flex_neighbor_radius > 0:
                        joint_flex_neighbors = find_flexible_neighbors(
                            model, joint_residues, radius=flex_neighbor_radius
                        )
                    nearby_joint_flex: List[str] = []
                    for neighbour in joint_flex_neighbors:
                        neighbour_id = f"{neighbour.get_parent().id}:{neighbour.id[1]}"
                        if neighbour_id not in nearby_joint_flex:
                            nearby_joint_flex.append(neighbour_id)
                    try:
                        min_info, site_energies = openmm_local_minimize_model(
                            model=model,
                            modified_residues=joint_residues,
                            residue_plddt=residue_plddt,
                            reference_geometry=geometry_cache if geometry_cache else None,
                            flexible_neighbors=joint_flex_neighbors if joint_flex_neighbors else None,
                            platform_name=openmm_platform,
                            max_iterations=minimize_max_iterations,
                            tolerance_kj_per_mol_nm=minimize_tolerance,
                            repulsion_scale=repulsion_scale,
                            repulsion_k=repulsion_k,
                            restraint_k_rigid=restraint_k_rigid,
                            restraint_k_backbone=restraint_k_backbone,
                            restraint_k_sidechain=restraint_k_sidechain,
                            restraint_k_neighbor=restraint_k_neighbor,
                            bond_k=bond_k,
                            angle_k=angle_k,
                        )
                        for residue in joint_residues:
                            joint_row = pose_per_residue.get(residue.full_id)
                            if joint_row is None:
                                continue
                            joint_row.update(min_info)
                            joint_row["flex_neighbors"] = ";".join(nearby_joint_flex)
                            joint_row["n_flex_neighbors"] = len(nearby_joint_flex)
                            joint_row["clashes_after_relax"] = count_new_atom_clashes(
                                structure,
                                residue,
                                [atom.get_name().strip() for atom in residue.get_atoms()
                                 if atom.get_name().strip() in {"P", "O1P", "O2P", "O3P", "H1P", "H2P", "H3P"}],
                                threshold=clash_threshold,
                            )
                            joint_row.update(site_energies.get(residue.full_id, {}))
                    except Exception as exc:
                        for residue in joint_residues:
                            joint_row = pose_per_residue.get(residue.full_id)
                            if joint_row is None:
                                continue
                            joint_row["relax_status"] = "FAILED"
                            joint_row["message"] = (
                                str(joint_row.get("message", ""))
                                + ("; " if joint_row.get("message") else "")
                                + f"joint minimization failed: {exc}"
                            ).strip()

            for residue in modified_residues_for_model:
                row = pose_per_residue.get(residue.full_id)
                if row is None or row.get("status") == "FAILED":
                    continue

                # Joint mode was handled once for all residues above.
                if relax_mode == "openmm_local" and minimization_coupling == "joint":
                    continue

                if relax_mode == "openmm_local":
                    this_site_residues = [residue]
                    flex_neighbors: List = []
                    if flex_neighbor_radius > 0:
                        flex_neighbors = find_flexible_neighbors(model, this_site_residues, radius=flex_neighbor_radius)
                        nearby_flex: List[str] = []
                        for neighbour in flex_neighbors:
                            neighbour_id = f"{neighbour.get_parent().id}:{neighbour.id[1]}"
                            if neighbour_id not in nearby_flex:
                                nearby_flex.append(neighbour_id)
                        row["flex_neighbors"] = ";".join(nearby_flex)
                        row["n_flex_neighbors"] = len(nearby_flex)

                    try:
                        min_info, site_energies = openmm_local_minimize_model(
                            model=model,
                            modified_residues=this_site_residues,
                            residue_plddt=residue_plddt,
                            reference_geometry=geometry_cache if geometry_cache else None,
                            flexible_neighbors=flex_neighbors if flex_neighbors else None,
                            platform_name=openmm_platform,
                            max_iterations=minimize_max_iterations,
                            tolerance_kj_per_mol_nm=minimize_tolerance,
                            repulsion_scale=repulsion_scale,
                            repulsion_k=repulsion_k,
                            restraint_k_rigid=restraint_k_rigid,
                            restraint_k_backbone=restraint_k_backbone,
                            restraint_k_sidechain=restraint_k_sidechain,
                            restraint_k_neighbor=restraint_k_neighbor,
                            bond_k=bond_k,
                            angle_k=angle_k,
                        )
                        row.update(min_info)
                        row["clashes_after_relax"] = count_new_atom_clashes(
                            structure,
                            residue,
                            [atom.get_name().strip() for atom in residue.get_atoms()
                             if atom.get_name().strip() in {"P", "O1P", "O2P", "O3P", "H1P", "H2P", "H3P"}],
                            threshold=clash_threshold,
                        )
                        row.update(site_energies.get(residue.full_id, {}))
                    except Exception as exc:
                        row["relax_status"] = "FAILED"
                        row["message"] = (
                            str(row.get("message", ""))
                            + ("; " if row.get("message") else "")
                            + f"minimization failed: {exc}"
                        ).strip()
                elif relax_mode == "torsion_scan":
                    row.update(relax_residue_by_torsion_scan(
                        structure, residue, residue.get_resname().strip().upper(), torsion_step_deg))

            pose_path = output_structure if actual_n_poses == 1 else _pose_output_path(output_structure, pose_rank + 1)
            save_structure(structure, output_path=pose_path,
                           residue_altloc_map=residue_altloc_map, prune_altlocs=prune_altlocs)
            pose_output_paths.append(pose_path)
            pose_report_rows_all.append(pose_rows)

    fieldnames = [
        "pose_rank", "model_index", "model_id", "chain_id", "resseq", "icode",
        "original_resname", "new_resname", "status", "message",
        "altloc_mode", "selected_altloc",
        "plddt_mean", "plddt_min", "plddt_category", "plddt_qc",
        "protonation", "n_shared_atoms", "n_added_atoms", "shared_atom_names", "added_atom_names",
        "clashes_below_threshold", "clashes_after_relax",
        "prescan_context", "site_order", "minimization_coupling", "applied_initial_rotations",
        "prescan_status", "prescan_score_before", "prescan_score_after", "prescan_rotations",
        "prescan_score_min", "prescan_score_max", "prescan_score_mean",
        "prescan_n_zero_clash", "prescan_n_total_steps", "prescan_zero_clash_arc",
        "n_flex_neighbors", "flex_neighbors",
        "relax_backend", "relax_status", "relax_energy_before", "relax_energy_after", "relax_rotations",
        "site_repulsion_energy_before", "site_repulsion_energy_after", "site_repulsion_energy_delta",
        "site_n_atoms", "site_repulsion_energy_delta_per_atom",
    ]
    for pose_path, pose_rows in zip(pose_output_paths, pose_report_rows_all):
        if actual_n_poses == 1:
            pose_report = report_path or str(Path(output_structure).with_suffix(Path(output_structure).suffix + ".report.tsv"))
        else:
            pose_report = str(Path(pose_path).with_suffix(".report.tsv"))
        with open(pose_report, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(pose_rows)

    import json as _json
    if report_path is None:
        report_path = str(Path(output_structure).with_suffix(Path(output_structure).suffix + ".report.tsv"))
    frames_path = str(Path(report_path).with_suffix(".prescan_frames.json"))
    frames_data = {}
    for row in report_rows:
        frames_json = row.get("prescan_frames_json", "")
        if frames_json:
            site_key = f"{row.get('chain_id', '')}:{row.get('resseq', '')}{str(row.get('icode', '')).strip()}"
            try:
                frames_data[site_key] = _json.loads(frames_json)
            except Exception:
                pass
    if frames_data:
        with open(frames_path, "w") as fh:
            fh.write(_json.dumps(frames_data, separators=(",", ":")))

    print(f"Generated {len(pose_output_paths)} independently minimized pose(s).")
    for pose_path in pose_output_paths:
        print(f"  {pose_path}")


if __name__ == "__main__":
    args = parse_args()
    run_grafting(
        input_structure=args.input_structure,
        output_structure=args.output_structure,
        site_tokens=args.sites,
        report_path=args.report,
        cache_dir=args.cache_dir,
        clash_threshold=args.clash_threshold,
        model_indices_text=args.model_indices,
        altloc_mode=args.altloc_mode,
        prune_altlocs=args.prune_altlocs,
        plddt_source=args.plddt_source,
        plddt_warn_below=args.plddt_warn_below,
        protonation=args.protonation,
        relax_mode=args.relax_mode,
        torsion_step_deg=args.torsion_step_deg,
        prescan_rotamer=args.prescan_rotamer,
        prescan_step_deg=args.prescan_step_deg,
        prescan_contact_weight=args.prescan_contact_weight,
        prescan_basic_weight=args.prescan_basic_weight,
        prescan_hbond_weight=args.prescan_hbond_weight,
        prescan_polar_clash_scale=args.prescan_polar_clash_scale,
        prescan_pack_weight=args.prescan_pack_weight,
        prescan_geom_prior_weight=args.prescan_geom_prior_weight,
        prescan_tpo_torsion_sigma=args.prescan_tpo_torsion_sigma,
        prescan_keep_kabsch_if_good=args.prescan_keep_kabsch_if_good,
        flex_neighbor_radius=args.flex_neighbor_radius,
        restraint_k_neighbor=args.restraint_k_neighbor,
        openmm_platform=args.openmm_platform,
        minimize_max_iterations=args.minimize_max_iterations,
        minimize_tolerance=args.minimize_tolerance,
        repulsion_scale=args.repulsion_scale,
        repulsion_k=args.repulsion_k,
        restraint_k_rigid=args.restraint_k_rigid,
        restraint_k_backbone=args.restraint_k_backbone,
        restraint_k_sidechain=args.restraint_k_sidechain,
        bond_k=args.bond_k,
        angle_k=args.angle_k,
        n_poses=args.n_poses,
        site_order=args.site_order,
        minimization_coupling=args.minimization_coupling,
    )
