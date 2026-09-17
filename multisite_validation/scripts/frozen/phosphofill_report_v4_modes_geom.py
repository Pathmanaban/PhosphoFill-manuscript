
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.SASA import ShrakeRupley

PHOSPHO_RESNAMES = {"SEP", "TPO", "PTR"}
PHOSPHO_ATOMS = {
    "SEP": ["P", "O1P", "O2P", "O3P", "H1P", "H2P", "H3P"],
    "TPO": ["P", "O1P", "O2P", "O3P", "H1P", "H2P", "H3P"],
    "PTR": ["P", "O1P", "O2P", "O3P", "H1P", "H2P", "H3P"],
}
PHOSPHATE_OXYGENS = {"O1P", "O2P", "O3P"}
BASIC_RES = {"LYS", "ARG", "HIS"}
ACIDIC_RES = {"ASP", "GLU"}
BASIC_ATOMS = {"NZ", "NE", "NH1", "NH2", "ND1", "NE2"}
ACIDIC_ATOMS = {"OD1", "OD2", "OE1", "OE2"}
POLAR_ATOM_INITIALS = {"N", "O", "S"}
VDW_RADII = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build per-site QC reports for one protein using unmodified and phospho-grafted structures."
    )
    p.add_argument("unmodified_structure")
    p.add_argument("modified_structure")
    p.add_argument("--sites", nargs="+", required=True, help="Examples: A:206 A:350 or 206 350 for single-chain files")
    p.add_argument("--output-prefix", required=True, help="Prefix for .json and .tsv")
    p.add_argument("--model-index", type=int, default=0)
    p.add_argument("--window-size", type=int, default=5, help="Residues on each side for local pLDDT mean and local CA RMSD")
    p.add_argument("--contact-radius", type=float, default=4.0)
    p.add_argument("--hb-radius", type=float, default=3.2, help="Heavy-atom donor/acceptor cutoff for hydrogen bonds")
    p.add_argument("--hb-min-angle", type=float, default=90.0, help="Minimum heavy-atom angle for hydrogen bonds")
    p.add_argument("--clash-threshold", type=float, default=2.0)
    p.add_argument("--sasa-probe-radius", type=float, default=1.4, help="Probe radius in Å for Shrake-Rupley SASA")
    p.add_argument("--sasa-n-points", type=int, default=960, help="Number of sphere points for Shrake-Rupley SASA")
    p.add_argument("--graft-report", default=None, help="Optional grafting TSV with site_repulsion_energy_* columns")
    p.add_argument("--graft-mode", choices=["joint", "single"], default="joint", help="Interpretation label mode for the graft run")
    return p.parse_args()


def parse_structure(path: str):
    suffix = Path(path).suffix.lower()
    parser = MMCIFParser(QUIET=True) if suffix in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    return parser.get_structure(Path(path).stem, path)


def parse_site_token(token: str) -> Tuple[Optional[str], int, str]:
    token = token.strip()
    chain_id = None
    site_part = token
    if ":" in token:
        chain_id, site_part = token.split(":", 1)
        chain_id = chain_id.strip() or None
    i = len(site_part)
    while i > 0 and site_part[i - 1].isalpha():
        i -= 1
    num_part = site_part[:i]
    icode = site_part[i:] if i < len(site_part) else " "
    if not num_part or not num_part.lstrip("-").isdigit():
        raise ValueError(f"Invalid site token: {token}")
    return chain_id, int(num_part), icode or " "


def pick_model(structure, model_index: int):
    models = list(structure.get_models())
    if model_index < 0 or model_index >= len(models):
        raise IndexError(f"Model index {model_index} out of range; found {len(models)} model(s)")
    return models[model_index]


def find_residue(model, chain_id: Optional[str], resseq: int, icode: str):
    chains = list(model.get_chains())
    if chain_id is None:
        if len(chains) != 1:
            raise ValueError(f"Site {resseq}{icode.strip()} needs an explicit chain ID because the file has {len(chains)} chains")
        chain = chains[0]
    else:
        chain = model[chain_id]
    for residue in chain.get_residues():
        hetflag, rnum, ins = residue.id
        if hetflag.strip() not in {"", "W"}:
            continue
        if rnum == resseq and (ins or " ") == (icode or " "):
            return residue
    raise KeyError(f"Residue not found: {chain.id}:{resseq}{icode.strip()}")


def get_atom_if_present(residue, atom_name: str):
    if atom_name in residue:
        atom = residue[atom_name]
        if atom.is_disordered():
            atom = atom.selected_child
        return atom
    return None


def get_all_model_atoms(model):
    atoms = []
    for atom in model.get_atoms():
        atoms.append(atom.selected_child if atom.is_disordered() else atom)
    return atoms


def plddt_class(v: Optional[float]) -> str:
    if v is None:
        return "unknown"
    if v >= 90:
        return "very_high"
    if v >= 70:
        return "confident"
    if v >= 50:
        return "low"
    return "very_low"


def residue_plddt(residue) -> Optional[float]:
    vals = [a.get_bfactor() for a in residue.get_atoms() if a.element != "H"]
    return float(np.mean(vals)) if vals else None


def local_window_plddt(chain, center_resseq: int, window: int) -> Optional[float]:
    vals = []
    for residue in chain.get_residues():
        hetflag, rnum, _ins = residue.id
        if hetflag.strip() not in {"", "W"}:
            continue
        if abs(rnum - center_resseq) <= window:
            rv = residue_plddt(residue)
            if rv is not None:
                vals.append(rv)
    return float(np.mean(vals)) if vals else None


def phosphate_anchor_atoms_for_unmodified(residue):
    mapping = {"SER": ["OG"], "THR": ["OG1"], "TYR": ["OH"]}
    out = []
    for name in mapping.get(residue.get_resname(), []):
        atom = get_atom_if_present(residue, name)
        if atom is not None:
            out.append(atom)
    return out


def phosphate_atoms(residue):
    names = PHOSPHO_ATOMS.get(residue.get_resname(), [])
    out = []
    for name in names:
        atom = get_atom_if_present(residue, name)
        if atom is not None:
            out.append(atom)
    return out


def residue_heavy_atoms(residue) -> List:
    """Return all non-hydrogen atoms for a residue."""
    out = []
    for atom in residue.get_atoms():
        real = atom.selected_child if atom.is_disordered() else atom
        elem = (real.element or real.get_name()[0]).strip().upper()
        if elem != "H":
            out.append(real)
    return out


def get_cb_or_ca(residue):
    """Return CB atom, falling back to CA for glycine."""
    atom = get_atom_if_present(residue, "CB")
    if atom is not None:
        return atom
    return get_atom_if_present(residue, "CA")


def cb_cb_contact_distances(
    site_residue,
    contact_residue_keys: Sequence[Tuple[str, int, str]],
    model,
) -> List[Dict]:
    """Compute CB–CB (or CA for GLY) distances between the site and each
    contact residue.

    Returns a list of dicts with chain, resseq, icode, resname, cb_dist.
    """
    site_cb = get_cb_or_ca(site_residue)
    if site_cb is None:
        return []
    results = []
    for chain_id, resseq, icode in contact_residue_keys:
        try:
            chain = model[chain_id]
        except KeyError:
            continue
        for residue in chain.get_residues():
            hetflag, rnum, ins = residue.id
            if hetflag.strip() not in {"", "W"}:
                continue
            if rnum == resseq and (ins or " ") == (icode or " "):
                partner_cb = get_cb_or_ca(residue)
                if partner_cb is not None:
                    d = float(np.linalg.norm(site_cb.coord - partner_cb.coord))
                    results.append({
                        "chain": chain_id,
                        "resseq": resseq,
                        "icode": icode,
                        "resname": residue.get_resname(),
                        "cb_distance": round(d, 3),
                    })
                break
    results.sort(key=lambda x: x["cb_distance"])
    return results


def cb_neighborhood_count(site_residue, model, radius: float = 8.0) -> int:
    """Count residues whose CB (or CA) is within *radius* of the site's CB."""
    return len(cb_neighborhood_residues(site_residue, model, radius))


def cb_neighborhood_residues(site_residue, model, radius: float = 8.0) -> List[Dict]:
    """Return all residues whose CB (or CA) is within *radius* of the site's CB,
    with their CB distance and residue label."""
    site_cb = get_cb_or_ca(site_residue)
    if site_cb is None:
        return []
    site_id = id(site_residue)
    results = []
    for chain in model.get_chains():
        for residue in chain.get_residues():
            if id(residue) == site_id:
                continue
            hetflag = residue.id[0].strip()
            if hetflag not in {"", "W"}:
                continue
            partner_cb = get_cb_or_ca(residue)
            if partner_cb is not None:
                d = float(np.linalg.norm(site_cb.coord - partner_cb.coord))
                if d <= radius:
                    results.append({
                        "resname": residue.get_resname(),
                        "chain": chain.id,
                        "resseq": residue.id[1],
                        "icode": (residue.id[2] or " ").strip() or None,
                        "cb_distance": round(d, 3),
                        "label": f"{residue.get_resname()} {chain.id}:{residue.id[1]}{(residue.id[2] or '').strip()}".strip(),
                    })
    results.sort(key=lambda x: x["cb_distance"])
    return results


def count_clashes_for_atoms(target_atoms, all_atoms, threshold: float) -> int:
    if not target_atoms or not all_atoms:
        return 0
    parent_residues = {id(a.get_parent()) for a in target_atoms}
    search_atoms = [a for a in all_atoms if id(a.get_parent()) not in parent_residues]
    if not search_atoms:
        return 0
    ns = NeighborSearch(search_atoms)
    clashing = set()
    for a in target_atoms:
        for b in ns.search(a.coord, threshold, level="A"):
            clashing.add(id(b))
    return len(clashing)


def nearest_contact_distance(target_atoms, all_atoms) -> Optional[float]:
    if not target_atoms or not all_atoms:
        return None
    parent_residues = {id(a.get_parent()) for a in target_atoms}
    search_atoms = [a for a in all_atoms if id(a.get_parent()) not in parent_residues]
    if not search_atoms:
        return None
    ns = NeighborSearch(search_atoms)
    best = None
    for a in target_atoms:
        nearby = ns.search(a.coord, 10.0, level="A")
        for b in nearby:
            d = float(np.linalg.norm(a.coord - b.coord))
            if best is None or d < best:
                best = d
    return best




# ---- Local anchor geometry and same-residue self-contact helpers ----
ANCHOR_GEOM = {
    "SER": {"anchor": "OG", "base": "CB", "torsions": [("torsion_N_CA_CB_OG", ["N", "CA", "CB", "OG"])]},
    "THR": {"anchor": "OG1", "base": "CB", "torsions": [("torsion_N_CA_CB_OG1", ["N", "CA", "CB", "OG1"])]},
    "TYR": {"anchor": "OH", "base": "CZ", "torsions": [("torsion_CD1_CE1_CZ_OH", ["CD1", "CE1", "CZ", "OH"])]},
    "SEP": {"anchor": "OG", "base": "CB", "torsions": [("torsion_N_CA_CB_OG", ["N", "CA", "CB", "OG"]), ("torsion_CA_CB_OG_P", ["CA", "CB", "OG", "P"])]},
    "TPO": {"anchor": "OG1", "base": "CB", "torsions": [("torsion_N_CA_CB_OG1", ["N", "CA", "CB", "OG1"]), ("torsion_CA_CB_OG1_P", ["CA", "CB", "OG1", "P"]), ("torsion_CG2_CB_OG1_P", ["CG2", "CB", "OG1", "P"])]},
    "PTR": {"anchor": "OH", "base": "CZ", "torsions": [("torsion_CD1_CE1_CZ_OH", ["CD1", "CE1", "CZ", "OH"]), ("torsion_CE1_CZ_OH_P", ["CE1", "CZ", "OH", "P"]), ("torsion_CE2_CZ_OH_P", ["CE2", "CZ", "OH", "P"])]},
}


def _distance(a, b) -> float:
    return float(np.linalg.norm(np.asarray(a, dtype=float) - np.asarray(b, dtype=float)))


def _angle(a, b, c) -> Optional[float]:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    c = np.asarray(c, dtype=float)
    v1 = a - b
    v2 = c - b
    denom = (np.linalg.norm(v1) * np.linalg.norm(v2)) + 1e-9
    if denom <= 1e-9:
        return None
    cosang = float(np.dot(v1, v2) / denom)
    return float(np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))))


def _dihedral(a, b, c, d) -> Optional[float]:
    p0 = np.asarray(a, dtype=float)
    p1 = np.asarray(b, dtype=float)
    p2 = np.asarray(c, dtype=float)
    p3 = np.asarray(d, dtype=float)
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
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


def compute_anchor_geometry(residue) -> Dict[str, Optional[float]]:
    resname = residue.get_resname().strip().upper()
    spec = ANCHOR_GEOM.get(resname)
    out: Dict[str, Optional[float]] = {
        "anchor_atom": None,
        "base_atom": None,
        "d_anchor_p": None,
        "angle_base_anchor_p": None,
    }
    if spec is None:
        return out
    anchor_name = spec["anchor"]
    base_name = spec["base"]
    out["anchor_atom"] = anchor_name
    out["base_atom"] = base_name
    anchor = get_atom_if_present(residue, anchor_name)
    base = get_atom_if_present(residue, base_name)
    p_atom = get_atom_if_present(residue, "P")
    if anchor is not None and p_atom is not None:
        out["d_anchor_p"] = round(_distance(anchor.coord, p_atom.coord), 4)
    if anchor is not None and base is not None and p_atom is not None:
        ang = _angle(base.coord, anchor.coord, p_atom.coord)
        out["angle_base_anchor_p"] = None if ang is None else round(ang, 4)
    for label, atom_names in spec.get("torsions", []):
        atoms = [get_atom_if_present(residue, nm) for nm in atom_names]
        if all(a is not None for a in atoms):
            dih = _dihedral(atoms[0].coord, atoms[1].coord, atoms[2].coord, atoms[3].coord)
            out[label] = None if dih is None else round(dih, 4)
        else:
            out[label] = None
    return out


def collect_same_residue_contacts(modified_residue, radius: float) -> List[Dict]:
    pho = phosphate_atoms(modified_residue)
    if not pho:
        return []
    resname = modified_residue.get_resname().strip().upper()
    spec = ANCHOR_GEOM.get(resname, {})
    anchor_name = spec.get("anchor")
    contacts: List[Dict] = []
    seen = set()
    for pa in pho:
        for atom in residue_heavy_atoms(modified_residue):
            an = atom.get_name().strip()
            if an in PHOSPHO_ATOMS.get(resname, []):
                continue
            if anchor_name and an == anchor_name:
                continue
            dist = float(np.linalg.norm(pa.coord - atom.coord))
            if dist <= radius:
                key = (pa.get_name(), an)
                if key in seen:
                    continue
                seen.add(key)
                contacts.append({
                    "phosphate_atom": pa.get_name(),
                    "self_atom": an,
                    "self_element": (atom.element or atom.get_name()[0]).strip().upper(),
                    "distance": round(dist, 3),
                })
    contacts.sort(key=lambda x: x["distance"])
    return contacts

def collect_contacts(target_atoms, all_atoms, radius: float) -> List[Dict]:
    if not target_atoms or not all_atoms:
        return []
    ns = NeighborSearch(all_atoms)
    seen = set()
    contacts = []
    for a in target_atoms:
        nearby = ns.search(a.coord, radius, level="A")
        for b in nearby:
            if b.get_parent() == a.get_parent():
                continue
            key = (a.get_full_id(), b.get_full_id())
            if key in seen:
                continue
            seen.add(key)
            res = b.get_parent()
            dist = float(np.linalg.norm(a.coord - b.coord))
            contacts.append(
                {
                    "phosphate_atom": a.get_name(),
                    "partner_atom": b.get_name(),
                    "partner_element": (b.element or b.get_name()[0]).strip().upper(),
                    "partner_resname": res.get_resname(),
                    "partner_chain": res.get_parent().id,
                    "partner_resseq": int(res.id[1]),
                    "partner_icode": (res.id[2] or " ").strip() or None,
                    "distance": round(dist, 3),
                }
            )
    contacts.sort(key=lambda x: x["distance"])
    return contacts


def residue_key_from_contact(contact: Dict) -> Tuple[str, int, str]:
    return (
        str(contact["partner_chain"]),
        int(contact["partner_resseq"]),
        str(contact["partner_icode"] or " "),
    )


def residue_label_from_contact(contact: Dict) -> str:
    return f"{contact['partner_resname']} {contact['partner_chain']}:{contact['partner_resseq']}{contact['partner_icode'] or ''}".strip()


def electrostatic_score_weighted(contacts: Sequence[Dict]) -> float:
    score = 0.0
    for c in contacts:
        d = c["distance"]
        if d < 0.1:
            continue
        weight = 1.0 / d
        if c["partner_resname"] in BASIC_RES and c["partner_atom"] in BASIC_ATOMS:
            score += weight
        if c["partner_resname"] in ACIDIC_RES and c["partner_atom"] in ACIDIC_ATOMS:
            score -= weight
    return round(score, 3)


def detect_salt_bridges(contacts: Sequence[Dict], max_dist: float = 4.0) -> List[Dict]:
    bridges = []
    for c in contacts:
        if (
            c["phosphate_atom"] in PHOSPHATE_OXYGENS
            and c["partner_resname"] in BASIC_RES
            and c["partner_atom"] in BASIC_ATOMS
            and c["distance"] <= max_dist
        ):
            bridges.append(c)
    return bridges


def contact_summary(contacts: Sequence[Dict]) -> Tuple[int, int, List[str], int, int]:
    basic_atom_contacts = 0
    acidic_atom_contacts = 0
    residue_labels = []
    residue_seen = set()
    basic_residues = set()
    acidic_residues = set()
    for c in contacts:
        label = residue_label_from_contact(c)
        if label not in residue_seen:
            residue_seen.add(label)
            residue_labels.append(label)
        if c["partner_resname"] in BASIC_RES and c["partner_atom"] in BASIC_ATOMS:
            basic_atom_contacts += 1
            basic_residues.add(label)
        if c["partner_resname"] in ACIDIC_RES and c["partner_atom"] in ACIDIC_ATOMS:
            acidic_atom_contacts += 1
            acidic_residues.add(label)
    return basic_atom_contacts, acidic_atom_contacts, residue_labels, len(basic_residues), len(acidic_residues)


def get_residue_by_contact(model, contact: Dict):
    chain = model[contact["partner_chain"]]
    target_icode = (contact["partner_icode"] or " ")
    for residue in chain.get_residues():
        hetflag, rnum, ins = residue.id
        if hetflag.strip() not in {"", "W"}:
            continue
        if rnum == int(contact["partner_resseq"]) and (ins or " ") == target_icode:
            return residue
    return None


def antecedent_atom_name(resname: str, atom_name: str) -> Optional[str]:
    mapping = {
        ("GLY", "N"): "CA",
        ("ALA", "N"): "CA",
        ("VAL", "N"): "CA",
        ("LEU", "N"): "CA",
        ("ILE", "N"): "CA",
        ("SER", "N"): "CA",
        ("THR", "N"): "CA",
        ("TYR", "N"): "CA",
        ("ASP", "N"): "CA",
        ("GLU", "N"): "CA",
        ("ASN", "N"): "CA",
        ("GLN", "N"): "CA",
        ("HIS", "N"): "CA",
        ("LYS", "N"): "CA",
        ("ARG", "N"): "CA",
        ("CYS", "N"): "CA",
        ("MET", "N"): "CA",
        ("PHE", "N"): "CA",
        ("TRP", "N"): "CA",
        ("PRO", "N"): "CA",
        ("O", "O"): None,
        ("SER", "OG"): "CB",
        ("THR", "OG1"): "CB",
        ("TYR", "OH"): "CZ",
        ("LYS", "NZ"): "CE",
        ("ARG", "NE"): "CD",
        ("ARG", "NH1"): "CZ",
        ("ARG", "NH2"): "CZ",
        ("HIS", "ND1"): "CG",
        ("HIS", "NE2"): "CD2",
        ("ASP", "OD1"): "CG",
        ("ASP", "OD2"): "CG",
        ("GLU", "OE1"): "CD",
        ("GLU", "OE2"): "CD",
        ("ASN", "OD1"): "CG",
        ("GLN", "OE1"): "CD",
        ("ASN", "ND2"): "CG",
        ("GLN", "NE2"): "CD",
        ("PTR", "O1P"): "P",
        ("PTR", "O2P"): "P",
        ("PTR", "O3P"): "P",
        ("SEP", "O1P"): "P",
        ("SEP", "O2P"): "P",
        ("SEP", "O3P"): "P",
        ("TPO", "O1P"): "P",
        ("TPO", "O2P"): "P",
        ("TPO", "O3P"): "P",
    }
    if atom_name == "O":
        return "C"
    return mapping.get((resname, atom_name))


def is_polar_contact(contact: Dict, hb_radius: float) -> bool:
    return contact["distance"] <= hb_radius and (contact["partner_element"][:1] in POLAR_ATOM_INITIALS)


def is_geometry_hbond(contact: Dict, model, modified_residue, max_dist: float = 3.2, min_angle: float = 90.0) -> bool:
    if contact["distance"] > max_dist:
        return False
    if contact["partner_element"][:1] not in POLAR_ATOM_INITIALS:
        return False
    partner_residue = get_residue_by_contact(model, contact)
    if partner_residue is None or contact["partner_atom"] not in partner_residue:
        return False
    partner_atom = get_atom_if_present(partner_residue, contact["partner_atom"])
    phos_atom = get_atom_if_present(modified_residue, contact["phosphate_atom"])
    if partner_atom is None or phos_atom is None:
        return False
    ant_name = antecedent_atom_name(partner_residue.get_resname(), partner_atom.get_name())
    if ant_name:
        ant_atom = get_atom_if_present(partner_residue, ant_name)
        if ant_atom is not None:
            v1 = phos_atom.coord - partner_atom.coord
            v2 = ant_atom.coord - partner_atom.coord
            denom = (np.linalg.norm(v1) * np.linalg.norm(v2)) + 1e-9
            cos_angle = float(np.dot(v1, v2) / denom)
            angle = float(np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0))))
            return angle >= min_angle
    return True


def compute_residue_sasa_map(model, probe_radius: float, n_points: int) -> Dict[Tuple[str, Tuple], float]:
    sr = ShrakeRupley(probe_radius=probe_radius, n_points=n_points)
    sr.compute(model, level="R")
    sasa_map = {}
    for chain in model.get_chains():
        for residue in chain.get_residues():
            sasa_map[(chain.id, residue.id)] = float(getattr(residue, "sasa", 0.0))
    return sasa_map


def compute_atom_sasa_map(model, probe_radius: float, n_points: int) -> Dict[int, float]:
    sr = ShrakeRupley(probe_radius=probe_radius, n_points=n_points)
    sr.compute(model, level="A")
    out: Dict[int, float] = {}
    for atom in model.get_atoms():
        real_atom = atom.selected_child if atom.is_disordered() else atom
        out[id(real_atom)] = float(getattr(real_atom, "sasa", 0.0))
    return out


def residue_sasa(sasa_map: Dict[Tuple[str, Tuple], float], residue) -> Optional[float]:
    return round(float(sasa_map.get((residue.get_parent().id, residue.id), 0.0)), 4)


def residue_shared_atom_sasa(residue, atom_sasa_map: Dict[int, float], shared_atom_names: Sequence[str]) -> float:
    total = 0.0
    for name in shared_atom_names:
        atom = get_atom_if_present(residue, name)
        if atom is not None:
            total += atom_sasa_map.get(id(atom), 0.0)
    return round(total, 4)


def neighbor_sasa_delta(
    sasa_map_before: Dict[Tuple[str, Tuple], float],
    sasa_map_after: Dict[Tuple[str, Tuple], float],
    contact_residue_ids: Sequence[Tuple[str, Tuple]],
) -> float:
    before = sum(sasa_map_before.get(rid, 0.0) for rid in contact_residue_ids)
    after = sum(sasa_map_after.get(rid, 0.0) for rid in contact_residue_ids)
    return round(after - before, 4)


def exposure_class_from_sasa(value: Optional[float]) -> str:
    if value is None:
        return "unknown"
    if value >= 50.0:
        return "exposed"
    if value >= 20.0:
        return "partially_exposed"
    return "buried"


def local_ca_rmsd(chain_un, chain_mod, center_resseq: int, window: int) -> Optional[float]:
    pairs = []
    residues_un = {r.id[1]: r for r in chain_un.get_residues() if r.id[0].strip() in {"", "W"}}
    residues_mod = {r.id[1]: r for r in chain_mod.get_residues() if r.id[0].strip() in {"", "W"}}
    for rnum, res_u in residues_un.items():
        if abs(rnum - center_resseq) > window:
            continue
        res_m = residues_mod.get(rnum)
        if res_m is None:
            continue
        ca_un = get_atom_if_present(res_u, "CA")
        ca_mod = get_atom_if_present(res_m, "CA")
        if ca_un is not None and ca_mod is not None:
            pairs.append((ca_un.coord, ca_mod.coord))
    if not pairs:
        return None
    diffs = np.array([a - b for a, b in pairs], dtype=float)
    return round(float(np.sqrt(np.mean(np.sum(diffs ** 2, axis=1)))), 4)


def phospho_quality_label(
    clash_after: Optional[int],
    strain_delta_per_atom: Optional[float],
    neighbor_delta: Optional[float],
    local_rmsd: Optional[float],
    prescan_n_zero_clash: Optional[int] = None,
    prescan_n_total_steps: Optional[int] = None,
    electrostatic_after: Optional[float] = None,
    plddt_site: Optional[float] = None,
) -> str:
    """Improved traffic-light label incorporating prescan freedom and
    electrostatic environment on top of steric/strain metrics."""
    if clash_after is None:
        return "unknown"

    strain_good = strain_delta_per_atom is None or strain_delta_per_atom < 0
    strain_bad = strain_delta_per_atom is not None and strain_delta_per_atom > 0
    neighbor_small = neighbor_delta is None or abs(neighbor_delta) <= 10.0
    rmsd_ok = local_rmsd is None or local_rmsd < 0.5
    rmsd_bad = local_rmsd is not None and local_rmsd > 0.8

    # Prescan: what fraction of orientations are clash-free?
    if prescan_n_zero_clash is not None and prescan_n_total_steps and prescan_n_total_steps > 0:
        freedom = prescan_n_zero_clash / prescan_n_total_steps
    else:
        freedom = None  # unknown

    # Electrostatic bonus/penalty
    electro_favorable = electrostatic_after is not None and electrostatic_after > 0.3

    # Low pLDDT reduces confidence — cap at amber even if steric fit is good
    low_confidence = plddt_site is not None and plddt_site < 50

    # Green: good fit, no clashes, ample rotational freedom
    if (clash_after == 0 and strain_good and neighbor_small and rmsd_ok
            and not low_confidence):
        if freedom is not None and freedom < 0.15:
            return "amber"  # narrow escape — downgrade
        return "green"

    # Red: severe problems
    if clash_after >= 3 or rmsd_bad:
        return "red"
    if strain_bad and (freedom is not None and freedom < 0.1):
        return "red"

    # Amber: everything else, but electrostatics can rescue a marginal case
    return "amber"


def phosphofill_confidence(
    plddt_site: Optional[float],
    plddt_window: Optional[float],
    clash_after: Optional[int],
    strain_delta_per_atom: Optional[float],
    prescan_n_zero_clash: Optional[int],
    prescan_n_total_steps: Optional[int],
    electrostatic_after: Optional[float],
    salt_bridges_after: int,
    local_rmsd: Optional[float],
) -> Optional[int]:
    """Compute a 0–100 PhosphoFill confidence score for a grafted site.

    Components (each 0–1, weighted):
      pLDDT trust      (0.20): input structure reliability
      Steric fit       (0.25): post-relax clash count
      Energy strain    (0.20): per-atom repulsion change
      Rotational free  (0.10): fraction of clash-free prescan orientations
      Electrostatics   (0.15): charge compatibility + salt bridges
      Structural stab  (0.10): local Cα RMSD
    """
    if clash_after is None:
        return None

    # 1. pLDDT trust
    plddt_min = min(
        plddt_site if plddt_site is not None else 100.0,
        plddt_window if plddt_window is not None else 100.0,
    )
    c_plddt = max(0.0, min(1.0, plddt_min / 100.0))

    # 2. Steric fit
    if clash_after == 0:
        c_steric = 1.0
    elif clash_after == 1:
        c_steric = 0.7
    elif clash_after == 2:
        c_steric = 0.4
    else:
        c_steric = 0.0

    # 3. Energy strain
    if strain_delta_per_atom is None or strain_delta_per_atom <= 0:
        c_strain = 1.0
    elif strain_delta_per_atom < 2.0:
        c_strain = max(0.0, 1.0 - strain_delta_per_atom / 2.0)
    else:
        c_strain = 0.0

    # 4. Rotational freedom
    if prescan_n_zero_clash is not None and prescan_n_total_steps and prescan_n_total_steps > 0:
        c_freedom = prescan_n_zero_clash / prescan_n_total_steps
    else:
        c_freedom = 0.5  # unknown — neutral

    # 5. Electrostatic compatibility (sigmoid-like)
    if electrostatic_after is None:
        c_electro = 0.5
    elif electrostatic_after > 1.0:
        c_electro = 1.0
    elif electrostatic_after > 0:
        c_electro = 0.5 + 0.5 * min(1.0, electrostatic_after / 1.0)
    elif electrostatic_after > -0.5:
        c_electro = 0.5 + electrostatic_after  # linearly from 0.5 down to 0
    else:
        c_electro = 0.0
    # Salt bridge bonus
    if salt_bridges_after > 0:
        c_electro = min(1.0, c_electro + 0.15 * salt_bridges_after)

    # 6. Structural stability
    if local_rmsd is None or local_rmsd < 0.2:
        c_stability = 1.0
    elif local_rmsd < 1.0:
        c_stability = max(0.0, 1.0 - (local_rmsd - 0.2) / 0.8)
    else:
        c_stability = 0.0

    score = (
        0.20 * c_plddt +
        0.25 * c_steric +
        0.20 * c_strain +
        0.10 * c_freedom +
        0.15 * c_electro +
        0.10 * c_stability
    )
    return round(100 * score)


def confidence_tier(score: Optional[int]) -> str:
    """Map a 0-100 confidence score to a tier label."""
    if score is None:
        return "unknown"
    if score >= 70:
        return "high"
    if score >= 40:
        return "moderate"
    return "low"


def interpretation_flags(
    plddt_site: Optional[float],
    plddt_window: Optional[float],
    clash_after: Optional[int],
    electrostatic_score_after: Optional[float],
    nearest_dist: Optional[float],
    quality_label: str,
    local_rmsd: Optional[float],
    neighbor_delta: Optional[float],
    strain_delta_per_atom: Optional[float],
) -> List[str]:
    flags: List[str] = []
    if plddt_site is not None and plddt_site < 70:
        flags.append("Low-confidence AlphaFold site")
    if plddt_window is not None and plddt_window < 70:
        flags.append("Low-confidence local neighborhood")
    if clash_after is not None:
        if clash_after == 0:
            flags.append("No residual phosphate clashes")
        elif clash_after <= 2:
            flags.append("Minor residual phosphate clashes")
        else:
            flags.append("Marked residual phosphate clashes")
    if strain_delta_per_atom is not None:
        if strain_delta_per_atom < 0:
            flags.append("Local steric repulsion improved after grafting/minimization")
        elif strain_delta_per_atom > 0:
            flags.append("Local steric repulsion worsened after grafting/minimization")
    if electrostatic_score_after is not None:
        if electrostatic_score_after > 0.6:
            flags.append("Nearby basic environment may support phosphate")
        elif electrostatic_score_after < -0.3:
            flags.append("Nearby acidic environment may disfavor phosphate")
    if nearest_dist is not None and nearest_dist < 2.5:
        flags.append("Very close phosphate contact after phosphorylation")
    if neighbor_delta is not None and abs(neighbor_delta) > 10.0:
        flags.append("Neighbor SASA changed substantially")
    if local_rmsd is not None:
        if local_rmsd < 0.3:
            flags.append("Minimal local backbone change after phosphorylation")
        elif local_rmsd <= 0.8:
            flags.append("Moderate local backbone adjustment after phosphorylation")
        else:
            flags.append("Substantial local backbone shift after phosphorylation")
    if quality_label == "green":
        flags.append("Overall phosphate accommodation is favorable")
    elif quality_label == "amber":
        flags.append("Overall phosphate accommodation is moderate")
    elif quality_label == "red":
        flags.append("Overall phosphate accommodation is poor")
    return flags


def intrinsic_site_fit_label(
    prescan_status: str,
    prescan_n_zero_clash: Optional[int],
    prescan_n_total_steps: Optional[int],
    plddt_site: Optional[float],
) -> str:
    """Label the site's intrinsic phosphate compatibility from the independent
    prescan performed on the unmodified input model."""
    if not prescan_status or prescan_status.upper() != "OK":
        return "unknown"
    freedom = None
    if prescan_n_zero_clash is not None and prescan_n_total_steps:
        if prescan_n_total_steps > 0:
            freedom = prescan_n_zero_clash / prescan_n_total_steps
    if freedom is None:
        return "unknown"
    low_confidence = plddt_site is not None and plddt_site < 50.0
    if freedom >= 0.40 and not low_confidence:
        return "favorable"
    if freedom >= 0.15:
        return "moderate"
    return "restricted"


def contextual_site_fit_label(quality_label: str) -> str:
    return {
        "green": "favorable",
        "amber": "moderate",
        "red": "poor",
    }.get(str(quality_label or "").strip().lower(), "unknown")


def overall_interpretation_label(
    plddt_site: Optional[float],
    plddt_window: Optional[float],
    quality_label: str,
    electrostatic_score_after: Optional[float],
    salt_bridges_after: int,
    prescan_context: str = "",
    intrinsic_fit: str = "unknown",
    graft_mode: str = "joint",
) -> str:
    low_conf = ((plddt_site is not None and plddt_site < 70) or (plddt_window is not None and plddt_window < 70))
    quality_text = {"green": "well accommodated", "amber": "partially accommodated", "red": "poorly accommodated"}.get(
        quality_label, "of uncertain accommodation"
    )
    if prescan_context == "independent_unmodified_model":
        prefix = f"Independent prescan suggests {intrinsic_fit} intrinsic site fit; "
    else:
        prefix = ""
    mode_text = "after single-site contextual minimization" if str(graft_mode).strip().lower() == "single" else "after joint contextual minimization"
    if low_conf:
        return prefix + f"low-confidence site with phosphate {quality_text} {mode_text}"
    if salt_bridges_after > 0:
        return prefix + f"confident site with phosphate {quality_text} {mode_text} and potential salt-bridge stabilization"
    if electrostatic_score_after is not None and electrostatic_score_after > 0.6:
        return prefix + f"confident site with phosphate {quality_text} {mode_text} in a favorable local charge environment"
    if quality_label == "red":
        return prefix + f"confident site with poor local phosphate accommodation {mode_text}"
    if quality_label == "amber":
        return prefix + f"confident site with moderate local phosphate accommodation {mode_text}"
    return prefix + f"confident site with favorable local phosphate accommodation {mode_text}"


def normalize_icode(value) -> str:
    text = "" if value is None else str(value).strip()
    return text if text else " "


def parse_optional_float(value) -> Optional[float]:
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    return float(text)


def parse_optional_int(value) -> Optional[int]:
    if value is None:
        return None
    text = str(value).strip()
    if text == "":
        return None
    return int(float(text))


def load_graft_report(path: Optional[str]) -> Dict[Tuple[str, int, str], Dict]:
    if not path:
        return {}
    out: Dict[Tuple[str, int, str], Dict] = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            chain_id = str(row.get("chain_id", "")).strip()
            if not chain_id:
                continue
            resseq = parse_optional_int(row.get("resseq"))
            if resseq is None:
                continue
            icode = normalize_icode(row.get("icode"))
            key = (chain_id, resseq, icode)
            out[key] = {
                "graft_status": row.get("status"),
                "graft_message": row.get("message"),
                "graft_clashes_after_relax": parse_optional_int(row.get("clashes_after_relax")),
                "site_repulsion_energy_before": parse_optional_float(row.get("site_repulsion_energy_before")),
                "site_repulsion_energy_after": parse_optional_float(row.get("site_repulsion_energy_after")),
                "site_repulsion_energy_delta": parse_optional_float(row.get("site_repulsion_energy_delta")),
                "site_n_atoms": parse_optional_int(row.get("site_n_atoms")),
                "site_repulsion_energy_delta_per_atom": parse_optional_float(row.get("site_repulsion_energy_delta_per_atom")),
                "prescan_context": row.get("prescan_context", ""),
                "applied_initial_rotations": row.get("applied_initial_rotations", ""),
                "prescan_status": row.get("prescan_status", ""),
                "prescan_score_before": parse_optional_float(row.get("prescan_score_before")),
                "prescan_score_after": parse_optional_float(row.get("prescan_score_after")),
                "prescan_rotations": row.get("prescan_rotations", ""),
                "prescan_score_min": parse_optional_float(row.get("prescan_score_min")),
                "prescan_score_max": parse_optional_float(row.get("prescan_score_max")),
                "prescan_score_mean": parse_optional_float(row.get("prescan_score_mean")),
                "prescan_n_zero_clash": parse_optional_int(row.get("prescan_n_zero_clash")),
                "prescan_n_total_steps": parse_optional_int(row.get("prescan_n_total_steps")),
                "prescan_zero_clash_arc": parse_optional_float(row.get("prescan_zero_clash_arc")),
                "n_flex_neighbors": parse_optional_int(row.get("n_flex_neighbors")),
                "flex_neighbors": row.get("flex_neighbors", ""),
            }
    return out


def write_reports(prefix: str, records: Sequence[Dict]) -> Tuple[Path, Path]:
    json_path = Path(prefix + ".json")
    tsv_path = Path(prefix + ".tsv")
    json_path.write_text(json.dumps(records, indent=2), encoding="utf-8")
    flat_fields = [
        "site_label", "chain_id", "residue_number", "insertion_code",
        "residue_name_original", "residue_name_modified",
        "anchor_atom_original", "base_atom_original", "d_anchor_p_original", "angle_base_anchor_p_original",
        "anchor_atom_modified", "base_atom_modified", "d_anchor_p_modified", "angle_base_anchor_p_modified",
        "torsion_N_CA_CB_OG", "torsion_CA_CB_OG_P",
        "torsion_N_CA_CB_OG1", "torsion_CA_CB_OG1_P", "torsion_CG2_CB_OG1_P",
        "torsion_CD1_CE1_CZ_OH", "torsion_CE1_CZ_OH_P", "torsion_CE2_CZ_OH_P",
        "same_residue_contact_count", "nearest_same_residue_contact_distance", "same_residue_contacts",
        "phosphofill_confidence", "confidence_tier",
        "plddt_site", "plddt_window_mean", "plddt_class",
        "clash_count_before", "clash_count_after",
        "nearest_contact_distance_before", "nearest_contact_distance_after",
        "nearby_basic_count_before", "nearby_basic_count_after",
        "nearby_acidic_count_before", "nearby_acidic_count_after",
        "electrostatic_score_weighted_before", "electrostatic_score_weighted_after",
        "salt_bridge_count_before", "salt_bridge_count_after",
        "salt_bridge_details_before", "salt_bridge_details_after",
        "polar_contact_count_before", "polar_contact_count_after",
        "hydrogen_bond_count_before", "hydrogen_bond_count_after",
        "sasa_unmodified", "sasa_modified", "sasa_delta",
        "sasa_shared_atoms_before", "sasa_shared_atoms_after",
        "neighbor_sasa_delta", "local_ca_rmsd",
        "cb_neighborhood_before", "cb_neighborhood_after",
        "cb_neighborhood_residues_before", "cb_neighborhood_residues_after",
        "cb_contact_distances_before", "cb_contact_distances_after",
        "site_repulsion_energy_before", "site_repulsion_energy_after",
        "site_repulsion_energy_delta", "site_repulsion_energy_delta_per_atom",
        "graft_mode", "prescan_context", "applied_initial_rotations", "optimization_strategy",
        "intrinsic_site_fit", "contextual_site_fit",
        "prescan_status", "prescan_score_before", "prescan_score_after", "prescan_rotations",
        "prescan_score_min", "prescan_score_max", "prescan_score_mean",
        "prescan_n_zero_clash", "prescan_n_total_steps", "prescan_zero_clash_arc",
        "n_flex_neighbors", "flex_neighbors",
        "pts_clash_penalty", "pts_energy_term", "pts_sasa_term",
        "phospho_quality_label",
        "exposure_class_unmodified", "exposure_class_modified",
        "overall_interpretation", "interpretation_flags",
        "nearby_contact_residues_before", "nearby_contact_residues_after",
        # Visualization-compatible aliases (counts under the names the vis JS expects)
        "polar_contacts_before", "polar_contacts_after",
        "hydrogen_bonds_before", "hydrogen_bonds_after",
        "salt_bridges_before", "salt_bridges_after",
    ]
    with tsv_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=flat_fields, delimiter="\t")
        writer.writeheader()
        for rec in records:
            row = {k: rec.get(k) for k in flat_fields}
            row["interpretation_flags"] = ";".join(rec.get("interpretation_flags", []))
            row["nearby_contact_residues_before"] = ";".join(rec.get("nearby_contact_residues_before", []))
            row["nearby_contact_residues_after"] = ";".join(rec.get("nearby_contact_residues_after", []))
            row["same_residue_contacts"] = ";".join(
                f"{c['phosphate_atom']}-{c['self_atom']}:{c['distance']}A" for c in rec.get("same_residue_contacts", [])
            )
            # Populate vis-compatible aliases with the count values
            row["polar_contacts_before"] = rec.get("polar_contact_count_before", 0)
            row["polar_contacts_after"] = rec.get("polar_contact_count_after", 0)
            row["hydrogen_bonds_before"] = rec.get("hydrogen_bond_count_before", 0)
            row["hydrogen_bonds_after"] = rec.get("hydrogen_bond_count_after", 0)
            row["salt_bridges_before"] = rec.get("salt_bridge_count_before", 0)
            row["salt_bridges_after"] = rec.get("salt_bridge_count_after", 0)
            writer.writerow(row)
    return json_path, tsv_path


def main():
    args = parse_args()
    unmod = parse_structure(args.unmodified_structure)
    mod = parse_structure(args.modified_structure)
    model_un = pick_model(unmod, args.model_index)
    model_mod = pick_model(mod, args.model_index)
    graft_map = load_graft_report(args.graft_report)

    all_atoms_un = get_all_model_atoms(model_un)
    all_atoms_mod = get_all_model_atoms(model_mod)
    sasa_map_un = compute_residue_sasa_map(model_un, args.sasa_probe_radius, args.sasa_n_points)
    sasa_map_mod = compute_residue_sasa_map(model_mod, args.sasa_probe_radius, args.sasa_n_points)
    atom_sasa_un = compute_atom_sasa_map(model_un, args.sasa_probe_radius, args.sasa_n_points)
    atom_sasa_mod = compute_atom_sasa_map(model_mod, args.sasa_probe_radius, args.sasa_n_points)

    records = []
    for token in args.sites:
        chain_id, resseq, icode = parse_site_token(token)
        res_un = find_residue(model_un, chain_id, resseq, icode)
        res_mod = find_residue(model_mod, chain_id, resseq, icode)
        pho_atoms = phosphate_atoms(res_mod)
        if not pho_atoms:
            raise ValueError(f"Modified residue at {token} does not contain phosphate atoms; found {res_mod.get_resname()}")
        chain_un = res_un.get_parent()
        chain_mod = res_mod.get_parent()

        geom_un = compute_anchor_geometry(res_un)
        geom_mod = compute_anchor_geometry(res_mod)
        self_contacts = collect_same_residue_contacts(res_mod, args.contact_radius)
        nearest_self_contact = self_contacts[0]["distance"] if self_contacts else None

        plddt_site_un = residue_plddt(res_un)
        plddt_window_un = local_window_plddt(chain_un, resseq, args.window_size)

        anchor_atoms_un = phosphate_anchor_atoms_for_unmodified(res_un)
        clash_before = count_clashes_for_atoms(anchor_atoms_un, all_atoms_un, args.clash_threshold) if anchor_atoms_un else None
        clash_after = count_clashes_for_atoms(pho_atoms, all_atoms_mod, args.clash_threshold)
        nearest_dist_before = nearest_contact_distance(anchor_atoms_un, all_atoms_un) if anchor_atoms_un else None
        nearest_dist_after = nearest_contact_distance(pho_atoms, all_atoms_mod)

        # Contact detection: use ALL heavy atoms of the residue in both
        # states for a fair before/after comparison.  The unmodified SER
        # has ~6 heavy atoms, the modified SEP has ~10 — the shared atoms
        # search the same space, and any difference comes only from the
        # added phosphate atoms.
        site_heavy_un = residue_heavy_atoms(res_un)
        site_heavy_mod = residue_heavy_atoms(res_mod)
        contacts_before = collect_contacts(site_heavy_un, all_atoms_un, args.contact_radius)
        contacts_after = collect_contacts(site_heavy_mod, all_atoms_mod, args.contact_radius)
        basic_count_before, acidic_count_before, contact_residues_before, _basic_res_before, _acidic_res_before = contact_summary(contacts_before)
        basic_count_after, acidic_count_after, contact_residues_after, _basic_res_after, _acidic_res_after = contact_summary(contacts_after)

        polar_before = [c for c in contacts_before if is_polar_contact(c, args.hb_radius)]
        polar_after = [c for c in contacts_after if is_polar_contact(c, args.hb_radius)]
        hbonds_before = [c for c in polar_before if is_geometry_hbond(c, model_un, res_un, args.hb_radius, args.hb_min_angle)]
        hbonds_after = [c for c in polar_after if is_geometry_hbond(c, model_mod, res_mod, args.hb_radius, args.hb_min_angle)]
        salt_bridges_before = detect_salt_bridges(contacts_before)
        salt_bridges_after = detect_salt_bridges(contacts_after)

        electro_before = electrostatic_score_weighted(contacts_before)
        electro_after = electrostatic_score_weighted(contacts_after)

        sasa_un = residue_sasa(sasa_map_un, res_un)
        sasa_mod = residue_sasa(sasa_map_mod, res_mod)
        sasa_delta = None if sasa_un is None or sasa_mod is None else round(sasa_mod - sasa_un, 4)

        shared_atom_names = sorted(set(atom.get_name() for atom in res_un.get_atoms()) & set(atom.get_name() for atom in res_mod.get_atoms()))
        sasa_shared_before = residue_shared_atom_sasa(res_un, atom_sasa_un, shared_atom_names)
        sasa_shared_after = residue_shared_atom_sasa(res_mod, atom_sasa_mod, shared_atom_names)

        contact_ids = []
        seen_ids = set()
        for c in contacts_before + contacts_after:
            rid = residue_key_from_contact(c)
            key = (rid[0], (" ", rid[1], rid[2]))
            if key not in seen_ids:
                seen_ids.add(key)
                contact_ids.append(key)
        neigh_delta = neighbor_sasa_delta(sasa_map_un, sasa_map_mod, contact_ids)
        local_rmsd = local_ca_rmsd(chain_un, chain_mod, resseq, args.window_size)

        # CB-level metrics: backbone/sidechain-root distances independent
        # of rotamer state.  Gives a stable reference frame to distinguish
        # "sidechain rearrangement" from "backbone shift".
        all_contact_keys = []
        for c in contacts_before + contacts_after:
            rid = residue_key_from_contact(c)
            k = (rid[0], rid[1], rid[2])
            if k not in [x for x in all_contact_keys]:
                all_contact_keys.append(k)
        # Deduplicate
        _seen_cb_keys: set = set()
        unique_contact_keys = []
        for k in all_contact_keys:
            if k not in _seen_cb_keys:
                _seen_cb_keys.add(k)
                unique_contact_keys.append(k)

        cb_dists_before = cb_cb_contact_distances(res_un, unique_contact_keys, model_un)
        cb_dists_after = cb_cb_contact_distances(res_mod, unique_contact_keys, model_mod)
        cb_neigh_list_before = cb_neighborhood_residues(res_un, model_un, radius=8.0)
        cb_neigh_list_after = cb_neighborhood_residues(res_mod, model_mod, radius=8.0)
        cb_neighborhood_before = len(cb_neigh_list_before)
        cb_neighborhood_after = len(cb_neigh_list_after)

        graft_key = (res_un.get_parent().id, resseq, normalize_icode(icode))
        graft_info = graft_map.get(graft_key, {})
        site_repulsion_before = graft_info.get("site_repulsion_energy_before")
        site_repulsion_after = graft_info.get("site_repulsion_energy_after")
        site_repulsion_delta = graft_info.get("site_repulsion_energy_delta")
        site_repulsion_delta_per_atom = graft_info.get("site_repulsion_energy_delta_per_atom")
        if graft_info.get("graft_clashes_after_relax") is not None:
            clash_after = graft_info["graft_clashes_after_relax"]

        quality = phospho_quality_label(
            clash_after, site_repulsion_delta_per_atom, neigh_delta, local_rmsd,
            prescan_n_zero_clash=graft_info.get("prescan_n_zero_clash"),
            prescan_n_total_steps=graft_info.get("prescan_n_total_steps"),
            electrostatic_after=electro_after,
            plddt_site=plddt_site_un,
        )
        intrinsic_fit = intrinsic_site_fit_label(
            graft_info.get("prescan_status", ""),
            graft_info.get("prescan_n_zero_clash"),
            graft_info.get("prescan_n_total_steps"),
            plddt_site_un,
        )
        contextual_fit = contextual_site_fit_label(quality)

        confidence = phosphofill_confidence(
            plddt_site=plddt_site_un,
            plddt_window=plddt_window_un,
            clash_after=clash_after,
            strain_delta_per_atom=site_repulsion_delta_per_atom,
            prescan_n_zero_clash=graft_info.get("prescan_n_zero_clash"),
            prescan_n_total_steps=graft_info.get("prescan_n_total_steps"),
            electrostatic_after=electro_after,
            salt_bridges_after=len(salt_bridges_after),
            local_rmsd=local_rmsd,
        )

        rec = {
            "site_label": f"{res_un.get_parent().id}:{resseq}{icode.strip()}",
            "chain_id": res_un.get_parent().id,
            "residue_number": resseq,
            "insertion_code": icode.strip() or None,
            "residue_name_original": res_un.get_resname(),
            "residue_name_modified": res_mod.get_resname(),
            "anchor_atom_original": geom_un.get("anchor_atom"),
            "base_atom_original": geom_un.get("base_atom"),
            "d_anchor_p_original": geom_un.get("d_anchor_p"),
            "angle_base_anchor_p_original": geom_un.get("angle_base_anchor_p"),
            "anchor_atom_modified": geom_mod.get("anchor_atom"),
            "base_atom_modified": geom_mod.get("base_atom"),
            "d_anchor_p_modified": geom_mod.get("d_anchor_p"),
            "angle_base_anchor_p_modified": geom_mod.get("angle_base_anchor_p"),
            "torsion_N_CA_CB_OG": geom_mod.get("torsion_N_CA_CB_OG"),
            "torsion_CA_CB_OG_P": geom_mod.get("torsion_CA_CB_OG_P"),
            "torsion_N_CA_CB_OG1": geom_mod.get("torsion_N_CA_CB_OG1"),
            "torsion_CA_CB_OG1_P": geom_mod.get("torsion_CA_CB_OG1_P"),
            "torsion_CG2_CB_OG1_P": geom_mod.get("torsion_CG2_CB_OG1_P"),
            "torsion_CD1_CE1_CZ_OH": geom_mod.get("torsion_CD1_CE1_CZ_OH"),
            "torsion_CE1_CZ_OH_P": geom_mod.get("torsion_CE1_CZ_OH_P"),
            "torsion_CE2_CZ_OH_P": geom_mod.get("torsion_CE2_CZ_OH_P"),
            "same_residue_contact_count": len(self_contacts),
            "nearest_same_residue_contact_distance": nearest_self_contact,
            "same_residue_contacts": self_contacts,
            "plddt_site": None if plddt_site_un is None else round(plddt_site_un, 2),
            "plddt_window_mean": None if plddt_window_un is None else round(plddt_window_un, 2),
            "plddt_class": plddt_class(plddt_site_un),
            "clash_count_before": clash_before,
            "clash_count_after": clash_after,
            "nearest_contact_distance_before": None if nearest_dist_before is None else round(nearest_dist_before, 3),
            "nearest_contact_distance_after": None if nearest_dist_after is None else round(nearest_dist_after, 3),
            "nearby_basic_count_before": basic_count_before,
            "nearby_basic_count_after": basic_count_after,
            "nearby_acidic_count_before": acidic_count_before,
            "nearby_acidic_count_after": acidic_count_after,
            "nearby_contact_residues_before": contact_residues_before,
            "nearby_contact_residues_after": contact_residues_after,
            "electrostatic_score_weighted_before": electro_before,
            "electrostatic_score_weighted_after": electro_after,
            "salt_bridges_before": salt_bridges_before,
            "salt_bridges_after": salt_bridges_after,
            "salt_bridge_count_before": len(salt_bridges_before),
            "salt_bridge_count_after": len(salt_bridges_after),
            "salt_bridge_details_before": ";".join(
                f"{c['partner_resname']} {c['partner_chain']}:{c['partner_resseq']} {c['partner_atom']} {c['distance']:.2f}A"
                for c in salt_bridges_before
            ),
            "salt_bridge_details_after": ";".join(
                f"{c['partner_resname']} {c['partner_chain']}:{c['partner_resseq']} {c['partner_atom']} {c['distance']:.2f}A"
                for c in salt_bridges_after
            ),
            "polar_contacts_before": polar_before,
            "polar_contacts_after": polar_after,
            "polar_contact_count_before": len(polar_before),
            "polar_contact_count_after": len(polar_after),
            "hydrogen_bonds_before": hbonds_before,
            "hydrogen_bonds_after": hbonds_after,
            "hydrogen_bond_count_before": len(hbonds_before),
            "hydrogen_bond_count_after": len(hbonds_after),
            "contacts_before": contacts_before,
            "contacts_after": contacts_after,
            "sasa_unmodified": sasa_un,
            "sasa_modified": sasa_mod,
            "sasa_delta": sasa_delta,
            "sasa_shared_atoms_before": sasa_shared_before,
            "sasa_shared_atoms_after": sasa_shared_after,
            "shared_atom_names": shared_atom_names,
            "neighbor_sasa_delta": neigh_delta,
            "local_ca_rmsd": local_rmsd,
            "cb_neighborhood_before": cb_neighborhood_before,
            "cb_neighborhood_after": cb_neighborhood_after,
            "cb_neighborhood_residues_before": ";".join(
                f"{d['label']} {d['cb_distance']}A"
                for d in cb_neigh_list_before
            ),
            "cb_neighborhood_residues_after": ";".join(
                f"{d['label']} {d['cb_distance']}A"
                for d in cb_neigh_list_after
            ),
            "cb_contact_distances_before": ";".join(
                f"{d['resname']} {d['chain']}:{d['resseq']} {d['cb_distance']}A"
                for d in cb_dists_before
            ),
            "cb_contact_distances_after": ";".join(
                f"{d['resname']} {d['chain']}:{d['resseq']} {d['cb_distance']}A"
                for d in cb_dists_after
            ),
            "cb_dists_before_raw": cb_dists_before,
            "cb_dists_after_raw": cb_dists_after,
            "site_repulsion_energy_before": site_repulsion_before,
            "site_repulsion_energy_after": site_repulsion_after,
            "site_repulsion_energy_delta": site_repulsion_delta,
            "site_repulsion_energy_delta_per_atom": site_repulsion_delta_per_atom,
            "graft_mode": args.graft_mode,
            "prescan_context": graft_info.get("prescan_context", ""),
            "applied_initial_rotations": graft_info.get("applied_initial_rotations", ""),
            "optimization_strategy": "sequential_grafted_prescan_then_sequential_site_local_minimization_per_pose",
            "prescan_status": graft_info.get("prescan_status", ""),
            "prescan_score_before": graft_info.get("prescan_score_before"),
            "prescan_score_after": graft_info.get("prescan_score_after"),
            "prescan_rotations": graft_info.get("prescan_rotations", ""),
            "intrinsic_site_fit": intrinsic_fit,
            "contextual_site_fit": contextual_fit,
            "n_flex_neighbors": graft_info.get("n_flex_neighbors"),
            "flex_neighbors": graft_info.get("flex_neighbors", ""),
            "pts_clash_penalty": clash_after,
            "pts_energy_term": site_repulsion_delta_per_atom,
            "pts_sasa_term": neigh_delta,
            "phospho_quality_label": quality,
            "phosphofill_confidence": confidence,
            "confidence_tier": confidence_tier(confidence),
            "prescan_score_min": graft_info.get("prescan_score_min"),
            "prescan_score_max": graft_info.get("prescan_score_max"),
            "prescan_score_mean": graft_info.get("prescan_score_mean"),
            "prescan_n_zero_clash": graft_info.get("prescan_n_zero_clash"),
            "prescan_n_total_steps": graft_info.get("prescan_n_total_steps"),
            "prescan_zero_clash_arc": graft_info.get("prescan_zero_clash_arc"),
            "exposure_class_unmodified": exposure_class_from_sasa(sasa_un),
            "exposure_class_modified": exposure_class_from_sasa(sasa_mod),
            "graft_info_found": bool(graft_info),
        }
        rec["interpretation_flags"] = interpretation_flags(
            rec["plddt_site"],
            rec["plddt_window_mean"],
            rec["clash_count_after"],
            rec["electrostatic_score_weighted_after"],
            rec["nearest_contact_distance_after"],
            rec["phospho_quality_label"],
            rec["local_ca_rmsd"],
            rec["neighbor_sasa_delta"],
            rec["site_repulsion_energy_delta_per_atom"],
        )
        rec["overall_interpretation"] = overall_interpretation_label(
            rec["plddt_site"],
            rec["plddt_window_mean"],
            rec["phospho_quality_label"],
            rec["electrostatic_score_weighted_after"],
            rec["salt_bridge_count_after"],
            prescan_context=rec.get("prescan_context", ""),
            intrinsic_fit=rec.get("intrinsic_site_fit", "unknown"),
            graft_mode=args.graft_mode,
        )
        records.append(rec)

    json_path, tsv_path = write_reports(args.output_prefix, records)
    print(f"Wrote JSON report: {json_path}")
    print(f"Wrote TSV report: {tsv_path}")


if __name__ == "__main__":
    main()
