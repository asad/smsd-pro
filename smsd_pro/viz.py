# SPDX-License-Identifier: Apache-2.0
# © 2025 BioInception PVT LTD.
# smsd_pro/viz.py – side-by-side PNG comparison of SMSD vs RDKit (with diff coloring)
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict, Iterable, Set
from io import BytesIO

from rdkit import Chem
from rdkit.Chem import Draw

from .engines import SMSD, SubstructureOptions, MCSOptions, MatchResult
from .chem import ChemOptions, SmartsPattern

@dataclass
class VizConfig:
    save_png: Optional[str] = None
    subimg_size: Tuple[int, int] = (500, 400)
    title_smsd: str = "SMSD"
    title_rdkit: str = "RDKit"
    # highlight colours (r,g,b) in 0..1
    common_color: Tuple[float,float,float] = (0.95, 0.4, 0.4)  # pink/red
    diff_color:   Tuple[float,float,float] = (0.2,  0.6, 1.0)  # blue

def _bond_indices_from_mapping(q: Chem.Mol, t: Chem.Mol, m: Dict[int,int]) -> List[int]:
    idxs: List[int] = []
    for b in q.GetBonds():
        i, j = int(b.GetBeginAtomIdx()), int(b.GetEndAtomIdx())
        if i in m and j in m:
            tb = t.GetBondBetweenAtoms(int(m[i]), int(m[j]))
            if tb is not None:
                idxs.append(int(tb.GetIdx()))
    return idxs

def _bond_indices_from_atomset(m: Chem.Mol, atoms: Iterable[int]) -> List[int]:
    aset: Set[int] = set(int(a) for a in atoms)
    idxs: List[int] = []
    for b in m.GetBonds():
        i, j = int(b.GetBeginAtomIdx()), int(b.GetEndAtomIdx())
        if i in aset and j in aset:
            idxs.append(int(b.GetIdx()))
    return idxs

def _to_pil(image_or_svg):
    try:
        from PIL import Image
    except Exception as e:
        raise RuntimeError("Pillow is required for saving PNGs.") from e
    if hasattr(image_or_svg, "save"):
        return image_or_svg
    if isinstance(image_or_svg, (bytes, bytearray)):
        return Image.open(BytesIO(image_or_svg))
    if hasattr(image_or_svg, "data"):
        return Image.open(BytesIO(image_or_svg.data))
    raise TypeError("Unsupported image object from RDKit.")

def _draw_pair_grid(mols: List[Chem.Mol],
                    legends: List[str],
                    highlightAtoms: List[List[int]],
                    highlightBonds: List[List[int]],
                    atomColorMaps: Optional[List[Dict[int,Tuple[float,float,float]]]] = None,
                    bondColorMaps: Optional[List[Dict[int,Tuple[float,float,float]]]] = None,
                    size: Tuple[int, int] = (500, 400)):
    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=2,
        subImgSize=size,
        legends=legends,
        highlightAtomLists=highlightAtoms,
        highlightBondLists=highlightBonds,
        highlightAtomColors=atomColorMaps,
        highlightBondColors=bondColorMaps,
        useSVG=False
    )
    return _to_pil(img)

# -------- Substructure (target vs target panels) --------
def viz_compare_substructure_png(query: str|Chem.Mol,
                                 target: str|Chem.Mol,
                                 *,
                                 smsd_chem: ChemOptions = ChemOptions(),
                                 smsd_opt: SubstructureOptions = SubstructureOptions(),
                                 mapping_index: int = 0,
                                 cfg: VizConfig = VizConfig()):
    """Render two panels of the *target* molecule:
    - Left: SMSD mapping
    - Right: RDKit mapping
    We colour atoms/bonds common to both mappings in cfg.common_color and
    atoms/bonds that are unique to the panel in cfg.diff_color.
    If no mapping exists for a panel, it is left unhighlighted.
    """
    smsd = SMSD(query, target, chem=smsd_chem)
    maps = smsd.substructure_all(smsd_opt)
    m_query = smsd.q
    m_target = smsd.t

    # SMSD mapping on target
    sms_map = maps[min(mapping_index, len(maps)-1)].mapping if maps else {}
    sms_atoms_t = sorted(sms_map.values())
    sms_bonds_t = _bond_indices_from_mapping(m_query, m_target, sms_map)

    # RDKit reference mapping (strip $name so RDKit can parse SMARTS)
    if isinstance(query, str):
        q_try = Chem.MolFromSmiles(query)
        if q_try is None:
            q_try = Chem.MolFromSmarts(SmartsPattern.strip_named_predicates(query))
    else:
        q_try = query
    t_try = Chem.MolFromSmiles(target) if isinstance(target, str) else target
    rd_matches = list(t_try.GetSubstructMatches(q_try, uniquify=True)) if q_try and t_try else []
    rdk_atoms_t = list(rd_matches[0]) if rd_matches else []
    rdk_bonds_t = _bond_indices_from_atomset(t_try, rdk_atoms_t) if rd_matches else []

    # Colour maps
    inter_atoms = set(sms_atoms_t).intersection(rdk_atoms_t)
    sa_only = set(sms_atoms_t) - inter_atoms
    rk_only = set(rdk_atoms_t) - inter_atoms

    inter_bonds = set(sms_bonds_t).intersection(rdk_bonds_t)
    sb_only = set(sms_bonds_t) - inter_bonds
    rb_only = set(rdk_bonds_t) - inter_bonds

    left_atom_colors = {int(i): cfg.common_color for i in inter_atoms}
    left_atom_colors.update({int(i): cfg.diff_color for i in sa_only})
    right_atom_colors = {int(i): cfg.common_color for i in inter_atoms}
    right_atom_colors.update({int(i): cfg.diff_color for i in rk_only})

    left_bond_colors = {int(i): cfg.common_color for i in inter_bonds}
    left_bond_colors.update({int(i): cfg.diff_color for i in sb_only})
    right_bond_colors = {int(i): cfg.common_color for i in inter_bonds}
    right_bond_colors.update({int(i): cfg.diff_color for i in rb_only})

    legends = [f"{cfg.title_smsd}: atoms={len(sms_atoms_t)}; bonds={len(sms_bonds_t)}",
               f"{cfg.title_rdkit}: atoms={len(rdk_atoms_t)}; bonds={len(rdk_bonds_t)}"]
    img = _draw_pair_grid(
        [m_target, m_target],
        legends,
        [sms_atoms_t, rdk_atoms_t],
        [sms_bonds_t, rdk_bonds_t],
        [left_atom_colors, right_atom_colors],
        [left_bond_colors, right_bond_colors],
        cfg.subimg_size
    )
    if cfg.save_png: img.save(cfg.save_png)
    return img

# -------- MCS (mol1 vs mol2 panels) --------
def viz_compare_mcs_png(mol1: str|Chem.Mol,
                        mol2: str|Chem.Mol,
                        *,
                        smsd_chem: ChemOptions = ChemOptions(),
                        mcs_opt: MCSOptions = MCSOptions(),
                        cfg: VizConfig = VizConfig()):
    """Render two panels:
    - Left: Mol1 with SMSD MCS highlights (vs Mol2)
    - Right: Mol2 with RDKit FMCS highlights (vs Mol1)
    Common atoms/bonds across the two references are coloured with cfg.common_color,
    unique ones in cfg.diff_color.
    """
    smsd = SMSD(mol1, mol2, chem=smsd_chem)
    mr = smsd.mcs_max(mcs_opt)
    m1 = smsd.q; m2 = smsd.t

    # SMSD mapping atoms/bonds on m1 and m2
    if mr:
        a1 = sorted(mr.mapping.keys())
        a2 = sorted(mr.mapping.values())
        b1 = _bond_indices_from_atomset(m1, a1)
        b2 = _bond_indices_from_mapping(m1, m2, mr.mapping)
    else:
        a1 = a2 = []; b1 = b2 = []

    # RDKit FMCS pattern and matches
    from rdkit.Chem import rdFMCS
    p = rdFMCS.MCSParameters(); p.MaximizeBonds = True; p.Timeout = 5
    res = rdFMCS.FindMCS([m1, m2], p)
    r1 = r2 = []; rb1 = rb2 = []
    if not res.canceled and res.smartsString:
        patt = Chem.MolFromSmarts(res.smartsString)
        if patt:
            m1m = list(m1.GetSubstructMatches(patt, uniquify=True))
            m2m = list(m2.GetSubstructMatches(patt, uniquify=True))
            r1 = list(m1m[0]) if m1m else []
            r2 = list(m2m[0]) if m2m else []
            rb1 = _bond_indices_from_atomset(m1, r1)
            rb2 = _bond_indices_from_atomset(m2, r2)

    # Intersections and diffs per side
    i1 = set(a1).intersection(r1); i2 = set(a2).intersection(r2)
    d1 = set(a1) - i1; d2 = set(a2) - i2
    ib1 = set(b1).intersection(rb1); ib2 = set(b2).intersection(rb2)
    db1 = set(b1) - ib1; db2 = set(b2) - ib2

    left_atom_colors = {int(i): cfg.common_color for i in i1}
    left_atom_colors.update({int(i): cfg.diff_color for i in d1})
    right_atom_colors = {int(i): cfg.common_color for i in i2}
    right_atom_colors.update({int(i): cfg.diff_color for i in d2})

    left_bond_colors = {int(i): cfg.common_color for i in ib1}
    left_bond_colors.update({int(i): cfg.diff_color for i in db1})
    right_bond_colors = {int(i): cfg.common_color for i in ib2}
    right_bond_colors.update({int(i): cfg.diff_color for i in db2})

    # legend sizes
    sms_atoms, sms_bonds = len(a1), len(b2)
    fm_atoms, fm_bonds = (len(r2), len(rb2)) if r2 else (0, 0)

    legends = [f"{cfg.title_smsd}: atoms={sms_atoms}; bonds={sms_bonds}",
               f"{cfg.title_rdkit}: atoms={fm_atoms}; bonds={fm_bonds}"]
    img = _draw_pair_grid(
        [m1, m2],
        legends,
        [a1, r2],
        [b1, rb2],
        [left_atom_colors, right_atom_colors],
        [left_bond_colors, right_bond_colors],
        cfg.subimg_size
    )
    if cfg.save_png: img.save(cfg.save_png)
    return img
