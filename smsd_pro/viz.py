# SPDX-License-Identifier: Apache-2.0
# © 2025 BioInception PVT LTD.
\
# smsd_pro/viz.py – side-by-side PNG comparison of SMSD vs RDKit
from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple, Optional, Dict
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

def _bond_indices_from_mapping(q: Chem.Mol, t: Chem.Mol, m: Dict[int,int]) -> List[int]:
    idxs: List[int] = []
    for b in q.GetBonds():
        i, j = int(b.GetBeginAtomIdx()), int(b.GetEndAtomIdx())
        if i in m and j in m:
            tb = t.GetBondBetweenAtoms(int(m[i]), int(m[j]))
            if tb is not None:
                idxs.append(int(tb.GetIdx()))
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
                    highlights: List[Tuple[List[int], List[int]]],
                    size: Tuple[int, int]):
    ha = [a for a,_ in highlights]
    hb = [b for _,b in highlights]
    img = Draw.MolsToGridImage(
        mols,
        molsPerRow=2,
        subImgSize=size,
        legends=legends,
        highlightAtomLists=ha,
        highlightBondLists=hb,
        useSVG=False
    )
    return _to_pil(img)

# -------- Substructure --------
def viz_compare_substructure_png(query: str|Chem.Mol,
                                 target: str|Chem.Mol,
                                 *,
                                 smsd_chem: ChemOptions = ChemOptions(),
                                 smsd_opt: SubstructureOptions = SubstructureOptions(),
                                 mapping_index: int = 0,
                                 cfg: VizConfig = VizConfig()):
    smsd = SMSD(query, target, chem=smsd_chem)
    maps = smsd.substructure_all(smsd_opt)
    q = smsd.q; t = smsd.t

    if not maps:
        img = _draw_pair_grid([q,t], ["Query (no match)", "Target"], [([],[]),([],[])], cfg.subimg_size)
        if cfg.save_png: img.save(cfg.save_png)
        return img

    mr: MatchResult = maps[min(mapping_index, len(maps)-1)]
    q_atoms = sorted(mr.mapping.keys())
    t_atoms = sorted(mr.mapping.values())
    t_bonds = _bond_indices_from_mapping(q, t, mr.mapping)

    # RDKit reference mapping (strip $name so RDKit can parse)
    if isinstance(query, str):
        q_try = Chem.MolFromSmiles(query)
        if q_try is None:
            q_try = Chem.MolFromSmarts(SmartsPattern.strip_named_predicates(query))
    else:
        q_try = query
    t_try = Chem.MolFromSmiles(target) if isinstance(target, str) else target
    rdkit_matches = list(t_try.GetSubstructMatches(q_try, uniquify=True)) if q_try and t_try else []
    atoms_rdkit = list(rdkit_matches[0]) if rdkit_matches else []

    legends = [f"{cfg.title_smsd}: atoms={len(q_atoms)}",
               f"{cfg.title_rdkit}: atoms={len(atoms_rdkit)}"]
    img = _draw_pair_grid(
        [q,t],
        legends,
        highlights=[(q_atoms, []), (t_atoms, t_bonds)],
        size=cfg.subimg_size
    )
    if cfg.save_png: img.save(cfg.save_png)
    return img

# -------- MCS --------
def viz_compare_mcs_png(mol1: str|Chem.Mol,
                        mol2: str|Chem.Mol,
                        *,
                        smsd_chem: ChemOptions = ChemOptions(),
                        mcs_opt: MCSOptions = MCSOptions(),
                        cfg: VizConfig = VizConfig()):
    smsd = SMSD(mol1, mol2, chem=smsd_chem)
    mr = smsd.mcs_max(mcs_opt)
    m1 = smsd.q; m2 = smsd.t
    if not mr:
        img = _draw_pair_grid([m1,m2], ["Mol1 (no MCS)", "Mol2"], [([],[]),([],[])], cfg.subimg_size)
        if cfg.save_png: img.save(cfg.save_png)
        return img

    a1 = sorted(mr.mapping.keys())
    a2 = sorted(mr.mapping.values())
    b2 = _bond_indices_from_mapping(m1, m2, mr.mapping)

    # RDKit FMCS reference (just to show sizes)
    from rdkit.Chem import rdFMCS
    p = rdFMCS.MCSParameters()
    p.MaximizeBonds = True; p.Timeout = 5
    res = rdFMCS.FindMCS([m1, m2], p)
    a_ref = b_ref = 0
    if not res.canceled and res.smartsString:
        patt = Chem.MolFromSmarts(res.smartsString)
        if patt:
            a_ref, b_ref = patt.GetNumAtoms(), patt.GetNumBonds()

    legends = [f"{cfg.title_smsd}: atoms={len(a1)}; bonds={len(b2)}",
               f"{cfg.title_rdkit}: atoms={a_ref}; bonds={b_ref}"]
    img = _draw_pair_grid([m1, m2], legends,
                          highlights=[(a1, []), (a2, b2)],
                          size=cfg.subimg_size)
    if cfg.save_png: img.save(cfg.save_png)
    return img
