from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

from rdkit import Chem
from rdkit.Chem import Draw
from .chem import ChemOptions, Standardiser
from .vf import SubstructureOptions, vf2pp_search
from .mcs import (
    MCSOptions, build_modular_product, bbmc_max_cliques, mask_to_mapping,
    largest_connected_in_query, mcgregor_extend
)

@dataclass(frozen=True)
class MatchResult:
    mapping: Dict[int,int]          # query_idx -> target_idx
    size: int
    algorithm: str
    tanimoto_atoms: Optional[float] = None
    tanimoto_bonds: Optional[float] = None

class SMSD:
    """
    Main façade combining substructure and MCS engines with chemistry options.
    """
    def __init__(self, query: Any, target: Any,
                 *, chem: ChemOptions = ChemOptions(),
                 standardise: bool = True,
                 std: Optional[Dict[str,Any]] = None):
        self.chem = chem
        self.std_opts = std or dict(largest_fragment=True, normalise=True, reionise=True,
                                    uncharge=True, canonical_tautomer=True, kekulise=False)
        stdzr = Standardiser()

        def _to_mol(x) -> Tuple[Chem.Mol, bool]:
            if isinstance(x, Chem.Mol):
                return (stdzr.run(x, **self.std_opts) if standardise else Chem.Mol(x)), False
            if isinstance(x, str):
                is_smarts = any(ch in x for ch in "*[]$;@/\\!?")
                m = Chem.MolFromSmarts(x) if is_smarts else Chem.MolFromSmiles(x)
                if m is None:
                    raise ValueError(f"Cannot parse input: {x}")
                if is_smarts:
                    return m, True
                else:
                    m = stdzr.run(m, **self.std_opts) if standardise else m
                    return m, False
            raise TypeError("Query/Target must be RDKit Mol or SMILES/SMARTS string.")

        self.q, self.query_is_smarts = _to_mol(query)
        self.t, _ = _to_mol(target)

        if not self.query_is_smarts:
            Chem.SanitizeMol(self.q)
        Chem.SanitizeMol(self.t)

    # ---- internals

    def _tan(self, m_atoms: int, m_bonds: int) -> Tuple[float,float]:
        a = m_atoms / max(1, (self.q.GetNumAtoms() + self.t.GetNumAtoms() - m_atoms))
        b = m_bonds / max(1, (self.q.GetNumBonds() + self.t.GetNumBonds() - m_bonds))
        return a, b

    def _dedupe(self, maps: List[Dict[int,int]], mode: str) -> List[Dict[int,int]]:
        if mode == "mapping":
            seen=set(); out=[]
            for m in maps:
                key = tuple(sorted(m.items()))
                if key in seen: continue
                seen.add(key); out.append(m)
            return out
        seen=set(); out=[]
        for m in maps:
            key = tuple(sorted(m.values()))
            if key in seen: continue
            seen.add(key); out.append(m)
        return out

    def subgraph_possible(self) -> bool:
        if self.query_is_smarts:
            return True
        from collections import Counter
        hq = Counter(a.GetAtomicNum() for a in self.q.GetAtoms())
        ht = Counter(a.GetAtomicNum() for a in self.t.GetAtoms())
        for z, k in hq.items():
            if ht.get(z, 0) < k: return False
        return True

    # ---- Substructure API

    def substructure_exists(self, opt: SubstructureOptions = SubstructureOptions()) -> bool:
        if not self.subgraph_possible(): return False
        return bool(vf2pp_search(self.q, self.t, self.chem, opt))

    def substructure_first(self, opt: SubstructureOptions = SubstructureOptions()) -> Optional[MatchResult]:
        if not self.subgraph_possible(): return None
        opt1 = SubstructureOptions(**{**opt.__dict__, "max_matches": 1})
        maps = vf2pp_search(self.q, self.t, self.chem, opt1)
        if not maps: return None
        m = maps[0]
        bonds_common = sum(1 for b in self.q.GetBonds()
                           if b.GetBeginAtomIdx() in m and b.GetEndAtomIdx() in m and
                           self.t.GetBondBetweenAtoms(int(m[b.GetBeginAtomIdx()]), int(m[b.GetEndAtomIdx()])) is not None)
        tana, tanb = self._tan(len(m), bonds_common)
        return MatchResult(m, len(m), "SUBSTRUCTURE_VF2PP", tana, tanb)

    def substructure_all(self, opt: SubstructureOptions = SubstructureOptions()) -> List[MatchResult]:
        if not self.subgraph_possible(): return []
        maps = vf2pp_search(self.q, self.t, self.chem, opt)
        maps = self._dedupe(maps, opt.uniquify_mode)
        out: List[MatchResult] = []
        for m in maps:
            bonds_common = sum(1 for b in self.q.GetBonds()
                               if b.GetBeginAtomIdx() in m and b.GetEndAtomIdx() in m and
                               self.t.GetBondBetweenAtoms(int(m[b.GetBeginAtomIdx()]), int(m[b.GetEndAtomIdx()])) is not None)
            tana, tanb = self._tan(len(m), bonds_common)
            out.append(MatchResult(m, len(m), "SUBSTRUCTURE_VF2PP", tana, tanb))
        return out

    # ---- MCS API

    def mcs_max(self, opt: MCSOptions = MCSOptions()) -> Optional[MatchResult]:
        nodes, adj = build_modular_product(self.q, self.t, self.chem)
        if not nodes:
            return None
        best_size, cliques = bbmc_max_cliques(adj, time_limit_s=opt.time_limit_s)
        if best_size == 0 or not cliques:
            return None
        mapping = mask_to_mapping(nodes, cliques[0])
        if opt.mcs_type == "MCCS":
            mapping = largest_connected_in_query(self.q, mapping)
        if opt.use_mcgregor_extend and mapping:
            mapping = mcgregor_extend(self.q, self.t, mapping, self.chem, time_limit_s=opt.extend_time_limit_s)
        bonds_common = sum(1 for b in self.q.GetBonds()
                           if b.GetBeginAtomIdx() in mapping and b.GetEndAtomIdx() in mapping and
                           self.t.GetBondBetweenAtoms(int(mapping[b.GetBeginAtomIdx()]), int(mapping[b.GetEndAtomIdx()])) is not None)
        tana, tanb = self._tan(len(mapping), bonds_common)
        algo = f"{opt.mcs_type}_CLIQUE_BBMC" + ("+MCG" if opt.use_mcgregor_extend else "")
        return MatchResult(mapping, len(mapping), algo, tana, tanb)

    # ---- Visual helpers (simple grid drawing)

    def _bond_indices_from_mapping(self, mapping: Dict[int,int]) -> List[int]:
        idxs = []
        for b in self.q.GetBonds():
            i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if i in mapping and j in mapping:
                tb = self.t.GetBondBetweenAtoms(int(mapping[i]), int(mapping[j]))
                if tb is not None:
                    idxs.append(tb.GetIdx())
        return idxs

    def grid_image_substructure(self, mr: Optional[MatchResult]) -> "PIL.Image.Image":
        q = self.q; t = self.t
        if mr is None:
            return Draw.MolsToGridImage([q,t], molsPerRow=2, subImgSize=(500,400),
                                        legends=["Query (no match)", "Target"])
        q_atoms = sorted(mr.mapping.keys())
        t_atoms = sorted(mr.mapping.values())
        t_bonds = self._bond_indices_from_mapping(mr.mapping)
        return Draw.MolsToGridImage([q,t], molsPerRow=2, subImgSize=(500,400),
                                    legends=[f"SMSD: atoms={len(q_atoms)}", f"Target: atoms={len(t_atoms)}"],
                                    highlightAtomLists=[q_atoms, t_atoms],
                                    highlightBondLists=[[], t_bonds])

# ---- CLI ----

def main():
    import argparse, sys, time
    p = argparse.ArgumentParser(description="SMSD Pro (BioInception) – substructure & MCS")
    p.add_argument("--query", help="Query SMILES/SMARTS for substructure")
    p.add_argument("--target", help="Target SMILES for substructure")
    p.add_argument("--mol1", help="Mol1 SMILES (MCS)")
    p.add_argument("--mol2", help="Mol2 SMILES (MCS)")
    p.add_argument("--mode", choices=["substructure","mcs"], required=True)
    p.add_argument("--mcs-type", default="MCIS", choices=["MCIS","MCCS"])
    p.add_argument("--timeout", type=float, default=5.0)
    p.add_argument("--extend", action="store_true", help="Use McGregor extend after clique")
    args = p.parse_args()

    if args.mode == "substructure":
        if not args.query or not args.target:
            print("Please provide --query and --target", file=sys.stderr); sys.exit(2)
        smsd = SMSD(args.query, args.target)
        t0 = time.perf_counter()
        mr = smsd.substructure_first()
        dt = time.perf_counter() - t0
        if mr:
            print(f"Substructure: True; size={mr.size}; tan_atoms={mr.tanimoto_atoms:.3f}; time={dt:.4f}s")
        else:
            print(f"Substructure: False; time={dt:.4f}s")
    else:
        if not args.mol1 or not args.mol2:
            print("Please provide --mol1 and --mol2", file=sys.stderr); sys.exit(2)
        smsd = SMSD(args.mol1, args.mol2)
        t0 = time.perf_counter()
        mr = smsd.mcs_max(MCSOptions(mcs_type=args.mcs_type, time_limit_s=args.timeout, use_mcgregor_extend=args.extend))
        dt = time.perf_counter() - t0
        if mr:
            print(f"MCS: {mr.algorithm}; size={mr.size}; tan_atoms={mr.tanimoto_atoms:.3f}; tan_bonds={mr.tanimoto_bonds:.3f}; time={dt:.4f}s")
        else:
            print(f"MCS: size=0; time={dt:.4f}s")
