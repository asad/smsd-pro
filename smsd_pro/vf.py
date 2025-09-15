from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
from rdkit import Chem
from .chem import ChemOptions, RingCache, atoms_compatible, bonds_compatible

@dataclass(frozen=True)
class SubstructureOptions:
    engine: str = "AUTO"              # kept for future selection
    connected_only: bool = True
    induced: bool = False
    wl_rounds: int = 0
    time_limit_s: Optional[float] = 5.0
    max_matches: Optional[int] = None
    seed: Optional[int] = 0
    uniquify_mode: str = "target_set" # "mapping" | "target_set"

class _VF2PPState:
    """State for VF2++-style mapping (frontier-based ordering, degree/ring look-ahead)."""
    def __init__(self, q: Chem.Mol, t: Chem.Mol, C: ChemOptions):
        self.q, self.t, self.C = q, t, C
        self.Nq, self.Nt = q.GetNumAtoms(), t.GetNumAtoms()
        self.rc_q, self.rc_t = RingCache(q), RingCache(t)

        self.q2t = np.full(self.Nq, -1, dtype=np.int32)
        self.t2q = np.full(self.Nt, -1, dtype=np.int32)
        self.depth = 0
        self.stack: List[Tuple[int,int]] = []

        self.qterm = np.zeros(self.Nq, dtype=np.int32)
        self.tterm = np.zeros(self.Nt, dtype=np.int32)

        self.QN = [[n.GetIdx() for n in q.GetAtomWithIdx(i).GetNeighbors()] for i in range(self.Nq)]
        self.TN = [[n.GetIdx() for n in t.GetAtomWithIdx(i).GetNeighbors()] for i in range(self.Nt)]
        self.qdeg = np.array([len(self.QN[i]) for i in range(self.Nq)], dtype=np.int32)
        self.tdeg = np.array([len(self.TN[i]) for i in range(self.Nt)], dtype=np.int32)

        self.compat = np.zeros((self.Nq, self.Nt), dtype=np.bool_)
        for i in range(self.Nq):
            qa = q.GetAtomWithIdx(i)
            for j in range(self.Nt):
                ta = t.GetAtomWithIdx(j)
                if atoms_compatible(qa, ta, self.rc_q, self.rc_t, self.C):
                    if self.qdeg[i] <= self.tdeg[j] + self.C.degree_slack:
                        self.compat[i, j] = True

    def finished(self) -> bool:
        return int(self.depth) == int(self.Nq)

    def add(self, i: int, j: int):
        self.depth += 1
        self.stack.append((i, j))
        self.q2t[i] = j
        self.t2q[j] = i
        for u in self.QN[i]:
            if self.q2t[u] == -1 and self.qterm[u] == 0:
                self.qterm[u] = self.depth
        for v in self.TN[j]:
            if self.t2q[v] == -1 and self.tterm[v] == 0:
                self.tterm[v] = self.depth

    def backtrack(self):
        i, j = self.stack.pop()
        for u in self.QN[i]:
            if self.qterm[u] == self.depth: self.qterm[u] = 0
        for v in self.TN[j]:
            if self.tterm[v] == self.depth: self.tterm[v] = 0
        self.q2t[i] = -1; self.t2q[j] = -1
        self.depth -= 1

    def _cand_targets(self, i: int) -> Iterable[int]:
        for j in np.flatnonzero(self.compat[i]):
            if self.t2q[j] != -1: continue
            ok = True
            for iqn in self.QN[i]:
                tj = self.q2t[iqn]
                if tj != -1:
                    qb = self.q.GetBondBetweenAtoms(int(i), int(iqn))
                    tb = self.t.GetBondBetweenAtoms(int(j), int(tj))
                    if not bonds_compatible(qb, tb, self.rc_q, self.rc_t, self.C):
                        ok = False; break
            if not ok: continue
            q_term = sum(1 for u in self.QN[i] if self.q2t[u]==-1 and self.qterm[u]>0)
            q_new  = sum(1 for u in self.QN[i] if self.q2t[u]==-1 and self.qterm[u]==0)
            t_term = sum(1 for v in self.TN[j] if self.t2q[v]==-1 and self.tterm[v]>0)
            t_new  = sum(1 for v in self.TN[j] if self.t2q[v]==-1 and self.tterm[v]==0)
            if q_term <= t_term and q_new <= t_new:
                yield int(j)

    def pick_next_q(self, connected_only: bool) -> Optional[int]:
        cand = [i for i in range(self.Nq) if self.q2t[i]==-1 and (self.qterm[i]>0 or not connected_only or self.depth==0)]
        if not cand: return None
        best, best_cnt, best_key = None, 1<<30, None
        for i in cand:
            cnt = 0
            for _ in self._cand_targets(i):
                cnt += 1
                if cnt >= best_cnt: break
            qa = self.q.GetAtomWithIdx(i)
            key = (-qa.GetDegree(), -int(qa.IsInRing()), -int(qa.GetIsAromatic()), i)
            if cnt < best_cnt or (cnt == best_cnt and (best_key is None or key < best_key)):
                best, best_cnt, best_key = i, cnt, key
        return best

def mapping_has_complete_rings(q: Chem.Mol, t: Chem.Mol, m: Dict[int,int]) -> bool:
    rc_q, rc_t = RingCache(q), RingCache(t)
    mapped_pairs = {(qi, qj): (m[qi], m[qj]) for b in q.GetBonds()
                    for qi,qj in [(b.GetBeginAtomIdx(), b.GetEndAtomIdx())]
                    if qi in m and qj in m}
    t_bond_set = set()
    for (qi,qj),(ti,tj) in mapped_pairs.items():
        b = t.GetBondBetweenAtoms(int(ti), int(tj))
        if b is None: return False
        t_bond_set.add(b.GetIdx())
    for ring in rc_q.bond_rings:
        mapped_bonds = []
        for bi in ring:
            bq = q.GetBondWithIdx(bi)
            qi, qj = bq.GetBeginAtomIdx(), bq.GetEndAtomIdx()
            if qi in m and qj in m:
                tb = t.GetBondBetweenAtoms(int(m[qi]), int(m[qj]))
                if tb is None: return False
                mapped_bonds.append(tb.GetIdx())
        if mapped_bonds and not any(set(mapped_bonds).issubset(R) for R in rc_t.bond_rings):
            return False
    return True

def vf2pp_search(q: Chem.Mol, t: Chem.Mol, C: ChemOptions, opt: SubstructureOptions) -> List[Dict[int,int]]:
    import time
    st = _VF2PPState(q, t, C)
    results: List[Dict[int,int]] = []
    seen: Set[Tuple[int,...]] = set()
    start = time.time()

    def rec() -> bool:
        if opt.time_limit_s is not None and (time.time() - start) > opt.time_limit_s:
            return False
        if st.finished():
            m = {i:int(st.q2t[i]) for i in range(st.Nq)}
            if C.complete_rings_only and not mapping_has_complete_rings(q, t, m):
                return False
            key = tuple(sorted(m.values())) if opt.uniquify_mode=="target_set" else tuple(sorted(m.items()))
            if key not in seen:
                seen.add(key); results.append(m)
            return True
        i = st.pick_next_q(connected_only=opt.connected_only)
        if i is None:
            return False
        any_hit = False
        for j in st._cand_targets(i):
            if opt.induced:
                ok_induced = True
                for iq in range(st.Nq):
                    tj = st.q2t[iq]
                    if tj != -1 and iq != i:
                        qb = q.GetBondBetweenAtoms(int(i), int(iq))
                        tb = t.GetBondBetweenAtoms(int(j), int(tj))
                        if (qb is None) != (tb is None):
                            ok_induced = False; break
                if not ok_induced: continue
            st.add(int(i), int(j))
            if rec():
                any_hit = True
                if opt.max_matches and len(results) >= opt.max_matches:
                    st.backtrack(); return True
            st.backtrack()
        return any_hit

    rec()
    return results
