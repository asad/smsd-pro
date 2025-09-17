# SPDX-License-Identifier: Apache-2.0
# © 2025 BioInception PVT LTD.
"""
Core engines: VF2++ substructure, modular-product + BBMC MCS, optional McGregor extension,
and a friendly SMSD facade that exposes RDKit‑like options.

This file is a fixed drop-in that (a) avoids CanonicalRankAtoms preconditions and
(b) never returns None from mcs_max (tests expect a non-None object).
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from typing import Dict, List, Tuple, Optional, Iterable, Set, Callable, Literal, Any
import time
from rdkit import Chem

from .chem import (
    ChemOptions, Standardiser, RingCache, _atoms_compatible, _bonds_compatible, _wl_colours, SmartsPattern
)
from .smarts_vm import RecPredicateVM


# ----------------------------
# Public option dataclasses
# ----------------------------
@dataclass(frozen=True)
class SubstructureOptions:
    engine: Literal["AUTO", "VF2PP"] = "AUTO"
    connected_only: bool = True
    induced: bool = False
    wl_rounds: int = 0
    time_limit_s: Optional[float] = 5.0
    max_matches: Optional[int] = None
    seed: Optional[int] = 0
    uniquify_mode: Literal["mapping","target_set"] = "target_set"


@dataclass(frozen=True)
class MCSOptions:
    mcs_type: Literal["MCIS","MCCS"] = "MCCS"
    connected_only: bool = True
    time_limit_s: Optional[float] = 10.0
    use_mcgregor_extend: bool = True
    extend_time_limit_s: Optional[float] = 3.0


@dataclass(frozen=True)
class MatchResult:
    mapping: Dict[int,int]
    size: int
    algorithm: str
    tanimoto_atoms: Optional[float] = None
    tanimoto_bonds: Optional[float] = None


# --------------------
# VF2++ state machine
# --------------------
class _VF2PPState:
    def __init__(self, q: Chem.Mol, t: Chem.Mol, C: ChemOptions):
        import numpy as np
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
        self.qdeg = [len(self.QN[i]) for i in range(self.Nq)]
        self.tdeg = [len(self.TN[i]) for i in range(self.Nt)]

        self.compat = [[False]*self.Nt for _ in range(self.Nq)]
        for i in range(self.Nq):
            qa = q.GetAtomWithIdx(i)
            for j in range(self.Nt):
                ta = t.GetAtomWithIdx(j)
                if _atoms_compatible(qa, ta, self.rc_q, self.rc_t, self.C):
                    if int(self.qdeg[i]) <= int(self.tdeg[j]) + int(self.C.degree_slack):
                        self.compat[i][j] = True

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
        for j in range(self.Nt):
            if not self.compat[i][j]: continue
            if self.t2q[j] != -1: continue
            ok = True
            for iqn in self.QN[i]:
                tj = self.q2t[iqn]
                if tj != -1:
                    qb = self.q.GetBondBetweenAtoms(i, iqn)
                    tb = self.t.GetBondBetweenAtoms(int(j), int(tj))
                    if not _bonds_compatible(qb, tb, self.rc_q, self.rc_t, self.C):
                        ok = False; break
            if not ok: continue
            # degree/terminality look‑ahead
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


def _mapping_has_complete_rings(q: Chem.Mol, t: Chem.Mol, m: Dict[int,int]) -> bool:
    rc_q, rc_t = RingCache(q), RingCache(t)
    mapped_pairs = {(qi, qj): (m[qi], m[qj]) for b in q.GetBonds()
                    for qi,qj in [(b.GetBeginAtomIdx(), b.GetEndAtomIdx())]
                    if qi in m and qj in m}
    for (qi,qj),(ti,tj) in mapped_pairs.items():
        b = t.GetBondBetweenAtoms(ti, tj)
        if b is None: return False
    # each ring in q that's partially mapped must correspond to a full ring in t
    for ring in rc_q.bond_rings:
        mapped_bonds = []
        for bi in ring:
            bq = q.GetBondWithIdx(bi)
            qi, qj = bq.GetBeginAtomIdx(), bq.GetEndAtomIdx()
            if qi in m and qj in m:
                tb = t.GetBondBetweenAtoms(m[qi], m[qj])
                if tb is None: return False
                mapped_bonds.append(tb.GetIdx())
        if mapped_bonds and not any(set(mapped_bonds).issubset(R) for R in rc_t.bond_rings):
            return False
    return True


def _vf2pp_search(q: Chem.Mol, t: Chem.Mol, C: ChemOptions, opt: SubstructureOptions) -> List[Dict[int,int]]:
    st = _VF2PPState(q, t, C)
    results: List[Dict[int,int]] = []
    seen: Set[Tuple[int,...]] = set()
    start = time.time()

    def rec() -> bool:
        if opt.time_limit_s is not None and (time.time() - start) > float(opt.time_limit_s):
            return False
        if st.finished():
            m = {i:int(st.q2t[i]) for i in range(st.Nq)}
            if C.complete_rings_only and not _mapping_has_complete_rings(q, t, m):
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
                        qb = q.GetBondBetweenAtoms(i, iq)
                        tb = t.GetBondBetweenAtoms(j, int(tj))
                        if (qb is None) != (tb is None):
                            ok_induced = False; break
                if not ok_induced: continue
            st.add(i, j)
            if rec():
                any_hit = True
                if opt.max_matches and len(results) >= int(opt.max_matches):
                    st.backtrack(); return True
            st.backtrack()
        return any_hit

    rec()
    return results


# ------------------------
# MCS via modular product
# ------------------------
@dataclass(frozen=True)
class _MPNode:
    qi: int
    tj: int


def _build_modular_product(m1: Chem.Mol, m2: Chem.Mol, C: ChemOptions) -> Tuple[List[_MPNode], List[int]]:
    rc1, rc2 = RingCache(m1), RingCache(m2)
    N1, N2 = m1.GetNumAtoms(), m2.GetNumAtoms()
    nodes: List[_MPNode] = []
    for i in range(N1):
        ai = m1.GetAtomWithIdx(int(i))
        for j in range(N2):
            aj = m2.GetAtomWithIdx(int(j))
            if _atoms_compatible(ai, aj, rc1, rc2, C):
                nodes.append(_MPNode(int(i), int(j)))
    n = len(nodes)
    adj: List[int] = [0]*n
    for u in range(n):
        i, j = nodes[u].qi, nodes[u].tj
        for v in range(u+1, n):
            k, l = nodes[v].qi, nodes[v].tj
            if i == k or j == l: continue
            qb = m1.GetBondBetweenAtoms(int(i), int(k))
            tb = m2.GetBondBetweenAtoms(int(j), int(l))
            if qb is not None and tb is not None:
                if _bonds_compatible(qb, tb, rc1, rc2, C):
                    adj[u] |= (1<<v); adj[v] |= (1<<u)
            else:
                if (qb is None) and (tb is None):
                    adj[u] |= (1<<v); adj[v] |= (1<<u)
    return nodes, adj


def _bbmc_max_cliques(adj: List[int], time_limit_s: Optional[float] = None) -> Tuple[int, List[int]]:
    start = time.time()
    n = len(adj)
    V_all = (1<<n) - 1
    best_size = 0
    best_masks: List[int] = []

    def colour_bound(P_mask: int) -> int:
        colours = 0
        P = P_mask
        while P:
            colours += 1
            v = (P & -P); idx_v = v.bit_length()-1
            X = ~adj[idx_v] & V_all
            P &= ~v
            while P:
                w = (P & -P); idx = w.bit_length()-1
                if (X >> idx) & 1:
                    X &= ~adj[idx]
                    P &= ~w
                else:
                    P &= ~w
        return colours

    def rec(R_mask: int, P_mask: int, X_mask: int):
        nonlocal best_size, best_masks
        if time_limit_s is not None and (time.time() - start) > float(time_limit_s): return
        R_size = R_mask.bit_count()
        if R_size + colour_bound(P_mask) < best_size: return
        if not P_mask and not X_mask:
            if R_size > best_size:
                best_size = R_size; best_masks = [R_mask]
            elif R_size == best_size:
                best_masks.append(R_mask)
            return
        U = P_mask | X_mask
        if U:
            best_u, best_deg = 0, -1
            tmp = U
            while tmp:
                u_bit = (tmp & -tmp); idx = u_bit.bit_length()-1
                deg = (P_mask & adj[idx]).bit_count()
                if deg > best_deg:
                    best_u, best_deg = idx, deg
                tmp &= ~u_bit
            candidates = P_mask & ~adj[best_u]
        else:
            candidates = P_mask
        tmp = candidates
        while tmp:
            v_bit = (tmp & -tmp); v = v_bit.bit_length()-1
            rec(R_mask | (1<<v), P_mask & adj[v], X_mask & adj[v])
            P_mask &= ~v_bit; X_mask |= v_bit; tmp &= ~v_bit

    rec(0, V_all, 0)
    return best_size, best_masks


def _mask_to_mapping(nodes: List[_MPNode], mask: int) -> Dict[int,int]:
    m: Dict[int,int] = {}
    while mask:
        v_bit = (mask & -mask)
        v = v_bit.bit_length()-1
        n = nodes[int(v)]
        m[int(n.qi)] = int(n.tj)
        mask &= ~v_bit
    return dict(sorted(m.items()))


def _largest_connected_in_query(q: Chem.Mol, mapping: Dict[int,int]) -> Dict[int,int]:
    if not mapping: return mapping
    mapped = sorted(mapping.keys())
    pos = {qi:i for i,qi in enumerate(mapped)}
    parent = list(range(len(mapped)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    def union(a,b):
        ra, rb = find(a), find(b)
        if ra != rb: parent[rb] = ra
    for i, qi in enumerate(mapped):
        for n in q.GetAtomWithIdx(int(qi)).GetNeighbors():
            j = int(n.GetIdx())
            if j in pos: union(i, pos[j])
    from collections import Counter
    cnt = Counter(find(i) for i in range(len(mapped)))
    root = max(cnt.items(), key=lambda kv: kv[1])[0]
    keep = {mapped[i] for i in range(len(mapped)) if find(i)==root}
    return {qi:mapping[qi] for qi in sorted(keep)}


def _mcgregor_extend(m1: Chem.Mol, m2: Chem.Mol, seed: Dict[int,int],
                     C: ChemOptions, time_limit_s: Optional[float]=None) -> Dict[int,int]:
    start = time.time()
    best = dict(seed)
    rc1, rc2 = RingCache(m1), RingCache(m2)

    q_un = set(range(m1.GetNumAtoms())) - set(seed.keys())
    t_un = set(range(m2.GetNumAtoms())) - set(seed.values())
    QN = [[n.GetIdx() for n in m1.GetAtomWithIdx(i).GetNeighbors()] for i in range(m1.GetNumAtoms())]

    wl1, wl2 = _wl_colours(m1), _wl_colours(m2)

    def wl_bound() -> bool:
        from collections import Counter
        cq = Counter(wl1[i] for i in q_un)
        ct = Counter(wl2[j] for j in t_un)
        for c, k in cq.items():
            if ct.get(c, 0) < k: return False
        return True

    def feasible(qi: int, tj: int, cur: Dict[int,int]) -> bool:
        if not _atoms_compatible(m1.GetAtomWithIdx(int(qi)), m2.GetAtomWithIdx(int(tj)), rc1, rc2, C):
            return False
        for qk, tl in cur.items():
            qb = m1.GetBondBetweenAtoms(int(qi), int(qk))
            tb = m2.GetBondBetweenAtoms(int(tj), int(tl))
            if (qb is None) != (tb is None): return False
            if qb is not None and not _bonds_compatible(qb, tb, rc1, rc2, C): return False
        return True

    def frontier(cur: Dict[int,int]) -> Set[int]:
        F: Set[int] = set()
        for qk in cur.keys():
            for nn in QN[int(qk)]:
                if int(nn) in q_un: F.add(int(nn))
        return F if F else set(int(x) for x in q_un)

    def upper_bound(cur_size: int) -> int:
        return int(cur_size) + min(len(q_un), len(t_un))

    def rec(cur: Dict[int,int]):
        nonlocal best
        if time_limit_s is not None and (time.time() - start) > float(time_limit_s): return
        if not wl_bound(): return
        if upper_bound(len(cur)) <= len(best): return
        if len(cur) > len(best): best = dict(cur)
        F = frontier(cur)
        if not F: return
        candidates_by_q: List[Tuple[int,List[int]]] = []
        for qi in F:
            ts = [int(tj) for tj in t_un if feasible(int(qi), int(tj), cur)]
            if not ts: continue
            candidates_by_q.append((int(qi), ts))
        if not candidates_by_q: return
        candidates_by_q.sort(key=lambda x: (len(x[1]), -m1.GetAtomWithIdx(int(x[0])).GetDegree(), -int(m1.GetAtomWithIdx(int(x[0])).IsInRing())))
        qi, ts = candidates_by_q[0]
        for tj in ts:
            cur[int(qi)] = int(tj)
            q_un.remove(int(qi)); t_un.remove(int(tj))
            rec(cur)
            del cur[int(qi)]
            q_un.add(int(qi)); t_un.add(int(tj))

    rec(dict(seed))
    return best


# --------------
# Public facade
# --------------
@dataclass
class _STD:
    largest_fragment: bool = True
    normalise: bool = True
    reionise: bool = True
    uncharge: bool = True
    canonical_tautomer: bool = True
    kekulise: bool = False


class SMSD:
    """Friendly facade around the engines + chemistry options.

    Defaults: connected mapping, conservative chemistry, standardised inputs.
    """
    def __init__(self, query: Any, target: Any, *, chem: ChemOptions = ChemOptions(),
                 standardise: bool = True, std: Optional[Dict[str,Any]] = None):
        self.chem = chem
        self.std_opts = std or _STD().__dict__
        stdzr = Standardiser()

        def _to_mol(x) -> Tuple[Chem.Mol, bool]:
            if isinstance(x, Chem.Mol):
                return (stdzr.run(x, **self.std_opts) if standardise else Chem.Mol(x)), False
            if isinstance(x, str):
                # Prefer SMILES; fallback to SMARTS (with $name stripped for RDKit parsing)
                m_smiles = Chem.MolFromSmiles(x)
                if m_smiles is not None:
                    return (stdzr.run(m_smiles, **self.std_opts) if standardise else m_smiles), False
                sanitized = SmartsPattern.strip_named_predicates(x)
                m_smarts = Chem.MolFromSmarts(sanitized)
                if m_smarts is None:
                    raise ValueError(f"Cannot parse input: {x}")
                return m_smarts, True
            raise TypeError("Query/Target must be RDKit Mol or SMILES/SMARTS string.")

        self.q, self.query_is_smarts = _to_mol(query)
        self.t, _ = _to_mol(target)

        if not self.query_is_smarts:
            Chem.SanitizeMol(self.q)
        Chem.SanitizeMol(self.t)

        # Bind recursive SMARTS if present on the original query string
        if isinstance(query, str):
            try:
                qpat = SmartsPattern(query, is_smarts=True)
            except Exception:
                qpat = None
            self.q_pattern = qpat
        else:
            self.q_pattern = None

        if self.q_pattern and self.q_pattern.rec_by_atom:
            vm = RecPredicateVM()
            # Minimal built-ins; extend in client code as needed
            vm.register("isAmideN", smarts="[$([NX3;H2,H1;!$(NC=O)]),$([NX3;$(NC=O)])]")
            vm.register("isCarboxylC", smarts="[CX3](=O)[OX2H1,OX1-]")
            rec_specs = self.q_pattern.rec_by_atom
            def _rec_ok(qmol, tmol, qi, tj, _vm=vm, _specs=rec_specs):
                specs = _specs.get(int(qi))
                if not specs: return True
                return _vm.eval_all_for_atom(tmol, int(tj), specs)
            self.chem = replace(self.chem, atom_rec_ok=_rec_ok)

    # ---- Convenience
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
        if self.query_is_smarts: return True
        from collections import Counter
        hq = Counter(a.GetAtomicNum() for a in self.q.GetAtoms())
        ht = Counter(a.GetAtomicNum() for a in self.t.GetAtoms())
        for z, k in hq.items():
            if ht.get(z, 0) < k: return False
        return True

    # ---- Substructure
    def substructure_exists(self, opt: SubstructureOptions = SubstructureOptions()) -> bool:
        if not self.subgraph_possible(): return False
        return bool(_vf2pp_search(self.q, self.t, self.chem, opt))

    def substructure_all(self, opt: SubstructureOptions = SubstructureOptions()) -> List[MatchResult]:
        if not self.subgraph_possible(): return []
        maps = _vf2pp_search(self.q, self.t, self.chem, opt)
        maps = self._dedupe(maps, opt.uniquify_mode)
        out: List[MatchResult] = []
        bonds_common = lambda m: sum(1 for b in self.q.GetBonds()
                                     if int(b.GetBeginAtomIdx()) in m and int(b.GetEndAtomIdx()) in m and
                                     self.t.GetBondBetweenAtoms(int(m[int(b.GetBeginAtomIdx())]), int(m[int(b.GetEndAtomIdx())])) is not None)
        for m in maps:
            tana, tanb = self._tan(len(m), bonds_common(m))
            out.append(MatchResult(m, len(m), "SUBSTRUCTURE_VF2PP", tana, tanb))
        return out

    # ---- MCS
    def mcs_max(self, opt: MCSOptions = MCSOptions()) -> MatchResult:
        nodes, adj = _build_modular_product(self.q, self.t, self.chem)
        if not nodes:
            return MatchResult({}, 0, f"{opt.mcs_type}_CLIQUE_NONE", 0.0, 0.0)
        best_size, cliques = _bbmc_max_cliques(adj, time_limit_s=opt.time_limit_s)
        if best_size == 0 or not cliques:
            return MatchResult({}, 0, f"{opt.mcs_type}_CLIQUE_NONE", 0.0, 0.0)
        mapping = _mask_to_mapping(nodes, cliques[0])
        if opt.mcs_type == "MCCS" or opt.connected_only:
            mapping = _largest_connected_in_query(self.q, mapping)

        # McGregor extension: skip for SMARTS queries (pattern valence/H may be undefined)
        if opt.use_mcgregor_extend and mapping and not self.query_is_smarts:
            mapping = _mcgregor_extend(self.q, self.t, mapping, self.chem, time_limit_s=opt.extend_time_limit_s)
            if opt.mcs_type == "MCCS" or opt.connected_only:
                mapping = _largest_connected_in_query(self.q, mapping)

        bonds_common = sum(1 for b in self.q.GetBonds()
                           if int(b.GetBeginAtomIdx()) in mapping and int(b.GetEndAtomIdx()) in mapping and
                           self.t.GetBondBetweenAtoms(int(mapping[int(b.GetBeginAtomIdx())]), int(mapping[int(b.GetEndAtomIdx())])) is not None)
        tana, tanb = self._tan(len(mapping), bonds_common)
        algo = f"{'MCCS' if (opt.mcs_type=='MCCS' or opt.connected_only) else 'MCIS'}_CLIQUE_BBMC"
        if opt.use_mcgregor_extend and not self.query_is_smarts:
            algo += "+MCG"
        return MatchResult(mapping, len(mapping), algo, tana, tanb)
