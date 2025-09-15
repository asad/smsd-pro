from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from rdkit import Chem
from .chem import ChemOptions, RingCache, atoms_compatible, bonds_compatible, wl_colours

@dataclass(frozen=True)
class MCSOptions:
    mcs_type: str = "MCIS"                 # "MCIS" | "MCCS"
    time_limit_s: Optional[float] = 10.0
    use_mcgregor_extend: bool = True
    extend_time_limit_s: Optional[float] = 3.0

@dataclass(frozen=True)
class _MPNode:
    qi: int
    tj: int

def build_modular_product(m1: Chem.Mol, m2: Chem.Mol, C: ChemOptions) -> Tuple[List[_MPNode], List[int]]:
    rc1, rc2 = RingCache(m1), RingCache(m2)
    N1, N2 = m1.GetNumAtoms(), m2.GetNumAtoms()
    nodes: List[_MPNode] = []
    for i in range(N1):
        ai = m1.GetAtomWithIdx(i)
        for j in range(N2):
            aj = m2.GetAtomWithIdx(j)
            if atoms_compatible(ai, aj, rc1, rc2, C):
                nodes.append(_MPNode(i, j))
    n = len(nodes)
    adj: List[int] = [0]*n
    for u in range(n):
        i, j = nodes[u].qi, nodes[u].tj
        for v in range(u+1, n):
            k, l = nodes[v].qi, nodes[v].tj
            if i == k or j == l:
                continue
            qb = m1.GetBondBetweenAtoms(int(i), int(k))
            tb = m2.GetBondBetweenAtoms(int(j), int(l))
            if qb is not None and tb is not None:
                if bonds_compatible(qb, tb, rc1, rc2, C):
                    adj[u] |= (1<<v); adj[v] |= (1<<u)
            else:
                if (qb is None) and (tb is None):
                    adj[u] |= (1<<v); adj[v] |= (1<<u)
    return nodes, adj

def bbmc_max_cliques(adj: List[int], time_limit_s: Optional[float] = None) -> Tuple[int, List[int]]:
    import time
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
            v = (P & -P)
            idx_v = v.bit_length()-1
            X = ~adj[idx_v] & V_all
            P &= ~v
            while P:
                w = (P & -P)
                idx = w.bit_length()-1
                if (X >> idx) & 1:
                    X &= ~adj[idx]
                    P &= ~w
                else:
                    P &= ~w
        return colours

    def rec(R_mask: int, P_mask: int, X_mask: int):
        nonlocal best_size, best_masks
        if time_limit_s is not None and (time.time() - start) > time_limit_s:
            return
        R_size = R_mask.bit_count()
        if R_size + colour_bound(P_mask) < best_size:
            return
        if not P_mask and not X_mask:
            if R_size > best_size:
                best_size = R_size
                best_masks = [R_mask]
            elif R_size == best_size:
                best_masks.append(R_mask)
            return
        U = P_mask | X_mask
        if U:
            best_u, best_deg = 0, -1
            tmp = U
            while tmp:
                u_bit = (tmp & -tmp)
                idx = u_bit.bit_length()-1
                deg = (P_mask & adj[idx]).bit_count()
                if deg > best_deg:
                    best_u, best_deg = idx, deg
                tmp &= ~u_bit
            candidates = P_mask & ~adj[best_u]
        else:
            candidates = P_mask

        tmp = candidates
        while tmp:
            v_bit = (tmp & -tmp)
            v = v_bit.bit_length()-1
            rec(R_mask | (1<<v), P_mask & adj[v], X_mask & adj[v])
            P_mask &= ~v_bit
            X_mask |= v_bit
            tmp &= ~v_bit

    rec(0, V_all, 0)
    return best_size, best_masks

def mask_to_mapping(nodes: List[_MPNode], mask: int) -> Dict[int,int]:
    m: Dict[int,int] = {}
    while mask:
        v_bit = (mask & -mask)
        v = v_bit.bit_length()-1
        n = nodes[v]
        m[n.qi] = n.tj
        mask &= ~v_bit
    return dict(sorted(m.items()))

def largest_connected_in_query(q: Chem.Mol, mapping: Dict[int,int]) -> Dict[int,int]:
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
        for n in q.GetAtomWithIdx(qi).GetNeighbors():
            j = n.GetIdx()
            if j in pos: union(i, pos[j])
    from collections import Counter
    cnt = Counter(find(i) for i in range(len(mapped)))
    root = max(cnt.items(), key=lambda kv: kv[1])[0]
    keep = {mapped[i] for i in range(len(mapped)) if find(i)==root}
    return {qi:mapping[qi] for qi in sorted(keep)}

# -------- McGregor extend --------

def mcgregor_extend(m1: Chem.Mol, m2: Chem.Mol, seed: Dict[int,int],
                    C: ChemOptions, time_limit_s: Optional[float]=None) -> Dict[int,int]:
    import time
    start = time.time()
    best = dict(seed)
    rc1, rc2 = RingCache(m1), RingCache(m2)

    q_un = set(range(m1.GetNumAtoms())) - set(seed.keys())
    t_un = set(range(m2.GetNumAtoms())) - set(seed.values())
    QN = [[n.GetIdx() for n in m1.GetAtomWithIdx(i).GetNeighbors()] for i in range(m1.GetNumAtoms())]

    wl1, wl2 = wl_colours(m1), wl_colours(m2)

    def wl_bound() -> bool:
        from collections import Counter
        cq = Counter(wl1[i] for i in q_un)
        ct = Counter(wl2[j] for j in t_un)
        for c, k in cq.items():
            if ct.get(c, 0) < k: return False
        return True

    def feasible(qi: int, tj: int, cur: Dict[int,int]) -> bool:
        if not atoms_compatible(m1.GetAtomWithIdx(qi), m2.GetAtomWithIdx(tj), rc1, rc2, C):
            return False
        for qk, tl in cur.items():
            qb = m1.GetBondBetweenAtoms(int(qi), int(qk))
            tb = m2.GetBondBetweenAtoms(int(tj), int(tl))
            if (qb is None) != (tb is None):
                return False
            if qb is not None and not bonds_compatible(qb, tb, rc1, rc2, C):
                return False
        return True

    def frontier(cur: Dict[int,int]) -> Set[int]:
        F: Set[int] = set()
        for qk in cur.keys():
            for nn in QN[qk]:
                if nn in q_un: F.add(nn)
        return F if F else set(q_un)

    def upper_bound(cur_size: int) -> int:
        return cur_size + min(len(q_un), len(t_un))

    def rec(cur: Dict[int,int]):
        nonlocal best
        if time_limit_s is not None and (time.time() - start) > time_limit_s:
            return
        if not wl_bound():
            return
        if upper_bound(len(cur)) <= len(best):
            return
        if len(cur) > len(best):
            best = dict(cur)
        F = frontier(cur)
        if not F:
            return
        candidates_by_q: List[Tuple[int,List[int]]] = []
        for qi in F:
            ts = [tj for tj in t_un if feasible(qi, tj, cur)]
            if not ts: continue
            candidates_by_q.append((qi, ts))
        if not candidates_by_q:
            return
        candidates_by_q.sort(key=lambda x: (len(x[1]),
                                            -m1.GetAtomWithIdx(x[0]).GetDegree(),
                                            -int(m1.GetAtomWithIdx(x[0]).IsInRing())))
        qi, ts = candidates_by_q[0]
        for tj in ts:
            cur[qi] = tj
            q_un.remove(qi); t_un.remove(tj)
            rec(cur)
            del cur[qi]
            q_un.add(qi); t_un.add(tj)

    rec(dict(seed))
    return best
