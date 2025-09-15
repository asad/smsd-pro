from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple
import hashlib

try:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
except Exception:
    from rdkit import Chem
    from rdkit.Chem import rdMolStandardize

# ------------------
# Standardiser
# ------------------

class Standardiser:
    """Conservative, predictable standardisation for medicinal chemistry."""
    def __init__(self):
        self._lfc = rdMolStandardize.LargestFragmentChooser()
        self._norm = rdMolStandardize.Normalizer()
        self._reion = rdMolStandardize.Reionizer()
        self._unchg = rdMolStandardize.Uncharger()
        self._taut  = rdMolStandardize.TautomerEnumerator()

    def run(self, mol: Chem.Mol, *,
            largest_fragment: bool = True,
            normalise: bool = True,
            reionise: bool = True,
            uncharge: bool = True,
            canonical_tautomer: bool = True,
            kekulise: bool = False,
            add_hs: bool = False,
            remove_hs: bool = True,
            clear_stereo: bool = False) -> Chem.Mol:
        m = Chem.Mol(mol)
        Chem.SanitizeMol(m)
        if largest_fragment:
            m = self._lfc.choose(m)
        if normalise:
            m = self._norm.normalize(m)
        if reionise:
            m = self._reion.reionize(m)
        if uncharge:
            m = self._unchg.uncharge(m)
        if clear_stereo:
            Chem.RemoveStereochemistry(m)
        if add_hs:
            m = Chem.AddHs(m)
        if kekulise:
            try:
                Chem.Kekulize(m, clearAromaticFlags=True)
            except Exception:
                pass
        if remove_hs:
            m = Chem.RemoveHs(m)
        if canonical_tautomer:
            try:
                m = self._taut.Canonicalize(m)
            except Exception:
                pass
        Chem.SanitizeMol(m)
        return m

# ------------------
# Options
# ------------------

@dataclass(frozen=True)
class ChemOptions:
    # atom-level
    match_atom_type: bool = True
    match_isotope: bool = False
    match_formal_charge: bool = True
    match_valence: bool = False
    aromaticity_mode: str = "strict"  # "strict" | "flexible"
    ring_matches_ring_only: bool = True
    ring_fusion_strict: bool = False
    ring_size_mode: str = "subset"    # "exact" | "subset" | "ignore"
    ring_size_tolerance: int = 0
    # bond-level
    match_bond_order: str = "strict"  # "strict" | "loose"
    bond_stereo: str = "defined"      # "off" | "defined" | "exact"
    # stereo tetrahedral
    use_chirality: bool = False
    # global
    degree_slack: int = 0
    # ring completeness (final mapping validation)
    complete_rings_only: bool = False

# ------------------
# WL colours & rings
# ------------------

def wl_colours(m: Chem.Mol, rounds: int = 2, include_chirality: bool = False) -> List[int]:
    N = m.GetNumAtoms()
    base = []
    for i in range(N):
        a = m.GetAtomWithIdx(i)
        inv = (
            a.GetAtomicNum(),
            a.GetFormalCharge(),
            a.GetTotalValence(),
            a.GetTotalNumHs(includeNeighbors=True),
            int(a.GetIsAromatic()),
            int(a.IsInRing()),
            int(a.GetChiralTag()) if include_chirality else 0,
        )
        base.append(int.from_bytes(hashlib.blake2b("|".join(map(str, inv)).encode(),
                                                  digest_size=4).digest(), "little"))
    adj = [[] for _ in range(N)]
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        adj[i].append(j); adj[j].append(i)
    cur = base
    for _ in range(rounds):
        nxt = []
        for i in range(N):
            neigh = sorted(cur[j] for j in adj[i])
            nxt.append(int.from_bytes(
                hashlib.blake2b((str(cur[i])+";"+",".join(map(str, neigh))).encode(),
                                digest_size=4).digest(), "little"))
        cur = nxt
    return cur

class RingCache:
    """Pre-computed ring data for atoms and bonds, plus ring-size look-ups."""
    def __init__(self, m: Chem.Mol):
        self.m = m
        self.N = m.GetNumAtoms()
        self.E = m.GetNumBonds()
        ri = m.GetRingInfo()
        self.atom_ring_count = [0]*self.N
        self.bond_rings: List[Set[int]] = []
        self.bond_ring_sizes: Dict[int, Set[int]] = {i:set() for i in range(self.E)}
        try:
            for ring in ri.BondRings():
                s = set(ring)
                self.bond_rings.append(s)
                ring_size = len(ring)
                atoms_in_ring = set()
                for bi in ring:
                    self.bond_ring_sizes[bi].add(ring_size)
                    b = m.GetBondWithIdx(bi)
                    atoms_in_ring.add(b.GetBeginAtomIdx()); atoms_in_ring.add(b.GetEndAtomIdx())
                for ai in atoms_in_ring: self.atom_ring_count[ai] += 1
        except Exception:
            for ring in ri.AtomRings():
                atoms = set(ring)
                for ai in atoms: self.atom_ring_count[ai] += 1

# ------------------
# Comparators
# ------------------

def atoms_compatible(q: Chem.Atom, t: Chem.Atom, rc_q: RingCache, rc_t: RingCache, C: ChemOptions) -> bool:
    if C.match_atom_type and q.GetAtomicNum() != t.GetAtomicNum():
        return False
    if C.match_isotope and q.GetIsotope() != t.GetIsotope():
        return False
    if C.match_formal_charge and q.GetFormalCharge() != t.GetFormalCharge():
        return False
    if C.aromaticity_mode == "strict":
        if q.GetIsAromatic() != t.GetIsAromatic():
            return False
    if C.ring_matches_ring_only and q.IsInRing() and not t.IsInRing():
        return False
    if C.ring_fusion_strict:
        if rc_q.atom_ring_count[q.GetIdx()] != rc_t.atom_ring_count[t.GetIdx()]:
            return False
    if C.match_valence and q.GetTotalValence() != t.GetTotalValence():
        return False
    if C.use_chirality:
        if q.GetChiralTag() != t.GetChiralTag():
            return False
    return True

def ring_size_ok_bond(qb: Chem.Bond, tb: Chem.Bond, rc_q: RingCache, rc_t: RingCache, C: ChemOptions) -> bool:
    if C.ring_size_mode == "ignore":
        return True
    q_sizes = rc_q.bond_ring_sizes.get(qb.GetIdx(), set())
    t_sizes = rc_t.bond_ring_sizes.get(tb.GetIdx(), set())
    if not q_sizes or not t_sizes:
        return C.ring_size_mode != "exact"
    if C.ring_size_mode == "exact":
        return any(sz in t_sizes for sz in q_sizes)
    for sq in q_sizes:
        for st in t_sizes:
            if abs(sq - st) <= C.ring_size_tolerance:
                return True
    return False

def bond_stereo_ok(qb: Chem.Bond, tb: Chem.Bond, mode: str) -> bool:
    if mode == "off":
        return True
    if qb.GetBondType() != Chem.BondType.DOUBLE or tb.GetBondType() != Chem.BondType.DOUBLE:
        return True
    sQ, sT = qb.GetStereo(), tb.GetStereo()
    if mode == "defined":
        return sQ != Chem.BondStereo.STEREONONE and sT != Chem.BondStereo.STEREONONE
    if mode == "exact":
        if sQ == Chem.BondStereo.STEREOANY:
            return sT != Chem.BondStereo.STEREONONE
        return sQ == sT and sQ != Chem.BondStereo.STEREONONE
    return True

def bonds_compatible(qb: Optional[Chem.Bond], tb: Optional[Chem.Bond],
                     rc_q: RingCache, rc_t: RingCache, C: ChemOptions) -> bool:
    if (qb is None) != (tb is None):
        return False
    if qb is None and tb is None:
        return True
    if C.ring_matches_ring_only and qb.IsInRing() and not tb.IsInRing():
        return False
    if C.match_bond_order == "strict":
        if qb.GetBondType() != tb.GetBondType():
            return False
    else:
        if qb.GetBondType() != tb.GetBondType():
            aromish = (qb.GetIsAromatic() and tb.GetIsAromatic())
            if not aromish:
                return False
    if qb.IsInRing() and tb.IsInRing():
        if not ring_size_ok_bond(qb, tb, rc_q, rc_t, C):
            return False
    if not bond_stereo_ok(qb, tb, C.bond_stereo):
        return False
    return True
