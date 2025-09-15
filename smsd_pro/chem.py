\
# smsd_pro/chem.py – chemistry utilities, invariants, comparators and SMARTS extraction
from __future__ import annotations
from dataclasses import dataclass, replace
from typing import Dict, List, Tuple, Optional, Iterable, Set, Callable, Literal, Any
import hashlib, functools, re

from rdkit import Chem
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
except Exception:
    from rdkit.Chem import rdMolStandardize  # type: ignore

# ----------------
# Standardisation
# ----------------

class Standardiser:
    """Conservative standardisation for med‑chem use. Safe on most inputs."""
    def __init__(self):
        self._lfc = rdMolStandardize.LargestFragmentChooser()
        self._norm = rdMolStandardize.Normalizer()
        self._reion = rdMolStandardize.Reionizer()
        self._unchg = rdMolStandardize.Uncharger()
        self._taut  = rdMolStandardize.TautomerEnumerator()

    def run(self, mol: Chem.Mol, *, largest_fragment=True, normalise=True, reionise=True,
            uncharge=True, canonical_tautomer=True, kekulise=False, add_hs=False,
            remove_hs=True, clear_stereo=False) -> Chem.Mol:
        m = Chem.Mol(mol)
        Chem.SanitizeMol(m)
        if largest_fragment: m = self._lfc.choose(m)
        if normalise:        m = self._norm.normalize(m)
        if reionise:         m = self._reion.reionize(m)
        if uncharge:         m = self._unchg.uncharge(m)
        if clear_stereo:     Chem.RemoveStereochemistry(m)
        if add_hs:           m = Chem.AddHs(m)
        if kekulise:
            try: Chem.Kekulize(m, clearAromaticFlags=True)
            except Exception: pass
        if remove_hs:        m = Chem.RemoveHs(m)
        if canonical_tautomer:
            try: m = self._taut.Canonicalize(m)
            except Exception: pass
        Chem.SanitizeMol(m)
        return m

# -----------------------------
# WL colours and ring caches
# -----------------------------

def _wl_colours(m: Chem.Mol, rounds: int = 2, include_chirality: bool = False) -> List[int]:
    N = m.GetNumAtoms()
    base: List[int] = []
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
    """Rings per atom and per bond, plus ring sizes, cached once per molecule."""
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

# -------------------------
# Chemistry match options
# -------------------------

@dataclass(frozen=True)
class ChemOptions:
    # atom‑level
    match_atom_type: bool = True
    match_isotope: bool = False
    match_formal_charge: bool = True
    match_valence: bool = False
    aromaticity_mode: Literal["strict","flexible"] = "strict"
    ring_matches_ring_only: bool = True
    ring_fusion_strict: bool = False
    ring_size_mode: Literal["exact","subset","ignore"] = "subset"
    ring_size_tolerance: int = 0
    # bond‑level
    match_bond_order: Literal["strict","loose"] = "strict"
    bond_stereo: Literal["off","defined","exact"] = "defined"
    # stereo tetrahedral
    use_chirality: bool = False
    # global
    degree_slack: int = 0
    # ring completeness (final mapping validation)
    complete_rings_only: bool = False
    # NEW: optional hook for recursive SMARTS
    atom_rec_ok: Optional[Callable[[Chem.Mol, Chem.Mol, int, int], bool]] = None

# ----------------
# SMARTS extractor
# ----------------

@dataclass(frozen=True)
class RecSpec:
    kind: str  # "smarts" or "name"
    expr: str

class SmartsPattern:
    """Extract constraints from SMARTS/SMILES and detect $() per atom."""
    def __init__(self, smarts_or_smiles: str, *, is_smarts: bool = True):
        self.src = smarts_or_smiles
        self.is_smarts = is_smarts
        self.m = Chem.MolFromSmarts(smarts_or_smiles) if is_smarts else Chem.MolFromSmiles(smarts_or_smiles)
        if self.m is None:
            raise ValueError(f"Cannot parse: {smarts_or_smiles}")
        self.N = self.m.GetNumAtoms()
        self.E = self.m.GetNumBonds()
        self.rec_by_atom: Dict[int, List[RecSpec]] = self._extract_recursive_specs(smarts_or_smiles) if is_smarts else {}

    @staticmethod
    def _extract_recursive_specs(smarts: str) -> Dict[int, List[RecSpec]]:
        two_letter = {"Cl","Br","Si","Se","Na","Li","Mg","Al","Ca","Fe","Zn","Sn","Ag","Pt","Au","Hg","Pb","Mn","Co","Ni","Cu"}
        i, n = 0, len(smarts)
        atom_idx = 0
        out: Dict[int, List[RecSpec]] = {}

        def parse_paren_block(s: str, start: int) -> Tuple[str, int]:
            depth, j = 1, start + 1
            while j < len(s) and depth > 0:
                c = s[j]
                if c == '(':
                    depth += 1
                elif c == ')':
                    depth -= 1
                j += 1
            return s[start+1:j-1], j

        while i < n:
            c = smarts[i]
            if c == '[':
                depth, j = 1, i+1
                while j < n and depth > 0:
                    if smarts[j] == '[': depth += 1
                    elif smarts[j] == ']': depth -= 1
                    j += 1
                content = smarts[i+1:j-1]
                k = 0
                while k < len(content):
                    pos = content.find('$(', k)
                    if pos == -1: break
                    block, nx = parse_paren_block(content, pos+1)
                    expr = block.strip()
                    if re.fullmatch(r'[A-Za-z_]\w*', expr):
                        out.setdefault(atom_idx, []).append(RecSpec("name", expr))
                    else:
                        out.setdefault(atom_idx, []).append(RecSpec("smarts", expr))
                    k = nx
                atom_idx += 1
                i = j
                continue

            if c.isalpha():
                if i+1 < n and smarts[i:i+2] in two_letter:
                    i += 2
                else:
                    i += 1
                atom_idx += 1
                continue

            i += 1
        return out

# ----------------------------
# Atom and bond compat checks
# ----------------------------

def _atoms_compatible(q: Chem.Atom, t: Chem.Atom, rc_q: RingCache, rc_t: RingCache, C: ChemOptions) -> bool:
    if C.match_atom_type and q.GetAtomicNum() != t.GetAtomicNum(): return False
    if C.match_isotope and q.GetIsotope() != t.GetIsotope(): return False
    if C.match_formal_charge and q.GetFormalCharge() != t.GetFormalCharge(): return False
    if C.aromaticity_mode == "strict":
        if q.GetIsAromatic() != t.GetIsAromatic(): return False
    if C.ring_matches_ring_only and q.IsInRing() and not t.IsInRing(): return False
    if C.ring_fusion_strict:
        if rc_q.atom_ring_count[q.GetIdx()] != rc_t.atom_ring_count[t.GetIdx()]: return False
    if C.match_valence and q.GetTotalValence() != t.GetTotalValence(): return False
    if C.use_chirality and q.GetChiralTag() != t.GetChiralTag(): return False
    # Recursive SMARTS hook (if present)
    if C.atom_rec_ok is not None:
        qmol = rc_q.m if hasattr(rc_q, "m") else q.GetOwningMol()
        tmol = rc_t.m if hasattr(rc_t, "m") else t.GetOwningMol()
        if not C.atom_rec_ok(qmol, tmol, int(q.GetIdx()), int(t.GetIdx())):
            return False
    return True

def _ring_size_ok_bond(qb: Chem.Bond, tb: Chem.Bond, rc_q: RingCache, rc_t: RingCache, C: ChemOptions) -> bool:
    if C.ring_size_mode == "ignore": return True
    q_sizes = rc_q.bond_ring_sizes.get(int(qb.GetIdx()), set())
    t_sizes = rc_t.bond_ring_sizes.get(int(tb.GetIdx()), set())
    if not q_sizes or not t_sizes:
        return C.ring_size_mode != "exact"
    if C.ring_size_mode == "exact":
        return any(sz in t_sizes for sz in q_sizes)
    for sq in q_sizes:
        for st in t_sizes:
            if abs(sq - st) <= C.ring_size_tolerance:
                return True
    return False

def _bond_stereo_ok(qb: Chem.Bond, tb: Chem.Bond, mode: Literal["off","defined","exact"]) -> bool:
    if mode == "off": return True
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

def _bonds_compatible(qb: Optional[Chem.Bond], tb: Optional[Chem.Bond],
                      rc_q: RingCache, rc_t: RingCache, C: ChemOptions) -> bool:
    if (qb is None) != (tb is None): return False
    if qb is None: return True
    if C.ring_matches_ring_only and qb.IsInRing() and not (tb and tb.IsInRing()): return False
    if C.match_bond_order == "strict":
        if qb.GetBondType() != tb.GetBondType(): return False
    else:
        if qb.GetBondType() != tb.GetBondType():
            if not (qb.GetIsAromatic() and tb.GetIsAromatic()): return False
    if qb.IsInRing() and tb.IsInRing():
        if not _ring_size_ok_bond(qb, tb, rc_q, rc_t, C): return False
    if not _bond_stereo_ok(qb, tb, C.bond_stereo): return False
    return True
