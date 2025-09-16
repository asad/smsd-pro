# SPDX-License-Identifier: Apache-2.0
# © 2025 BioInception PVT LTD.

"""Evaluator for recursive SMARTS predicates: both $name and $(SMARTS) blocks.

We anchor $(SMARTS) so the first atom in the pattern must map to the candidate atom.

Named predicates can be registered either as functions or as SMARTS strings.

"""

from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Dict, Optional, Tuple, List
from functools import lru_cache
from rdkit import Chem


@dataclass(frozen=True)
class RecSpec:
    kind: str  # "smarts" or "name"
    expr: str


def _anchor_match_smartspattern(smarts: str, mol: Chem.Mol, atom_idx: int) -> bool:
    """Return True iff the first atom of SMARTS `smarts` maps to `atom_idx` in `mol`."""
    patt = Chem.MolFromSmarts(smarts)
    if not patt:
        return False
    for match in mol.GetSubstructMatches(patt, uniquify=True):
        if match and match[0] == atom_idx:
            return True
    return False


class RecPredicateVM:
    """Cached evaluator for recursive SMARTS $() and named $predicates."""
    def __init__(self, anchor_matcher: Optional[Callable[[str, Chem.Mol, int], bool]] = None):
        self._named: Dict[str, Callable[[Chem.Mol, int], bool]] = {}
        self._anchor = anchor_matcher or _anchor_match_smartspattern

    def register(self, name: str, *, fn: Optional[Callable[[Chem.Mol, int], bool]] = None,
                 smarts: Optional[str] = None) -> None:
        if fn is None and smarts is None:
            raise ValueError("Provide either fn= or smarts= for a named predicate.")
        if smarts is not None:
            def _fn(m: Chem.Mol, i: int, _s=smarts):
                return self.eval_smarts(m, i, _s)
            self._named[name] = lru_cache(maxsize=50000)(_fn)
        else:
            self._named[name] = lru_cache(maxsize=50000)(fn)  # type: ignore[arg-type]

    @lru_cache(maxsize=100000)
    def _cache_named(self, mol_id: int, atom_idx: int, name: str) -> bool:
        if name not in self._named:
            raise KeyError(f"Unknown recursive predicate: ${name}")
        return self._named[name](self._mol_lookup[mol_id], atom_idx)

    @lru_cache(maxsize=100000)
    def _cache_smarts(self, mol_id: int, atom_idx: int, smarts: str) -> bool:
        return self._anchor(smarts, self._mol_lookup[mol_id], atom_idx)

    _mol_lookup: Dict[int, Chem.Mol] = {}

    @classmethod
    def _remember_mol(cls, mol: Chem.Mol) -> int:
        mid = id(mol)
        if mid not in cls._mol_lookup:
            cls._mol_lookup[mid] = mol
        return mid

    def eval_named(self, mol: Chem.Mol, atom_idx: int, name: str) -> bool:
        mid = self._remember_mol(mol)
        return self._cache_named(mid, atom_idx, name)

    def eval_smarts(self, mol: Chem.Mol, atom_idx: int, smarts: str) -> bool:
        mid = self._remember_mol(mol)
        return self._cache_smarts(mid, atom_idx, smarts)

    def eval_all_for_atom(self, mol: Chem.Mol, atom_idx: int, specs: List[RecSpec]) -> bool:
        for spec in specs:
            if spec.kind == "name":
                if not self.eval_named(mol, atom_idx, spec.expr):
                    return False
            else:
                if not self.eval_smarts(mol, atom_idx, spec.expr):
                    return False
        return True
