# SPDX-License-Identifier: Apache-2.0
# © 2025 BioInception PVT LTD.

"""
Benchmark runner for SMSD vs RDKit reference (substructure and FMCS).

- Writes all outputs under test_output/
- CSV: bench_results.csv
- JSON: bench_results.json
- Cases: bench_cases.csv and bench_cases.json
- PNGs: selected visual comparisons

Design goals:
- Emphasise connected mappings (chemist-friendly). Disconnected is off unless explicitly enabled.
- Offer RDKit-like chemistry profiles for like-for-like comparison.
- Include timings (ms) and Tanimoto-like scores on atoms and bonds.
"""
import os, time, csv, json, math, re
from dataclasses import dataclass
from typing import List, Tuple

from rdkit import Chem, RDLogger
from rdkit.Chem import rdFMCS

from smsd_pro.engines import SMSD, MCSOptions, SubstructureOptions
from smsd_pro.chem import ChemOptions, rdkit_substructure_profile, rdkit_fmcs_profile
from smsd_pro.viz import viz_compare_substructure_png, viz_compare_mcs_png, VizConfig

RDLogger.DisableLog('rdApp.error')

ROOT = os.path.dirname(os.path.dirname(__file__))
OUT = os.path.join(ROOT, "test_output")
os.makedirs(OUT, exist_ok=True)

@dataclass
class Case:
    name: str
    query: str
    target: str
    chem: ChemOptions
    notes: str = ""


# ----------------------
# Chemistry presets
# ----------------------
CHEM_DEFAULT = ChemOptions()
CHEM_FLEX    = rdkit_substructure_profile()
CHEM_STEREO  = ChemOptions(bond_stereo="exact", use_chirality=True)
CHEM_NOCHARGE = ChemOptions(match_formal_charge=False)


def _tanimoto(a_common: int, a1: int, a2: int) -> float:
    denom = max(1, a1 + a2 - a_common)
    return a_common / denom


def _bonds_common_from_mapping(q: Chem.Mol, t: Chem.Mol, mp: dict) -> int:
    c = 0
    for b in q.GetBonds():
        i, j = int(b.GetBeginAtomIdx()), int(b.GetEndAtomIdx())
        if i in mp and j in mp and t.GetBondBetweenAtoms(int(mp[i]), int(mp[j])) is not None:
            c += 1
    return c


def build_cases() -> List[Case]:
    cases: List[Case] = []

    # 1) Chains of varying sizes
    for k in range(2, 12):
        q = "C"*k
        t = "C"*(k+2)
        cases.append(Case(f"Chain {k} in {k+2}", q, t, CHEM_DEFAULT))

    # 2) Benzene in naphthalene (strict vs flexible)
    cases.append(Case("Benzene in Naphthalene (strict)", "c1ccccc1", "c1ccc2ccccc2c1", CHEM_DEFAULT))
    cases.append(Case("Benzene in Naphthalene (flexible)", "c1ccccc1", "C1=CC=CC=C1c2ccccc2", CHEM_FLEX))

    # 3) SMARTS stereo: trans/cis (parity with RDKit)
    cases.append(Case("Trans-2-butene SMARTS vs trans", "C/C=C\\C", "C/C=C\\C", CHEM_STEREO))
    cases.append(Case("Trans-2-butene SMARTS vs cis",   "C/C=C\\C", "C/C=C/C",  CHEM_STEREO))

    # 4) Disconnected query (C.C) in ethane (dedup by target set → one)
    cases.append(Case("C.C in Ethane", "C.C", "CC", CHEM_DEFAULT))

    # 5) Ring sizes and ring-to-ring constraints (5–8 member) – MCSS MCS like RDKit
    for n in range(5, 9):
        cases.append(Case(f"Ring size {n} in polycycle", f"C1" + "C"*(n-1) + "1", "c1ccc2ccccc2c1", CHEM_DEFAULT))

    # 6) Functional groups (amide & charge-flexible carboxylate/acid)
    cases.append(Case("Amide in peptide", "NC(=O)", "NCC(=O)NCCC", CHEM_DEFAULT))
    cases.append(Case("Carboxylate vs acid", "[O-]C=O", "OC=O", ChemOptions(**{**CHEM_FLEX.__dict__, "match_formal_charge": False})))

    # 7) Fusion strictness
    CHEM_FUSED = ChemOptions(ring_fusion_strict=True, ring_matches_ring_only=True)
    cases.append(Case("Ring fusion strictness", "c1ccc2ccccc2c1", "c1ccc2c(c1)cccc2", CHEM_FUSED))

    # 8) Chirality required – enforce in both engines
    cases.append(Case("Chiral centre required", "[C@](F)(Cl)Br", "C(F)(Cl)Br", ChemOptions(**{**CHEM_DEFAULT.__dict__, "use_chirality": True})))

    # 9) Aromatic vs Kekulé MCS (flexible, MCSS like RDKit)
    cases.append(Case("Aromatic vs Kekulé (flexible)", "c1ccccc1", "C1=CC=CC=C1", CHEM_FLEX))

    # 10) Drug-like pair (realistic)
    cases.append(Case(
        "Drug-like pair",
        "NC(=O)c1[nH]c2ccccc2c1S(=O)(=O)N1CCOC(C(=O)N2CCc3c(Br)cccc3C2)C1",
        "NC(=O)c1[nH]c2ccccc2c1S(=O)(=O)N1CCOC(C(=O)NCCOc2ccccc2Br)C1",
        CHEM_DEFAULT
    ))

    # 11) Recursive SMARTS example
    cases.append(Case("Recursive SMARTS carbox", "[C;$(C(=O)O)]", "CC(=O)O", CHEM_DEFAULT,
                      notes="SMARTS with recursive $() anchoring"))

    return cases


def run() -> str:
    cases = build_cases()

    # Persist the input cases
    with open(os.path.join(OUT, "bench_cases.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["case","query","target","chem_profile","notes"])
        for c in cases:
            w.writerow([c.name, c.query, c.target, _chem_label(c.chem), c.notes])
    with open(os.path.join(OUT, "bench_cases.json"), "w") as f:
        json.dump([c.__dict__ | {"chem": _chem_label(c.chem)} for c in cases], f, indent=2)

    rows = []
    for c in cases:
        # --- SMSD substructure ---
        smsd = SMSD(c.query, c.target, chem=c.chem)
        t0 = time.time()
        ss_exists = smsd.substructure_exists(SubstructureOptions(connected_only=True))
        ss_smsd_ms = (time.time() - t0)*1000.0

        # RDKit reference substructure (strip $name so RDKit can parse the query)
        q_try = Chem.MolFromSmiles(c.query) or Chem.MolFromSmarts(_strip_named_preds(c.query))
        t_try = Chem.MolFromSmiles(c.target) or Chem.MolFromSmarts(_strip_named_preds(c.target))
        t0 = time.time()
        ss_rdkit = False
        if q_try and t_try:
            ss_rdkit = bool(t_try.HasSubstructMatch(q_try))
        ss_rdkit_ms = (time.time() - t0)*1000.0

        # --- SMSD MCS (connected common subgraph by default) ---
        t0 = time.time()
        mr = smsd.mcs_max(MCSOptions(mcs_type="MCCS", connected_only=True, use_mcgregor_extend=True))
        mcs_smsd_ms = (time.time() - t0)*1000.0
        if mr:
            m_atoms = len(mr.mapping)
            m_bonds = _bonds_common_from_mapping(smsd.q, smsd.t, mr.mapping)
            tan_a = _tanimoto(m_atoms, smsd.q.GetNumAtoms(), smsd.t.GetNumAtoms())
            tan_b = _tanimoto(m_bonds, smsd.q.GetNumBonds(), smsd.t.GetNumBonds())
            m_algo = mr.algorithm
        else:
            m_atoms = m_bonds = 0
            tan_a = tan_b = 0.0
            m_algo = ""

        # --- RDKit FMCS ---
        fmcs_atoms = fmcs_bonds = 0
        fmcs_tan_a = fmcs_tan_b = 0.0
        fmcs_ms = 0.0
        if q_try and t_try:
            p = rdFMCS.MCSParameters()
            p.MaximizeBonds = True; p.Timeout = 5; p.BondCompare = rdFMCS.BondCompare.CompareOrder
            t0 = time.time()
            res = rdFMCS.FindMCS([q_try, t_try], p)
            fmcs_ms = (time.time()-t0)*1000.0
            if not res.canceled and res.smartsString:
                patt = Chem.MolFromSmarts(res.smartsString)
                if patt:
                    fmcs_atoms, fmcs_bonds = patt.GetNumAtoms(), patt.GetNumBonds()
                    fmcs_tan_a = _tanimoto(fmcs_atoms, q_try.GetNumAtoms(), t_try.GetNumAtoms())
                    fmcs_tan_b = _tanimoto(fmcs_bonds, q_try.GetNumBonds(), t_try.GetNumBonds())

        rows.append({
            "case": c.name, "query": c.query, "target": c.target, "chem_profile": _chem_label(c.chem),
            "ss_smsd": int(ss_exists), "ss_smsd_ms": round(ss_smsd_ms, 3),
            "ss_rdkit": int(ss_rdkit), "ss_rdkit_ms": round(ss_rdkit_ms, 3),
            "mcs_smsd_atoms": m_atoms, "mcs_smsd_bonds": m_bonds,
            "mcs_smsd_tan_atoms": _round4(tan_a), "mcs_smsd_tan_bonds": _round4(tan_b),
            "mcs_smsd_algo": m_algo, "mcs_smsd_ms": round(mcs_smsd_ms, 3),
            "fmcs_atoms": fmcs_atoms, "fmcs_bonds": fmcs_bonds,
            "fmcs_tan_atoms": _round4(fmcs_tan_a), "fmcs_tan_bonds": _round4(fmcs_tan_b),
            "fmcs_ms": round(fmcs_ms, 3), "notes": c.notes
        })

    # Persist results
    csv_path = os.path.join(OUT, "bench_results.csv")
    json_path = os.path.join(OUT, "bench_results.json")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    with open(json_path, "w") as f:
        json.dump(rows, f, indent=2)

    # A few PNGs
    viz_compare_substructure_png("c1ccccc1","c1ccc2ccccc2c1",
                                 cfg=VizConfig(save_png=os.path.join(OUT, "benzene_in_naphthalene_sub.png")))
    viz_compare_mcs_png("c1ccccc1","c1ccc2ccccc2c1",
                        cfg=VizConfig(save_png=os.path.join(OUT, "benzene_in_naphthalene_mcs.png")))
    viz_compare_mcs_png(
        "NC(=O)c1[nH]c2ccccc2c1S(=O)(=O)N1CCOC(C(=O)N2CCc3c(Br)cccc3C2)C1",
        "NC(=O)c1[nH]c2ccccc2c1S(=O)(=O)N1CCOC(C(=O)NCCOc2ccccc2Br)C1",
        cfg=VizConfig(save_png=os.path.join(OUT, "Drug-like_pair_mcs.png"))
    )

    print(f"Wrote {csv_path}")
    print(f"Wrote {json_path}")
    print(f"PNG examples are in {OUT}")
    return csv_path


def _strip_named_preds(s: str) -> str:
    # remove $name but don't touch $()
    out = []
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == '[':
            depth, j = 1, i+1
            while j < n and depth > 0:
                if s[j] == '[': depth += 1
                elif s[j] == ']': depth -= 1
                j += 1
            content = s[i+1:j-1]
            content = re.sub(r"\$[A-Za-z_]\w*(?!\s*\()", "", content)
            content = re.sub(r";{2,}", ";", content).strip(";")
            out.append("[" + content + "]")
            i = j
        else:
            out.append(c); i += 1
    return "".join(out)


def _chem_label(c: ChemOptions) -> str:
    # compact description for logs
    return ",".join([
        f"arom={'flex' if c.aromaticity_mode=='flexible' else 'strict'}",
        f"bond={'loose' if c.match_bond_order=='loose' else 'strict'}",
        f"ring={'subset±%d'%c.ring_size_tolerance if c.ring_size_mode=='subset' else c.ring_size_mode}",
        f"fusion={'strict' if c.ring_fusion_strict else 'off'}",
        f"charge={'on' if c.match_formal_charge else 'off'}",
    ])


def _round4(x: float) -> float:
    return float(f"{x:.4f}")


if __name__ == "__main__":
    run()
