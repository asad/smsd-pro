# SPDX-License-Identifier: Apache-2.0
# © 2025 BioInception PVT LTD.
#!/usr/bin/env python
# scripts/bench.py — curated benchmark producing CSV/JSON/Markdown + PNGs under test_output/
import os, csv, time, json, argparse, shutil, hashlib
from dataclasses import dataclass, asdict
from typing import List, Tuple, Dict, Any, Optional

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdFMCS

from smsd_pro.engines import SMSD, SubstructureOptions, MCSOptions
from smsd_pro.chem import ChemOptions
from smsd_pro.viz import viz_compare_substructure_png, viz_compare_mcs_png, VizConfig

# hide RDKit parse noise for SMARTS edge-cases
RDLogger.DisableLog('rdApp.error')

ROOT = os.path.dirname(os.path.dirname(__file__))
OUT  = os.path.join(ROOT, "test_output")
os.makedirs(OUT, exist_ok=True)

@dataclass
class BenchRow:
    case: str
    query: str
    target: str
    chem_profile: str
    # Substructure (exists + timings for SMSD and RDKit)
    ss_smsd: int
    ss_smsd_ms: int
    ss_rdkit: int
    ss_rdkit_ms: int
    # MCS sizes/algos + timing
    mcs_smsd_atoms: int
    mcs_smsd_bonds: int
    mcs_smsd_tan_atoms: float
    mcs_smsd_tan_bonds: float
    mcs_smsd_algo: str
    mcs_smsd_ms: int
    fmcs_atoms: int
    fmcs_bonds: int
    fmcs_tan_atoms: float
    fmcs_tan_bonds: float
    fmcs_ms: int
    notes: str = ""

def _slug(x: str) -> str:
    return hashlib.sha1(x.encode()).hexdigest()[:8]

def _fmcs(m1: Chem.Mol, m2: Chem.Mol, timeout: int = 5):
    p = rdFMCS.MCSParameters()
    p.MaximizeBonds = True
    p.Timeout = timeout
    res = rdFMCS.FindMCS([m1, m2], p)
    patt = Chem.MolFromSmarts(res.smartsString) if (not res.canceled and res.smartsString) else None
    return patt

def _tan_atoms_bonds(q: Chem.Mol, t: Chem.Mol, m_atoms: int, m_bonds: int) -> Tuple[float,float]:
    a = m_atoms / max(1, (q.GetNumAtoms() + t.GetNumAtoms() - m_atoms))
    b = m_bonds / max(1, (q.GetNumBonds() + t.GetNumBonds() - m_bonds))
    return round(a, 4), round(b, 4)

def _bond_indices_from_atomset(m: Chem.Mol, atoms):
    aset = set(atoms)
    out = []
    for b in m.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in aset and j in aset:
            out.append(int(b.GetIdx()))
    return out

# ---- Case description ----
@dataclass(frozen=True)
class Case:
    name: str
    query: str
    target: str
    chem: ChemOptions = ChemOptions()
    notes: str = ""

def _chem_name(c: ChemOptions) -> str:
    parts = []
    if c.aromaticity_mode != 'strict': parts.append(f"arom={c.aromaticity_mode}")
    if c.match_bond_order != 'strict': parts.append(f"bond={c.match_bond_order}")
    if c.ring_size_mode != 'subset' or c.ring_size_tolerance != 0:
        parts.append(f"ring={c.ring_size_mode}±{c.ring_size_tolerance}")
    if c.bond_stereo != 'defined': parts.append(f"stereo={c.bond_stereo}")
    if c.use_chirality: parts.append("chirality")
    if c.ring_fusion_strict: parts.append("fusion=strict")
    if c.match_formal_charge is False: parts.append("nocharge")
    return ",".join(parts) or "default"

# ---- Runner helpers ----
def run_substruct(rowset: List[BenchRow], case: Case, limit_png: Optional[int], png_counter: List[int]):
    # SMSD
    smsd = SMSD(case.query, case.target, chem=case.chem)
    t0 = time.time()
    hit = smsd.substructure_exists(SubstructureOptions())
    ss_smsd_ms = int(1000*(time.time() - t0))

    # RDKit reference
    if isinstance(case.query, str):
        q_try = Chem.MolFromSmiles(case.query) or Chem.MolFromSmarts(case.query)
    else:
        q_try = case.query
    t_try = Chem.MolFromSmiles(case.target) if isinstance(case.target, str) else case.target

    t1 = time.time()
    rd_hits = list(t_try.GetSubstructMatches(q_try, uniquify=True)) if (q_try and t_try) else []
    ss_rdkit_ms = int(1000*(time.time() - t1))

    # MCS metrics (SMSD + RDKit FMCS)
    mr = smsd.mcs_max(MCSOptions(mcs_type="MCIS"))
    sms_atoms = int(getattr(mr, "size", 0) or 0)
    # common bonds for SMSD mapping on target
    sms_bonds = 0
    if mr:
        sms_bonds = sum(1 for b in smsd.q.GetBonds()
                        if int(b.GetBeginAtomIdx()) in mr.mapping and int(b.GetEndAtomIdx()) in mr.mapping and
                        smsd.t.GetBondBetweenAtoms(int(mr.mapping[int(b.GetBeginAtomIdx())]), int(mr.mapping[int(b.GetEndAtomIdx())])) is not None)
    tana, tanb = _tan_atoms_bonds(smsd.q, smsd.t, sms_atoms, sms_bonds)
    algo = getattr(mr, "algorithm", "") or ""
    # timing for SMSD MCS (fresh run for timing isolation)
    t2 = time.time(); _ = smsd.mcs_max(MCSOptions(mcs_type="MCIS")); mcs_smsd_ms = int(1000*(time.time()-t2))

    patt = _fmcs(smsd.q, smsd.t)
    if patt:
        fa, fb = patt.GetNumAtoms(), patt.GetNumBonds()
        t3 = time.time(); _ = rdFMCS.FindMCS([smsd.q, smsd.t], rdFMCS.MCSParameters()); fmcs_ms = int(1000*(time.time()-t3))
        fta, ftb = _tan_atoms_bonds(smsd.q, smsd.t, fa, fb)
    else:
        fa = fb = fmcs_ms = 0
        fta = ftb = 0.0

    # Save visuals (limit if requested)
    if (limit_png is None) or (png_counter[0] < int(limit_png)):
        viz_compare_substructure_png(case.query, case.target,
                                     cfg=VizConfig(save_png=os.path.join(OUT, f"{case.name}_sub.png")))
        viz_compare_mcs_png(case.query, case.target,
                            cfg=VizConfig(save_png=os.path.join(OUT, f"{case.name}_mcs.png")))
        png_counter[0] += 2

    rowset.append(BenchRow(
        case=case.name, query=case.query, target=case.target, chem_profile=_chem_name(case.chem),
        ss_smsd=int(bool(hit)), ss_smsd_ms=ss_smsd_ms,
        ss_rdkit=int(bool(rd_hits)), ss_rdkit_ms=ss_rdkit_ms,
        mcs_smsd_atoms=sms_atoms, mcs_smsd_bonds=sms_bonds,
        mcs_smsd_tan_atoms=tana, mcs_smsd_tan_bonds=tanb,
        mcs_smsd_algo=algo, mcs_smsd_ms=mcs_smsd_ms,
        fmcs_atoms=int(fa), fmcs_bonds=int(fb), fmcs_tan_atoms=fta, fmcs_tan_bonds=ftb, fmcs_ms=int(fmcs_ms),
        notes=case.notes
    ))

def main(only: List[str] | None = None, limit_png: int | None = None) -> None:
    rows: List[BenchRow] = []
    png_counter = [0]

    # Chemistry profiles
    CHEM_DEFAULT = ChemOptions()
    CHEM_FLEX    = ChemOptions(aromaticity_mode='flexible', match_bond_order='loose',
                               ring_size_mode='subset', ring_size_tolerance=1)
    CHEM_STEREO  = ChemOptions(bond_stereo='exact')
    CHEM_FUSED   = ChemOptions(ring_fusion_strict=True, ring_matches_ring_only=True)

    def C(name, q, t, chem=CHEM_DEFAULT, notes=""): return Case(name, q, t, chem, notes)

    # --- Curated + your requested cases ---
    cases: List[Case] = []
    # 1) Chains
    for k in range(2, 12):
        cases.append(C(f"Chain {k} in {k+2}", "C"*k, "C"*(k+2)))
    # 2) Benzene in naphthalene (strict vs flexible)
    cases.append(C("Benzene in Naphthalene (strict)", "c1ccccc1", "c1ccc2ccccc2c1", CHEM_DEFAULT))
    cases.append(C("Benzene in Naphthalene (flexible)", "c1ccccc1", "C1=CC=CC=C1c2ccccc2", CHEM_FLEX))
    # 3) SMARTS stereo: trans vs cis
    cases.append(C("Trans-2-butene SMARTS vs trans", "C/C=C\\C", "C/C=C\\C", CHEM_STEREO))
    cases.append(C("Trans-2-butene SMARTS vs cis",   "C/C=C\\C", "C/C=C/C",  CHEM_STEREO))
    # 4) Disconnected query
    cases.append(C("C.C in Ethane", "C.C", "CC"))
    # 5) Ring sizes & MCS checks
    for n in range(5, 9):
        cases.append(C(f"Ring size {n} in polycycle", f"C1" + "C"*(n-1) + "1", "c1ccc2ccccc2c1", CHEM_DEFAULT))
    # 6) Functional groups
    cases.append(C("Amide in peptide", "NC(=O)", "NCC(=O)NCCC", CHEM_DEFAULT))
    cases.append(C("Carboxylate vs acid", "[O-]C=O", "OC=O", ChemOptions(**{**CHEM_FLEX.__dict__, "match_formal_charge": False})))
    # 7) Fusion strictness
    cases.append(C("Ring fusion strictness", "c1ccc2ccccc2c1", "c1ccc2c(c1)cccc2", CHEM_FUSED))
    # 8) Chirality
    cases.append(C("Chiral centre required", "[C@](F)(Cl)Br", "C(F)(Cl)Br", ChemOptions(**{**CHEM_DEFAULT.__dict__, "use_chirality": True})))
    # 9) Aromatic vs Kekulé (flexible)
    cases.append(C("Aromatic vs Kekulé (flexible)", "c1ccccc1", "C1=CC=CC=C1", CHEM_FLEX))
    # 10) Drug-like pair
    cases.append(C("Drug-like pair",
                   "NC(=O)c1[nH]c2ccccc2c1S(=O)(=O)N1CCOC(C(=O)N2CCc3c(Br)cccc3C2)C1",
                   "NC(=O)c1[nH]c2ccccc2c1S(=O)(=O)N1CCOC(C(=O)NCCOc2ccccc2Br)C1"))

    # Extra: recursive SMARTS anchor case
    cases.append(C("Recursive SMARTS carbox", "[C;$(C(=O)O)]", "CC(=O)O", notes="SMARTS with recursive $() anchoring"))

    # Filter
    if only:
        names = set(only)
        cases = [c for c in cases if c.name in names]

    # Run
    for case in cases:
        run_substruct(rows, case, limit_png, png_counter)

    # Write CSV
    csv_path = os.path.join(OUT, "bench_cases.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([*BenchRow.__annotations__.keys()])
        for r in rows:
            w.writerow([r.case, r.query, r.target, r.chem_profile,
                        r.ss_smsd, r.ss_smsd_ms, r.ss_rdkit, r.ss_rdkit_ms,
                        r.mcs_smsd_atoms, r.mcs_smsd_bonds, r.mcs_smsd_tan_atoms, r.mcs_smsd_tan_bonds,
                        r.mcs_smsd_algo, r.mcs_smsd_ms, r.fmcs_atoms, r.fmcs_bonds, r.fmcs_tan_atoms, r.fmcs_tan_bonds, r.fmcs_ms, r.notes])

    # Write JSON
    json_path = os.path.join(OUT, "bench_cases.json")
    with open(json_path, "w") as f:
        json.dump([asdict(r) for r in rows], f, indent=2)

    # Index markdown with images
    md_path = os.path.join(OUT, "cases_index.md")
    with open(md_path, "w") as f:
        f.write("# SMSD-Pro Bench: Cases & Images\n\n")
        f.write(f"_Generated: {time.ctime()}_\n\n")
        for r in rows:
            sub = f"{r.case.replace(' ', '_')}_sub.png"
            mcs = f"{r.case.replace(' ', '_')}_mcs.png"
            ssub = os.path.join(OUT, sub)
            smcs = os.path.join(OUT, mcs)
            f.write(f"## {r.case}\n\n**Query:** `{r.query}`  \\n**Target:** `{r.target}`  \\n**Chem:** `{r.chem_profile}`\n\n")
            if os.path.exists(ssub):
                f.write(f"![{r.case} substructure]({sub})\n\n")
            if os.path.exists(smcs):
                f.write(f"![{r.case} mcs]({mcs})\n\n")

    # Snapshot this script for reproducibility (avoids confusion with the source file)
    try:
        shutil.copy2(__file__, os.path.join(OUT, "bench_snapshot.py"))
    except Exception:
        pass

    print(f"Wrote {len(rows)} rows to {csv_path}")
    if limit_png:
        print(f"Images written under {OUT} (capped to {limit_png} cases)")
    else:
        print(f"Images written under {OUT}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Run curated SMSD-Pro benchmarks; outputs live in test_output/")
    ap.add_argument("--only", nargs="*", help="Run only these case names (exact match)")
    ap.add_argument("--limit-png", type=int, default=None, help="Max number of cases to emit PNGs for (each case makes two PNGs)")
    args = ap.parse_args()
    main(args.only, args.limit_png)
