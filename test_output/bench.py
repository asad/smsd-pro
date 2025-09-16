# SPDX-License-Identifier: Apache-2.0
# © 2025 BioInception PVT LTD.
#!/usr/bin/env python
# scripts/bench.py — curated benchmark producing CSV + PNGs under test_output/
import os, csv, time, json, shutil, argparse, hashlib
from dataclasses import dataclass, asdict
from typing import List, Tuple, Dict, Any

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import rdFMCS

from smsd_pro.engines import SMSD, SubstructureOptions, MCSOptions
from smsd_pro.chem import ChemOptions
from smsd_pro.viz import viz_compare_substructure_png, viz_compare_mcs_png, VizConfig

# quiet RDKit parse chatter
RDLogger.DisableLog('rdApp.error')

ROOT = os.path.dirname(os.path.dirname(__file__))
OUT  = os.path.join(ROOT, "test_output")
os.makedirs(OUT, exist_ok=True)

@dataclass
class BenchRow:
    case: str
    query: str
    target: str
    ss_exists: int
    ss_time_ms: int
    mcs_size_atoms: int
    mcs_algo: str
    fmcs_atoms: int
    fmcs_bonds: int
    notes: str = ""

def _slug(x: str) -> str:
    return hashlib.sha1(x.encode()).hexdigest()[:8]

def _fmcs_sizes(m1: Chem.Mol, m2: Chem.Mol) -> Tuple[int,int]:
    p = rdFMCS.MCSParameters()
    p.MaximizeBonds = True
    p.Timeout = 5
    res = rdFMCS.FindMCS([m1, m2], p)
    if res.canceled or not res.smartsString:
        return 0, 0
    patt = Chem.MolFromSmarts(res.smartsString)
    if not patt:
        return 0, 0
    return int(patt.GetNumAtoms()), int(patt.GetNumBonds())

def main(selected: List[str] | None = None, limit_png: int | None = None) -> None:
    # Curated cases: include SMARTS, RDKit FMCS toughies, and blog/discussion examples.
    cases: List[Tuple[str,str,str,str]] = [
        # (name, query, target, notes)
        ("smartscarbox", "[C;$(C(=O)O)]", "CC(=O)O", "SMARTS with recursive $() anchoring"),
        ("benzene_in_naphthalene", "c1ccccc1", "c1ccc2ccccc2c1", "Aromatics; ring subset"),
        ("stereo_butene_exact", "C/C=C\\C", "C/C=C\\C", "E/Z equality"),
        # RDKit blog: para vs ortho linkage
        ("rdkit_blog_degenerate", "Nc1ccc(cc1)C-Cc1c(N)cccc1", "Nc1ccc(cc1)C=Cc1c(N)cccc1", "RDKit MCS degeneracy"),
        # RDKit GitHub #1585: slow pair
        ("gh_1585_slowA", "c1cc(c(c(c1)Cl)N2c3cc(cc(c3CNC2=O)c4ccc(cc4F)F)N5CCNCC5)Cl", "CCNc1cc(c2c(c1)N(C(=O)NC2)c3ccc(cc3)n4ccc-5ncnc5c4)c6ccnnc6", "GitHub #1585 example A"),
        # ring options discussion
        ("mlist_ring_opts", "C(NC1=NCCN1)c1cccs1", "OCCNC(=O)c1cccs1", "rdkit-discuss ringMatches example"),
        # Java tests (representative)
        ("java_like_amide", "CC(=O)NCC", "CC(=O)NCCO", "amide/amine growth"),
    ]

    # Add size‑graded alkane chain tests
    for k in range(2, 12):
        cases.append((f"alk_{k}", "C"*k, "C"*(k+2), "alkane growth"))

    # Filtering by selection name(s)
    if selected:
        names = set(selected)
        cases = [c for c in cases if c[0] in names]

    rows: List[BenchRow] = []
    png_count = 0

    for name, q, t, notes in cases:
        smsd = SMSD(q, t, chem=ChemOptions())
        t0 = time.time()
        ok = smsd.substructure_exists(SubstructureOptions())
        ss_ms = int(1000*(time.time() - t0))

        mr = smsd.mcs_max(MCSOptions(mcs_type="MCIS"))
        mcs_atoms = int(getattr(mr, "size", 0) or 0)
        mcs_algo  = getattr(mr, "algorithm", "") or ""

        # ref FMCS sizes
        m1 = smsd.q; m2 = smsd.t
        fa, fb = _fmcs_sizes(m1, m2)

        rows.append(BenchRow(
            case=name, query=q, target=t, ss_exists=int(bool(ok)), ss_time_ms=ss_ms,
            mcs_size_atoms=mcs_atoms, mcs_algo=mcs_algo, fmcs_atoms=fa, fmcs_bonds=fb, notes=notes
        ))

        # Save SMSD vs RDKit side-by-side PNGs (cap if limit_png is set)
        if limit_png is None or png_count < int(limit_png):
            viz_compare_substructure_png(q, t, cfg=VizConfig(save_png=os.path.join(OUT, f"{name}_sub.png")))
            viz_compare_mcs_png(q, t, cfg=VizConfig(save_png=os.path.join(OUT, f"{name}_mcs.png")))
            png_count += 2

    # Write CSV
    csv_path = os.path.join(OUT, "bench_cases.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([*BenchRow.__annotations__.keys()])
        for r in rows:
            w.writerow([r.case, r.query, r.target, r.ss_exists, r.ss_time_ms, r.mcs_size_atoms, r.mcs_algo, r.fmcs_atoms, r.fmcs_bonds, r.notes])

    # Write JSON
    json_path = os.path.join(OUT, "bench_cases.json")
    with open(json_path, "w") as f:
        json.dump([asdict(r) for r in rows], f, indent=2)

    # Write an index with image links
    md_path = os.path.join(OUT, "cases_index.md")
    with open(md_path, "w") as f:
        f.write("# SMSD-Pro Bench: Cases & Images\n\n")
        f.write(f"_Generated: {time.ctime()}_\n\n")
        for r in rows:
            sub = f"{r.case}_sub.png"
            mcs = f"{r.case}_mcs.png"
            f.write(f"## {r.case}\n\n**Query:** `{r.query}`  \n**Target:** `{r.target}`  \n\n")
            if os.path.exists(os.path.join(OUT, sub)):
                f.write(f"![{r.case} substructure]({sub})\n\n")
            if os.path.exists(os.path.join(OUT, mcs)):
                f.write(f"![{r.case} mcs]({mcs})\n\n")

    # Copy this script to test_output for reproducibility, as requested
    try:
        shutil.copy2(__file__, os.path.join(OUT, "bench.py"))
    except Exception:
        pass

    print(f"Wrote {len(rows)} rows to {csv_path}")
    print(f"Images written under {OUT} (capped to {limit_png} cases)" if limit_png else f"Images written under {OUT}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Run curated SMSD-Pro benchmarks; outputs live in test_output/")
    ap.add_argument("--only", nargs="*", help="Run only these case names")
    ap.add_argument("--limit-png", type=int, default=None, help="Max number of cases to emit PNGs for (each case makes two PNGs)")
    args = ap.parse_args()
    main(args.only, args.limit_png)
