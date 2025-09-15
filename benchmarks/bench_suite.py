import time, os
from dataclasses import dataclass
from typing import Tuple
from rdkit import Chem
from rdkit.Chem import rdFMCS
from smsd_pro.engines import SMSD
from smsd_pro.mcs import MCSOptions
from smsd_pro.vf import SubstructureOptions
from smsd_pro.chem import ChemOptions
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(__file__), "out")
os.makedirs(OUT, exist_ok=True)

CHEM_DEFAULT = ChemOptions()

@dataclass
class Case:
    name: str
    q: str
    t: str

CASES = [Case(f"Chain {k} in {k+2}", "C"*k, "C"*(k+2)) for k in range(2, 18)]
CASES += [Case("Benzene in naph", "c1ccccc1", "c1ccc2ccccc2c1")]

def rdkit_fmcs(m1, m2, timeout=5):
    p = rdFMCS.MCSParameters()
    p.MaximizeBonds = True
    p.Timeout = int(timeout)
    p.AtomCompare = rdFMCS.AtomCompare.CompareElements
    p.BondCompare = rdFMCS.BondCompare.CompareOrder
    res = rdFMCS.FindMCS([m1, m2], p)
    if res.canceled: return (0,0)
    patt = Chem.MolFromSmarts(res.smartsString)
    return (patt.GetNumAtoms(), patt.GetNumBonds()) if patt else (0,0)

def main():
    xs = []; ys_smsd = []; ys_fmcs = []
    for cs in CASES:
        smsd = SMSD(cs.q, cs.t, chem=CHEM_DEFAULT)
        t0 = time.perf_counter()
        mr = smsd.mcs_max(MCSOptions(mcs_type="MCIS", time_limit_s=3.0))
        dt_smsd = time.perf_counter() - t0
        m1, m2 = Chem.MolFromSmiles(cs.q), Chem.MolFromSmiles(cs.t)
        t0 = time.perf_counter()
        a,b = rdkit_fmcs(m1, m2, timeout=3)
        dt_rd = time.perf_counter() - t0
        xs.append(cs.name)
        ys_smsd.append(dt_smsd)
        ys_fmcs.append(dt_rd)
        print(f"{cs.name:>20}: SMSD {dt_smsd:.4f}s | FMCS {dt_rd:.4f}s (FMCS atoms={a})")

    plt.figure(figsize=(10,4))
    plt.plot(range(len(xs)), ys_smsd, label="SMSD Pro (MCIS)")
    plt.plot(range(len(xs)), ys_fmcs, label="RDKit FMCS")
    plt.xticks(range(len(xs)), xs, rotation=70, fontsize=8)
    plt.ylabel("Seconds")
    plt.legend()
    plt.tight_layout()
    fp = os.path.join(OUT, "benchmark_times.png")
    plt.savefig(fp, dpi=160)
    print("Saved:", fp)

if __name__ == "__main__":
    main()
