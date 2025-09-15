\
# scripts/bench.py – small benchmark producing a CSV and two PNGs
from __future__ import annotations
import time, csv, os
from statistics import median
from rdkit import Chem
from smsd_pro.engines import SMSD, SubstructureOptions, MCSOptions
from smsd_pro.chem import ChemOptions
from smsd_pro.viz import viz_compare_substructure_png, viz_compare_mcs_png, VizConfig
import matplotlib.pyplot as plt

CASES = [
    ("Chain 6 in 8", "CCCCCC", "CCCCCCCC"),
    ("Benzene in Naphthalene", "c1ccccc1", "c1ccc2ccccc2c1"),
    ("Amide in peptide", "NC(=O)", "NCC(=O)NCCC"),
    ("Carboxylate vs acid", "[O-]C=O", "OC=O"),
]

def run_once(name, q, t):
    chem = ChemOptions()
    smsd = SMSD(q, t, chem=chem)
    t0 = time.perf_counter()
    ok = smsd.substructure_exists(SubstructureOptions(connected_only=True))
    dt1 = time.perf_counter() - t0

    t0 = time.perf_counter()
    mr = smsd.mcs_max(MCSOptions(mcs_type="MCIS"))
    dt2 = time.perf_counter() - t0
    return {"name": name, "sub_ms": 1000*dt1, "mcs_ms": 1000*dt2, "sub_ok": bool(ok), "mcs_size": (mr.size if mr else 0)}

def main(out_dir="test_output"):
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    for name, q, t in CASES:
        for _ in range(5):
            rows.append(run_once(name, q, t))
    csv_path = os.path.join(out_dir, "bench.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name","sub_ms","mcs_ms","sub_ok","mcs_size"])
        w.writeheader(); w.writerows(rows)

    # summary bar chart
    labels = [name for (name,_,_) in CASES]
    sub_med = []
    mcs_med = []
    for name in labels:
        xs = [r["sub_ms"] for r in rows if r["name"]==name]
        ys = [r["mcs_ms"] for r in rows if r["name"]==name]
        sub_med.append(median(xs)); mcs_med.append(median(ys))

    plt.figure(figsize=(7,4))
    plt.title("SMSD Pro – median runtimes (ms)")
    plt.bar(labels, sub_med)
    plt.xticks(rotation=20, ha="right")
    plt.ylabel("ms")
    plt.tight_layout()
    png1 = os.path.join(out_dir, "bench_sub_ms.png")
    plt.savefig(png1); plt.close()

    plt.figure(figsize=(7,4))
    plt.title("SMSD Pro – median MCS runtimes (ms)")
    plt.bar(labels, mcs_med)
    plt.xticks(rotation=20, ha="right")
    plt.ylabel("ms")
    plt.tight_layout()
    png2 = os.path.join(out_dir, "bench_mcs_ms.png")
    plt.savefig(png2); plt.close()

    print("Wrote:", csv_path, png1, png2)

if __name__ == "__main__":
    main()
