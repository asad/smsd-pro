import os, time, csv
from rdkit import Chem
from smsd_pro.engines import SMSD, MCSOptions, SubstructureOptions
from smsd_pro.chem import ChemOptions
from smsd_pro.viz import viz_compare_substructure_png, viz_compare_mcs_png, VizConfig

OUT = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_output")
os.makedirs(OUT, exist_ok=True)

def run():
    rows = []
    pairs = [
        ("C"*k, "C"*(k+2)) for k in range(2, 12)
    ] + [
        ("c1ccccc1", "c1ccc2ccccc2c1"),
        ("C/C=C\\C", "C/C=C\\C"),
        ("[C;$(C(=O)O)]","CC(=O)O"),
    ]
    for q,t in pairs:
        smsd = SMSD(q, t, chem=ChemOptions())
        t0 = time.time(); ok = smsd.substructure_exists(SubstructureOptions()); dt = time.time()-t0
        mr, tm = smsd.mcs_max(MCSOptions()), 0.0
        rows.append([q,t, int(ok), getattr(mr,'size',0), round(dt,4)])
    with open(os.path.join(OUT, "bench.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["query","target","ss_exists","mcs_size","ss_time"])
        w.writerows(rows)
    viz_compare_substructure_png("c1ccccc1","c1ccc2ccccc2c1", cfg=VizConfig(save_png=os.path.join(OUT,"viz_sub.png")))
    viz_compare_mcs_png("c1ccccc1","c1ccc2ccccc2c1", cfg=VizConfig(save_png=os.path.join(OUT,"viz_mcs.png")))

if __name__ == "__main__":
    run()
