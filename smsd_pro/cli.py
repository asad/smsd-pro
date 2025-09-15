\
# smsd_pro/cli.py – small CLI to run substructure/MCS and export PNGs
from __future__ import annotations
import argparse, json, time, sys
from rdkit import Chem
from .engines import SMSD, SubstructureOptions, MCSOptions
from .chem import ChemOptions
from .viz import viz_compare_substructure_png, viz_compare_mcs_png, VizConfig

def main(argv=None):
    p = argparse.ArgumentParser(prog="smsd-pro", description="SMSD Pro CLI")
    sub = p.add_subparsers(dest="cmd", required=True)

    ps = sub.add_parser("substructure", help="Run substructure search")
    ps.add_argument("query")
    ps.add_argument("target")
    ps.add_argument("--smiles-target", action="store_true", help="Treat target as SMILES even if SMARTS-y")
    ps.add_argument("--timeout", type=float, default=5.0)
    ps.add_argument("--png", type=str, default=None, help="Save side-by-side PNG")
    ps.add_argument("--flex", action="store_true", help="Flexible aromatic/bond order")
    ps.add_argument("--connected-only", action="store_true", default=True)

    pm = sub.add_parser("mcs", help="Run MCS (MCIS/MCCS)")
    pm.add_argument("mol1")
    pm.add_argument("mol2")
    pm.add_argument("--type", choices=["MCIS","MCCS"], default="MCIS")
    pm.add_argument("--timeout", type=float, default=10.0)
    pm.add_argument("--extend", type=float, default=3.0, help="McGregor extend seconds (0 to disable)")
    pm.add_argument("--png", type=str, default=None, help="Save side-by-side PNG")

    args = p.parse_args(argv)

    chem = ChemOptions()
    if getattr(args, "flex", False):
        chem = ChemOptions(aromaticity_mode="flexible", match_bond_order="loose", ring_matches_ring_only=False)

    if args.cmd == "substructure":
        smsd = SMSD(args.query, args.target, chem=chem)
        t0 = time.time()
        ok = smsd.substructure_exists(SubstructureOptions(time_limit_s=args.timeout,
                                                          connected_only=args.connected_only))
        dt = time.time() - t0
        print(json.dumps({"exists": bool(ok), "time_s": round(dt, 4)}))
        if args.png:
            viz_compare_substructure_png(args.query, args.target, smsd_chem=chem,
                                         cfg=VizConfig(save_png=args.png))
            print(f"PNG saved to {args.png}")
    else:
        smsd = SMSD(args.mol1, args.mol2, chem=chem)
        t0 = time.time()
        mr = smsd.mcs_max(MCSOptions(mcs_type=args.type, time_limit_s=args.timeout,
                                     use_mcgregor_extend=(args.extend>0), extend_time_limit_s=max(0.0,args.extend)))
        dt = time.time() - t0
        if mr:
            out = {"size": mr.size, "algo": mr.algorithm, "tanimoto_atoms": round(mr.tanimoto_atoms or 0,3),
                   "tanimoto_bonds": round(mr.tanimoto_bonds or 0,3), "time_s": round(dt, 4)}
        else:
            out = {"size": 0, "algo": args.type, "time_s": round(dt,4)}
        print(json.dumps(out))
        if args.png:
            viz_compare_mcs_png(args.mol1, args.mol2, smsd_chem=chem,
                                cfg=VizConfig(save_png=args.png))
            print(f"PNG saved to {args.png}")

if __name__ == "__main__":
    main()
