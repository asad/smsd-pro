# SPDX-License-Identifier: Apache-2.0
# © 2025 BioInception PVT LTD.

"""
SMSD Pro CLI
------------
Find substructure(s) or MCS between two molecules, report mapping(s), Jaccard/Tanimoto
(similarity on atoms/bonds), timings, and save a PNG with highlights.

The PNG and any tabular/JSON outputs are written to ./test_output by default.
"""

from __future__ import annotations
import argparse, os, json, time, sys
from dataclasses import asdict
from typing import Dict, Tuple, Optional, List

from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.error')  # keep CLI clean

from smsd_pro.engines import SMSD, SubstructureOptions, MCSOptions
from smsd_pro.chem import ChemOptions
from smsd_pro.viz import viz_compare_substructure_png, viz_compare_mcs_png, VizConfig

# ----------------------
# Helpers / calculations
# ----------------------
def bonds_common(q: Chem.Mol, t: Chem.Mol, m: Dict[int,int]) -> int:
    cnt = 0
    for b in q.GetBonds():
        i, j = int(b.GetBeginAtomIdx()), int(b.GetEndAtomIdx())
        if i in m and j in m and t.GetBondBetweenAtoms(int(m[i]), int(m[j])) is not None:
            cnt += 1
    return cnt

def tanimoto(n_common: int, n1: int, n2: int) -> float:
    denom = max(1, n1 + n2 - n_common)
    return float(n_common) / float(denom)

# ----------------------
# Profiles to mimic RDKit expectations
# ----------------------
def profile(name: str) -> Tuple[str, ChemOptions]:
    n = name.lower()
    if n in ("default", "strict"):
        return "default", ChemOptions()
    if n in ("flex", "flexible"):
        return "flexible", ChemOptions(
            aromaticity_mode="flexible",
            match_bond_order="loose",
            ring_size_mode="subset",
            ring_size_tolerance=1,
            ring_matches_ring_only=False,
            bond_stereo="defined",
            use_chirality=False,
        )
    if n in ("rdkit-substruct", "rdkit"):
        return "rdkit-substruct", ChemOptions(
            aromaticity_mode="strict",
            match_bond_order="strict",
            ring_matches_ring_only=True,
            bond_stereo="defined",
            use_chirality=False,
        )
    if n in ("rdkit-fmcs", "fmcs"):
        return "rdkit-fmcs", ChemOptions(
            aromaticity_mode="flexible",
            match_bond_order="loose",
            ring_size_mode="subset",
            ring_size_tolerance=1,
            ring_matches_ring_only=True,
            bond_stereo="off",
            use_chirality=False,
        )
    raise SystemExit(f"Unknown chem profile: {name}")

# ----------------------
# Main execution
# ----------------------
def run(args) -> int:
    out_dir = args.output_dir or os.path.join(os.getcwd(), "test_output")
    os.makedirs(out_dir, exist_ok=True)

    prof_name, chem = profile(args.profile)

    subopt = SubstructureOptions(
        connected_only=not args.allow_disconnected,
        induced=args.induced,
        wl_rounds=0,
        time_limit_s=args.timeout,
        max_matches=args.max_matches,
        uniquify_mode=args.uniquify,
    )

    mcsopt = MCSOptions(
        mcs_type=args.mcs_type.upper(),
        time_limit_s=args.timeout,
        use_mcgregor_extend=not args.no_mcgregor,
        extend_time_limit_s=max(0.0, args.extend_time),
    )

    smsd = SMSD(args.query, args.target, chem=chem)

    result = {
        "query": args.query,
        "target": args.target,
        "chem_profile": prof_name,
        "mode": args.mode,
        "substructure": None,
        "mcs": None,
    }

    if args.mode in ("substructure", "both"):
        t0 = time.perf_counter()
        maps = smsd.substructure_all(subopt)
        ss_ms = (time.perf_counter() - t0) * 1000.0
        if maps:
            m = maps[0].mapping
            q_atoms = len(m)
            q_bonds = bonds_common(smsd.q, smsd.t, m)
            tana = tanimoto(q_atoms, smsd.q.GetNumAtoms(), smsd.t.GetNumAtoms())
            tanb = tanimoto(q_bonds, smsd.q.GetNumBonds(), smsd.t.GetNumBonds())
            result["substructure"] = {
                "exists": True,
                "mappings_found": len(maps),
                "mapping": m,
                "atoms_mapped": q_atoms,
                "bonds_mapped": q_bonds,
                "tanimoto_atoms": round(tana, 6),
                "tanimoto_bonds": round(tanb, 6),
                "ms": round(ss_ms, 3),
                "algorithm": "SUBSTRUCTURE_VF2PP"
            }
        else:
            result["substructure"] = {"exists": False, "mappings_found": 0, "ms": round(ss_ms, 3)}

        if args.png:
            png_path = os.path.join(out_dir, args.png if args.png.lower().endswith(".png") else args.png + ".png")
            viz_compare_substructure_png(args.query, args.target,
                                         smsd_chem=chem, smsd_opt=subopt,
                                         cfg=VizConfig(save_png=png_path))
            result.setdefault("artifacts", {})["substructure_png"] = png_path

    if args.mode in ("mcs", "both"):
        t0 = time.perf_counter()
        mr = smsd.mcs_max(mcsopt)
        mcs_ms = (time.perf_counter() - t0) * 1000.0
        if mr:
            m = mr.mapping
            a = len(m)
            b = bonds_common(smsd.q, smsd.t, m)
            tana = tanimoto(a, smsd.q.GetNumAtoms(), smsd.t.GetNumAtoms())
            tanb = tanimoto(b, smsd.q.GetNumBonds(), smsd.t.GetNumBonds())
            result["mcs"] = {
                "mapping": m,
                "atoms": a,
                "bonds": b,
                "tanimoto_atoms": round(tana, 6),
                "tanimoto_bonds": round(tanb, 6),
                "algorithm": mr.algorithm,
                "ms": round(mcs_ms, 3)
            }
        else:
            result["mcs"] = {"mapping": {}, "atoms": 0, "bonds": 0, "algorithm": "NONE", "ms": round(mcs_ms, 3)}

        if args.png:
            png_path = os.path.join(out_dir, args.png if args.png.lower().endswith(".png") else args.png + ".png")
            viz_compare_mcs_png(args.query, args.target,
                                smsd_chem=chem, mcs_opt=mcsopt,
                                cfg=VizConfig(save_png=png_path))
            result.setdefault("artifacts", {})["mcs_png"] = png_path

    # Save JSON/CSV if requested
    if args.json:
        json_path = os.path.join(out_dir, args.json if args.json.lower().endswith(".json") else args.json + ".json")
        with open(json_path, "w") as f:
            json.dump(result, f, indent=2)
        result.setdefault("artifacts", {})["json"] = json_path

    if args.csv:
        import csv
        csv_path = os.path.join(out_dir, args.csv if args.csv.lower().endswith(".csv") else args.csv + ".csv")
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["mode","chem_profile","atoms_mapped","bonds_mapped","tan_atoms","tan_bonds","ms","algorithm"])
            if args.mode in ("substructure", "both") and result["substructure"]:
                ss = result["substructure"]
                w.writerow(["substructure", prof_name, ss.get("atoms_mapped",0), ss.get("bonds_mapped",0),
                            ss.get("tanimoto_atoms",0.0), ss.get("tanimoto_bonds",0.0), ss.get("ms",0.0),
                            ss.get("algorithm","")])
            if args.mode in ("mcs", "both") and result["mcs"]:
                mm = result["mcs"]
                w.writerow(["mcs", prof_name, mm.get("atoms",0), mm.get("bonds",0),
                            mm.get("tanimoto_atoms",0.0), mm.get("tanimoto_bonds",0.0), mm.get("ms",0.0),
                            mm.get("algorithm","")])
        result.setdefault("artifacts", {})["csv"] = csv_path

    # Console summary
    def _fmt(d: Optional[dict]) -> str:
        if not d: return "—"
        if "exists" in d:
            if not d["exists"]:
                return f"no match (ms={d['ms']})"
            return f"atoms={d.get('atoms_mapped',0)} bonds={d.get('bonds_mapped',0)} tanA={d.get('tanimoto_atoms',0.0)} tanB={d.get('tanimoto_bonds',0.0)} ms={d['ms']}"
        return f"atoms={d.get('atoms',0)} bonds={d.get('bonds',0)} tanA={d.get('tanimoto_atoms',0.0)} tanB={d.get('tanimoto_bonds',0.0)} ms={d['ms']} algo={d.get('algorithm','')}"

    print(f"[SMSD] mode={args.mode} profile={prof_name}")
    print(f"  Query : {args.query}")
    print(f"  Target: {args.target}")
    if args.mode in ('substructure', 'both'):
        print("  Substructure:", _fmt(result['substructure']))
    if args.mode in ('mcs', 'both'):
        print("  MCS        :", _fmt(result['mcs']))
    if "artifacts" in result:
        for k,v in result["artifacts"].items():
            print(f"  {k}: {v}")

    return 0

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="smsd.py",
        description="SMSD Pro CLI: substructure and MCS with rich chemistry options; PNG overlays and metrics."
    )
    p.add_argument("--mode", choices=["substructure","mcs","both"], default="both",
                   help="Which operation to run.")
    p.add_argument("--q","--query", dest="query", required=True, help="Query SMILES/SMARTS (or single-mol file path).")
    p.add_argument("--t","--target", dest="target", required=True, help="Target SMILES (or single-mol file path).")

    # Chemistry profile
    p.add_argument("--profile", default="default",
                   choices=["default","flexible","flex","rdkit-substruct","rdkit","rdkit-fmcs","fmcs"],
                   help="Convenience profiles to mimic RDKit or use flexible matching.")

    # Substructure options
    p.add_argument("--allow-disconnected", action="store_true",
                   help="Permit disconnected mappings for substructure (default: connected only).")
    p.add_argument("--induced", action="store_true", help="Require induced subgraph for substructure.")
    p.add_argument("--max-matches", type=int, default=None, help="Cap #substructure mappings to enumerate.")
    p.add_argument("--uniquify", choices=["mapping","target_set"], default="target_set",
                   help="Deduplicate mappings by exact mapping or by target-atom set.")

    # MCS options
    p.add_argument("--mcs-type", choices=["MCIS","MCCS"], default="MCCS",
                   help="MCIS=induced; MCCS=connected common subgraph in query (default).")
    p.add_argument("--no-mcgregor", action="store_true",
                   help="Disable McGregor extension step (seed grow).")
    p.add_argument("--extend-time", type=float, default=3.0, help="McGregor extension time limit (s).")

    # Common
    p.add_argument("--timeout", type=float, default=10.0, help="Time limit (s) for search.")

    # Output
    p.add_argument("--png", default=None, help="Filename for PNG (saved under test_output).")
    p.add_argument("--json", default=None, help="Filename for JSON summary (saved under test_output).")
    p.add_argument("--csv", default=None, help="Filename for CSV summary (saved under test_output).")
    p.add_argument("--output-dir", default=None, help="Where to write outputs (default: ./test_output).")
    return p

def main(argv=None):
    argv = argv or sys.argv[1:]
    args = build_parser().parse_args(argv)
    return run(args)

if __name__ == "__main__":
    raise SystemExit(main())
