# Bench: What It Runs and What You Get

Run:
```bash
python scripts/bench.py --limit-png 20
```

Outputs (all under `test_output/`):

- `bench_cases.csv` and `bench_cases.json` — metrics per pair (substructure existence time, MCS size/algorithm)  
- PNGs for each case: `*_sub.png` and `*_mcs.png`  
- `cases_index.md` — a simple gallery of all generated images

Included examples cover aromatics, E/Z stereo, chain growth, recursive SMARTS like `[C;$(C(=O)O)]` vs `CC(=O)O`, and a few “slow” community cases.

**Background:** This bench continues the lineage of the **Small Molecule Subgraph Detector (SMSD)** by Rahman *et al.* (J. Cheminf. 2009, 1:12; doi:10.1186/1758-2946-1-12). SMSD Pro modernises the algorithms and adds a modular‑product/BBMC MCS core with optional McGregor extension.
