# Benchmark Guide

Run:
```bash
python scripts/bench.py
```

Outputs to `test_output/`:
- `bench_results.csv` and `bench_results.json` with per‑case metrics:
  - substructure: existence, time (ms)
  - MCS: atoms, bonds, Tanimoto (atoms/bonds), algorithm, time (ms)
  - RDKit substructure and FMCS reference sizes/times for context
- PNGs: selected side‑by‑side comparisons

Includes cases from chains, rings, aromatic vs Kekulé, stereo, functional groups, recursive SMARTS, and one drug‑like pair.

Please cite the original **SMSD** paper: Rahman *et al.* J. Cheminformatics (2009) 1:12. DOI:10.1186/1758-2946-1-12.
