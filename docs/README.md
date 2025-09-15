# SMSD Pro

**SMSD Pro** is a compact, fast set of engines for **substructure search** and **maximum common substructure (MCS)** with rich, medicinal‑chemistry aware matching.
It is independent of RDKit’s matching/FMCS implementations (we only use RDKit for molecule IO and chemistry metadata).

> Highlights
>
> - Substructure: **VF2++** style backtracker with terminal look‑ahead, degree/ring pruning, stereo checks.
> - MCS (MCIS/MCCS): **BBMC** (Tomita‑like) maximum clique on the modular product + **McGregor** seed→extend.
> - **Recursive SMARTS `$()`** via a tiny cached predicate VM, *without touching the engines*.
> - Strict/flexible aromaticity, ring‑only matching, ring size comparators, bond order modes, chirality and E/Z.
> - Visualiser to export side‑by‑side PNGs comparing **SMSD vs RDKit** mappings.
> - Benchmarks + 100 auto‑generated tests covering simple to very difficult cases.

## Quick start

```bash
pip install -e .[dev]
pytest -q
python -m smsd_pro.cli --help
```

See `docs/WHITEPAPER.md` for an algorithmic overview.
