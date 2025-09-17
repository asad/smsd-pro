# SMSD Pro

**Small Molecule Subgraph Detector – Pro Edition (v1.0.0)**  
Fast, configurable substructure and MCS engines with chemistry‑aware constraints and clear RDKit parity modes.

> Heritage: This work builds on the ideas in the original **SMSD** paper:  
> S. A. Rahman, M. Bashton, G. L. Holliday, R. Schrader, and J. M. Thornton.  
> *Small Molecule Subgraph Detector (SMSD) toolkit.* **Journal of Cheminformatics** (2009) 1:12. DOI:10.1186/1758-2946-1-12.

**Repository**: https://github.com/asad/smsd-pro  
**Author**: Syed Asad Rahman (BioInception PVT LTD)

---

## Why SMSD Pro

- **Chemistry‑aware**: atom/bond comparators for aromaticity, charge, valence, ring size/tightness, stereo.
- **Modern cores**: VF2++ for substructure; modular product + bit‑parallel BBMC for maximum clique; optional McGregor grow step.
- **Profiles**: out‑of‑the‑box configurations to mimic RDKit behaviour (`rdkit-strict`, `rdkit-flexible`) or to explore very loose matching.
- **SMARTS‑safe**: recursive SMARTS `$()` are supported; pattern queries avoid unsafe valence calls during extension.
- **CLI & PNGs**: command‑line tool `smsd-pro` prints sizes, Jaccard/Tanimoto, mappings, and saves comparison PNGs.

---

## Install

```bash
# (recommended) new env with RDKit
conda create -n smsdpro -c conda-forge python=3.11 rdkit -y
conda activate smsdpro

# from source
pip install -e .
```

> Packaging follows PEP 621; license uses SPDX string. No legacy license classifiers are used.

---

## Quick start (CLI)

```bash
# Substructure (first hit by default)
smsd-pro ss --query "c1ccccc1" --target "c1ccc2ccccc2c1" --profile rdkit-strict --png test_output/ss.png

# Substructure – enumerate all unique target placements
smsd-pro ss --query "C.C" --target "CC" --all-matches --uniq target_set

# MCS (connected by default; like FMCS)
smsd-pro mcs --mol1 "NC(=O)c1[nH]c2ccccc2c1S(=O)(=O)N1CCOC(C(=O)N2CCc3c(Br)cccc3C2)C1"              --mol2 "NC(=O)c1[nH]c2ccccc2c1S(=O)(=O)N1CCOC(C(=O)NCCOc2ccccc2Br)C1"              --profile rdkit-strict --png test_output/mcs.png
```

### Selected flags

- `--profile` one of: `default`, `rdkit-strict`, `rdkit-flexible`, `maximal-loose`.
- `--connected-only/--no-connected-only` (substructure).
- `--mcs-type` one of: `MCCS` (connected; default) or `MCIS` (induced).
- `--all-matches` (substructure) or `--max-matches N`.
- `--timeout-s` time limit for search; `--extend-timeout-s` for McGregor.
- `--png` save a side‑by‑side PNG with highlights into `test_output/`.

Run `smsd-pro --help` for the complete reference.

---

## Python API

```python
from smsd_pro import SMSD, ChemOptions, SubstructureOptions, MCSOptions, profiles

smsd = SMSD("c1ccccc1", "c1ccc2ccccc2c1", chem=profiles.rdkit_strict())
hit = smsd.substructure_exists(SubstructureOptions())  # True

mr = smsd.mcs_max(MCSOptions(mcs_type="MCCS"))
print(mr.size, mr.tanimoto_atoms, mr.tanimoto_bonds)
```

---

## Tests & Bench

```bash
pytest -q
python scripts/bench.py
```

Outputs are written to **`test_output/`**:
- `bench_results.csv` & `bench_results.json` – per‑case metrics (speed, atoms, bonds, Tanimoto).
- PNGs – selected SMSD vs RDKit/FMSC visuals with matched atoms/bonds highlighted.

---

## Citation & Attribution

If you use SMSD Pro in academic or industrial work, please cite both:
1. **SMSD Pro (this software)** – see `CITATION.cff` at the repository root.
2. The original **SMSD** paper: Rahman *et al.* J. Cheminf. (2009) 1:12. DOI:10.1186/1758-2946-1-12.

---

## License & Contribution

Licensed under **Apache-2.0**. See `LICENSE` for terms, warranty disclaimer, and limitation of liability.  
Contributions are welcome – please read `CONTRIBUTING.md` and the `CODE_OF_CONDUCT.md`.

> Public release supported by **BioInception PVT LTD** to enable community reuse in medicinal chemistry and graph algorithms research.
