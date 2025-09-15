# SMSD Pro — BioInception PVT LTD

**SMSD Pro** provides fast, RDKit-independent engines for:

- **Substructure search** (VF2++-style with terminal look-ahead, degree/ring pruning, options for chirality, bond stereo, ring and aromaticity controls).
- **Maximum Common Substructure (MCS)** via an **exact** induced approach:
  - Modular product + **bit-parallel clique (BBMC)** (Tomita-like pivot and sequential colouring bound).
  - Optional **McGregor seed→extend** to further enlarge mappings under chemistry constraints.
- SMARTS support (including E/Z handling through directional single propagation).

### Install

```bash
# Recommended: create an environment and install RDKit from conda-forge
conda create -n smsd-pro -c conda-forge python=3.11 rdkit
conda activate smsd-pro

# Then install SMSD Pro
pip install .
```

If you prefer pip-only and your platform has RDKit wheels:
```bash
pip install "rdkit-pypi>=2023.9.1"
pip install .
```

### Quick start

```python
from smsd_pro.engines import SMSD, ChemOptions, SubstructureOptions, MCSOptions

smsd = SMSD("c1ccccc1", "c1ccc2ccccc2c1")
print("Substructure?", smsd.substructure_exists())

m = smsd.mcs_max(MCSOptions(mcs_type="MCCS"))
print("MCS size:", m.size if m else 0)
```

Command line:

```bash
smsd-pro --query "c1ccccc1" --target "c1ccc2ccccc2c1" --mode substructure
smsd-pro --mol1 "CCCCCC" --mol2 "CCCCCCCC" --mode mcs --mcs-type MCIS
```

See **`docs/WHITEPAPER.md`** for algorithm details and citations, and **`tests/`** for 100+ checks.
