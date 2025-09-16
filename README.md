# SMSD Pro 1.0.0

Fast, chemistry‑aware **substructure** and **MCS** (Maximum Common Subgraph) engines with visualisation. Built for med‑chem search workloads, powered by RDKit for IO and chemistry and our own exact VF2++ subgraph isomorphism and bit‑parallel maximum‑clique core.

**Heritage.** SMSD Pro continues the ideas introduced in the original **Small Molecule Subgraph Detector (SMSD)** toolkit by S. A. Rahman *et al.* (J. Cheminformatics, 2009) — see the citation below — modernised for today’s Python/RDKit stacks and extended with modular‑product MCS, McGregor‑style extension, recursive SMARTS, and benchmarkable visuals.

- **Repository:** <https://github.com/asad/smsd-pro>  
- **Author:** Syed Asad Rahman (BioInception PVT LTD)

---

## Install

Python **3.9–3.13** supported.

### Quick start (local editable install)

```bash
# Optional: create env
python -m venv .venv && source .venv/bin/activate

# RDKit (PyPI build) + project
pip install "rdkit-pypi>=2025.03.1"
pip install -e .
```

> If you prefer extras: `pip install -e .[rdkit]` also pulls `rdkit-pypi` for you.

**Why you saw the previous packaging error.** We fixed `pyproject.toml` to use a proper **SPDX license expression** (`Apache-2.0`) and **removed license classifiers** (deprecated by setuptools/PEP 639). If you still have an older file in your working tree, replace it with the one in this release.

---

## Test

```bash
pytest -q
```

What to expect:
- **Unit tests:** fast checks over 100+ SMILES/SMARTS pairs for both substructure and MCS.  
- **Visual tests:** PNGs written to **`test_output/`** (created automatically).

Outputs:
- `test_output/benzene_in_naphthalene.png`
- `test_output/benzene_mcs.png`

---

## Bench

```bash
python scripts/bench.py --limit-png 20
```

This writes everything under **`test_output/`**:
- `bench_cases.csv` – per‑pair metrics (existence time, MCS size, algorithm flags)
- `bench_cases.json` – same in JSON
- per‑case PNGs: `*_sub.png`, `*_mcs.png`
- `cases_index.md` – quick index with image links

The curated set includes **SMARTS with recursive `$()`**, classic aromatic/ring cases, stereo examples, a few RDKit “slow pairs”, and chain growth series.

---

## Usage (API)

```python
from smsd_pro import SMSD, SubstructureOptions, MCSOptions, ChemOptions

smiles_q = "c1ccccc1"
smiles_t = "c1ccc2ccccc2c1"

smsd = SMSD(smiles_q, smiles_t, chem=ChemOptions())
exists = smsd.substructure_exists(SubstructureOptions())

mcs = smsd.mcs_max(MCSOptions(mcs_type="MCIS"))  # or "MCCS"
print(mcs.size, mcs.algorithm, mcs.tanimoto_atoms, mcs.tanimoto_bonds)
```

See `smsd_pro/viz.py` for side‑by‑side plotting helpers.

---

## How it works (very short)

- **Substructure:** exact VF2++‑style backtracking with chemistry‑aware pruning, ring size options, bond stereo modes, and optional recursive SMARTS `$()` anchored at the current atom.
- **MCS:** modular product + bit‑parallel maximum clique (BBMC colouring bound). Optional **McGregor extension** grows the seed greedily when the query is a concrete molecule (skipped for SMARTS to avoid RDKit preconditions).
- **Standardisation:** conservative RDKit `MolStandardize` pipeline (largest fragment, normalise, reionise, uncharge, canonical tautomer).

For a full description, see the whitepaper in `docs/WHITEPAPER.md`.

---

## Citation

If you use SMSD Pro in scientific work, please cite:

- **SMSD (2009):** S. A. Rahman, M. Bashton, G. L. Holliday, R. Schrader and J. M. Thornton. *Small Molecule Subgraph Detector (SMSD) toolkit.* **Journal of Cheminformatics** 2009, 1:12. doi:10.1186/1758-2946-1-12.

A software citation is provided in **CITATION.cff**.

---

## License & Notice

Licensed under the **Apache License 2.0**.  
© 2025 BioInception PVT LTD.

This software is provided **“as is”**, without warranties or conditions of any kind. You are responsible for compliance with all laws, regulations, and third‑party rights in your jurisdiction. By using the software you agree that neither the authors nor BioInception shall be liable for any claim, damages, or other liability arising from its use.

We ask — as a goodwill gesture — that derivative works and publications **acknowledge SMSD Pro** and include a reference to the 2009 SMSD paper above.

---

## Acknowledgements

Built with love on RDKit. Our thanks to the community for decades of algorithms on subgraph isomorphism, maximum clique, and cheminformatics.

