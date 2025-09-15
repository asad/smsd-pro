# smsd-pro


conda create -n smsd-pro -c conda-forge python=3.11 rdkit
conda activate smsd-pro

pip install -e .
pytest -q           # 100+ tests, creates PNGs in test_output
python scripts/bench.py


# Substructure (benzene in naphthalene) + PNG
smsd-pro substructure "c1ccccc1" "c1ccc2ccccc2c1" --png benzene_in_naph.png

# MCS (MCIS) + PNG
smsd-pro mcs "NC(=O)CC" "O=C(N)CCC" --type MCIS --png mcs_example.png

# Flexible matching demo
smsd-pro substructure "c1ccccc1" "C1=CC=CC=C1" --flex


ChemOptions(
  match_atom_type=True,
  match_formal_charge=True,
  aromaticity_mode="strict",     # "flexible" allows ring aromatics ↔ Kekulé
  ring_matches_ring_only=True,
  ring_size_mode="subset",       # "exact" | "subset" | "ignore"
  match_bond_order="strict",     # "loose" treats aromatics as compatible with single/double in rings
  bond_stereo="defined",         # "off" | "defined" | "exact"
  use_chirality=False,
  complete_rings_only=False,
)


flex = ChemOptions(aromaticity_mode="flexible", match_bond_order="loose", ring_matches_ring_only=False)



from smsd_pro import SMSD, ChemOptions

# Carboxylic acid carbon via $()
smsd = SMSD("[C;$(C(=O)O)]", "CC(=O)O", chem=ChemOptions())
print("Carboxyl C present?", smsd.substructure_exists())

# Amide nitrogen via named predicate ($isAmideN)
smsd = SMSD("[N;$isAmideN]", "CC(=O)NCC")
print("Amide N present?", smsd.substructure_exists())



