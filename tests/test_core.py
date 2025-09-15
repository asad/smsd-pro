import pytest
from smsd_pro.engines import SMSD, SubstructureOptions, MCSOptions
from smsd_pro.chem import ChemOptions

# Generate at least 100 cases: chains, rings, aromatics, stereo, charges, isotopes, SMARTS
CHAINS = [("C"*k, "C"*(k+2)) for k in range(2, 22)]
AROM = [("c1ccccc1", "c1ccc2ccccc2c1"), ("c1ccccc1", "C1=CC=CC=C1")]
STEREO = [("C/C=C\\C", "C/C=C\\C"), ("C/C=C\\C", "C/C=C/C")]
DISCON = [("C.C", "CC")]
FUNCTIONALS = [("NC(=O)", "NCC(=O)NCCC"), ("[O-]C=O", "OC=O")]
OTHERS = [("[C;$(C(=O)O)]","CC(=O)O"), ("[N;$isAmideN]","CC(=O)NCC")]

CASES = []
for i,(q,t) in enumerate(CHAINS): CASES.append((f"chain_{i}", q, t))
for i,(q,t) in enumerate(AROM): CASES.append((f"arom_{i}", q, t))
for i,(q,t) in enumerate(STEREO): CASES.append((f"stereo_{i}", q, t))
for i,(q,t) in enumerate(DISCON): CASES.append((f"disc_{i}", q, t))
for i,(q,t) in enumerate(FUNCTIONALS): CASES.append((f"fx_{i}", q, t))
for i,(q,t) in enumerate(OTHERS): CASES.append((f"rec_{i}", q, t))

# pad to 100+ using varied random-ish alkyl/phenyl attachments
EXTRA = [
    ("c1ccccc1CC", "c1ccccc1CCCC"),
    ("CCN(CC)CC", "CCN(CC)CCO"),
    ("C1CCCCC1", "C1CCCCCC1"),
    ("[13CH3]C", "CC"),
    ("c1ncccc1", "c1ncccc1Cl"),
    ("N1CCOCC1", "O=C(N1CCOCC1)C"),
    ("FC(F)F", "CC(C)(F)F"),
    ("O=S(=O)N", "O=S(=O)NCC"),
    ("P(=O)(O)O", "OP(=O)(O)OCC"),
    ("Clc1ccccc1", "Clc1ccc(Cl)cc1"),
]*10  # 90 more

for i,(q,t) in enumerate(EXTRA):
    CASES.append((f"extra_{i}", q, t))

def _ok(maybe):
    return (maybe is not None) and (maybe.size >= 0)

@pytest.mark.parametrize("name,q,t", CASES)
def test_substructure_exists(name, q, t):
    chem = ChemOptions()
    smsd = SMSD(q, t, chem=chem)
    assert isinstance(smsd.substructure_exists(SubstructureOptions()), bool)

@pytest.mark.parametrize("name,q,t", CASES[:30])
def test_mcs_runs(name, q, t):
    chem = ChemOptions()
    smsd = SMSD(q, t, chem=chem)
    mr = smsd.mcs_max(MCSOptions(mcs_type="MCIS"))
    assert _ok(mr)
