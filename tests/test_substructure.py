import pytest, time
from rdkit import Chem
from smsd_pro.engines import SMSD
from smsd_pro.vf import SubstructureOptions
from smsd_pro.chem import ChemOptions

def run(q, t, chem=ChemOptions(), opt=SubstructureOptions()):
    smsd = SMSD(q, t, chem=chem, standardise=True)
    return smsd.substructure_exists(opt), smsd.substructure_all(SubstructureOptions(**{**opt.__dict__, "max_matches": 64}))

def test_benzene_in_naph():
    ok, maps = run("c1ccccc1", "c1ccc2ccccc2c1")
    assert ok
    assert len(maps) >= 2

def test_trans_cis_smarts():
    chem = ChemOptions(bond_stereo="exact")
    ok1, _ = run("C/C=C\\C", "C/C=C\\C", chem=chem)
    ok2, _ = run("C/C=C\\C", "C/C=C/C", chem=chem)
    assert ok1 and not ok2

@pytest.mark.parametrize("k", list(range(2, 22)))
def test_chain_in_chain(k):
    q = "C"*k; t = "C"*(k+2)
    ok, maps = run(q, t)
    assert ok
    assert len(maps) >= 1

def test_disconnected_query():
    ok, _ = run("C.C", "CC")
    assert ok

def test_ring_strict():
    ok, _ = run("C1CCCC1", "C1CCCCC1")
    assert not ok
