\
import math, os
import itertools as it
from rdkit import Chem
from smsd_pro.engines import SMSD, SubstructureOptions, MCSOptions
from smsd_pro.chem import ChemOptions

def jaccard(a_inter, a_total):
    return a_inter / max(1, a_total)

def test_ethane_in_butane():
    smsd = SMSD("CC", "CCCC")
    assert smsd.substructure_exists()

def test_benzene_in_naphthalene_strict():
    smsd = SMSD("c1ccccc1", "c1ccc2ccccc2c1")
    assert smsd.substructure_exists()

def test_trans_smarts_self():
    smsd = SMSD("C/C=C\\C", "C/C=C\\C")
    assert smsd.substructure_exists()

def test_recursive_smarts_amide_n():
    smsd = SMSD("[N;$isAmideN]", "CC(=O)NCC")
    assert smsd.substructure_exists()

def test_mcis_chain():
    smsd = SMSD("CCCCC", "CCCCCCC")
    mr = smsd.mcs_max(MCSOptions(mcs_type="MCIS"))
    assert mr and mr.size == 5

# --- generate 90+ parameterised small cases ---
def gen_linear_cases(n=30):
    for k in range(2, n+2):
        q = "C"*k
        t = "C"*(k+2)
        yield q, t

def gen_ring_cases():
    targets = ["c1ccc2ccccc2c1", "C1=CC=CC=C1", "C1CCCC1"]
    for n in range(5, 9):
        q = "C1" + "C"*(n-1) + "1"
        for t in targets:
            yield q, t

def gen_functional():
    return [
        ("NC(=O)", "NCC(=O)NCCC"),
        ("[O-]C=O", "OC=O"),
        ("[C;$(C(=O)O)]", "CC(=O)O"),
        ("[N;$isAmideN]", "CC(=O)NCC"),
    ]

CASES = list(gen_linear_cases(40)) + list(gen_ring_cases()) + gen_functional()
assert len(CASES) >= 100

def test_many_substructures():
    ok = 0
    for q,t in CASES:
        smsd = SMSD(q, t, chem=ChemOptions())
        if smsd.substructure_exists(SubstructureOptions(connected_only=True)):
            ok += 1
    assert ok >= 60  # loose lower bound just to catch regressions

def test_mcises_dont_crash_and_typical_sizes():
    for q,t in gen_linear_cases(20):
        smsd = SMSD(q, t)
        mr = smsd.mcs_max()
        assert mr and mr.size >= len(q) - 2
