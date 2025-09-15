import pytest
from smsd_pro.engines import SMSD
from smsd_pro.mcs import MCSOptions
from smsd_pro.chem import ChemOptions

@pytest.mark.parametrize("a,b,expected", [
    ("CC", "CCCC", 2),
    ("CCC", "CCCCC", 3),
    ("c1ccccc1", "c1ccc2ccccc2c1", 6),
])
def test_mcis_sizes(a,b,expected):
    smsd = SMSD(a, b)
    mr = smsd.mcs_max(MCSOptions(mcs_type="MCIS", time_limit_s=2.0))
    assert mr and mr.size >= expected

def test_mccs_connected():
    smsd = SMSD("CC.CC", "CCCC")
    mr = smsd.mcs_max(MCSOptions(mcs_type="MCCS", time_limit_s=2.0))
    assert mr and mr.size >= 2
