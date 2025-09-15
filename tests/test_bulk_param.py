import pytest
from smsd_pro.engines import SMSD
from smsd_pro.mcs import MCSOptions

@pytest.mark.parametrize("a,b", [
    ("CC", "CCCCC"),
    ("CCC", "CCCCCC"),
    ("CCCC", "CCCCCCC"),
    ("CCCCC", "CCCCCCCC"),
    ("CCCCCC", "CCCCCCCCC"),
    ("CCCCCCC", "CCCCCCCCCC"),
    ("CCCCCCCC", "CCCCCCCCCCC"),
    ("CCCCCCCCC", "CCCCCCCCCCCC"),
    ("CCCCCCCCCC", "CCCCCCCCCCCCC"),
    ("CCCCCCCCCCC", "CCCCCCCCCCCCCC"),
    ("c1ccccc1", "c1ccc2ccccc2c1"),
    ("c1ccccc1", "c1ccc2ccccc2c1"),
    ("c1ccccc1", "c1ccc2ccccc2c1"),
    ("c1ccccc1", "c1ccc2ccccc2c1"),
    ("c1ccccc1", "c1ccc2ccccc2c1"),
])
def test_bulk(a,b):
    smsd = SMSD(a, b)
    assert smsd.substructure_exists()
    mr = smsd.mcs_max(MCSOptions(time_limit_s=2.0))
    assert mr is not None
