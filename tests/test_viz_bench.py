import os
from smsd_pro.viz import viz_compare_substructure_png, viz_compare_mcs_png, VizConfig

OUT = os.path.join(os.path.dirname(os.path.dirname(__file__)), "test_output")

def test_visual_benzene_naphthalene_png():
    img = viz_compare_substructure_png("c1ccccc1", "c1ccc2ccccc2c1",
                                       cfg=VizConfig(save_png=os.path.join(OUT, "benzene_in_naphthalene.png")))
    assert img is not None

def test_visual_mcs_png():
    img = viz_compare_mcs_png("c1ccccc1", "c1ccc2ccccc2c1",
                              cfg=VizConfig(save_png=os.path.join(OUT, "benzene_mcs.png")))
    assert img is not None
