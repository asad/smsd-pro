.PHONY: install test bench viz

install:
\tpip install -e .

test:
\tpytest -q

bench:
\tpython -m benchmarks.bench_suite

viz:
\tpython - <<'PY'
from smsd_pro.engines import SMSD
from smsd_pro.mcs import MCSOptions
smsd = SMSD("c1ccccc1","c1ccc2ccccc2c1")
mr = smsd.mcs_max(MCSOptions())
img = smsd.grid_image_substructure(mr)
img.save("viz_demo.png")
print("Saved viz_demo.png")
PY
