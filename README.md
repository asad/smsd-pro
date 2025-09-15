# smsd-pro


conda create -n smsd-pro -c conda-forge python=3.11 rdkit
conda activate smsd-pro

#Install SMSD Pro

pip install smsd-pro         # or: cd smsd-pro && pip install .
# for dev mode:
pip install -e smsd-pro


#Run tests (100+)
pytest -q

#Benchmarks & graphs
python -m benchmarks.bench_suite
# saves a plot at: benchmarks/out/benchmark_times.png


#CLI examples

# Substructure
smsd-pro --mode substructure --query "c1ccccc1" --target "c1ccc2ccccc2c1"

# MCS (MCIS)
smsd-pro --mode mcs --mol1 "CCCCCC" --mol2 "CCCCCCCC" --mcs-type MCIS --timeout 5 --extend



# Python code
from smsd_pro.engines import SMSD, MCSOptions

smsd = SMSD("c1ccccc1","c1ccc2ccccc2c1")
mr = smsd.mcs_max(MCSOptions())
img = smsd.grid_image_substructure(mr)
img.save("viz_demo.png")



# Substructure
smsd-pro --mode substructure --query "c1ccccc1" --target "c1ccc2ccccc2c1"

# MCS (MCIS)
smsd-pro --mode mcs --mol1 "CCCCCC" --mol2 "CCCCCCCC" --mcs-type MCIS --timeout 5 --extend

