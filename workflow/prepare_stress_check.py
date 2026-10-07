"""Prepare +/- isotropic strain checks on CO2; retain every real CP2K output."""
import argparse,json
from pathlib import Path
from ase.io import read,write
from local_compute import cp2k_input
p=argparse.ArgumentParser();p.add_argument('runs',type=Path);a=p.parse_args()
source=read(a.runs/'smoke-001/co2.relaxed.extxyz')
for name,strain in [('minus',-.001),('zero',0),('plus',.001)]:
    out=a.runs/'stress-check'/name;out.mkdir(parents=True,exist_ok=True)
    atoms=source.copy();atoms.set_cell(atoms.cell*(1+strain),scale_atoms=True)
    write(out/'co2.relaxed.extxyz',atoms);cp2k_input(out,'co2')
    path=out/'co2.inp';s=path.read_text().replace('EPS_SCF 1.0E-6','EPS_SCF 1.0E-8');path.write_text(s)
(a.runs/'stress-check/settings.json').write_text(json.dumps({'isotropic_strain':.001,'volume_A3':source.get_volume(),'purpose':'finite-difference stress sign/units integration test, not material convergence'}))
