import argparse
from pathlib import Path
from ase.io import read,write
p=argparse.ArgumentParser();p.add_argument('runs',type=Path);p.add_argument('--system',choices=['co2','adsorbed'],default='co2');p.add_argument('--dataset',type=Path);a=p.parse_args()
d=a.dataset or a.runs/f'{a.system}-dataset'; atoms=read(a.runs/f'smoke-001/{a.system}.relaxed.extxyz')
types=['C','O'] if a.system=='co2' else ['Ca','P','O','H','C']
write(d/'co2.data',atoms,format='lammps-data',atom_style='atomic',specorder=types,masses=True)
header='''units metal
atom_style atomic
boundary p p p
read_data co2.data
pair_style deepmd co2-smoke.pth
pair_coeff * * C O
neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
timestep 0.00025
velocity all create 100.0 791 mom yes rot yes dist gaussian
fix integ all nve
thermo 1
thermo_style custom step temp pe ke etotal
'''
header=header.replace('co2-smoke.pth',f'{a.system}-smoke.pth').replace('pair_coeff * * C O','pair_coeff * * '+' '.join(types))
(d/'lammps-first.in').write_text(header+'''run 5
write_restart co2.restart
write_dump all custom first-final.dump id type xu yu zu vx vy vz modify sort id format float %.16g
''')
(d/'lammps-continuous.in').write_text(header+'''run 10
write_dump all custom continuous-final.dump id type xu yu zu vx vy vz modify sort id format float %.16g
''')
(d/'lammps-restart.in').write_text(('''units metal
atom_style atomic
read_restart co2.restart
pair_style deepmd co2-smoke.pth
pair_coeff * * C O
neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
fix integ all nve
thermo 1
thermo_style custom step temp pe ke etotal
run 5
write_dump all custom restart-final.dump id type xu yu zu vx vy vz modify sort id format float %.16g
''').replace('co2-smoke.pth',f'{a.system}-smoke.pth').replace('pair_coeff * * C O','pair_coeff * * '+' '.join(types)))
