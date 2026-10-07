"""Audit same-geometry numerical repeatability at retained CP2K restart boundaries."""
import argparse,json,re
from pathlib import Path
import numpy as np
from ase.io import read
from aimd_to_deepmd import frames
from refine_pipeline import forces_from_log,HARTREE,digest
def audit(run):
    folders=sorted(run.glob('dft-relax-attempt*'))+[run/'dft-relax'];checks=[]
    for previous,current in zip(folders,folders[1:]):
        saved=read(previous/'opt-1.restart',format='cp2k-restart');start=read(current/'opt.relaxed.extxyz')
        assert saved.get_chemical_symbols()==start.get_chemical_symbols()
        coordinate=float(abs(saved.positions-start.positions).max());cell=float(abs(saved.cell.array-start.cell.array).max())
        xyz=frames(previous/'opt-pos-1.xyz')[-1];old_energy=xyz[1]*HARTREE
        energies=re.findall(r'ENERGY\| Total FORCE_EVAL.*?energy \[hartree\]\s+([-+\d.Ee]+)',(current/'opt.log').read_text())
        new_energy=float(energies[0])*HARTREE
        old_force=forces_from_log(previous/'opt.log')[-1];new_force=forces_from_log(current/'opt.log')[0]
        assert old_force.shape==new_force.shape==(len(saved),3)
        force=float(abs(old_force-new_force).max());energy=abs(old_energy-new_energy)
        checks.append({'from_attempt':previous.name,'to_attempt':current.name,'coordinate_max_difference_A':coordinate,'cell_max_difference_A':cell,
            'energy_difference_eV':energy,'force_component_max_difference_eV_A':force,
            'tolerances':{'coordinate_A':1e-7,'cell_A':1e-10,'energy_eV':1e-6,'force_eV_A':1e-4},
            'passed':coordinate<1e-7 and cell<1e-10 and energy<1e-6 and force<1e-4,
            'sources':{str(p.relative_to(run)):digest(p) for p in [previous/'opt-1.restart',previous/'opt-pos-1.xyz',previous/'opt.log',current/'opt.relaxed.extxyz',current/'opt.log']}})
    return {'provisional':json.loads((run/'dft-relax/job.json').read_text()).get('status')=='running','checks':checks,'restart_repeatability_gate':all(x['passed'] for x in checks) if checks else None,'scope':'Same saved state, independently restarted SCF. This is not proof of production accuracy or bitwise identical later optimizer trajectories.','manuscript_eligible':False}
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();q=audit(a.run)
    (a.run/'cp2k_restart_audit.json').write_text(json.dumps(q,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in x.items() if k not in ['sources','tolerances']} for x in q['checks']],indent=2))
if __name__=='__main__':main()
