"""Explicit structure preparation and matched-reference screening tools."""
import argparse, json, shutil
from pathlib import Path
import numpy as np

def main(command):
    from ase.io import read, write
    p=argparse.ArgumentParser()
    if command=='slabs':
        p.add_argument('--bulk',required=True,type=Path); p.add_argument('--out',required=True,type=Path)
        p.add_argument('--hkl',nargs=3,type=int,default=[0,0,1]); p.add_argument('--layers',type=int,default=1)
        p.add_argument('--vacuum',type=float,default=12)
        a=p.parse_args()
        from pymatgen.core import Structure
        from pymatgen.io.ase import AseAtomsAdaptor
        from ase.build import surface
        s=Structure.from_file(a.bulk)
        if not s.is_ordered: raise ValueError('Explicit disorder resolution required before slab generation')
        if a.layers<1 or a.vacuum<=0 or not any(a.hkl): raise ValueError('Invalid slab parameters')
        slab=surface(AseAtomsAdaptor.get_atoms(s),a.hkl,a.layers,vacuum=a.vacuum,periodic=True)
        a.out.parent.mkdir(parents=True,exist_ok=True); write(a.out,slab)
    elif command=='place':
        from ase.build import molecule
        p.add_argument('--slab',required=True,type=Path); p.add_argument('--out',required=True,type=Path)
        p.add_argument('--height',type=float,default=3.0); p.add_argument('--xy',nargs=2,type=float,required=True)
        a=p.parse_args(); slab=read(a.slab)
        if a.height<=0: raise ValueError('Height must be positive')
        gas=molecule('CO2'); gas.rotate('z','x'); gas.translate([*a.xy,slab.positions[:,2].max()+a.height]-gas.positions[0])
        combined=slab+gas; combined.set_constraint(slab.constraints)
        a.out.parent.mkdir(parents=True,exist_ok=True); write(a.out,combined)
        a.out.with_suffix('.groups.json').write_text(json.dumps({'co2_ids':list(range(len(slab)+1,len(combined)+1)),
            'co2_oxygen_ids':[len(slab)+i+1 for i,x in enumerate(gas) if x.symbol=='O']}))
    elif command=='energy':
        from local_compute import adsorption_energy
        p.add_argument('references',type=Path); a=p.parse_args()
        print(json.dumps({'adsorption_energy_eV':adsorption_energy(json.loads(a.references.read_text()))}))
    elif command in ['rank','export']:
        p.add_argument('results',nargs='+',type=Path); p.add_argument('--out',required=True,type=Path)
        p.add_argument('--allow-smoke',action='store_true'); p.add_argument('--top',type=int,default=5)
        a=p.parse_args(); rows=[]
        if a.top<1: raise ValueError('top must be positive')
        for path in a.results:
            r=json.loads(path.read_text())
            energy=r['adsorption_energy_eV']
            if not np.isfinite(energy): raise ValueError('Nonfinite adsorption energy')
            if not a.allow_smoke and (not r['all_relaxations_converged'] or not r.get('manuscript_eligible',False)):
                raise ValueError('Unreviewed/unconverged result; only --allow-smoke can inspect it')
            refs=r['references']
            from local_compute import adsorption_energy
            if abs(adsorption_energy(refs)-energy)>1e-8: raise ValueError('Reference energy mismatch')
            rows.append({'source':str(path.resolve()),'adsorption_energy_eV':energy,'manuscript_eligible':False,
                         'model_sha256':refs['slab']['model_sha256']})
        if len({x['model_sha256'] for x in rows})!=1: raise ValueError('Mixed model rankings forbidden')
        rows.sort(key=lambda r:r['adsorption_energy_eV'])
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(json.dumps(rows[:a.top] if command=='export' else rows,indent=2))
