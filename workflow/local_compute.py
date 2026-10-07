"""Bounded HAP/CO2 functional smoke calculations, never production evidence."""
from __future__ import annotations
import argparse, hashlib, json, platform, subprocess, sys, time
from datetime import datetime, timezone
from collections import Counter
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')

def prepare(out):
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor
    from ase.build import surface, molecule
    from ase.constraints import FixAtoms
    from ase.io import write
    source = ROOT / '3.ML_AIMD/Step4_lammps_scaleup/HAP_ConventionalCell.cif'
    s = Structure.from_file(source)
    keep = []
    for site in s:
        if site.is_ordered:
            keep.append(site)
        else:
            element = next(iter(site.species)).symbol
            allowed = {'H': [0.0608, 0.5608], 'O': [0.1978, 0.6978]}
            if element not in allowed:
                raise ValueError('Unexpected disorder; explicit ordering required')
            if any(abs(site.frac_coords[2] - z) < 1e-4 for z in allowed[element]):
                keep.append(site)
    ordered = Structure(s.lattice, [next(iter(x.species)).symbol for x in keep], [x.frac_coords for x in keep])
    bulk = AseAtomsAdaptor.get_atoms(ordered)
    assert len(bulk) == 44 and bulk.get_chemical_formula() == 'H2Ca10O26P6'
    distance = bulk.get_all_distances(mic=True)
    np.fill_diagonal(distance, np.inf)
    if distance.min() < 0.7:
        raise ValueError('Unphysical bulk contact')
    slab = surface(bulk, (0, 0, 1), 1, vacuum=12, periodic=True)
    slab.set_constraint(FixAtoms(mask=slab.positions[:, 2] < np.median(slab.positions[:, 2])))
    gas = molecule('CO2')
    gas.rotate('z', 'x')
    gas.set_cell(slab.cell)
    gas.set_pbc(True)
    gas.center()
    top = max((i for i, a in enumerate(slab) if a.symbol == 'Ca'), key=lambda i: slab.positions[i, 2])
    adsorbate = gas.copy()
    adsorbate.translate(slab.positions[top] + [0, 0, 2.8] - adsorbate.positions[0])
    combined = slab + adsorbate
    combined.set_constraint(slab.constraints)
    groups = {'slab_ids': list(range(1, len(slab)+1)), 'co2_ids': list(range(len(slab)+1, len(combined)+1)),
              'co2_oxygen_ids': [len(slab)+i+1 for i, a in enumerate(adsorbate) if a.symbol == 'O']}
    for name, atoms in [('bulk', bulk), ('slab', slab), ('co2', gas), ('adsorbed', combined)]:
        write(out / f'{name}.extxyz', atoms)
    save(out/'groups.json', groups)
    save(out/'structure_review.json', {'source': str(source.relative_to(ROOT)), 'source_sha256': sha(source),
         'formula': bulk.get_chemical_formula(), 'atoms': len(bulk), 'minimum_bulk_distance_A': float(distance.min()),
         'OH_ordering_fractional_z': {'O': [0.1978,0.6978], 'H': [0.0608,0.5608]},
         'facet': [0,0,1], 'bulk_repetitions': 1, 'vacuum_each_side_A':12,
         'review': 'computational sanity checks only; human ordering/termination review pending',
         'manuscript_eligible': False})

def adsorption_energy(records):
    if set(records) != {'slab','co2','adsorbed'}:
        raise ValueError('Three matched references required')
    for r in records.values():
        if not np.isfinite(r['energy_eV']):
            raise ValueError('Nonfinite energy')
    if len({r['model_sha256'] for r in records.values()}) != 1:
        raise ValueError('Mixed model references')
    if all('composition' in r for r in records.values()):
        if Counter(records['slab']['composition']) + Counter(records['co2']['composition']) != Counter(records['adsorbed']['composition']):
            raise ValueError('References do not conserve chemical composition')
        if not all(np.allclose(r['cell_A'], records['slab']['cell_A']) for r in records.values()):
            raise ValueError('Reference cells differ')
    return records['adsorbed']['energy_eV']-records['slab']['energy_eV']-records['co2']['energy_eV']

def mace_run(out, model, steps):
    import torch
    from ase.io import read, write
    from ase.optimize import FIRE
    from mace.calculators import MACECalculator
    if not torch.cuda.is_available():
        raise RuntimeError('This profile requires CUDA; no silent CPU fallback')
    torch.set_num_threads(4)
    torch.cuda.set_per_process_memory_fraction(0.35)
    calc = MACECalculator(model_paths=str(model), device='cuda', default_dtype='float64')
    results = {}
    for name in ['slab','co2','adsorbed']:
        atoms = read(out/f'{name}.extxyz')
        atoms.calc = calc
        initial = float(atoms.get_potential_energy())
        opt = FIRE(atoms, logfile=str(out/f'{name}.opt.log'), trajectory=str(out/f'{name}.traj'), maxstep=0.04)
        converged = opt.run(fmax=0.05, steps=steps)
        energy = float(atoms.get_potential_energy())
        forces = atoms.get_forces()
        if not np.isfinite(forces).all():
            raise RuntimeError('Nonfinite forces')
        results[name] = {'energy_eV': energy, 'initial_energy_eV':initial,
                         'composition':dict(Counter(atoms.get_chemical_symbols())), 'cell_A':atoms.cell.array.tolist(),
                         'max_free_force_eV_A':float(np.linalg.norm(forces,axis=1).max()),
                         'converged':bool(converged), 'steps':opt.nsteps, 'model_sha256':sha(model)}
        write(out/f'{name}.relaxed.extxyz',atoms)
        # Independent reload must reproduce the saved geometry's single-point energy.
        restarted=read(out/f'{name}.relaxed.extxyz'); restarted.calc=calc
        results[name]['reload_energy_difference_eV']=float(restarted.get_potential_energy()-energy)
        if abs(results[name]['reload_energy_difference_eV']) > 1e-6:
            raise RuntimeError('Reload energy mismatch')
    save(out/'mace_results.json', {'references':results, 'adsorption_energy_eV':adsorption_energy(results),
        'all_relaxations_converged':all(r['converged'] for r in results.values()),
        'peak_cuda_allocated_bytes':torch.cuda.max_memory_allocated(), 'manuscript_eligible':False,
        'interpretation':'functional smoke only; unconverged slab size, vacuum, ordering, sites and ML uncertainty'})

def cp2k_input(out, name):
    from ase.io import read
    atoms = read(out/f'{name}.relaxed.extxyz')
    cell='\n'.join('      '+axis+' '+' '.join(f'{v:.12f}' for v in vector) for axis,vector in zip('ABC',atoms.cell))
    coord='\n'.join(f'      {a.symbol} '+' '.join(f'{v:.12f}' for v in a.position) for a in atoms)
    q={'Ca':10,'P':5,'O':6,'H':1,'C':4}
    kinds='\n'.join(f'    &KIND {el}\n      BASIS_SET DZVP-MOLOPT-SR-GTH\n      POTENTIAL GTH-PBE-q{q[el]}\n    &END KIND' for el in sorted(set(atoms.get_chemical_symbols())))
    text=f'''&GLOBAL
  PROJECT {name}
  RUN_TYPE ENERGY_FORCE
  PRINT_LEVEL MEDIUM
&END GLOBAL
&FORCE_EVAL
  METHOD QUICKSTEP
  STRESS_TENSOR ANALYTICAL
  &DFT
    BASIS_SET_FILE_NAME /opt/cp2k/data/BASIS_MOLOPT
    POTENTIAL_FILE_NAME /opt/cp2k/data/GTH_POTENTIALS
    &MGRID
      CUTOFF 300
      REL_CUTOFF 40
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-10
    &END QS
    &SCF
      EPS_SCF 1.0E-6
      MAX_SCF 80
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        EPS_SCF 1.0E-6
        MAX_SCF 10
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
{cell}
      PERIODIC XYZ
    &END CELL
    &COORD
{coord}
    &END COORD
{kinds}
  &END SUBSYS
  &PRINT
    &STRESS_TENSOR ON
      STRESS_UNIT bar
    &END STRESS_TENSOR
    &FORCES ON
    &END FORCES
  &END PRINT
&END FORCE_EVAL
'''
    (out/f'{name}.inp').write_text(text)

def main():
    p=argparse.ArgumentParser(__doc__)
    p.add_argument('stage',choices=['prepare','mace','cp2k-inputs'])
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--model',type=Path)
    p.add_argument('--steps',type=int,default=12)
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    if a.steps < 1 or a.steps > 200: p.error('smoke steps must be 1..200')
    if (a.out/f'{a.stage}.manifest.json').exists():
        p.error('Preserving existing run manifest; choose a fresh output directory')
    start=time.time()
    manifest={'schema_version':1,'run_id':a.out.name,'stage':a.stage,'backend':'local','status':'running',
       'configuration':vars(a)|{'out':str(a.out),'model':str(a.model)},'environment':{'python':sys.version,'platform':platform.platform()},
       'validation':{'passed':False,'reviewer':None,'criteria':['finite energies and forces','saved geometry energy repeatability < 1e-6 eV']},
       'manuscript_eligible':False,'started_at':datetime.now(timezone.utc).isoformat(),
       'inputs':[{'path':str(f),'sha256':sha(f)} for f in sorted(a.out.glob('*.extxyz'))],
       'source_sha256':sha(__file__), 'command':sys.argv}
    if a.model and a.model.is_file(): manifest['inputs'].append({'path':str(a.model),'sha256':sha(a.model)})
    try:
        manifest['git']={'commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                         'branch':subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),
                         'dirty':bool(subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip())}
    except (OSError,subprocess.CalledProcessError): manifest['git']={'unavailable':True}
    try:
        if a.stage=='prepare': prepare(a.out)
        elif a.stage=='mace': mace_run(a.out,a.model,a.steps)
        else:
            for n in ['slab','co2','adsorbed']: cp2k_input(a.out,n)
        manifest['status']='completed'; manifest['validation']['passed']=True
    except Exception as e:
        manifest['status']='failed'; manifest['error']=repr(e); raise
    finally:
        manifest['elapsed_seconds']=time.time()-start
        manifest['finished_at']=datetime.now(timezone.utc).isoformat()
        manifest['outputs']=[{'path':str(f.name),'sha256':sha(f)} for f in sorted(a.out.iterdir()) if f.is_file() and not f.name.endswith('manifest.json')]
        save(a.out/f'{a.stage}.manifest.json',manifest)

if __name__=='__main__': main()
