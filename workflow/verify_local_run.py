"""Inspect real saved artifacts; functional status is separate from scientific validity."""
import argparse,hashlib,json,re
from pathlib import Path
import numpy as np
from local_compute import adsorption_energy
from cp2k_validate import parse
p=argparse.ArgumentParser();p.add_argument('runs',type=Path);p.add_argument('--output',type=Path);a=p.parse_args();r=a.runs
report={'manuscript_eligible':False,'checks':{},'limitations':['No production convergence','OH ordering/termination unreviewed','AIMD only 1.25 fs/replica','MLIP statistical/domain validation inadequate']}
def check(name,fn):
    try:report['checks'][name]={'functional_pass':True,'details':fn()}
    except Exception as e:report['checks'][name]={'functional_pass':False,'error':repr(e)}
def mace():
    data=json.loads((r/'smoke-final/mace_results.json').read_text());refs=data['references']
    assert data['all_relaxations_converged']
    assert abs(adsorption_energy(refs)-data['adsorption_energy_eV'])<1e-8
    assert all(abs(x['reload_energy_difference_eV'])<1e-6 for x in refs.values())
    return data
def cp2k():
    data={n:parse(r/f'smoke-final/{n}.cp2k.log') for n in ['slab','co2','adsorbed']}
    data['adsorption_energy_eV']=data['adsorbed']['energy_eV']-data['slab']['energy_eV']-data['co2']['energy_eV']
    return data
def stress():
    data=json.loads((r/'stress-check/comparison.json').read_text());assert data['passed'];return data
def mlip():
    d=r/'adsorbed-dataset-v2';seen=set();counts={}
    for split in ['train','valid','test']:
        arrays={n:np.load(d/split/'set.000'/f'{n}.npy') for n in ['coord','force','energy','box']}
        assert all(np.isfinite(x).all() for x in arrays.values())
        counts[split]=len(arrays['coord']);assert counts[split]==5
        assert all(len(x)==5 for x in arrays.values())
        for coords in arrays['coord']:
            key=hashlib.sha256(np.round(coords,8).tobytes()).hexdigest()
            assert key not in seen,'Duplicate configuration leakage';seen.add(key)
    text=(d/'test.log').read_text(); rmse=float(re.findall(r'Force\s+RMSE\s*:\s*([\d.eE+-]+)',text)[-1])
    assert np.isfinite(rmse) and (d/'adsorbed-smoke.pth').is_file()
    return {'frames':counts,'unique_geometries':len(seen),'heldout_force_rmse_eV_A':rmse,'provisional_force_gate_pass':rmse<.1,'virial_validation':'not performed; no labels'}
def lammps():
    d=r/'adsorbed-dataset-v2'
    x=np.loadtxt(d/'continuous-final.dump',skiprows=9);y=np.loadtxt(d/'restart-final.dump',skiprows=9)
    assert x.shape==y.shape==(47,8) and np.isfinite(x).all() and np.isfinite(y).all()
    error=float(abs(x-y).max());assert error<1e-8
    for name in ['first','restart','continuous']:
        text=(d/f'lammps-{name}.console.log').read_text();assert 'Total wall time:' in text and 'ERROR' not in text
    return {'atoms':47,'max_saved_state_restart_difference':error,'total_steps':10,'time_fs':2.5}
def mattergen():
    d=r/'mattergen-smoke';g=json.loads((d/'generate.manifest.json').read_text());e=json.loads((d/'evaluate.manifest.json').read_text())
    assert g['status']==e['status']=='completed'
    metrics=json.loads((d/'metrics.json').read_text())
    from ase.io import read
    atoms=read(d/'generated_crystals.extxyz'); g['actual_formula']=atoms.get_chemical_formula();g['hap_chemistry_eligible']=set(['Ca','P','O','H']).issubset(atoms.get_chemical_symbols())
    filter_summary=json.loads((d/'filter_summary.json').read_text())
    assert filter_summary['input_count']==1 and filter_summary['chemistry_rejected']==1 and filter_summary['saved_count']==0
    return {'generation':g,'evaluation':e,'metrics':metrics,'filter':filter_summary,'sample_size':1}
for name,fn in [('mace',mace),('cp2k',cp2k),('stress',stress),('dataset_and_mlip',mlip),('lammps_restart',lammps),('mattergen',mattergen)]:check(name,fn)
report['all_functional_stages_complete']=all(v['functional_pass'] for v in report['checks'].values())
out=a.output or r/'verification.json';out.write_text(json.dumps(report,indent=2,allow_nan=False));print(json.dumps({k:v['functional_pass'] for k,v in report['checks'].items()},indent=2))
if not report['all_functional_stages_complete']:raise SystemExit(1)
