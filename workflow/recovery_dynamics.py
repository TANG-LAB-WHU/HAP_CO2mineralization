"""Audit the separately retained recovery trajectory without replacing original gates."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from ase.io import read,write
from refine_pipeline import chemistry,TYPES

def thermo(path):
    text=path.read_text(encoding='utf-8')
    if 'Total wall time:' not in text or 'ERROR:' in text:raise ValueError('Incomplete or failed LAMMPS log')
    rows=[];active=False;finished=False
    for line in text.splitlines():
        if line.split()==['Step','Time','Temp','PotEng','KinEng','TotEng','Press']:
            active=True;continue
        if active:
            if line.lstrip().startswith('Loop time of'):
                finished=True;break
            fields=line.split()
            if len(fields)!=7:raise ValueError('Malformed thermodynamic row')
            try:rows.append([float(x) for x in fields])
            except ValueError as exc:raise ValueError('Non-numeric thermodynamic row') from exc
    array=np.array(rows)
    if not finished or not rows or not np.isfinite(array).all():raise ValueError('Missing or nonfinite thermodynamic block')
    if len(array)>1 and (np.any(np.diff(array[:,0])<=0) or np.any(np.diff(array[:,1])<=0)):
        raise ValueError('Thermodynamic steps/times are not strictly increasing')
    return array

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();r=a.run;d=r/'lammps-recovered'
    if (r/'recovered-dynamics.json').exists():raise FileExistsError('Preserve the prior recovery trajectory audit')
    t=thermo(d/'continuous.log');mapping=thermo(d/'mapping.log');assert t.shape==(101,7) and np.isfinite(t).all()
    assert np.array_equal(t[:,0],np.arange(0,1001,10)) and abs(t[-1,1]-.25)<1e-12
    frames=read(d/'trajectory.dump',format='lammps-dump-text',index=':',specorder=TYPES);assert len(frames)==101
    write(d/'trajectory.extxyz',frames);np.save(d/'thermo.npy',t)
    x=np.loadtxt(d/'mapping.dump',skiprows=9);pred=np.load(r/'wider-fullbatch/test_predictions.npz')
    assert np.array_equal(x[:,0],np.arange(1,48))
    force_error=float(abs(x[:,5:8]-pred['predicted_force_eV_A'][0]).max())
    energy_error=float(abs(mapping[0,3]-pred['predicted_energy_eV'][0]))
    records=[]
    for k,frame in enumerate(frames):
        row=chemistry(frame);row['time_fs']=float(t[k,1]*1000);records.append(row)
    drift=float(np.polyfit(t[:,1],(t[:,5]-t[0,5])/47*1000,1)[0])
    result={'saved_frames':len(frames),'time_fs':250,'timestep_fs':.25,'ensemble':'NVE','constraints':'all 47 atoms mobile',
        'engine_force_max_difference_eV_A':force_error,'engine_energy_difference_eV':energy_error,
        'engine_mapping_gate':force_error<1e-5 and energy_error<1e-5,
        'temperature_K':{'initial':float(t[0,2]),'final':float(t[-1,2]),'minimum':float(t[:,2].min()),'maximum':float(t[:,2].max())},
        'geometry_failed_frame_count':sum(not x['geometry_gate'] for x in records),'geometry_gate':all(x['geometry_gate'] for x in records),
        'first_geometry_failure_time_fs':next((x['time_fs'] for x in records if not x['geometry_gate']),None),
        'energy_drift_meV_atom_ps':drift,'drift_gate':abs(drift)<5,'maximum_energy_excursion_meV_atom':float(abs(t[:,5]-t[0,5]).max()/47*1000),
        'frames':records,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [d/'model.pth',d/'continuous.in',d/'continuous.log',d/'trajectory.dump',d/'mapping.dump']},
        'scope':'One retrospective 250 fs trajectory with a validation-selected optimizer recovery; no untouched prospective holdout, half-step recovery run or production stability claim.',
        'original_gates_changed':False,'production_ready':False,'manuscript_eligible':False}
    (r/'recovered-dynamics.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['frames','source_sha256']},indent=2))
if __name__=='__main__':main()
