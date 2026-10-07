"""Audit actual fixed-model restart, timestep and 1 ps runs."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from ase.io import read,write
from ase.geometry import find_mic
from refine_pipeline import chemistry,TYPES
from recovery_dynamics import thermo
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();r=a.run;d=r/'lammps-recovered'
    if (r/'numerical-followup-results.json').exists():raise FileExistsError('Preserve existing numerical follow-up')
    for name in ['first','restart','halfstep','long']:
        q=json.loads((d/(name+'.console.execution.json')).read_text());assert q['status']=='completed' and q['returncode']==0
    x=np.loadtxt(d/'continuous-final.dump',skiprows=9);y=np.loadtxt(d/'restart-final.dump',skiprows=9);half=np.loadtxt(d/'halfstep-final.dump',skiprows=9)
    assert np.array_equal(x[:,:2],y[:,:2]) and np.array_equal(x[:,:2],half[:,:2])
    position=float(abs(x[:,2:5]-y[:,2:5]).max());velocity=float(abs(x[:,5:8]-y[:,5:8]).max())
    initial=read(d/'initial.traj');delta=find_mic(x[:,2:5]-half[:,2:5],initial.cell,True)[0]
    t=thermo(d/'long.log');h=thermo(d/'halfstep.log')
    assert t.shape==(401,7) and np.array_equal(t[:,0],np.arange(0,4001,10)) and abs(t[-1,1]-1)<1e-12
    assert h.shape==(201,7) and abs(h[-1,1]-.25)<1e-12
    frames=read(d/'trajectory-1ps.dump',format='lammps-dump-text',index=':',specorder=TYPES);assert len(frames)==401
    write(d/'trajectory-1ps.extxyz',frames);np.save(d/'thermo-1ps.npy',t);np.save(d/'thermo-halfstep.npy',h)
    records=[]
    for k,frame in enumerate(frames):
        row=chemistry(frame);row.update(time_fs=float(t[k,1]*1000),CO2_bonds_A=[float(frame.get_distance(44,j,mic=True)) for j in [45,46]])
        records.append(row)
    drift=lambda a:float(np.polyfit(a[:,1],(a[:,5]-a[0,5])/47*1000,1)[0])
    result={'restart':{'position_max_difference_A':position,'velocity_max_difference_A_ps':velocity,'strict_original_gate_pass':max(position,velocity)<1e-8},
        'half_timestep':{'timestep_fs':.125,'duration_fs':250,'energy_drift_meV_atom_ps':drift(h),'drift_gate':abs(drift(h))<5,
            'endpoint_minimum_image_position_RMS_A':float(np.sqrt(np.mean(np.sum(delta**2,axis=1)))),
            'endpoint_velocity_max_difference_A_ps':float(abs(x[:,5:8]-half[:,5:8]).max()),'position_convergence_gate':'No post hoc tolerance assigned'},
        'one_ps':{'duration_fs':1000,'timestep_fs':.25,'saved_frames':401,'geometry_failed_frame_count':sum(not x['geometry_gate'] for x in records),
            'geometry_gate':all(x['geometry_gate'] for x in records),'first_geometry_failure_time_fs':next((x['time_fs'] for x in records if not x['geometry_gate']),None),
            'temperature_K':{'initial':float(t[0,2]),'final':float(t[-1,2]),'minimum':float(t[:,2].min()),'maximum':float(t[:,2].max())},
            'energy_drift_meV_atom_ps':drift(t),'drift_gate':abs(drift(t))<5,'maximum_energy_excursion_meV_atom':float(abs(t[:,5]-t[0,5]).max()/47*1000)},
        'frames':records,'model_sha256':hashlib.sha256((d/'model.pth').read_bytes()).hexdigest(),
        'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [d/'long.in',d/'long.log',d/'trajectory-1ps.dump',d/'restart.log',d/'restart-final.dump',d/'halfstep.log',d/'halfstep-final.dump']},
        'scope':'Fixed validation-selected recovery; all original thresholds retained. One ps is a bounded numerical/stability diagnostic, not equilibrium, statistical independence or production convergence. Retrospective force test still fails.',
        'original_gates_changed':False,'production_ready':False,'manuscript_eligible':False}
    (r/'numerical-followup-results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['frames','source_sha256']},indent=2))
if __name__=='__main__':main()
