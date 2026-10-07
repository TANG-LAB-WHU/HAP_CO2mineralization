"""Compare the actual failure geometry to DFT and audit a full-duration half step."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from ase.io import read,write
from refine_pipeline import chemistry,forces_from_log,TYPES
from recovery_dynamics import thermo
from cp2k_validate import parse
ROOT=Path(__file__).resolve().parents[1]
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();r=a.run;d=r/'lammps-recovered';fd=r/'first-failure-dft'
    if (r/'final-instability-results.json').exists():raise FileExistsError('Preserve prior final instability audit')
    plan=json.loads((r/'final-instability-plan.json').read_text());job=json.loads((fd/'job.json').read_text());probe={'status':job['status'],'selected_frame':plan['DFT_probe']}
    if job['status']=='completed':
        cp=parse(fd/'probe.log');assert cp['output_sha256']==job['result']['output_sha256']
        index=plan['DFT_probe']['source_frame_index'];lines=(d/'trajectory-1ps.dump').read_text().splitlines();start=index*56
        assert lines[start]=='ITEM: TIMESTEP' and int(lines[start+1])==int(plan['DFT_probe']['time_fs']/.25)
        raw=np.array([[float(x) for x in line.split()] for line in lines[start+9:start+56]])
        assert np.array_equal(raw[:,0],np.arange(1,48));model_f=raw[:,8:11];dft_f=forces_from_log(fd/'probe.log')[-1];error=model_f-dft_f
        np.savez(fd/'force-comparison.npz',model_force_eV_A=model_f,DFT_force_eV_A=dft_f,types=raw[:,1].astype(int))
        t=np.load(d/'thermo-1ps.npy');truth0=float(np.load(ROOT/'.research-runs/refine-001/dataset/test/set.000/energy.npy')[0])
        dft_delta=cp['energy_eV']-truth0;model_delta=float(t[index,3]-t[0,3])
        probe.update(DFT_energy_change_from_start_eV=dft_delta,model_energy_change_from_start_eV=model_delta,
            energy_change_error_eV=model_delta-dft_delta,force_component_RMSE_eV_A=float(np.sqrt(np.mean(error**2))),
            force_max_abs_error_eV_A=float(abs(error).max()),DFT_fmax_eV_A=float(np.linalg.norm(dft_f,axis=1).max()),
            per_element_force_RMSE_eV_A={el:float(np.sqrt(np.mean(error[raw[:,1]==i+1]**2))) for i,el in enumerate(TYPES)},
            DFT_source_sha256=cp['output_sha256'],scope='One post hoc first-failure snapshot. DFT and model energy changes each use their own same-method initial energy; no mixed-method reference subtraction. Not a prospective generalization estimate and not used for fitting.')
    else:probe.update(error=job.get('error'),scope='DFT diagnosis unavailable; no replacement values generated.')
    t=thermo(d/'long-half.log');assert t.shape==(801,7) and abs(t[-1,1]-1)<1e-12 and np.array_equal(t[:,0],np.arange(0,8001,10))
    frames=read(d/'trajectory-1ps-half.dump',format='lammps-dump-text',index=':',specorder=TYPES);assert len(frames)==401
    write(d/'trajectory-1ps-half.extxyz',frames);np.save(d/'thermo-1ps-half.npy',t)
    rows=[]
    for i,frame in enumerate(frames):
        q=chemistry(frame);q['time_fs']=float(t[2*i,1]*1000);rows.append(q)
    drift=float(np.polyfit(t[:,1],(t[:,5]-t[0,5])/47*1000,1)[0])
    result={'first_failure_DFT_probe':probe,'one_ps_halfstep':{'timestep_fs':.125,'duration_fs':1000,'saved_frames':401,
        'geometry_failed_frame_count':sum(not x['geometry_gate'] for x in rows),'geometry_gate':all(x['geometry_gate'] for x in rows),
        'first_geometry_failure_time_fs':next((x['time_fs'] for x in rows if not x['geometry_gate']),None),
        'temperature_K':{'initial':float(t[0,2]),'final':float(t[-1,2]),'minimum':float(t[:,2].min()),'maximum':float(t[:,2].max())},
        'energy_drift_meV_atom_ps':drift,'drift_gate':abs(drift)<5,'maximum_energy_excursion_meV_atom':float(abs(t[:,5]-t[0,5]).max()/47*1000)},
        'frames':rows,'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [d/'long-half.in',d/'long-half.log',d/'trajectory-1ps-half.dump',r/'final-instability-plan.json']},
        'original_gates_changed':False,'further_fitting_performed':False,'production_ready':False,'manuscript_eligible':False}
    (r/'final-instability-results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['frames','source_sha256']},indent=2))
if __name__=='__main__':main()
