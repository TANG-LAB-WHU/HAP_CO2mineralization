"""Export aligned CP2K XYZ labels for an isolated CO2 plumbing test only."""
import argparse,json,re,hashlib
from pathlib import Path
import numpy as np
from ase.io import read
from cp2k_validate import parse,HARTREE_EV
BOHR_A=0.529177210903
def frames(path):
    lines=path.read_text().splitlines(); out=[]; n=int(lines[0]); stride=n+2
    if len(lines)%stride: raise ValueError('Truncated XYZ')
    for j in range(0,len(lines),stride):
        if int(lines[j])!=n: raise ValueError('Atom count changed')
        head=lines[j+1]
        step=re.search(r'i\s*=\s*(\d+)',head); energy=re.search(r'E\s*=\s*([-+\d.Ee]+)',head)
        if not step: raise ValueError('Step identity missing')
        rows=[line.split() for line in lines[j+2:j+stride]]
        out.append((int(step[1]),None if not energy else float(energy[1]),[x[0] for x in rows],np.array([x[1:4] for x in rows],float)))
    return out
def main():
    p=argparse.ArgumentParser(); p.add_argument('runs',type=Path);p.add_argument('--system',choices=['co2','adsorbed'],default='co2');p.add_argument('--output',type=Path);a=p.parse_args(); root=a.runs
    out=a.output or root/f'{a.system}-dataset'
    if out.exists():raise FileExistsError('Preserving dataset: '+str(out))
    out.mkdir()
    type_map=['C','O'] if a.system=='co2' else ['Ca','P','O','H','C']
    split={}; seen=set()
    for i,role in enumerate(['train','valid','test']):
        d=root/f'aimd-{a.system}-{i}'; parse(d/'md.log')
        xyz=frames(d/f'{a.system}_{i}-pos-1.xyz'); frc=frames(d/f'{a.system}_{i}-frc-1.xyz')
        if len(xyz)!=len(frc) or len(xyz)!=6: raise ValueError('Expected six exactly aligned MD frames')
        if any(x[0]!=f[0] or x[2]!=f[2] for x,f in zip(xyz,frc)): raise ValueError('Frame or atom identities mismatch')
        cells=np.loadtxt(d/f'{a.system}_{i}-1.cell')
        if cells.shape!=(len(xyz),12) or not np.array_equal(cells[:,0],[x[0] for x in xyz]):raise ValueError('Cell step alignment mismatch')
        if (np.linalg.det(cells[:,2:11].reshape(-1,3,3))<=0).any():raise ValueError('Invalid cell')
        # Replicas share the same starting geometry. It must not leak between splits.
        xyz,frc,cells=xyz[1:],frc[1:],cells[1:]
        for x in xyz:
            key=hashlib.sha256(np.round(x[3],8).tobytes()).hexdigest()
            if key in seen:raise ValueError('Duplicate geometry across frames/splits')
            seen.add(key)
        dest=out/role; (dest/'set.000').mkdir(parents=True,exist_ok=True)
        arrays={'coord':np.array([x[3].flatten() for x in xyz]),'force':np.array([x[3].flatten() for x in frc])*HARTREE_EV/BOHR_A,
                'energy':np.array([x[1] for x in xyz])*HARTREE_EV,'box':cells[:,2:11]}
        for name,v in arrays.items():
            if not np.isfinite(v).all(): raise ValueError('Nonfinite label')
            np.save(dest/'set.000'/f'{name}.npy',v)
        (dest/'type_map.raw').write_text('\n'.join(type_map)+'\n')
        np.savetxt(dest/'type.raw',[type_map.index(x) for x in xyz[0][2]],fmt='%d')
        split[role]={'source':d.name,'steps':[x[0] for x in xyz],'frames':len(xyz),'input_sha256':{f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in d.iterdir() if f.suffix in ['.xyz','.cell','.inp','.log']}}
    (out/'provenance.json').write_text(json.dumps({'splits':split,'temperature_K':[100,200,300],'time_per_replica_fs':1.25,
      'units':{'energy':'eV','force':'eV/A','box':'A'},'shared_initial_frame_removed':True,'unique_geometry_count':len(seen),
      'scope':a.system+' software plumbing only; nearby replicas are still correlated and coverage is inadequate','manuscript_eligible':False},indent=2))
    conf={'model':{'type_map':type_map,'descriptor':{'type':'se_e2_a','sel':[8,16] if a.system=='co2' else [32,32,96,32,16],'rcut_smth':4.5,'rcut':5.0,'neuron':[8,16,16],'axis_neuron':4,'seed':715},
       'fitting_net':{'neuron':[16,16],'seed':716}},'learning_rate':{'type':'exp','start_lr':0.001,'stop_lr':0.0001,'decay_steps':10},
       'loss':{'type':'ener','start_pref_e':1,'limit_pref_e':1,'start_pref_f':1,'limit_pref_f':1,'start_pref_v':0,'limit_pref_v':0},
       'training':{'training_data':{'systems':['train'],'batch_size':1},'validation_data':{'systems':['valid'],'batch_size':1,'numb_btch':1},
        'numb_steps':20,'seed':717,'disp_file':'lcurve.out','disp_freq':1,'save_freq':10}}
    (out/'input.json').write_text(json.dumps(conf,indent=2))
if __name__=='__main__': main()
