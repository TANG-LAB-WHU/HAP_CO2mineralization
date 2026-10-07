"""Select one predeclared optimizer recovery by validation; retain test reuse caveat."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import torch
from deepmd.infer import DeepPot
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--baseline',type=Path,required=True);a=p.parse_args();r=a.run;b=a.baseline
    if (r/'fullbatch-evaluation.json').exists():raise FileExistsError('Preserve previous model selection')
    plan=load(r/'fullbatch-recovery-plan.json');fit=load(r/'wider-fullbatch/result.json')
    assert fit['training_frames']==16 and fit['test_set_read'] is False
    torch.set_num_threads(4);torch.cuda.set_per_process_memory_fraction(.35)
    model=DeepPot(str(r/'wider-fullbatch/model.pth'));results={}
    def evaluate(role):
        d=b/'dataset'/role;coords=np.load(d/'set.000/coord.npy');box=np.load(d/'set.000/box.npy');types=np.loadtxt(d/'type.raw',dtype=int)
        e,f,_=model.eval(coords,box,types);refe=np.load(d/'set.000/energy.npy').reshape(-1);reff=np.load(d/'set.000/force.npy').reshape(f.shape)
        assert np.isfinite(e).all() and np.isfinite(f).all()
        np.savez(r/'wider-fullbatch'/f'{role}_predictions.npz',reference_energy_eV=refe,predicted_energy_eV=e.reshape(-1),reference_force_eV_A=reff,predicted_force_eV_A=f)
        return {'frames':len(e),'force_rmse_eV_A':float(np.sqrt(np.mean((f-reff)**2))),'energy_rmse_eV_atom':float(np.sqrt(np.mean(((e.reshape(-1)-refe)/len(types))**2))),
            'per_element_force_rmse_eV_A':{s:float(np.sqrt(np.mean((f[:,types==i]-reff[:,types==i])**2))) for i,s in enumerate(model.get_type_map())}}
    for role in ['train','valid']:results[role]=evaluate(role)
    original=load(b/'model_evaluation.json')['models']['wider']['valid']['force_rmse_eV_A']
    selected=results['valid']['force_rmse_eV_A']<original
    selection={'selection_rule':plan['selection'],'original_validation_force_eV_A':original,'candidate_validation_force_eV_A':results['valid']['force_rmse_eV_A'],'recovery_selected':selected,'model_sha256':sha(r/'wider-fullbatch/model.pth')}
    # Save the decision before opening any test arrays.
    (r/'validation-selection.json').write_text(json.dumps(selection,indent=2)+'\n')
    if selected:results['test']=evaluate('test')
    out={**selection,'candidate':results,'training_fit':fit,'test_reuse_caveat':plan['old_test_caveat'],'test_fitted':False,'source_sha256':sha(Path(__file__)),
         'plan_sha256':sha(r/'fullbatch-recovery-plan.json'),'manuscript_eligible':False,'production_ready':False}
    if selected:out['retrospective_gates']={'force_rmse_below_0p1':results['test']['force_rmse_eV_A']<.1,'energy_rmse_below_10_meV_atom':results['test']['energy_rmse_eV_atom']<.01}
    (r/'fullbatch-evaluation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
