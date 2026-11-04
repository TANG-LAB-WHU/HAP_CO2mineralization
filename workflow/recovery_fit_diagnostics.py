"""Inspect actual DeepMD energy/force consistency and parameter gradients, without fitting."""
import argparse,copy,hashlib,json,sys
from pathlib import Path
import numpy as np
import torch
from ase.io import read
from deepmd.infer import DeepPot
from deepmd.pt.model.model import get_model
from deepmd.pt.train.wrapper import ModelWrapper

def main():
    p=argparse.ArgumentParser();p.add_argument('--baseline',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError('Preserve the existing diagnosis')
    torch.set_num_threads(4);torch.cuda.set_per_process_memory_fraction(.35)
    folder=a.baseline/'models/overfit';data=a.baseline/'dataset/overfit';conf=json.loads((folder/'out.json').read_text())
    state=torch.load(folder/'model.ckpt-1500.pt',map_location='cpu',weights_only=False)
    wrapper=ModelWrapper(get_model(copy.deepcopy(conf['model']))).to('cuda');wrapper.load_state_dict(state['model']);wrapper.train()
    coord=np.load(data/'set.000/coord.npy').reshape(1,-1,3);box=np.load(data/'set.000/box.npy').reshape(1,3,3);types=np.loadtxt(data/'type.raw',dtype=int)
    truth_f=np.load(data/'set.000/force.npy').reshape(1,-1,3);truth_e=np.load(data/'set.000/energy.npy').reshape(-1)
    ct=torch.tensor(coord,dtype=torch.float64,device='cuda');bt=torch.tensor(box,dtype=torch.float64,device='cuda');tt=torch.tensor(types[None,:],device='cuda')
    tf=torch.tensor(truth_f,dtype=torch.float64,device='cuda')
    def predict(c=ct):return wrapper(c,tt,box=bt,inference_only=True)[0]
    pred=predict();f=pred['force'];e=pred['energy'];loss=((f-tf)**2).mean()
    params=[(n,v) for n,v in wrapper.named_parameters() if v.requires_grad]
    gradients=torch.autograd.grad(loss,[v for _,v in params],allow_unused=True)
    groups={};gradient_rows=[]
    for (name,value),g in zip(params,gradients):
        group='descriptor' if 'descriptor' in name else 'fitting'
        groups.setdefault(group,{'parameters':0,'gradient_norm_squared':0.,'nonzero_gradient_tensors':0,'unused_tensors':0})
        groups[group]['parameters']+=value.numel()
        if g is None:groups[group]['unused_tensors']+=1;continue
        norm=float(torch.linalg.vector_norm(g));groups[group]['gradient_norm_squared']+=norm**2;groups[group]['nonzero_gradient_tensors']+=int(norm>0)
        gradient_rows.append({'parameter':name,'gradient_norm':norm,'parameter_norm':float(torch.linalg.vector_norm(value.detach()))})
    for group in groups.values():group['gradient_norm']=float(np.sqrt(group.pop('gradient_norm_squared')))
    totalnorm=float(np.sqrt(sum(x['gradient_norm']**2 for x in groups.values())))
    assert totalnorm>0 and np.isfinite(totalnorm),'Missing/nonfinite force-loss gradients'
    values=[]
    for sign in [-1,1]:
        wrapper.load_state_dict(state['model'])
        with torch.no_grad():
            for (_,parameter),g in zip(params,gradients):
                if g is not None:parameter.add_(g,alpha=sign*1e-5/totalnorm)
        values.append(float(((predict()['force']-tf)**2).mean().detach()))
    wrapper.load_state_dict(state['model'])
    fd_parameter=(values[1]-values[0])/(2e-5)
    frozen=DeepPot(str(folder/'model.pth'));ef,ff,_=frozen.eval(coord,box,types)
    error=ff-truth_f;flat=np.argsort(abs(error).ravel())[-8:][::-1];checks=[]
    selected=set(int(x) for x in flat[:4]);selected.add(int(np.argmax(abs(truth_f))))
    for index in sorted(selected):
        atom,axis=divmod(index,3)
        for delta in [1e-3,1e-4]:
            plus=coord.copy();minus=coord.copy();plus[0,atom,axis]+=delta;minus[0,atom,axis]-=delta
            ep=frozen.eval(plus,box,types)[0].item();em=frozen.eval(minus,box,types)[0].item();fd=-(ep-em)/(2*delta)
            checks.append({'atom_id':atom+1,'axis':'xyz'[axis],'delta_A':delta,'analytic_force_eV_A':float(ff[0,atom,axis]),'finite_difference_force_eV_A':fd,'absolute_difference_eV_A':abs(fd-ff[0,atom,axis])})
    species=conf['model']['type_map'];atoms=read(a.baseline/'dataset/train.traj',index=0)
    result={'scope':'Original one-frame fit diagnosis; no parameter update or test-set selection.',
        'one_frame_force_rmse_eV_A':float(np.sqrt(np.mean(error**2))),'one_frame_energy_error_eV_atom':float((ef.item()-truth_e[0])/len(types)),
        'frozen_vs_trainable_force_max_difference_eV_A':float(abs(ff-f.detach().cpu().numpy()).max()),'frozen_vs_trainable_energy_difference_eV':float(abs(ef.item()-e.detach().item())),
        'force_loss_gradient_groups':groups,'top_gradient_tensors':sorted(gradient_rows,key=lambda x:x['gradient_norm'],reverse=True)[:12],
        'parameter_direction_derivative':{'analytic':totalnorm,'finite_difference':fd_parameter,'relative_difference':abs(fd_parameter-totalnorm)/totalnorm,'parameter_step':1e-5},
        'energy_force_finite_differences':checks,
        'largest_force_residuals':[{'atom_id':int(i//3+1),'element':species[types[i//3]],'axis':'xyz'[i%3],'reference_eV_A':float(truth_f.ravel()[i]),'prediction_eV_A':float(ff.ravel()[i]),'error_eV_A':float(error.ravel()[i])} for i in flat],
        'fixed_atom_count':sum(len(x.index) for x in atoms.constraints),'fitted_atom_count':len(types),'reference_force_sum_eV_A':truth_f[0].sum(axis=0).tolist(),
        'training_force_scale_rms_eV_A':float(np.sqrt(np.mean(truth_f**2))),
        'checkpoint_sha256':hashlib.sha256((folder/'model.ckpt-1500.pt').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'manuscript_eligible':False}
    a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
