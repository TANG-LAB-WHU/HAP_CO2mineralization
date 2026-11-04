"""Bounded full-batch optimization of a fixed training set; never reads test labels."""
import argparse,copy,hashlib,json,time
from pathlib import Path
import numpy as np
import torch
from deepmd.infer import DeepPot
from deepmd.pt.model.model import get_model
from deepmd.pt.train.wrapper import ModelWrapper

def main():
    p=argparse.ArgumentParser();p.add_argument('--checkpoint',type=Path,required=True);p.add_argument('--resolved-config',type=Path,required=True);p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--max-outer',type=int,default=30);p.add_argument('--force-target',type=float,default=.01);p.add_argument('--max-closures',type=int,default=1200);p.add_argument('--max-seconds',type=float,default=1200);p.add_argument('--resume-lbfgs',action='store_true');a=p.parse_args()
    if min(a.max_outer,a.max_closures,a.max_seconds,a.force_target)<=0:raise ValueError('Positive optimization limits required')
    a.out.mkdir(parents=True,exist_ok=False);(a.out/'optimizer-source.py').write_bytes(Path(__file__).read_bytes());torch.set_num_threads(4);torch.cuda.set_per_process_memory_fraction(.35)
    config=json.loads(a.resolved_config.read_text());state=torch.load(a.checkpoint,map_location='cpu',weights_only=False)
    wrapper=ModelWrapper(get_model(copy.deepcopy(config['model']))).to('cuda');wrapper.load_state_dict(state['model']);wrapper.train()
    types=np.loadtxt(a.data/'type.raw',dtype=int);natoms=len(types)
    coord=np.load(a.data/'set.000/coord.npy').reshape(-1,natoms,3);box=np.load(a.data/'set.000/box.npy').reshape(-1,3,3);types=np.loadtxt(a.data/'type.raw',dtype=int)
    reference_f=np.load(a.data/'set.000/force.npy').reshape(-1,natoms,3);reference_e=np.load(a.data/'set.000/energy.npy').reshape(-1,1)
    ct=torch.tensor(coord,dtype=torch.float64,device='cuda');bt=torch.tensor(box,dtype=torch.float64,device='cuda');tt=torch.tensor(np.tile(types,(len(coord),1)),device='cuda')
    tf=torch.tensor(reference_f,dtype=torch.float64,device='cuda');te=torch.tensor(reference_e,dtype=torch.float64,device='cuda')
    optimizer=torch.optim.LBFGS(wrapper.parameters(),lr=1.,max_iter=20,max_eval=25,history_size=20,line_search_fn='strong_wolfe',tolerance_grad=1e-10,tolerance_change=1e-13)
    if a.resume_lbfgs:
        if state.get('optimizer_kind')!='torch.optim.LBFGS':raise ValueError('Cannot resume a different optimizer')
        optimizer.load_state_dict(state['optimizer'])
    class BudgetReached(Exception):pass
    calls=0;start=time.monotonic();history=[]
    def forward():
        pred=wrapper(ct,tt,box=bt,inference_only=True)[0]
        lf=((pred['force']-tf)**2).mean();le=(((pred['energy']-te)/natoms)**2).mean()
        return lf+le,lf,le
    def closure():
        nonlocal calls
        if calls>=a.max_closures or time.monotonic()-start>a.max_seconds:raise BudgetReached('Explicit optimization budget reached')
        optimizer.zero_grad();loss,_,_=forward()
        if not torch.isfinite(loss):raise RuntimeError('Nonfinite training loss')
        loss.backward();calls+=1
        return loss
    def record(outer):
        loss,lf,le=forward();row={'outer':outer,'closure_calls':calls,'loss':float(loss.detach()),'force_rmse_eV_A':float(lf.detach().sqrt()),'energy_rmse_eV_atom':float(le.detach().sqrt()),'elapsed_seconds':time.monotonic()-start}
        history.append(row);print(json.dumps(row),flush=True);(a.out/'history.json').write_text(json.dumps(history,indent=2),encoding='utf-8');return row
    result=record(0);status='budget_complete'
    for outer in range(1,a.max_outer+1):
        accepted=copy.deepcopy(wrapper.state_dict());accepted_optimizer=copy.deepcopy(optimizer.state_dict())
        try:optimizer.step(closure)
        except BudgetReached as exc:
            wrapper.load_state_dict(accepted);optimizer.load_state_dict(accepted_optimizer);status='bounded_stop';(a.out/'stop_reason.txt').write_text(str(exc));break
        result=record(outer)
        if result['force_rmse_eV_A']<a.force_target and result['energy_rmse_eV_atom']<.001:status='internal_target_reached';break
    torch.save({'model':wrapper.state_dict(),'optimizer':optimizer.state_dict(),'optimizer_kind':'torch.optim.LBFGS','closure_calls':calls},a.out/'lbfgs-state.pt')
    wrapper.eval();torch.jit.script(wrapper.model['Default']).save(str(a.out/'model.pth'))
    exported=DeepPot(str(a.out/'model.pth'));e,f,_=exported.eval(coord,box,types)
    final={'force_rmse_eV_A':float(np.sqrt(np.mean((f-reference_f)**2))),'energy_rmse_eV_atom':float(np.sqrt(np.mean(((e.reshape(-1,1)-reference_e)/natoms)**2)))}
    assert abs(final['force_rmse_eV_A']-result['force_rmse_eV_A'])<1e-8,'Frozen model does not match trainable model'
    report={'status':status,'training_frames':len(coord),'closure_calls':calls,'final':final,'original_force_gate_eV_A':.02 if len(coord)==1 else None,
        'original_force_gate_pass':final['force_rmse_eV_A']<.02 if len(coord)==1 else None,'test_set_read':False,
        'loss':'mean force-component squared error + mean per-atom energy squared error; original physical units retained',
        'optimizer':{'type':'torch.optim.LBFGS','lr':1.,'history_size':20,'max_outer':a.max_outer,'max_iter_per_outer':20,'max_eval_per_outer':25,'hard_closure_cap':a.max_closures,'hard_elapsed_cap_seconds':a.max_seconds,'resumed_optimizer_history':a.resume_lbfgs,'line_search':'strong_wolfe'},
        'source_checkpoint_sha256':hashlib.sha256(a.checkpoint.read_bytes()).hexdigest(),'resolved_config_sha256':hashlib.sha256(a.resolved_config.read_bytes()).hexdigest(),
        'dataset_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (a.data/'set.000').glob('*.npy')},'model_sha256':hashlib.sha256((a.out/'model.pth').read_bytes()).hexdigest(),
        'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'manuscript_eligible':False}
    (a.out/'result.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8');print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':main()
