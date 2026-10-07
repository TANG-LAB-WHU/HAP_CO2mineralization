"""Plot actual recovery predictions/dynamics and animate saved frames only."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from matplotlib.lines import Line2D
from ase.io import read
from plot_refinement import COLORS
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();r=a.run.resolve();d=r/'lammps-recovered'
    if (r/'reproducibility-manifest.json').exists():raise RuntimeError('Captured recovery is immutable')
    evaluation=load(r/'fullbatch-evaluation.json');numerics=load(r/'numerical-followup-results.json');old=ROOT/'.research-runs/refine-001'
    baseline=load(old/'model_evaluation.json')['models']['wider'];t=np.load(d/'thermo-1ps.npy');h=np.load(d/'thermo-1ps-half.npy');previous=np.load(old/'lammps/thermo.npy')
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(10,7.8));fig.suptitle('recovery-001 | Optimizer improvement; force gate still fails',x=.06,ha='left')
    xx=np.arange(3);roles=['train','valid','test']
    ax[0,0].bar(xx-.18,[baseline[k]['force_rmse_eV_A'] for k in roles],.36,label='Original wider',color='#c47662')
    ax[0,0].bar(xx+.18,[evaluation['candidate'][k]['force_rmse_eV_A'] for k in roles],.36,label='Recovery',color='#28968a')
    ax[0,0].axhline(.1,color='#a52323',ls='--',label='Test force gate');ax[0,0].set_xticks(xx,['Training','Validation','Reused test']);ax[0,0].set(ylabel='Force RMSE (eV/A)',title='Same labels and architecture');ax[0,0].legend(fontsize=8)
    ax[0,1].plot(previous[:,1]*1000,previous[:,2],color='#a52323',label='Original: 250 fs');ax[0,1].plot(t[:,1]*1000,t[:,2],color='#28968a',label='Recovery: 1 ps');ax[0,1].set(xlabel='Time (fs)',ylabel='NVE temperature (K)',title='Both trajectories are actual calculations');ax[0,1].legend(fontsize=8)
    times=np.array([x['time_fs'] for x in numerics['frames']]);oh=np.array([x['O_H_distances_A'] for x in numerics['frames']])
    for i in range(oh.shape[1]):ax[1,0].plot(times,oh[:,i],label=f'OH {i+1}')
    ax[1,0].axhspan(.8,1.2,color='#dbe9e4',alpha=.5,zorder=-1);ax[1,0].set(xlabel='Recovery time (fs)',ylabel='Nearest O-H distance (A)',title=f"Geometry failures: {numerics['one_ps']['geometry_failed_frame_count']}/401");ax[1,0].legend(fontsize=8)
    ax[1,1].plot(t[:,1]*1000,(t[:,5]-t[0,5])/47*1000,label='dt=0.25 fs; 1 ps')
    ax[1,1].plot(h[:,1]*1000,(h[:,5]-h[0,5])/47*1000,label='dt=0.125 fs; 1 ps')
    ax[1,1].set(xlabel='Time (fs)',ylabel='Total-energy change (meV/atom)',title='Same 1 ps duration; different timesteps');ax[1,1].legend(fontsize=8)
    fig.text(.06,.02,'Retrospective local diagnostics. One ps does not establish equilibrium, production accuracy or mineralization.',fontsize=8)
    fig.tight_layout(rect=[0,.04,1,.94]);fig.savefig(r/'10_recovered_dynamics.pdf');fig.savefig(r/'10_recovered_dynamics.png',dpi=170);plt.close(fig)
    final=load(r/'final-instability-results.json');probe=final['first_failure_DFT_probe']
    if probe['status']=='completed':
        fig,ax=plt.subplots(1,2,figsize=(10,4.6));fig.suptitle('Post hoc DFT diagnostic at the first failed configuration',x=.06,ha='left')
        ax[0].bar(['DFT PBE','Learned model'],[probe['DFT_energy_change_from_start_eV'],probe['model_energy_change_from_start_eV']],color=['#28968a','#c47662']);ax[0].axhline(0,color='#777777',lw=.8);ax[0].set(ylabel='Energy change from initial state (eV)',title=f"Same 47-atom cell; {probe['selected_frame']['time_fs']:.0f} fs")
        f=np.load(r/'first-failure-dft/force-comparison.npz');truth=f['DFT_force_eV_A'];prediction=f['model_force_eV_A'];types=f['types'];bound=max(abs(truth).max(),abs(prediction).max())*1.05
        for i,element in enumerate(['Ca','P','O','H','C']):
            mask=types==i+1;ax[1].scatter(truth[mask].ravel(),prediction[mask].ravel(),s=15,label=element,color=COLORS[element],alpha=.8)
        ax[1].plot([-bound,bound],[-bound,bound],color='#777777',ls='--');ax[1].set(xlabel='DFT force (eV/A)',ylabel='Model force (eV/A)',title=f"Force RMSE {probe['force_component_RMSE_eV_A']:.3f} eV/A");ax[1].legend(fontsize=8,ncol=3)
        fig.text(.06,.02,'One selected failure snapshot; no further fitting, no unbiased-test or material-property claim.',fontsize=8)
        fig.tight_layout(rect=[0,.06,1,.9]);fig.savefig(r/'12_failure_dft.pdf');fig.savefig(r/'12_failure_dft.png',dpi=170);plt.close(fig)
    path=d/'trajectory-1ps.extxyz';frames=read(path,index=':');indices=np.arange(0,len(frames),4);chosen=[frames[i] for i in indices]
    assert len(frames)==len(t)==401 and len(chosen)==101
    fig,axes=plt.subplots(1,2,figsize=(10,6));fig.subplots_adjust(top=.82,bottom=.14)
    allpos=np.concatenate([x.positions for x in chosen]);lo=allpos.min(axis=0)-1;hi=allpos.max(axis=0)+1;scat=[];colors=[COLORS[s] for s in frames[0].get_chemical_symbols()]
    for axis,view in zip(axes,[(0,2),(0,1)]):
        scat.append(axis.scatter([],[],s=45,edgecolors='white',linewidths=.3));axis.set(xlim=(lo[view[0]],hi[view[0]]),ylim=(lo[view[1]],hi[view[1]]),xlabel='xyz'[view[0]]+' (A)',ylabel='xyz'[view[1]]+' (A)');axis.set_aspect('equal')
    title=fig.suptitle('');fig.legend(handles=[Line2D([],[],marker='o',ls='',color=c,label=s) for s,c in COLORS.items()],ncol=5,loc='upper center',bbox_to_anchor=(.5,.91),fontsize=11)
    fig.text(.02,.03,'Actual saved frames; no interpolation. FAILED force / 1 ps geometry validation; no mineralization claim.',fontsize=9,color='#9c6416')
    def update(k):
        index=indices[k]
        for scatter,view in zip(scat,[(0,2),(0,1)]):scatter.set_offsets(chosen[k].positions[:,view]);scatter.set_color(colors)
        geometry='PASS' if numerics['frames'][index]['geometry_gate'] else 'FAIL'
        title.set_text(f"recovery-001 | NVE | t = {t[index,1]*1000:.0f} fs | dt = 0.25 fs\nT = {t[index,2]:.1f} K | geometry {geometry} | all atoms mobile")
        return scat+[title]
    animation=FuncAnimation(fig,update,frames=len(chosen),interval=120,blit=False);animation.save(r/'11_recovered_actual_1ps.gif',writer=PillowWriter(fps=8));plt.close(fig)
    sources=[r/'fullbatch-evaluation.json',r/'numerical-followup-results.json',r/'final-instability-results.json',r/'first-failure-dft/force-comparison.npz',d/'thermo-1ps.npy',d/'thermo-1ps-half.npy',path,old/'model_evaluation.json',old/'lammps/thermo.npy']
    record={'script_sha256':sha(Path(__file__)),'sources':[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p)} for p in sources],
        'gif':{'frames':101,'source_frame_indices':indices.tolist(),'source_time_fs':(t[indices,1]*1000).tolist(),'interpolated':False},'manuscript_eligible':False}
    (r/'recovery-figure-manifest.json').write_text(json.dumps(record,indent=2)+'\n');print('Rendered recovery comparison and 101 actual 1 ps trajectory frames')
if __name__=='__main__':main()
