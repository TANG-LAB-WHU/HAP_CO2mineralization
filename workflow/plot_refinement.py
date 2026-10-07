"""Render only actual saved calculations, with explicit missing-stage panels and hashes."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation,PillowWriter
from matplotlib.lines import Line2D
from ase.io import read
COLORS={'Ca':'#279d8f','P':'#dd9331','O':'#d34b4b','H':'#9ea8b5','C':'#25364b'}
plt.rcParams.update({'font.size':12,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white'})
def load(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def empty(ax,label):ax.text(.5,.5,'MISSING / NOT COMPLETED\n'+label,ha='center',va='center',transform=ax.transAxes,color='#b42318');ax.set_axis_off()
def layout(r,name,title,axes=(2,2),size=(10,8)):
    fig,axs=plt.subplots(*axes,figsize=size,squeeze=False);fig.suptitle(f'{r.name} | {title}',fontsize=14,x=.03,ha='left');return fig,axs.ravel()
def finish(fig,d,name,sources,top=.94):
    fig.text(.02,.01,'Local diagnostic / pilot evidence only. No production or mineralization claim.',fontsize=8,color='#666666')
    fig.tight_layout(rect=[0,.035,1,top]);fig.savefig(d/f'{name}.png',dpi=170);fig.savefig(d/f'{name}.pdf');plt.close(fig)
    return {'figure':name,'sources':[{'path':str(p),'sha256':sha(p)} for p in sources if p.is_file()]}
def geom(ax,a,axes,title):
    s=np.array(a.get_chemical_symbols());pos=a.positions
    for el,c in COLORS.items():
        sel=s==el;ax.scatter(pos[sel,axes[0]],pos[sel,axes[1]],c=c,s=35 if el!='H' else 15,label=el,edgecolors='white',linewidths=.3,zorder=3)
    for i in range(len(a)):
        for j in range(i+1,len(a)):
            if {s[i],s[j]} in [{'P','O'},{'H','O'},{'C','O'}] and np.linalg.norm(pos[i]-pos[j])< (1.9 if 'P' in [s[i],s[j]] else 1.5):
                ax.plot(pos[[i,j],axes[0]],pos[[i,j],axes[1]],color='#8b949e',lw=.7,zorder=1)
    ax.set(xlabel='xyz'[axes[0]]+' (A)',ylabel='xyz'[axes[1]]+' (A)',title=title);ax.set_aspect('equal',adjustable='datalim')
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--animate',action='store_true');a=p.parse_args();r=a.run;d=r/'figures';d.mkdir(exist_ok=True);manifest=[]
    if (r/'reproducibility-manifest.json').exists():raise RuntimeError('Captured evidence is immutable; copy the run to a separate review directory before re-rendering')
    fig,ax=layout(r,'structure','Structure and chemical identity')
    sources=[]
    for k,label in enumerate(['parallel','antiparallel']):
        path=r/'structures'/label/'slab.traj'
        if path.exists():geom(ax[k],read(path),(0,2),label+' OH: initial slab');sources.append(path)
        else:empty(ax[k],label)
    path=r/'dft-relax/final.traj'
    if not path.exists():path=r/'structures/antiparallel/adsorbed.refined.traj'
    if path.exists():
        atoms=read(path);geom(ax[2],atoms,(0,2),'Selected configuration: side');geom(ax[3],atoms,(0,1),'Selected configuration: top');sources.append(path)
        fig.legend(handles=[Line2D([],[],marker='o',ls='',color=c,label=s) for s,c in COLORS.items()],ncol=5,loc='upper right',bbox_to_anchor=(.98,.94),fontsize=12)
    else:empty(ax[2],'selected geometry');empty(ax[3],'selected geometry')
    manifest.append(finish(fig,d,'01_structure',sources))
    fig,ax=layout(r,'relax','MACE relaxation and DFT refinement');sources=[]
    for variant,style in [('parallel','-'),('antiparallel','--')]:
        for name,color in [('slab','#2674ac'),('adsorbed','#d78125')]:
            path=r/f'structures/{variant}/{name}.relax.json'
            if path.exists():
                data=load(path);step=[x['step'] for x in data];energy=np.array([x['energy_eV'] for x in data]);force=[x['fmax_eV_A'] for x in data]
                ax[0].plot(step,energy-energy[0],style,color=color,label=variant+'/'+name);ax[1].semilogy(step,force,style,color=color,label=variant+'/'+name);sources.append(path)
    path=r/'structures/antiparallel/adsorbed.refine.json'
    if path.exists():
        data=load(path);ax[1].semilogy(np.array([x['step'] for x in data])+400,[x['fmax_eV_A'] for x in data],color='#772f8c',label='LBFGS continuation');sources.append(path)
    ax[0].set(xlabel='Optimizer step',ylabel='E - initial E (eV)',title='Within-method energy descent');ax[0].legend(fontsize=12)
    ax[1].axhline(.03,color='black',ls=':',label='0.03 eV/A gate');ax[1].set(xlabel='Optimizer step',ylabel='Maximum force (eV/A)',title='All atoms mobile');handles,labels=ax[1].get_legend_handles_labels();ax[1].legend(handles[-2:],labels[-2:],fontsize=12)
    paths=sorted(r.glob('dft-relax-attempt*/opt-pos-1.xyz'))
    if (r/'dft-relax/opt-pos-1.xyz').exists():paths.append(r/'dft-relax/opt-pos-1.xyz')
    if paths:
        from aimd_to_deepmd import frames
        from refine_pipeline import forces_from_log,HARTREE
        energy_offset=None;force_offset=0;frame_offset=0
        for path in paths:
            data=frames(path);ener=np.array([x[1] for x in data])*HARTREE
            if energy_offset is None:energy_offset=ener[0]
            ax[2].plot(np.arange(len(ener))+frame_offset,ener-energy_offset,label=path.parent.name)
            fs=forces_from_log(path.parent/'opt.log')
            ax[3].semilogy(np.arange(len(fs))+force_offset,[np.linalg.norm(x,axis=1).max() for x in fs],label=path.parent.name)
            if path.parent.name=='dft-relax' and (r/'optimizer_trial_audit.json').exists():
                trial=load(r/'optimizer_trial_audit.json')
                rejected=[x for x in trial['evaluations'] if not x['matches_initial_or_accepted_saved_energy']]
                if rejected:
                    peak=max(rejected,key=lambda x:x['fmax_eV_A'])
                    ax[3].annotate('Line-search trial\nnot accepted',xy=(force_offset+peak['force_evaluation']-1,peak['fmax_eV_A']),
                        xytext=(.73,.72),textcoords='axes fraction',ha='center',fontsize=10,arrowprops={'arrowstyle':'->','color':'#555555'})
                sources.append(r/'optimizer_trial_audit.json')
            force_offset+=len(fs);frame_offset+=len(ener);sources.extend([path,path.parent/'opt.log'])
        ax[2].set(xlabel='Saved frame (attempts concatenated)',ylabel='E - first saved E (eV)',title='CP2K PBE, force-tested grid');ax[2].legend(fontsize=10)
        ax[3].axhline(.05,color='black',ls=':');ax[3].set(xlabel='Force evaluation (including restart)',ylabel='Maximum force (eV/A)',title='DFT force gate = 0.05 eV/A');ax[3].legend(fontsize=10)
    else:empty(ax[2],'DFT optimization');empty(ax[3],'DFT forces')
    manifest.append(finish(fig,d,'02_relaxation',sources))
    fig,ax=layout(r,'reference','Reference consistency and numerical sensitivity');path=r/'reference_audit.json';sources=[]
    if path.exists():
        data=load(path);names=list(data['settings']);values=[data['settings'][x]['adsorption_energy_eV'] for x in names]
        ax[0].bar(names,values,color=['#3975a5']*3+['#ba7334']);ax[0].axhline(0,color='black',lw=.7);ax[0].set(ylabel='Three-reference difference (eV)',title='MACE geometries; includes\nsubstrate reconstruction');ax[0].tick_params(axis='x',rotation=20)
        delta=np.array(values[:3])-values[2];ax[1].plot([300,500,700],delta,'o-');ax[1].axhspan(-.02,.02,color='#70b798',alpha=.2);ax[1].set(xlabel='Plane-wave cutoff (Ry)',ylabel='Difference from 700 Ry (eV)',title='Relative cutoff: 40 / 60 / 60 Ry')
        for name in ['slab','co2','adsorbed']:ax[2].plot(names,[data['settings'][x][name]['fmax_eV_A'] for x in names],'o-',label=name)
        ax[2].set(ylabel='CP2K maximum force (eV/A)',title='MACE geometry is not\na stationary DFT structure');ax[2].tick_params(axis='x',rotation=20);ax[2].legend()
        ax[3].text(.05,.85,'Cutoff gate: '+str(data['cutoff_gate'])+'\n|500-700| = %.4f eV'%data['cutoff_change_eV']+'\n\nBasis / vacuum / slab thickness:\nNOT converged\nDispersion: a method change\nNo cross-method absolute energies mixed',transform=ax[3].transAxes,va='top');ax[3].set_axis_off();sources=[path]
    else:
        for x in ax:empty(x,'completed CP2K reference audit')
    fragments=[r/'mace_fragments.json',r/'cp2k_fragments.json']
    if all(p.exists() for p in fragments):
        ax[3].clear();ax[3].set_axis_on();keys=['vertical_interaction_eV','slab_deformation_eV','co2_deformation_eV'];x=np.arange(3)
        for file,label,shift,color in zip(fragments,['MACE','CP2K PBE500'],[-.18,.18],['#3975a5','#ba7334']):
            data=load(file);ax[3].bar(x+shift,[data[k] for k in keys],width=.36,label=label,color=color)
        ax[3].set_xticks(x,['Interaction','Substrate','Molecule']);ax[3].tick_params(axis='x',rotation=15);ax[3].axhline(0,color='black',lw=.7)
        ax[3].set(ylabel='Same-method energy term (eV)',title='Decomposition of the total difference');ax[3].legend(fontsize=11);sources.extend(fragments)
    manifest.append(finish(fig,d,'03_cp2k_references',sources))
    fig,ax=layout(r,'scf','SCF behavior and force-grid sensitivity');sources=[]
    import re
    for label,path in [('500/60',r/'references/pbe-500/adsorbed/adsorbed.log'),('900/80',r/'force-grid/900/grid.log'),('1100/100',r/'force-grid/1100/grid.log')]:
        if path.exists():
            rows=[]
            for line in path.read_text(errors='replace').splitlines():
                if re.match(r'^\s*\d+\s+OT DIIS',line):
                    fields=line.split();rows.append([float(fields[5]),float(fields[6])])
            if rows:
                rows=np.array(rows);ax[0].semilogy(np.arange(1,len(rows)+1),rows[:,0],label=label);ax[1].semilogy(np.arange(1,len(rows)+1),np.maximum(abs(rows[:,1]-rows[-1,1])*27.211386245988,1e-12),label=label);sources.append(path)
    ax[0].axhline(1e-8,color='black',ls=':');ax[0].set(xlabel='Cumulative SCF update',ylabel='Reported SCF convergence measure',title='Fine-grid SCF target = 1e-8');ax[0].legend(fontsize=12)
    ax[1].set(xlabel='Cumulative SCF update',ylabel='|E - final SCF E| (eV)',title='Within each fixed-geometry calculation');ax[1].legend(fontsize=12)
    path=r/'force_grid_audit.json'
    if path.exists():
        fine=load(path);refs=load(r/'reference_audit.json')['settings'];delta=np.array(refs['pbe-500']['adsorbed']['force_eV_A'])-np.array(refs['pbe-700']['adsorbed']['force_eV_A']);coarse=float(np.sqrt(np.mean(delta**2)))
        ax[2].bar(['500/60 vs700/60','900/80 vs1100/100'],[coarse,fine['force_component_rmse_eV_A']],color=['#c87534','#3a977b']);ax[2].axhline(.02,color='black',ls=':',label='force-label gate');ax[2].set(ylabel='Force-component RMSE (eV/A)');ax[2].tick_params(axis='x',rotation=10);ax[2].legend(fontsize=12)
        ax[3].text(.05,.9,f"Fine-grid gate: {fine['force_gate']}\nSelected labels: 1100/100 Ry, SCF 1e-8\nMaximum component change: {fine['max_component_difference_eV_A']:.5f} eV/A\n\nOne-geometry sensitivity test.\nBasis / k-point / cell convergence:\nnot established.\nTwo numerical settings do not provide\nstatistical error bars.",transform=ax[3].transAxes,va='top');ax[3].set_axis_off();sources.extend([path,r/'reference_audit.json'])
    else:empty(ax[2],'force-grid audit');empty(ax[3],'force-grid audit')
    manifest.append(finish(fig,d,'03b_scf_force_grid',sources))
    fig,ax=layout(r,'aimd','Disjoint short CP2K replicas',axes=(3,2),size=(10,11));path=r/'aimd_audit.json';sources=[]
    if path.exists():
        for q in load(path)['replicas']:
            t=q['time_fs'];label=f"R{q['replica']} {q['split']}"
            ax[0].plot(t,q['temperature_K'],label=label);ax[1].plot(t,(np.array(q['conserved_eV'])-q['conserved_eV'][0])*1000/47);ax[2].plot(t,q['fmax_eV_A']);ax[3].plot(t,q['rms_displacement_A']);ax[4].plot(t,q['minimum_distance_A']);ax[5].plot(t,q['potential_acf'])
        for x,y,title in zip(ax,['Temperature (K)','Conserved-energy change (meV/atom)','Maximum force (eV/A)','RMS displacement (A)','Minimum pair distance (A)','Normalized potential autocorrelation'],['NVE: not temperature-controlled','Short-run drift diagnostic','Finite force labels','COM translation removed','Contact check','No decorrelation claim']):x.set(xlabel='Time / lag (fs)',ylabel=y,title=title)
        handles,labels=ax[0].get_legend_handles_labels();fig.legend(handles,labels,ncol=4,loc='upper center',bbox_to_anchor=(.52,.95),fontsize=10)
        ax[4].axhline(.7,color='black',ls=':');sources=[path]
    else:
        for x in ax:empty(x,'aligned AIMD labels')
    manifest.append(finish(fig,d,'04_aimd',sources,top=.90))
    fig,ax=layout(r,'splits','Dataset coverage and leakage audit');path=r/'dataset_audit.json';sources=[]
    if path.exists():
        q=load(path);index=q['index'];roles=['train','valid','test'];colors=['#3975a5','#d49230','#429177'];baseline=np.mean([x['energy_eV'] for x in index if x['role']=='train'])
        energy_bins=np.histogram_bin_edges([(x['energy_eV']-baseline)*1000/47 for x in index],bins=8);force_bins=np.histogram_bin_edges([x['fmax_eV_A'] for x in index],bins=8)
        for role,c in zip(roles,colors):
            subset=[x for x in index if x['role']==role];ax[0].hist([(x['energy_eV']-baseline)*1000/47 for x in subset],alpha=.5,label=role,color=c,bins=energy_bins);ax[1].hist([x['fmax_eV_A'] for x in subset],alpha=.5,label=role,color=c,bins=force_bins)
        ax[0].set(xlabel='(E - training mean) / atom (meV)',ylabel='Frames');ax[1].set(xlabel='Maximum force (eV/A)',ylabel='Frames');ax[0].legend();ax[1].legend()
        dist=np.load(r/'dataset/geometry_distance_matrix.npy');im=ax[2].imshow(dist,cmap='viridis');fig.colorbar(im,ax=ax[2],label='Pair-distance-vector RMS difference (A)');ax[2].set(xlabel='Frame index (train / valid / test)',ylabel='Frame index')
        for boundary in np.cumsum([q['frames_per_split'][x] for x in roles])[:-1]-.5:
            ax[2].axhline(boundary,color='white',lw=.8);ax[2].axvline(boundary,color='white',lw=.8)
        counts=' / '.join(str(q['frames_per_split'][x]) for x in ['train','valid','test'])
        ax[3].text(.05,.9,f"Train / valid / test: {counts}\nDisjoint replicas: {q['replica_disjoint']}\nCross-split minimum: {q['cross_split_min_pair_distance_rms_A']:.4f} A\nPredefined threshold: {q['threshold_A']:.3f} A\n\nWithin-replica frames remain correlated.\nSmall local basin; no broad transfer claim.\nNo virial labels or pressure accuracy claim.",va='top',transform=ax[3].transAxes);ax[3].set_axis_off();sources=[path,r/'dataset/geometry_distance_matrix.npy']
    else:
        for x in ax:empty(x,'dataset split audit')
    manifest.append(finish(fig,d,'05_dataset',sources))
    fig,ax=layout(r,'models','Model convergence and held-out prediction',axes=(3,2),size=(10,11));sources=[]
    for name in ['overfit','compact','wider']:
        path=r/'models'/name/'lcurve.out'
        if path.exists():
            q=np.genfromtxt(path,names=True);names=q.dtype.names
            for key,style in [('rmse_f_trn','-'),('rmse_f_val','--')]:
                if key in names:ax[0].semilogy(q[names[0]],q[key],style,label=name+('/train' if key.endswith('trn') else '/validation'))
            sources.append(path)
    ax[0].axhline(.1,color='black',ls=':');ax[0].set(xlabel='Training update',ylabel='Force RMSE (eV/A)',title='Overfit is a separate implementation check')
    if sources:ax[0].legend(fontsize=9)
    else:empty(ax[0],'DeepMD learning curves')
    path=r/'model_evaluation.json'
    if path.exists():
        q=load(path);selected=q['selected'];pred=r/'models'/selected/'test_predictions.npz';z=np.load(pred);ef=(z['predicted_energy_eV']-z['reference_energy_eV'])/47;ff=z['predicted_force_eV_A']-z['reference_force_eV_A'];truth=z['reference_force_eV_A'].ravel();guess=z['predicted_force_eV_A'].ravel()
        ax[1].scatter(truth,guess,s=5,alpha=.35);bounds=[min(truth.min(),guess.min()),max(truth.max(),guess.max())];ax[1].plot(bounds,bounds,'k:',lw=1);ax[1].set(xlabel='CP2K force component (eV/A)',ylabel='DeepMD force component (eV/A)',title=selected+': held-out test')
        ax[2].hist(ff.ravel(),bins=35,color='#3975a5');ax[2].set(xlabel='Force-component error (eV/A)',ylabel='Components')
        centered=(z['reference_energy_eV']-z['reference_energy_eV'].mean())*1000/47;predcenter=(z['predicted_energy_eV']-z['reference_energy_eV'].mean())*1000/47
        ax[3].scatter(centered,predcenter);bounds=[min(centered.min(),predcenter.min()),max(centered.max(),predcenter.max())];ax[3].plot(bounds,bounds,'k:');ax[3].set(xlabel='CP2K energy, centered (meV/atom)',ylabel='DeepMD energy, same center (meV/atom)')
        ax[4].hist(ef*1000,bins=8,color='#d49230');ax[4].axvline(10,color='black',ls=':');ax[4].axvline(-10,color='black',ls=':');ax[4].set(xlabel='Energy error (meV/atom)',ylabel='Frames')
        metrics=q['models'][selected]['test'];ax[5].text(.03,.9,f"Selected on validation: {selected}\nTest force RMSE: {metrics['force_rmse_eV_A']:.4f} eV/A\nTest energy RMSE: {1000*metrics['energy_rmse_eV_atom']:.3f} meV/atom\nForce gate (<0.1): {q['heldout_force_gate']}\nEnergy gate (<10 meV): {q['heldout_energy_gate']}\nSymmetry gate: {q['symmetry_gate']}\n\nNo independent reaction or\ntemperature transfer claim.",transform=ax[5].transAxes,va='top');ax[5].set_axis_off();sources.extend([path,pred])
    else:
        for x in ax[1:]:empty(x,'independent model evaluation')
    manifest.append(finish(fig,d,'06_models',sources))
    fig,ax=layout(r,'lammps','LAMMPS trajectory evolution and transfer');path=r/'lammps_audit.json';sources=[]
    if path.exists():
        q=load(path);t=np.load(r/'lammps/thermo.npy');ax[0].plot(t[:,1]*1000,t[:,2]);ax[0].set(xlabel='Time (fs)',ylabel='Temperature (K)',title='NVE; 200 K initial velocities')
        ax[1].plot(t[:,1]*1000,(t[:,5]-t[0,5])*1000/47,label='dt=0.25 fs');half=np.load(r/'lammps/halfstep-thermo.npy');ax[1].plot(half[:,1]*1000,(half[:,5]-half[0,5])*1000/47,label='dt=0.125 fs');ax[1].legend(fontsize=12);ax[1].set(xlabel='Time (fs)',ylabel='Total-energy change (meV/atom)',title=f"0.25 fs drift = {q['energy_drift_meV_atom_ps']:.3f} meV/atom/ps")
        frames=q['frames'];times=[x['time_fs'] for x in frames];ax[2].plot(times,[x['minimum_distance_A'] for x in frames],label='minimum contact');ax[2].plot(times,[x['rms_displacement_A'] for x in frames],label='RMS displacement');ax[2].set(xlabel='Time (fs)',ylabel='Distance (A)');ax[2].legend()
        ax[3].text(.03,.9,f"Engine force mismatch: {q['engine_force_max_difference_eV_A']:.2e} eV/A\nEngine energy mismatch: {q['engine_energy_difference_eV']:.2e} eV\nRestart position difference: {q['restart_position_max_difference_A']:.2e} A\nRestart velocity difference: {q['restart_velocity_max_difference_A_per_ps']:.2e} A/ps\nMapping gate: {q['mapping_gate']}\nRestart gate: {q['restart_gate']}\nGeometry gate: {q['chemistry_gate']}\n\n250 fs does not establish long-term stability.\nInterpret together with held-out model errors.",transform=ax[3].transAxes,va='top');ax[3].set_axis_off();sources=[path,r/'lammps/thermo.npy',r/'lammps/halfstep-thermo.npy']
    else:
        for x in ax:empty(x,'completed LAMMPS diagnostics')
    manifest.append(finish(fig,d,'07_lammps',sources))
    if a.animate:
        path=r/'lammps/trajectory.extxyz'
        if not path.exists():raise FileNotFoundError('No real trajectory available; refusing animation substitute')
        frames=read(path,index=':');chosen=frames[::2];fig,axes=plt.subplots(1,2,figsize=(10,6));fig.subplots_adjust(top=.82,bottom=.14);colors=[COLORS[s] for s in frames[0].get_chemical_symbols()];scat=[]
        allpos=np.concatenate([x.positions for x in chosen]);lo=allpos.min(axis=0)-1;hi=allpos.max(axis=0)+1
        for axis,view in zip(axes,[(0,2),(0,1)]):
            scatter=axis.scatter([],[],s=45,edgecolors='white',linewidths=.3);scat.append(scatter);axis.set(xlim=(lo[view[0]],hi[view[0]]),ylim=(lo[view[1]],hi[view[1]]),xlabel='xyz'[view[0]]+' (A)',ylabel='xyz'[view[1]]+' (A)');axis.set_aspect('equal')
        title=fig.suptitle('');fig.legend(handles=[Line2D([],[],marker='o',ls='',color=c,label=s) for s,c in COLORS.items()],ncol=5,loc='upper center',bbox_to_anchor=(.5,.91),fontsize=11)
        fig.text(.02,.03,'Actual saved frames; no interpolation. FAILED force/geometry validation; not mineralization evidence.',fontsize=9,color='#a52323')
        def update(i):
            for scatter,view in zip(scat,[(0,2),(0,1)]):scatter.set_offsets(chosen[i].positions[:,view]);scatter.set_color(colors)
            title.set_text(f'{r.name} | LAMMPS NVE | t = {i*5:.1f} fs | dt = 0.25 fs\nDry HAP + CO2; all atoms mobile; atom IDs preserved')
            return scat+[title]
        animation=FuncAnimation(fig,update,frames=len(chosen),interval=120,blit=False);animation.save(d/'08_lammps_actual.gif',writer=PillowWriter(fps=8));plt.close(fig)
        manifest.append({'figure':'08_lammps_actual.gif','frames':len(chosen),'saved_frame_interval_fs':5.,'source':str(path),'source_sha256':sha(path),'interpolation':False})
    (d/'manifest.json').write_text(json.dumps({'run':r.name,'script_sha256':sha(Path(__file__)),'figures':manifest,'manuscript_eligible':False},indent=2),encoding='utf-8')
if __name__=='__main__':main()
