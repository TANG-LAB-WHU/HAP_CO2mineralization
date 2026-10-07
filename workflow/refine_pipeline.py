"""Bounded, provenance-preserving local HAP model validation; never production evidence."""
from __future__ import annotations
import argparse, hashlib, json, re, shutil, sys, time
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.geometry import find_mic

ROOT=Path(__file__).resolve().parents[1]
HARTREE=27.211386245988
BOHR=0.529177210903
TYPES=['Ca','P','O','H','C']

def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,obj): Path(p).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def load(p): return json.loads(Path(p).read_text(encoding='utf-8'))

def chemistry(a):
    s=np.array(a.get_chemical_symbols());d=a.get_all_distances(mic=True);np.fill_diagonal(d,np.inf)
    po=[int(np.sum(d[i,s=='O']<1.9)) for i in np.where(s=='P')[0]]
    oh=[float(np.min(d[i,s=='O'])) for i in np.where(s=='H')[0]]
    charges=np.array([{'Ca':2,'P':5,'O':-2,'H':1,'C':4}[x] for x in s])
    carbon=[]
    for i in np.where(s=='C')[0]:
        nearest=sorted(np.where(s=='O')[0],key=lambda j:d[i,j])[:4]
        carbon.append({'carbon_id':int(i+1),'oxygen_ids':[int(j+1) for j in nearest],'distances_A':[float(d[i,j]) for j in nearest],'oxygen_count_within_1p6A':int(np.sum(d[i,s=='O']<1.6))})
    return {'carbon_oxygen_contacts':carbon,'formula':a.get_chemical_formula(),'natoms':len(a),'minimum_distance_A':float(d.min()),
        'P_O_coordination_under_1p9A':po,'O_H_distances_A':oh,'formal_charge_e':int(charges.sum()),
        'formal_charge_dipole_z_eA':float(charges@a.positions[:,2]),'cell_A':a.cell.array.tolist(),
        'pbc':a.pbc.tolist(),'fixed_atom_count':sum(len(c.index) for c in a.constraints),
        'geometry_gate':bool(charges.sum()==0 and d.min()>.7 and all(x==4 for x in po) and all(.8<x<1.2 for x in oh))}

def validate_sampling_profile(cfg):
    # The preassigned eight-start design and downstream time labels are deliberate.
    # Reject unsupported edits instead of silently simulating a different plan.
    sampling=cfg['sampling']
    if sampling['replicas']!=8 or sampling['timestep_fs']!=0.25:
        raise ValueError('This bounded design supports exactly 8 replicas at 0.25 fs; implement and review a new design before changing these fields')
    if sampling['steps_per_replica']<1 or sampling['frame_stride']<1:
        raise ValueError('Sampling step count and frame stride must be positive')
    if not 1<=sampling['discard_initial_steps']<=sampling['steps_per_replica']:
        raise ValueError('Initial-step exclusion must leave sampled frames')

def prepare(r,cfg):
    from ase.build import surface,molecule
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor
    source=ROOT/'3.ML_AIMD/Step4_lammps_scaleup/HAP_ConventionalCell.cif'
    s=Structure.from_file(source); records={}
    for label,allowed in [('parallel',{'O':[.1978,.6978],'H':[.0608,.5608]}),('antiparallel',{'O':[.1978,.8022],'H':[.0608,.9392]})]:
        dest=r/'structures'/label;dest.mkdir(parents=True)
        keep=[x for x in s if x.is_ordered or any(abs(x.frac_coords[2]-z)<1e-4 for z in allowed[next(iter(x.species)).symbol])]
        bulk=AseAtomsAdaptor.get_atoms(Structure(s.lattice,[next(iter(x.species)).symbol for x in keep],[x.frac_coords for x in keep]))
        assert len(bulk)==44 and bulk.get_chemical_formula()=='H2Ca10O26P6'
        slab=surface(bulk,(0,0,1),1,vacuum=12,periodic=True)
        gas=molecule('CO2');gas.rotate('z','x');gas.set_cell(slab.cell);gas.pbc=True;gas.center()
        top=max((i for i,a in enumerate(slab) if a.symbol=='Ca'),key=lambda i:slab.positions[i,2])
        ads=gas.copy();ads.translate(slab.positions[top]+[0,0,2.8]-ads.positions[0]);combined=slab+ads
        records[label]={'OH_ordering':allowed,'source_sha256':digest(source),'systems':{}}
        for name,a in [('bulk',bulk),('slab',slab),('co2',gas),('adsorbed',combined)]:
            a.set_constraint();write(dest/f'{name}.extxyz',a);write(dest/f'{name}.traj',a)
            records[label]['systems'][name]=chemistry(a)
            assert chemistry(a)['geometry_gate']
        save(dest/'groups.json',{'slab_ids':list(range(1,45)),'co2_ids':[45,46,47],'co2_oxygen_ids':[46,47],
             'constraint':'all atoms mobile in MACE, CP2K, AIMD and LAMMPS; no implicit XYZ constraints'})
    save(r/'structure_audit.json',{'variants':records,'selected':'antiparallel',
        'rationale':'Initial neutral slab has antiparallel OH and zero formal-charge z dipole; compare energetics without claiming stable ordering.',
        'limitations':['one-cell dry (001) slab','formal ionic dipole is not electronic polarization','no hydration or termination/thickness convergence'],
        'baseline_correction':'Baseline extxyz dropped FixAtoms; prior actual relaxations were unconstrained despite documentation.', 'manuscript_eligible':False})

def mace(r,cfg):
    import torch
    from mace.calculators import MACECalculator
    from ase.optimize import FIRE
    assert torch.cuda.is_available();torch.set_num_threads(4);torch.cuda.set_per_process_memory_fraction(.35)
    model=ROOT/cfg['mace']['model'];calc=MACECalculator(model_paths=str(model),device='cuda',default_dtype='float64')
    allresults={}
    for variant in ['parallel','antiparallel']:
        d=r/'structures'/variant;results={}
        for name in ['slab','co2','adsorbed']:
            a=read(d/f'{name}.traj');assert not a.constraints;a.calc=calc;trace=[]
            opt=FIRE(a,logfile=str(d/f'{name}.opt.log'),trajectory=str(d/f'{name}.relax.traj'),maxstep=.035)
            def record(): trace.append({'step':opt.nsteps,'energy_eV':float(a.get_potential_energy()),'fmax_eV_A':float(np.linalg.norm(a.get_forces(),axis=1).max())})
            opt.attach(record);ok=opt.run(fmax=cfg['mace']['fmax_eV_A'],steps=cfg['mace']['max_steps'])
            write(d/f'{name}.relaxed.extxyz',a);write(d/f'{name}.relaxed.traj',a)
            save(d/f'{name}.relax.json',trace)
            results[name]={'energy_eV':float(a.get_potential_energy()),'converged':bool(ok),'steps':opt.nsteps,
                'fmax_eV_A':float(np.linalg.norm(a.get_forces(),axis=1).max()),'chemistry':chemistry(a)}
        results['adsorption_energy_eV']=results['adsorbed']['energy_eV']-results['slab']['energy_eV']-results['co2']['energy_eV']
        allresults[variant]=results;save(r/'mace.json',{'variants':allresults,'model_sha256':digest(model),'peak_cuda_bytes':torch.cuda.max_memory_allocated(),'manuscript_eligible':False})
    assert all(x[n]['converged'] for x in allresults.values() for n in ['slab','co2','adsorbed']), 'Relaxation gate failed'

def cp_input(a,d,name,cutoff=500,relcut=60,eps='1e-7',dispersion=False,steps=0,temp=200,seed=1,wfn=None):
    import local_compute
    d.mkdir(parents=True,exist_ok=True);write(d/f'{name}.relaxed.extxyz',a)
    local_compute.cp2k_input(d,name);s=(d/f'{name}.inp').read_text()
    s=s.replace('CUTOFF 300',f'CUTOFF {cutoff}').replace('REL_CUTOFF 40',f'REL_CUTOFF {relcut}').replace('1.0E-6',eps)
    if dispersion:
        s=s.replace('    &END XC\n','''      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3
          PARAMETER_FILE_NAME /opt/cp2k/data/dftd3.dat
          REFERENCE_FUNCTIONAL PBE
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
''')
    if wfn:
        shutil.copyfile(wfn,d/'initial.wfn');s=s.replace('  &DFT','  &DFT\n    WFN_RESTART_FILE_NAME initial.wfn',1).replace('    &SCF','    &SCF\n      SCF_GUESS RESTART',1)
    if steps:
        s=s.replace('RUN_TYPE ENERGY_FORCE',f'RUN_TYPE MD\n  SEED {seed}')
        s+=f'''\n&MOTION
  &MD
    ENSEMBLE NVE
    STEPS {steps}
    TIMESTEP 0.25
    TEMPERATURE {temp}
  &END MD
  &PRINT
    &TRAJECTORY
      &EACH
        MD 1
      &END EACH
    &END TRAJECTORY
    &FORCES
      &EACH
        MD 1
      &END EACH
    &END FORCES
    &CELL
      &EACH
        MD 1
      &END EACH
    &END CELL
    &RESTART
      &EACH
        MD 4
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION
'''
    (d/f'{name}.inp').write_text(s,encoding='utf-8')

def cp_prepare(r,cfg):
    base=r/'structures'/'antiparallel';jobs=[]
    for setting in cfg['cp2k']['reference_settings']:
        for name in ['slab','co2','adsorbed']:
            d=r/'references'/setting['id']/name;a=read(base/(f'{name}.refined.traj' if (base/f'{name}.refined.traj').exists() else f'{name}.relaxed.traj'))
            cp_input(a,d,name,setting['cutoff'],setting['relcut'],setting['eps'],setting['dispersion'])
            jobs.append({'path':str(d.relative_to(r)),'input':name+'.inp','output':name+'.log','setting':setting})
    save(r/'reference_jobs.json',jobs)

def repair_d3(r,cfg):
    setting=next(x for x in cfg['cp2k']['reference_settings'] if x['id']=='pbe-d3-500')
    base=r/'structures'/'antiparallel'
    for name in ['slab','co2','adsorbed']:
        d=r/'references'/setting['id']/name
        if d.exists():raise FileExistsError('Archive the original failed attempt first')
        a=read(base/(f'{name}.refined.traj' if (base/f'{name}.refined.traj').exists() else f'{name}.relaxed.traj'))
        cp_input(a,d,name,setting['cutoff'],setting['relcut'],setting['eps'],setting['dispersion'])

def mace_continue(r,cfg):
    import torch
    from mace.calculators import MACECalculator
    from ase.optimize import LBFGS
    torch.set_num_threads(4);torch.cuda.set_per_process_memory_fraction(.35)
    d=r/'structures'/'antiparallel';a=read(d/'adsorbed.relaxed.traj');a.calc=MACECalculator(model_paths=str(ROOT/cfg['mace']['model']),device='cuda',default_dtype='float64')
    opt=LBFGS(a,logfile=str(d/'adsorbed.refine.log'),trajectory=str(d/'adsorbed.refine.traj'),maxstep=.04);trace=[]
    opt.attach(lambda:trace.append({'step':opt.nsteps,'energy_eV':float(a.get_potential_energy()),'fmax_eV_A':float(np.linalg.norm(a.get_forces(),axis=1).max())}))
    ok=opt.run(fmax=.03,steps=300);write(d/'adsorbed.refined.traj',a);write(d/'adsorbed.refined.extxyz',a);save(d/'adsorbed.refine.json',trace)
    save(r/'mace_refined.json',{'converged':bool(ok),'energy_eV':float(a.get_potential_energy()),'fmax_eV_A':trace[-1]['fmax_eV_A'],'chemistry':chemistry(a),'steps':opt.nsteps,'manuscript_eligible':False})
    assert ok and chemistry(a)['geometry_gate'],'Refined geometry failed'

def forces_from_log(path):
    text=Path(path).read_text(errors='replace');blocks=text.split('FORCES| Atomic forces [hartree/bohr]')[1:];out=[]
    for block in blocks:
        rows=re.findall(r'^\s*FORCES\|\s+\d+\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)',block.split('SUM OF ATOMIC FORCES')[0],re.M)
        if rows:out.append(np.array(rows,float)*HARTREE/BOHR)
    return out

def references_audit(r,cfg):
    from cp2k_validate import parse
    records={}
    for setting in cfg['cp2k']['reference_settings']:
        results={}
        for name in ['slab','co2','adsorbed']:
            d=r/'references'/setting['id']/name;results[name]=parse(d/f'{name}.log')
            f=forces_from_log(d/f'{name}.log')[-1];results[name]['fmax_eV_A']=float(np.linalg.norm(f,axis=1).max())
            results[name]['force_eV_A']=f.tolist()
        results['adsorption_energy_eV']=results['adsorbed']['energy_eV']-results['slab']['energy_eV']-results['co2']['energy_eV'];records[setting['id']]=results
    change=abs(records['pbe-700']['adsorption_energy_eV']-records['pbe-500']['adsorption_energy_eV'])
    save(r/'reference_audit.json',{'settings':records,'cutoff_change_eV':change,'cutoff_gate':change<cfg['gates']['cutoff_adsorption_change_eV'],
       'interpretation':'Single points on MACE geometries; different relaxed substrate basins can contaminate total adsorption differences. D3 comparison is a method change, not proof of accuracy.','manuscript_eligible':False})

def dft_prepare(r,cfg):
    # Energy convergence alone does not certify force labels. Run a separate grid gate.
    from run_refinement import cp_jobs
    from cp2k_validate import parse
    a=read(r/'structures/antiparallel/adsorbed.refined.traj');jobs=[]
    for cutoff,relcut in [(900,80),(1100,100)]:
        fd=r/'force-grid'/str(cutoff)
        cp_input(a,fd,'grid',cutoff=cutoff,relcut=relcut,eps='1e-8',wfn=r/'references/pbe-500/adsorbed/adsorbed-RESTART.wfn')
        jobs.append({'path':str(fd.relative_to(r)),'input':'grid.inp','output':'grid.log'})
    save(r/'force_grid_jobs.json',jobs);cp_jobs(r,jobs)
    f1=forces_from_log(r/'force-grid/900/grid.log')[-1];f2=forces_from_log(r/'force-grid/1100/grid.log')[-1];delta=f1-f2
    audit={'force_component_rmse_eV_A':float(np.sqrt(np.mean(delta**2))),'max_component_difference_eV_A':float(abs(delta).max()),
        'tested_settings':[[900,80],[1100,100]],'force_gate':bool(np.sqrt(np.mean(delta**2))<.02),'selected':{'cutoff':1100,'relcut':100,'eps':'1e-8'},
        'scope':'One adsorbed geometry force-grid sensitivity; not global numerical or basis convergence','manuscript_eligible':False}
    save(r/'force_grid_audit.json',audit)
    assert audit['force_gate'],'Force-grid convergence gate failed; preserve results and refine before labels'
    d=r/'dft-relax';a=read(r/'structures/antiparallel/adsorbed.refined.traj')
    cp_input(a,d,'opt',cutoff=1100,relcut=100,eps='1e-8',wfn=r/'force-grid/1100/grid-RESTART.wfn')
    p=d/'opt.inp';s=p.read_text().replace('RUN_TYPE ENERGY_FORCE','RUN_TYPE GEO_OPT')
    s+='''\n&MOTION
  &GEO_OPT
    OPTIMIZER BFGS
    MAX_ITER 60
    MAX_FORCE 0.000972345
    RMS_FORCE 0.00048617
    MAX_DR 0.003
    RMS_DR 0.0015
  &END GEO_OPT
  &PRINT
    &TRAJECTORY
      &EACH
        GEO_OPT 1
      &END EACH
    &END TRAJECTORY
    &FORCES
      &EACH
        GEO_OPT 1
      &END EACH
    &END FORCES
  &END PRINT
&END MOTION
''';p.write_text(s)
    save(r/'dft_jobs.json',[{'path':'dft-relax','input':'opt.inp','output':'opt.log'}])

def dft_resume(r,cfg):
    """Prepare a same-method restart after a preserved wall-time failure only."""
    old=r/'dft-relax-attempt1';d=r/'dft-relax'
    assert not d.exists(), 'Archive the failed folder before preparing its continuation'
    job=load(old/'job.json');execution=load(old/'console.execution.json')
    assert job['status']=='failed' and 'TimeoutExpired' in execution.get('error','')
    restart=old/'opt-1.restart';s=restart.read_text()
    completed=int(re.search(r'STEP_START_VAL\s+(\d+)',s)[1])
    assert 0<completed<60, 'The original 60-step bound must remain meaningful'
    assert '&BFGS' not in s
    s=s.replace('   &END GEO_OPT','''     &BFGS
       RESTART_HESSIAN T
       RESTART_FILE_NAME initial.hessian
       TRUST_RADIUS 0.05
     &END BFGS
   &END GEO_OPT''',1)
    a=read(restart,format='cp2k-restart');assert chemistry(a)['geometry_gate']
    d.mkdir();(d/'opt.inp').write_text(s,encoding='utf-8')
    write(d/'opt.relaxed.extxyz',a)
    shutil.copyfile(old/'opt-RESTART.wfn',d/'initial.wfn')
    shutil.copyfile(old/'opt-BFGS.Hessian',d/'initial.hessian')
    save(r/'dft_continuation.json',{'cause':'wall-time limit, not accepted as convergence','archived_attempt':str(old.relative_to(r)),
        'completed_step_counter':completed,'total_step_bound':60,'optimizer':'BFGS with saved Hessian','trust_radius_A':.05,'previous_trust_radius_A':.25,'criteria_unchanged':True,
        'electronic_method_unchanged':True,'restart_sha256':digest(restart),'wfn_sha256':digest(d/'initial.wfn'),
        'hessian_sha256':digest(d/'initial.hessian'),'manuscript_eligible':False})
    save(r/'dft_jobs.json',[{'path':'dft-relax','input':'opt.inp','output':'opt.log'}])

def dft_extend(r,cfg):
    """One bounded same-method extension after a recorded 60-step failure."""
    old=r/'dft-relax-attempt2';d=r/'dft-relax'
    assert not d.exists(), 'Archive the stopped attempt first'
    assert load(old/'job.json')['status']=='completed'
    log=(old/'opt.log').read_text()
    assert 'MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED' in log and 'GEOMETRY OPTIMIZATION COMPLETED' not in log
    audit=load(old/'audit-result.json');assert audit['chemistry']['geometry_gate']
    restart=old/'opt-1.restart';s=restart.read_text()
    counter=int(re.search(r'STEP_START_VAL\s+(\d+)',s)[1]);assert counter==60
    assert re.search(r'MAX_ITER\s+60\b',s)
    s=re.sub(r'MAX_ITER\s+60\b','MAX_ITER 90',s,count=1)
    s=re.sub(r'RESTART_FILE_NAME\s+initial.hessian', 'RESTART_FILE_NAME initial.hessian',s)
    assert 'RESTART_HESSIAN' in s and 'initial.hessian' in s
    a=read(restart,format='cp2k-restart');assert chemistry(a)['geometry_gate']
    d.mkdir();(d/'opt.inp').write_text(s,encoding='utf-8');write(d/'opt.relaxed.extxyz',a)
    shutil.copyfile(old/'opt-RESTART.wfn',d/'initial.wfn');shutil.copyfile(old/'opt-BFGS.Hessian',d/'initial.hessian')
    save(r/'dft_extension.json',{'cause':'60-step optimization limit reached without convergence; retained as failed gate',
        'archived_attempt':str(old.relative_to(r)),'completed_step_counter':counter,'new_total_bound':90,'additional_steps':30,
        'same_force_and_displacement_criteria':True,'same_electronic_method':True,'same_0p05A_trust_radius':True,
        'previous_final_fmax_eV_A':audit['fmax_eV_A'],'restart_sha256':digest(restart),'wfn_sha256':digest(d/'initial.wfn'),
        'hessian_sha256':digest(d/'initial.hessian'),'manuscript_eligible':False})
    save(r/'dft_jobs.json',[{'path':'dft-relax','input':'opt.inp','output':'opt.log'}])

def dft_diagnose(r,cfg):
    import struct,datetime
    d=r/'dft-relax';assert load(d/'job.json')['status']=='completed'
    aa=read(d/'opt-pos-1.xyz',index='-2:');base=read(d/'opt.relaxed.extxyz');step=aa[-1].positions-aa[-2].positions
    raw=(d/'opt-BFGS.Hessian').read_bytes();offset=0;rows=[]
    while offset<len(raw):
        size=struct.unpack_from('<i',raw,offset)[0];offset+=4
        assert size>0 and size%8==0
        rows.append(np.frombuffer(raw,dtype='<f8',count=size//8,offset=offset));offset+=size
        assert struct.unpack_from('<i',raw,offset)[0]==size;offset+=4
    h=np.array(rows);assert h.shape==(len(base)*3,len(base)*3)
    w,v=np.linalg.eigh((h+h.T)/2);basis=np.zeros((len(base)*3,3))
    for i in range(3):basis[i::3,i]=1/np.sqrt(len(base))
    log=(d/'opt.log').read_text();displacements=np.linalg.norm(step,axis=1);fs=forces_from_log(d/'opt.log')
    save(r/'optimizer_diagnosis-step90.json',{'recorded_at_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'result':load(r/'dft_relaxation.json'),'last_step_max_atom_displacement_A':float(displacements.max()),
        'mean_translation_A':step.mean(axis=0).tolist(),'mass_COM_translation_A':np.average(step,axis=0,weights=base.get_masses()).tolist(),
        'largest_displacement_atoms':[{'id':int(i+1),'element':base[i].symbol,'displacement_A':float(displacements[i])} for i in np.argsort(displacements)[-6:][::-1]],
        'smallest_Hessian_eigenvalues_CP2K_units':w[:8].tolist(),'largest_Hessian_eigenvalue_CP2K_units':float(w[-1]),
        'Hessian_positive':bool(w[0]>0),'Hessian_condition_number':float(w[-1]/w[0]),
        'soft_mode_rigid_translation_overlaps':np.sum((basis.T@v[:,:8])**2,axis=0).tolist(),
        'force_history_eV_A':[float(np.linalg.norm(x,axis=1).max()) for x in fs],
        'last_reported_criteria':[x.strip() for x in log.splitlines() if 'is converged' in x][-4:],
        'interpretation':'The force gate is satisfied, while BFGS displacement criteria fail. Positive soft internal modes and small COM drift do not establish a stalled or singular optimizer. A fresh limited-memory history is a bounded algorithm check, not a relaxed criterion or proof of global equilibrium.',
        'sources':{p.name:digest(p) for p in [d/'opt.log',d/'opt-pos-1.xyz',d/'opt-1.restart',d/'opt-BFGS.Hessian']},'manuscript_eligible':False})

def dft_lbfgs_prepare(r,cfg):
    old=r/'dft-relax-attempt3';d=r/'dft-relax';assert not d.exists()
    previous=load(old/'audit-result.json');assert previous['force_gate'] and previous['chemistry']['geometry_gate'] and not previous['optimizer_converged']
    assert load(old/'job.json')['status']=='completed'
    s=(old/'opt-1.restart').read_text();assert int(re.search(r'STEP_START_VAL\s+(\d+)',s)[1])==90
    s=re.sub(r'OPTIMIZER\s+BFGS\b','OPTIMIZER LBFGS',s,count=1)
    s=re.sub(r'MAX_ITER\s+90\b','MAX_ITER 110',s,count=1)
    s=re.sub(r'\n\s*&BFGS\b.*?&END BFGS','',s,count=1,flags=re.S)
    s=s.replace('   &END GEO_OPT','''     &LBFGS
       MAX_H_RANK 5
       MAX_F_PER_ITER 8
     &END LBFGS
   &END GEO_OPT''',1)
    assert '&BFGS' not in s and '&LBFGS' in s
    a=read(old/'opt-1.restart',format='cp2k-restart');assert chemistry(a)['geometry_gate']
    d.mkdir();(d/'opt.inp').write_text(s);write(d/'opt.relaxed.extxyz',a);shutil.copyfile(old/'opt-RESTART.wfn',d/'initial.wfn')
    save(r/'dft_lbfgs_recovery.json',{'previous_BFGS_steps':90,'new_total_step_bound':110,'additional_LBFGS_steps':20,
        'fresh_optimizer_history':True,'LBFGS_rank':5,'max_force_evaluations_per_iteration':8,
        'LBFGS_trust_radius':'CP2K default negative value: no coordinate bounds; line search controls updates',
        'same_force_displacement_and_electronic_criteria':True,'constraints':'all atoms mobile; no constraints added',
        'starting_fmax_eV_A':previous['fmax_eV_A'],'diagnosis':'optimizer_diagnosis-step90.json',
        'restart_sha256':digest(old/'opt-1.restart'),'wfn_sha256':digest(d/'initial.wfn'),'manuscript_eligible':False})
    save(r/'dft_jobs.json',[{'path':'dft-relax','input':'opt.inp','output':'opt.log'}])

def dft_audit(r,cfg):
    from cp2k_validate import parse
    d=r/'dft-relax';q=parse(d/'opt.log');a=read(d/'opt-pos-1.xyz',index=-1);base=read(d/'opt.relaxed.extxyz');a.cell=base.cell;a.pbc=True
    from aimd_to_deepmd import frames
    xyz=frames(d/'opt-pos-1.xyz')[-1]
    assert xyz[2]==base.get_chemical_symbols() and abs(xyz[1]*HARTREE-q['energy_eV'])<2e-6,'Final optimizer positions/energy mismatch'
    f=forces_from_log(d/'opt.log')[-1];assert f.shape==(len(a),3),'Incomplete final force record'
    q.update(fmax_eV_A=float(np.linalg.norm(f,axis=1).max()),chemistry=chemistry(a),optimizer_converged='GEOMETRY OPTIMIZATION COMPLETED' in (d/'opt.log').read_text())
    q['force_gate']=q['fmax_eV_A']<.05;q['manuscript_eligible']=False
    write(d/'final.traj',a);write(d/'final.extxyz',a);save(r/'dft_relaxation.json',q)
    from audit_cp2k_restarts import audit
    save(r/'cp2k_restart_audit.json',audit(r))
    assert q['chemistry']['geometry_gate'],'DFT geometry chemistry failed'
    assert q['force_gate'] and q['optimizer_converged'],'DFT optimization gate failed; do not accept this sampling origin'

def fragments_prepare(r,cfg):
    # Fixed-geometry fragments separate interaction from deformation and reconstruction.
    a=read(r/'structures/antiparallel/adsorbed.refined.traj');jobs=[]
    for name,x in [('slab',a[:44]),('co2',a[44:])]:
        d=r/'fragments'/name;cp_input(x,d,name);jobs.append({'path':str(d.relative_to(r)),'input':name+'.inp','output':name+'.log'})
    save(r/'fragment_jobs.json',jobs)

def fragment_mace(r,cfg):
    import torch
    from mace.calculators import MACECalculator
    torch.set_num_threads(4);torch.cuda.set_per_process_memory_fraction(.35)
    calc=MACECalculator(model_paths=str(ROOT/cfg['mace']['model']),device='cuda',default_dtype='float64');records={}
    a=read(r/'structures/antiparallel/adsorbed.refined.traj')
    for name,x in [('adsorbed',a),('slab_fragment',a[:44]),('co2_fragment',a[44:])]:x.calc=calc;records[name]=float(x.get_potential_energy())
    refs=load(r/'mace.json')['variants']['antiparallel'];records['vertical_interaction_eV']=records['adsorbed']-records['slab_fragment']-records['co2_fragment']
    records['slab_deformation_eV']=records['slab_fragment']-refs['slab']['energy_eV'];records['co2_deformation_eV']=records['co2_fragment']-refs['co2']['energy_eV']
    records['total_reference_difference_eV']=records['vertical_interaction_eV']+records['slab_deformation_eV']+records['co2_deformation_eV'];save(r/'mace_fragments.json',records)

def fragments_audit(r,cfg):
    from cp2k_validate import parse
    refs=load(r/'reference_audit.json')['settings']['pbe-500'];slab=parse(r/'fragments/slab/slab.log');gas=parse(r/'fragments/co2/co2.log')
    e=refs['adsorbed']['energy_eV']-slab['energy_eV']-gas['energy_eV'];deform=slab['energy_eV']-refs['slab']['energy_eV'];gasdeform=gas['energy_eV']-refs['co2']['energy_eV']
    save(r/'cp2k_fragments.json',{'vertical_interaction_eV':e,'slab_deformation_eV':deform,'co2_deformation_eV':gasdeform,
        'total_reference_difference_eV':e+deform+gasdeform,'interpretation':'The isolated MACE-relaxed slab reference has higher PBE single-point energy than the adsorbed-geometry slab fragment. These PBE single points do not establish DFT stationary basin identities; the total is not pure molecular binding.','basis_superposition_correction':False,'manuscript_eligible':False})

def aimd_prepare(r,cfg):
    base=read(r/'dft-relax/final.traj');jobs=[];plans=[]
    temperatures=[100,200,300,400,150,350,125,325];angles=[-60,-20,20,60,-40,40,-30,30];heights=[-.2,0,.2,.4,-.1,.3,-.05,.25]
    for i in range(8):
        role='train' if i<4 else ('valid' if i<6 else 'test');rng=np.random.default_rng(cfg['sampling']['seed']+i)
        a=base.copy();a.positions+=rng.normal(0,.025,(len(a),3));ads=a[44:];center=ads.positions[0].copy();ads.rotate(angles[i],'z',center=center);a.positions[44:]=ads.positions+[0,0,heights[i]]
        assert chemistry(a)['geometry_gate'];d=r/'aimd'/f'replica-{i:02d}'
        electronic=load(r/'force_grid_audit.json')['selected']
        cp_input(a,d,'md',**electronic,steps=cfg['sampling']['steps_per_replica'],temp=temperatures[i],seed=cfg['sampling']['seed']+100+i,wfn=r/'dft-relax/opt-RESTART.wfn')
        plan={'id':i,'split':role,'seed':cfg['sampling']['seed']+100+i,'temperature_K':temperatures[i],'rotation_z_deg':angles[i],'co2_height_offset_A':heights[i],'rattle_sigma_A':.025,'initial_chemistry':chemistry(a)}
        save(d/'plan.json',plan);write(d/'initial.traj',a);plans.append(plan);jobs.append({'path':str(d.relative_to(r)),'input':'md.inp','output':'md.log'})
    save(r/'sampling_plan.json',{'replicas':plans,'design':'Split assigned before labels/model fitting; entire independent initial condition/trajectory belongs to exactly one split. Temperatures/orientations span a small local basin, not material-wide coverage.'})
    save(r/'aimd_jobs.json',jobs)

def pair_features(a):
    d=a.get_all_distances(mic=True);return d[np.triu_indices(len(a),1)]

def export(r,cfg):
    from aimd_to_deepmd import frames
    from cp2k_validate import parse
    allframes=[];audits=[]
    for i in range(8):
        d=r/'aimd'/f'replica-{i:02d}';parse(d/'md.log');plan=load(d/'plan.json');base=read(d/'initial.traj')
        x=frames(d/'md-pos-1.xyz');f=frames(d/'md-frc-1.xyz');cells=np.loadtxt(d/'md-1.cell');ener=np.loadtxt(d/'md-1.ener')
        expected=cfg['sampling']['steps_per_replica']+1
        if len(x)!=len(f) or len(x)!=expected or len(cells)!=expected or len(ener)!=expected:raise ValueError('Incomplete trajectory')
        if not np.array_equal(cells[:,0],np.arange(expected)) or not np.array_equal(ener[:,0],np.arange(expected)):raise ValueError('Cell/energy step mismatch')
        logforces=forces_from_log(d/'md.log');assert len(logforces)==expected
        rmsds=[];mins=[];fmax=[];features=[];force_sums=[]
        for j,(xx,ff) in enumerate(zip(x,f)):
            if xx[0]!=j or ff[0]!=j or xx[2]!=ff[2] or xx[2]!=base.get_chemical_symbols():raise ValueError('Atom ordering/frame ID mismatch')
            a=base.copy();a.positions=xx[3];a.cell=cells[j,2:11].reshape(3,3)
            force=ff[3]*HARTREE/BOHR;energy=xx[1]*HARTREE
            if not np.isfinite(force).all() or not np.isfinite(energy) or np.linalg.det(a.cell)<=0:raise ValueError('Invalid label')
            assert np.allclose(force,logforces[j],atol=3e-6), 'Independent force conversion mismatch'
            assert abs(energy-ener[j,4]*HARTREE)<2e-6, 'Independent energy alignment mismatch'
            chem=chemistry(a);assert chem['geometry_gate'],'Trajectory chemistry failed'
            force_sums.append(force.sum(axis=0).tolist())
            disp=find_mic(a.positions-base.positions,a.cell,pbc=True)[0];disp-=disp.mean(axis=0);rmsds.append(float(np.sqrt((disp**2).sum(axis=1).mean())))
            mins.append(chem['minimum_distance_A']);fmax.append(float(np.linalg.norm(force,axis=1).max()));features.append(pair_features(a))
            if j>=cfg['sampling']['discard_initial_steps'] and j%cfg['sampling']['frame_stride']==0:
                allframes.append({'atoms':a,'force':force,'energy':energy,'replica':i,'step':j,'role':plan['split'],'feature':features[-1]})
        v=ener[:,4]-ener[:,4].mean();acf=np.correlate(v,v,'full')[len(v)-1:];acf=acf/acf[0] if acf[0]>0 else np.ones_like(acf)
        drift=np.polyfit(ener[:,1]/1000,ener[:,5]*HARTREE/len(base)*1000,1)[0]
        audits.append({'replica':i,'split':plan['split'],'steps':ener[:,0].tolist(),'time_fs':ener[:,1].tolist(),'temperature_K':ener[:,3].tolist(),
            'potential_eV':(ener[:,4]*HARTREE).tolist(),'conserved_eV':(ener[:,5]*HARTREE).tolist(),'energy_drift_meV_atom_ps':float(drift),
            'net_force_vectors_eV_A':force_sums,'rms_displacement_A':rmsds,'minimum_distance_A':mins,'fmax_eV_A':fmax,'potential_acf':acf.tolist(),'duration_fs':float(ener[-1,1]),
            'decorrelation_established':False,'raw_sha256':{p.name:digest(p) for p in d.iterdir() if p.suffix in ['.xyz','.cell','.ener','.inp','.log']}})
    features=np.array([x['feature'] for x in allframes]);dist=np.sqrt(((features[:,None,:]-features[None,:,:])**2).mean(axis=2));np.fill_diagonal(dist,np.inf)
    roles=np.array([x['role'] for x in allframes]);cross=roles[:,None]!=roles[None,:];minimum=float(dist[cross].min());threshold=cfg['sampling']['near_geometry_rmsd_A']
    # Never move heldout structures into training or silently discard hard labels.
    if minimum<threshold:raise ValueError(f'Near-geometry cross-split leakage {minimum} < {threshold}; redesign independent starts')
    dataset=r/'dataset';dataset.mkdir();index=[]
    for role in ['train','valid','test']:
        chosen=[x for x in allframes if x['role']==role];dest=dataset/role;(dest/'set.000').mkdir(parents=True)
        arrays={'coord':np.array([x['atoms'].positions.flatten() for x in chosen]),'box':np.array([x['atoms'].cell.array.flatten() for x in chosen]),
           'force':np.array([x['force'].flatten() for x in chosen]),'energy':np.array([x['energy'] for x in chosen])}
        for name,value in arrays.items():np.save(dest/'set.000'/f'{name}.npy',value)
        (dest/'type_map.raw').write_text('\n'.join(TYPES)+'\n');np.savetxt(dest/'type.raw',[TYPES.index(s) for s in chosen[0]['atoms'].get_chemical_symbols()],fmt='%d')
        for x in chosen:index.append({'role':role,'replica':x['replica'],'step':x['step'],'energy_eV':x['energy'],'fmax_eV_A':float(np.linalg.norm(x['force'],axis=1).max())})
        write(dataset/f'{role}.traj',[x['atoms'] for x in chosen])
    finite_dist=np.where(np.isfinite(dist),dist,0);np.save(dataset/'geometry_distance_matrix.npy',finite_dist)
    replica_ids=np.array([x['replica'] for x in allframes]);coverage=[]
    for i in range(8):
        own=replica_ids==i;other=~own
        within=dist[np.ix_(own,own)];finite=within[np.isfinite(within)]
        coverage.append({'replica':i,'split':str(roles[own][0]),'selected_frames':int(own.sum()),
            'maximum_within_replica_pair_distance_rms_A':float(finite.max()) if len(finite) else 0.,
            'minimum_other_replica_pair_distance_rms_A':float(dist[np.ix_(own,other)].min()),
            'minimum_other_split_pair_distance_rms_A':float(dist[np.ix_(own,roles!=roles[own][0])].min())})
    save(r/'sampling_coverage.json',{'metric':'RMS difference of all minimum-image pair distances with fixed atom identities',
        'replicas':coverage,'statistical_independence_established':False,'effective_sample_size_estimated':False,
        'interpretation':'Distinct seeded initial conditions and disjoint trajectories are local structural holdouts. Four fs cannot establish equilibrium, decorrelation, or an effective independent sample size. Passing the geometric-separation gate is not proof of statistical independence.',
        'manuscript_eligible':False})
    save(r/'dataset_audit.json',{'index':index,'frames_per_split':{x:int(sum(roles==x)) for x in ['train','valid','test']},'cross_split_min_pair_distance_rms_A':minimum,
        'threshold_A':threshold,'replica_disjoint':True,'near_geometry_gate':minimum>=threshold,'units':{'energy':'eV','force':'eV/A','cell':'A'},
        'net_force_max_norm_eV_A':float(max(np.linalg.norm(x['force'].sum(axis=0)) for x in allframes)),
        'translation_constraint_force_rmse_floor_eV_A':float(np.sqrt(np.mean(np.array([x['force'].sum(axis=0)/len(x['atoms']) for x in allframes])**2))),
        'type_map':TYPES,'virial_labels':False,'reason_no_virial':'Fixed-cell small sample targets energies/forces; no pressure generalization claimed.',
        'sampling_caveat':'Within each trajectory frames remain correlated. Entire trajectories are disjoint, but eight nearby starts do not establish broad generalization.','manuscript_eligible':False})
    save(r/'aimd_audit.json',{'replicas':audits,'manuscript_eligible':False})

def train_prepare(r,cfg):
    from ase.neighborlist import neighbor_list
    # Check descriptor neighbor capacity on all sampled structures, including PBC images.
    maxima=np.zeros(len(TYPES),dtype=int)
    for role in ['train','valid','test']:
        for atoms in read(r/f'dataset/{role}.traj',index=':'):
            ii,jj=neighbor_list('ij',atoms,5.0)
            for t,el in enumerate(TYPES):
                selected=np.array(atoms.get_chemical_symbols())[jj]==el
                maxima[t]=max(maxima[t],np.bincount(ii[selected],minlength=len(atoms)).max())
    assert np.all(maxima<=np.array([32,32,96,32,16])),'Insufficient descriptor neighbor capacity'
    save(r/'descriptor_audit.json',{'type_map':TYPES,'rcut_A':5.0,'neighbor_capacity':[32,32,96,32,16],'observed_max_neighbors':maxima.tolist(),
         'capacity_gate':True,'electrostatics':'No explicit long-range term; restricted to fixed composition/cell local test','manuscript_eligible':False})
    dataset=r/'dataset';over=dataset/'overfit';shutil.copytree(dataset/'train',over)
    for p in (over/'set.000').glob('*.npy'):np.save(p,np.load(p)[:1])
    for label,desc,fit,steps in [('overfit',[8,16,32],[32,32],cfg['training']['overfit_steps']),('compact',[8,16,32],[32,32],cfg['training']['compact_steps']),('wider',[16,32,64],[64,64,64],cfg['training']['wider_steps'])]:
        d=r/'models'/label;d.mkdir(parents=True);train='overfit' if label=='overfit' else 'train';valid='overfit' if label=='overfit' else 'valid'
        conf={'model':{'type_map':TYPES,'descriptor':{'type':'se_e2_a','sel':[32,32,96,32,16],'rcut_smth':4.5,'rcut':5.,'neuron':desc,'axis_neuron':4,'precision':'float64','seed':814},'fitting_net':{'neuron':fit,'precision':'float64','seed':815}},
          'learning_rate':{'type':'exp','start_lr':.001,'stop_lr':1e-5,'decay_steps':max(100,steps//20)},
          'loss':{'type':'ener','start_pref_e':.1,'limit_pref_e':1,'start_pref_f':100,'limit_pref_f':10,'start_pref_v':0,'limit_pref_v':0},
          'training':{'training_data':{'systems':['../../dataset/'+train],'batch_size':1 if label=='overfit' else 4},'validation_data':{'systems':['../../dataset/'+valid],'batch_size':1 if label=='overfit' else 4,'numb_btch':1 if label=='overfit' else 2},'numb_steps':steps,'seed':816,'disp_file':'lcurve.out','disp_freq':50,'save_freq':steps//2}}
        save(d/'input.json',conf)

def evaluate(r,cfg):
    from deepmd.infer import DeepPot
    import torch
    torch.set_num_threads(4);torch.cuda.set_per_process_memory_fraction(.35)
    results={};models={}
    def predict(model,role,folder):
        d=r/'dataset'/role;coords=np.load(d/'set.000/coord.npy');box=np.load(d/'set.000/box.npy');types=np.loadtxt(d/'type.raw',dtype=int)
        e,f,v=model.eval(coords,box,types);truth_e=np.load(d/'set.000/energy.npy');truth_f=np.load(d/'set.000/force.npy').reshape(f.shape)
        e=e.reshape(-1);diff_f=f-truth_f;diff_e=(e-truth_e)/len(types)
        np.savez(folder/f'{role}_predictions.npz',reference_energy_eV=truth_e,predicted_energy_eV=e,reference_force_eV_A=truth_f,predicted_force_eV_A=f)
        records=load(r/'dataset_audit.json')['index']
        selected_records=records[:1] if role=='overfit' else [x for x in records if x['role']==role]
        assert len(selected_records)==len(e),'Prediction/frame provenance length mismatch'
        replica_ids=np.array([x['replica'] for x in selected_records]);per_replica=[]
        for replica in np.unique(replica_ids):
            mask=replica_ids==replica
            per_replica.append({'replica':int(replica),'frames':int(mask.sum()),
                'energy_rmse_eV_atom':float(np.sqrt(np.mean(diff_e[mask]**2))),
                'force_rmse_eV_A':float(np.sqrt(np.mean(diff_f[mask]**2)))})
        return {'per_replica':per_replica,'frames':len(e),'energy_rmse_eV_atom':float(np.sqrt(np.mean(diff_e**2))),'force_rmse_eV_A':float(np.sqrt(np.mean(diff_f**2))),
          'force_p95_abs_error_eV_A':float(np.quantile(abs(diff_f),.95)),'force_max_abs_error_eV_A':float(abs(diff_f).max()),
          'per_element_force_rmse_eV_A':{el:float(np.sqrt(np.mean(diff_f[:,types==i,:]**2))) for i,el in enumerate(TYPES)}}
    for label in ['overfit','compact','wider']:
        d=r/'models'/label;model=DeepPot(str(d/'model.pth'));models[label]=model
        roles=['overfit'] if label=='overfit' else ['train','valid'];results[label]={role:predict(model,role,d) for role in roles}
        results[label]['model_sha256']=digest(d/'model.pth')
    selected=min(['compact','wider'],key=lambda n:results[n]['valid']['force_rmse_eV_A'])
    results[selected]['test']=predict(models[selected],'test',r/'models'/selected)
    model=models[selected];d=r/'dataset/test';coord=np.load(d/'set.000/coord.npy')[:1];box=np.load(d/'set.000/box.npy')[:1];types=np.loadtxt(d/'type.raw',dtype=int)
    e,f,v=model.eval(coord,box,types)
    translated=coord.reshape(1,-1,3)+np.array([.73,-.41,.28]);et,ft,_=model.eval(translated,box,types)
    perm=np.arange(len(types));rng=np.random.default_rng(920)
    for t in np.unique(types):ids=np.where(types==t)[0];perm[ids]=rng.permutation(ids)
    ep,fp,_=model.eval(coord.reshape(1,-1,3)[:,perm],box,types[perm])
    angle=.37;rot=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    er,fr,_=model.eval(coord.reshape(1,-1,3)@rot.T,box.reshape(1,3,3)@rot.T,types)
    sym={'translation_energy_eV':float(abs(et-e).max()),'translation_force_eV_A':float(abs(ft-f).max()),
         'permutation_energy_eV':float(abs(ep-e).max()),'permutation_force_eV_A':float(abs(fp-f[:,perm]).max()),
         'rotation_energy_eV':float(abs(er-e).max()),'rotation_force_eV_A':float(abs(fr-f@rot.T).max())}
    save(r/'model_evaluation.json',{'models':results,'selected':selected,'selection_rule':'lowest full validation-set force RMSE; test used only after architecture selection',
        'overfit_gate':results['overfit']['overfit']['force_rmse_eV_A']<cfg['gates']['overfit_force_rmse_eV_A'],
        'heldout_force_gate':results[selected]['test']['force_rmse_eV_A']<cfg['gates']['force_rmse_eV_A'],
        'heldout_energy_gate':results[selected]['test']['energy_rmse_eV_atom']<cfg['gates']['energy_rmse_eV_atom'],
        'symmetry':sym,'symmetry_gate':max(sym.values())<1e-5,'manuscript_eligible':False})

def lammps_prepare(r,cfg):
    selected=load(r/'model_evaluation.json')['selected'];model=r/'models'/selected/'model.pth';d=r/'lammps';d.mkdir()
    a=read(r/'dataset/test.traj',index=0);write(d/'system.data',a,format='lammps-data',atom_style='atomic',specorder=TYPES,masses=True);write(d/'initial.traj',a);shutil.copyfile(model,d/'model.pth')
    shared='''units metal
atom_style atomic
boundary p p p
read_data system.data
pair_style deepmd model.pth
pair_coeff * * Ca P O H C
neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
timestep 0.00025
group adsorbate id 45 46 47
group CO2_oxygen id 46 47
velocity all create 200.0 951 mom yes rot yes dist gaussian
fix integ all nve
thermo 10
thermo_style custom step time temp pe ke etotal press
thermo_modify format float %.16g
'''
    (d/'mapping.in').write_text(shared+'''run 0
write_dump all custom mapping.dump id type xu yu zu fx fy fz modify sort id format float %.16g
''')
    dump='''dump trajectory all custom 10 trajectory.dump id type xu yu zu vx vy vz fx fy fz
dump_modify trajectory sort id format float %.16g
'''
    (d/'continuous.in').write_text(shared+dump+'''run 1000
write_dump all custom continuous-final.dump id type xu yu zu vx vy vz modify sort id format float %.16g
''')
    (d/'halfstep.in').write_text(shared.replace('timestep 0.00025','timestep 0.000125')+'''run 2000
write_dump all custom halfstep-final.dump id type xu yu zu vx vy vz modify sort id format float %.16g
''')
    (d/'first.in').write_text(shared+'''run 500
write_restart state.restart
''')
    (d/'restart.in').write_text('''units metal
atom_style atomic
read_restart state.restart
pair_style deepmd model.pth
pair_coeff * * Ca P O H C
neighbor 1.0 bin
neigh_modify every 1 delay 0 check yes
fix integ all nve
thermo 10
thermo_style custom step time temp pe ke etotal press
thermo_modify format float %.16g
run 500
write_dump all custom restart-final.dump id type xu yu zu vx vy vz modify sort id format float %.16g
''')
    save(d/'provenance.json',{'selected_model':selected,'model_sha256':digest(model),'source_geometry':'dataset/test.traj frame 0',
        'types':TYPES,'atom_ids':'1-based original ordering; CO2=45,46,47; O=46,47','steps':1000,'timestep_fs':.25,'time_fs':250,
        'constraints':'all atoms mobile','ensemble':'NVE','temperature_initial_K':200,'manuscript_eligible':False,
        'warning':'Diagnostic trajectory from a tiny local model; failed heldout gates remain failed; not production dynamics.'})

def lammps_audit(r,cfg):
    from deepmd.infer import DeepPot
    from ase.io import read as ase_read
    d=r/'lammps';model=DeepPot(str(d/'model.pth'));a=read(d/'initial.traj');types=[TYPES.index(s) for s in a.get_chemical_symbols()]
    e,f,_=model.eval(a.positions.reshape(1,-1,3),a.cell.array.reshape(1,3,3),types);mapping=np.loadtxt(d/'mapping.dump',skiprows=9)
    assert np.array_equal(mapping[:,0],np.arange(1,48)) and np.array_equal(mapping[:,1],np.array(types)+1)
    force_error=float(abs(mapping[:,5:8]-f[0]).max())
    def thermo(path):
        s=path.read_text();assert 'Total wall time:' in s and 'ERROR' not in s
        lines=s.splitlines();rows=[];active=False
        for line in lines:
            if line.strip().startswith('Step') and 'TotEng' in line:active=True;continue
            if active:
                fields=line.split()
                if len(fields)==7:
                    try:rows.append([float(x) for x in fields]);continue
                    except ValueError:pass
                if rows:break
        return np.array(rows)
    mapping_t=thermo(d/'mapping.log');energy_error=float(abs(mapping_t[0,3]-e[0,0]));t=thermo(d/'continuous.log');np.save(d/'thermo.npy',t)
    x=np.loadtxt(d/'continuous-final.dump',skiprows=9);y=np.loadtxt(d/'restart-final.dump',skiprows=9);assert np.array_equal(x[:,:2],y[:,:2]),'Restart atom identity mismatch'
    position_difference=float(abs(x[:,2:5]-y[:,2:5]).max());velocity_difference=float(abs(x[:,5:8]-y[:,5:8]).max());restart=max(position_difference,velocity_difference)
    h=np.loadtxt(d/'halfstep-final.dump',skiprows=9)
    assert np.array_equal(x[:,:2],h[:,:2]) and np.isfinite(h).all(),'Half-timestep atom/state mismatch'
    half_position_delta=x[:,2:5]-h[:,2:5]
    half_state={'final_position_rms_difference_A':float(np.sqrt(np.mean(np.sum(half_position_delta**2,axis=1)))),
        'final_position_max_component_difference_A':float(abs(half_position_delta).max()),
        'final_velocity_max_component_difference_A_per_ps':float(abs(x[:,5:8]-h[:,5:8]).max()),
        'scope':'Same initial state and physical duration, different timesteps; diagnostic values without an after-the-fact acceptance threshold.'}
    traj=ase_read(d/'trajectory.dump',format='lammps-dump-text',index=':',specorder=TYPES)
    assert len(traj)==101 and np.isfinite(t).all();write(d/'trajectory.extxyz',traj)
    metrics=[];base=traj[0]
    for k,b in enumerate(traj):
        c=chemistry(b);disp=find_mic(b.positions-base.positions,b.cell,True)[0];disp-=disp.mean(axis=0)
        c.update(time_fs=k*10*.25,rms_displacement_A=float(np.sqrt((disp**2).sum(axis=1).mean())),co2_angle_deg=float(b.get_angle(45,44,46,mic=True)),co2_bonds_A=[float(b.get_distance(44,i,mic=True)) for i in [45,46]])
        metrics.append(c)
    drift=float(np.polyfit(t[:,1],(t[:,5]-t[0,5])/47*1000,1)[0]);half=thermo(d/'halfstep.log');np.save(d/'halfstep-thermo.npy',half);halfdrift=float(np.polyfit(half[:,1],(half[:,5]-half[0,5])/47*1000,1)[0])
    save(r/'lammps_audit.json',{'engine_energy_difference_eV':energy_error,'engine_force_max_difference_eV_A':force_error,
        'mapping_gate':force_error<cfg['gates']['engine_force_max_difference_eV_A'] and energy_error<1e-5,'restart_max_difference':restart,'restart_position_max_difference_A':position_difference,'restart_velocity_max_difference_A_per_ps':velocity_difference,'restart_ID_type_match':True,
        'restart_gate':restart<cfg['gates']['restart_max_difference'],'energy_drift_meV_atom_ps':drift,'drift_gate':abs(drift)<cfg['gates']['nve_drift_meV_atom_ps'],
        'half_timestep_state_comparison':half_state,'half_timestep_fs':.125,'half_timestep_drift_meV_atom_ps':halfdrift,'timestep_drift_gate':abs(halfdrift)<cfg['gates']['nve_drift_meV_atom_ps'] and abs(drift)<cfg['gates']['nve_drift_meV_atom_ps'],
        'chemistry_gate':all(m['geometry_gate'] for m in metrics),'frames':metrics,'time_fs':250,'manuscript_eligible':False})

def optimizer_trial_audit(r,cfg):
    from aimd_to_deepmd import frames
    d=r/'dft-relax';log=(d/'opt.log').read_text()
    energies=list(re.finditer(r'ENERGY\| Total FORCE_EVAL.*?energy \[hartree\]\s+([-+\d.Ee]+)',log))
    force_headers=list(re.finditer(re.escape('FORCES| Atomic forces [hartree/bohr]'),log))
    previous=r/'dft-relax-attempt3/opt-pos-1.xyz'
    # CP2K's current trajectory writes the accepted state, not necessarily the initial restart.
    # The preceding saved frame is independently linked by cp2k_restart_audit.json.
    saved=np.array([frames(previous)[-1][1]*HARTREE]+[x[1]*HARTREE for x in frames(d/'opt-pos-1.xyz')])
    forces=forces_from_log(d/'opt.log');evaluations=[]
    assert len(force_headers)==len(forces)
    for i,(header,force) in enumerate(zip(force_headers,forces)):
        preceding=[e for e in energies if e.end()<header.start()];assert preceding
        energy=float(preceding[-1][1])*HARTREE
        evaluations.append({'force_evaluation':i+1,'energy_eV':energy,'fmax_eV_A':float(np.linalg.norm(force,axis=1).max()),
            'matches_initial_or_accepted_saved_energy':bool(np.min(abs(saved-energy))<2e-6)})
    save(r/'optimizer_trial_audit.json',{'evaluations':evaluations,'saved_energies_eV':saved.tolist(),
        'interpretation':'Force callbacks include line-search trial points. Evaluations absent from the saved initial/accepted geometry energies are not labeled as accepted configurations. Only the final accepted state passes to AIMD.',
        'sources':{p.relative_to(r).as_posix():digest(p) for p in [d/'opt.log',d/'opt-pos-1.xyz',previous,r/'cp2k_restart_audit.json']},'manuscript_eligible':False})

def report(r,cfg):
    if (r/'dft_lbfgs_recovery.json').exists():optimizer_trial_audit(r,cfg)
    required=['structure_audit.json','mace.json','mace_refined.json','reference_audit.json','mace_fragments.json','cp2k_fragments.json','force_grid_audit.json','dft_relaxation.json','cp2k_restart_audit.json','aimd_audit.json','dataset_audit.json','sampling_coverage.json','model_evaluation.json','lammps_audit.json']
    result={'run_id':r.name,'scope':cfg['scope'],'manuscript_eligible':False,'sources':{},'gates':{},'missing':[]}
    for name in required:
        p=r/name
        if p.exists():result['sources'][name]={'sha256':digest(p),'data':load(p)}
        else:result['missing'].append(name)
    for name in ['optimizer_diagnosis-step90.json','dft_lbfgs_recovery.json','optimizer_trial_audit.json','reference_interpretation_note.json']:
        if (r/name).exists():result['sources'][name]={'sha256':digest(r/name),'data':load(r/name)}
    def gate(name,file,key,meaning):
        value=result['sources'].get(file,{}).get('data',{}).get(key)
        result['gates'][name]={'passed':value,'meaning':meaning,'source':file}
    gate('cutoff_reference','reference_audit.json','cutoff_gate','500 vs700 Ry adsorption difference <0.02eV; not basis/slab convergence')
    gate('mace_refined','mace_refined.json','converged','Selected MACE geometry satisfies0.03eV/A maximum-force threshold')
    gate('dft_optimizer','dft_relaxation.json','optimizer_converged','CP2K reports all optimizer convergence criteria satisfied')
    gate('dft_force','dft_relaxation.json','force_gate','Final DFT maximum force <0.05eV/A')
    gate('force_grid','force_grid_audit.json','force_gate','900/80 vs1100/100Ry force-component RMSE <0.02eV/A on selected geometry')
    if result['sources'].get('cp2k_restart_audit.json',{}).get('data',{}).get('checks'):
        gate('cp2k_restart','cp2k_restart_audit.json','restart_repeatability_gate','Same-state CP2K restart: E within1e-6eV and force components within1e-4eV/A')
    gate('split_leakage','dataset_audit.json','near_geometry_gate','Disjoint replicas and cross-split pair-distance RMS >=0.025A')
    gate('overfit_sanity','model_evaluation.json','overfit_gate','One-frame force fit <0.02eV/A; not generalization')
    gate('independent_force','model_evaluation.json','heldout_force_gate','Full independent test force RMSE <0.1eV/A')
    gate('independent_energy','model_evaluation.json','heldout_energy_gate','Full independent test energy RMSE <0.01eV/atom')
    gate('model_symmetry','model_evaluation.json','symmetry_gate','Rigid translation/rotation and equal-element permutation invariance/covariance')
    gate('engine_mapping','lammps_audit.json','mapping_gate','Python and LAMMPS agree on same geometry, types, energy and forces')
    gate('lammps_restart','lammps_audit.json','restart_gate','500+500 restart matches1000 steps: positions within1e-8A and velocities within1e-8A/ps')
    gate('short_nve_drift','lammps_audit.json','drift_gate','250fs NVE drift magnitude <5meV/atom/ps; short-duration only')
    gate('timestep_drift','lammps_audit.json','timestep_drift_gate','Both 0.25fs and0.125fs over250fs NVE drift magnitudes <5meV/atom/ps')
    gate('trajectory_geometry','lammps_audit.json','chemistry_gate','Charge, contact, PO4/OH geometry remain sane during diagnostic trajectory')
    result['functional_complete']=not result['missing'];result['all_pilot_gates_pass']=not result['missing'] and all(x['passed'] is True for x in result['gates'].values())
    result['production_ready']=False;result['unmet_production_requirements']=['hydrated environments and termination selection','slab thickness/lateral size/vacuum/basis/k-point convergence','electronic polarization and explicit long-range model validation','independent reaction-path and broad configuration coverage','long-time dynamics/timestep/ensemble convergence','review of chemical mechanisms and claim-specific publication eligibility']
    save(r/'verification.json',result);print(json.dumps({k:v['passed'] for k,v in result['gates'].items()},indent=2))

def main():
    p=argparse.ArgumentParser();p.add_argument('stage');p.add_argument('--run',type=Path,required=True);p.add_argument('--config',type=Path,default=ROOT/'workflow/profiles/rtx4090-refine.json');a=p.parse_args()
    r=a.run.resolve();cfg=load(a.config);validate_sampling_profile(cfg);r.mkdir(parents=True,exist_ok=True)
    marker=r/(a.stage+'.manifest.json')
    if marker.exists():raise FileExistsError(f'Preserving {marker}; use a fresh run or a new stage attempt')
    snapshots=r/'source-snapshots';snapshots.mkdir(exist_ok=True)
    for src in [Path(__file__),a.config,Path(__file__).with_name('local_compute.py'),Path(__file__).with_name('run_refinement.py')]:
        target=snapshots/(digest(src)+'-'+src.name)
        if not target.exists():shutil.copyfile(src,target)
    m={'stage':a.stage,'status':'running','source_sha256':digest(__file__),'config_sha256':digest(a.config),'command':sys.argv,'manuscript_eligible':False};start=time.time();save(marker,m)
    try:
        globals()[a.stage.replace('-','_')](r,cfg);m['status']='completed'
    except Exception as exc:m['status']='failed';m['error']=repr(exc);raise
    finally:m['elapsed_seconds']=time.time()-start;save(marker,m)

if __name__=='__main__':main()
