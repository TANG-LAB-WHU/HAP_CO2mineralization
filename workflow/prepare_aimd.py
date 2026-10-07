from pathlib import Path
import argparse
p=argparse.ArgumentParser(); p.add_argument('runs',type=Path);p.add_argument('--system',choices=['co2','adsorbed'],default='co2');p.add_argument('--steps',type=int,default=5)
a=p.parse_args(); r=a.runs.resolve(); system=a.system
source=(r/f'smoke-001/{system}.inp').read_text()
for i,temp in enumerate([100,200,300]):
    d=r/f'aimd-{system}-{i}'
    if d.exists(): raise FileExistsError(d)
    d.mkdir()
    s=source.replace(f'PROJECT {system}',f'PROJECT {system}_{i}').replace('RUN_TYPE ENERGY_FORCE','RUN_TYPE MD\n  SEED '+str(710+i))
    # Same geometry, cell, basis and electronic settings; reuse only the initial wavefunction.
    import shutil
    shutil.copyfile(r/f'smoke-001/{system}-RESTART.wfn',d/'initial.wfn')
    s=s.replace('  &DFT','  &DFT\n    WFN_RESTART_FILE_NAME initial.wfn',1)
    s=s.replace('    &SCF','    &SCF\n      SCF_GUESS RESTART',1)
    s += f'''
&MOTION
  &MD
    ENSEMBLE NVE
    STEPS {a.steps}
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
        MD 1
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION
'''
    (d/'md.inp').write_text(s)
