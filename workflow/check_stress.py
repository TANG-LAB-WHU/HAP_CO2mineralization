import argparse,json
from pathlib import Path
import numpy as np
from cp2k_validate import parse,virial
p=argparse.ArgumentParser();p.add_argument('directory',type=Path);a=p.parse_args()
records={n:parse(a.directory/n/'co2.log') for n in ['minus','zero','plus']}
settings=json.loads((a.directory/'settings.json').read_text())
stress=records['zero']['cp2k_stress'];v=virial(stress,settings['volume_A3'],records['zero']['stress_unit'])
numerical=(records['plus']['energy_eV']-records['minus']['energy_eV'])/(2*settings['isotropic_strain'])
analytical=-float(np.trace(v))
# Grid/basis Pulay effects require cutoff convergence; this bounded check tests sign and scale.
tolerance=max(.1,abs(analytical)*.3)
result={'finite_difference_dE_dstrain_eV':numerical,'minus_trace_virial_eV':analytical,
 'absolute_error_eV':abs(numerical-analytical),'tolerance_eV':tolerance,
 'passed':abs(numerical-analytical)<tolerance,'records':records,'manuscript_eligible':False}
(a.directory/'comparison.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if not result['passed']:raise SystemExit('Stress finite-difference mismatch; preserve result and investigate')
