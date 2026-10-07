"""Check CP2K completion/SCF and parse physical quantities without inventing labels."""
import argparse, json, re, hashlib
from pathlib import Path
import numpy as np
HARTREE_EV=27.211386245988

def parse(path):
    text=Path(path).read_text(errors='replace')
    if 'PROGRAM ENDED AT' not in text or 'SCF run converged' not in text or 'SCF run NOT converged' in text or '[ABORT]' in text:
        raise ValueError(f'Incomplete or unconverged CP2K output: {path}')
    energies=re.findall(r'ENERGY\| Total FORCE_EVAL.*?energy \[hartree\]\s+([-+\d.Ee]+)',text)
    if not energies: raise ValueError('Missing energy')
    result={'energy_eV':float(energies[-1])*HARTREE_EV,'scf_converged':True,
            'output_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest()}
    m=re.search(r'(?:STRESS\|\s*)?Analytical stress tensor \[(bar|GPa)\](.*?)(?:Trace|Determinant|Eigenvalues)',text,re.S)
    if m:
        rows=re.findall(r'(?:STRESS\|)?\s+[xyz]\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)\s+([-+\d.Ee]+)',m[2])
        if len(rows)!=3: raise ValueError('Malformed stress tensor')
        result['cp2k_stress']=np.array(rows,dtype=float).tolist(); result['stress_unit']=m[1]
    if not np.isfinite(result['energy_eV']): raise ValueError('Nonfinite energy')
    return result

def virial(stress, volume, unit):
    # CP2K compressive-positive stress; ASE has the opposite sign.
    factors={'bar':1e5/1.602176634e11,'GPa':1e9/1.602176634e11}
    if unit not in factors or not np.isfinite(volume) or volume<=0: raise ValueError('Invalid volume/unit')
    tensor=np.asarray(stress,dtype=float)
    if tensor.shape!=(3,3) or not np.isfinite(tensor).all(): raise ValueError('Invalid tensor')
    return tensor*volume*factors[unit]

def main():
    p=argparse.ArgumentParser(); p.add_argument('directory',type=Path); a=p.parse_args()
    records={name:parse(a.directory/f'{name}.cp2k.log') for name in ['slab','co2','adsorbed']}
    records['adsorption_energy_eV']=records['adsorbed']['energy_eV']-records['slab']['energy_eV']-records['co2']['energy_eV']
    records['manuscript_eligible']=False
    records['interpretation']='PBE single points on smoke MACE geometries; no DFT relaxation, basis/cutoff/slab/vacuum convergence or dispersion validation'
    (a.directory/'cp2k_results.json').write_text(json.dumps(records,indent=2))
    print(json.dumps(records,indent=2))
if __name__=='__main__':main()
