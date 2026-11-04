"""Fetch exactly the tested MACE-MP-0 small release; verify before use."""
import argparse,hashlib,json,urllib.request
from pathlib import Path
URL='https://github.com/ACEsuit/mace-mp/releases/download/mace_mp_0/2023-12-10-mace-128-L0_energy_epoch-249.model'
SHA='2ddb079cee0e131eaaf6912ba581b394551ead283e95c99cfe78c605d10b5736'
p=argparse.ArgumentParser();p.add_argument('--out',type=Path,default=Path('.research-runs/models'));p.add_argument('--mattersim',action='store_true');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
if a.mattersim:
    URL='https://raw.githubusercontent.com/microsoft/mattersim/main/pretrained_models/mattersim-v1.0.0-1M.pth'
    SHA='28b0b0b0f13efefee06b47ea4c9105a26bd3e2c8396da193430da96b3b49a8be'
dest=a.out/('mattersim-v1.0.0-1M.pth' if a.mattersim else 'mace-mp-small.model')
if not dest.exists():
    part=dest.with_suffix('.part')
    if part.exists():raise FileExistsError('Inspect previous incomplete download: '+str(part))
    with urllib.request.urlopen(URL,timeout=60) as source,part.open('wb') as target:
        total=0
        while chunk:=source.read(1024*1024):
            total+=len(chunk)
            if total>50*1024*1024:raise ValueError('Unexpected model size')
            target.write(chunk)
    if hashlib.sha256(part.read_bytes()).hexdigest()!=SHA:raise ValueError('Checkpoint checksum mismatch')
    part.rename(dest)
if hashlib.sha256(dest.read_bytes()).hexdigest()!=SHA:raise ValueError('Existing checkpoint checksum mismatch')
(a.out/('mattersim-provenance.json' if a.mattersim else 'provenance.json')).write_text(json.dumps({'url':URL,'sha256':SHA,'bytes':dest.stat().st_size,'license':'see upstream model release','purpose':'functional smoke; no HAP accuracy guarantee'},indent=2))
print(dest)
