"""Capture local reproducibility metadata and hashes; no credentials or remote actions."""
import argparse,datetime,hashlib,json,platform,subprocess,sys,shutil
from pathlib import Path
import psutil
p=argparse.ArgumentParser();p.add_argument('runs',type=Path);a=p.parse_args();r=a.runs.resolve();repo=r.parent
def command(args):
    q=subprocess.run(args,cwd=repo,capture_output=True,text=True,encoding='utf-8',errors='replace')
    return {'args':args,'returncode':q.returncode,'stdout':q.stdout,'stderr':q.stderr}
(r/'environment/final-source.patch').write_text(command(['git','diff','HEAD'])['stdout'],encoding='utf-8')
data={'recorded_at_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'computer':platform.node(),'python':sys.version,
 'ram':dict(psutil.virtual_memory()._asdict()),'disk_D':dict(shutil.disk_usage(str(repo))._asdict()),
 'git':command(['git','rev-parse','HEAD']), 'branch':command(['git','branch','--show-current']),
 'status':command(['git','status','--short']), 'nvidia':command(['nvidia-smi','--query-gpu=name,memory.total,driver_version,memory.used,utilization.gpu','--format=csv']),
 'docker_images':command(['docker','image','inspect','deepmodeling/deepmd-kit:3.1.0_cuda129','cp2k/cp2k:2025.1_mpich_generic_cuda_P100_psmp','--format','{{.Id}} {{json .RepoDigests}}']),
 'containers':command(['docker','stats','--no-stream','--format','{{.Name}} {{.MemUsage}} {{.CPUPerc}}']), 'artifacts':[]}
for directory in ['environment','models','entrypoint-smoke','smoke-001','relax-200','smoke-final','stress-check','co2-dataset','adsorbed-dataset','adsorbed-dataset-v2','mattergen-smoke']+[f'aimd-{s}-{i}' for s in ['co2','adsorbed'] for i in range(3)]:
    for f in sorted((r/directory).rglob('*')):
        if f.is_file():
            with f.open('rb') as handle:digest=hashlib.file_digest(handle,'sha256').hexdigest()
            data['artifacts'].append({'path':str(f.relative_to(r)),'bytes':f.stat().st_size,'sha256':digest,'manuscript_eligible':False})
(r/'reproducibility-manifest.json').write_text(json.dumps(data,indent=2))
print('Recorded',len(data['artifacts']),'files; source commit and dirty status retained explicitly')
