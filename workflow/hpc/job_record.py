"""Standard-library Slurm provenance writer; execution never grants scientific eligibility."""
import argparse,datetime,hashlib,json,os,platform,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
def sha(p):
    value=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):value.update(block)
    return value.hexdigest()
def ref(p):
    p=Path(p).resolve();return {'path':str(p),'sha256':sha(p)}
def validated_upstream(path):
    data=json.loads(Path(path).read_text())
    if data.get('status')!='completed' or data.get('validation',{}).get('passed') is not True:
        raise ValueError('Validated predecessor manifest required; afterok alone is insufficient')
    outputs=data.get('outputs',[])
    if not outputs:raise ValueError('Predecessor must identify actual output artifacts')
    for output in outputs:
        artifact=Path(output['path'])
        if not artifact.is_absolute():artifact=Path(path).resolve().parent/artifact
        if not artifact.is_file() or sha(artifact)!=output.get('sha256'):
            raise ValueError(f'Predecessor output missing or changed: {artifact}')
    return data

def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def main():
    p=argparse.ArgumentParser();p.add_argument('action',choices=['start','finish']);p.add_argument('--run',type=Path,required=True);p.add_argument('--stage',choices=['mlip_training','scaleup']);p.add_argument('--input');p.add_argument('--image');p.add_argument('--upstream');p.add_argument('--model');p.add_argument('--restart-file');p.add_argument('--returncode',type=int);p.add_argument('--output',action='append',default=[]);a=p.parse_args()
    run=a.run.resolve();job=os.environ.get('SLURM_JOB_ID','manual')
    if not job.replace('_','').replace('-','').isalnum():raise ValueError('Invalid job identifier')
    path=run/f'job-{job}.manifest.json'
    if a.action=='start':
        if path.exists():raise FileExistsError('Preserving previous job manifest')
        upstream=validated_upstream(a.upstream)
        status=subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)
        if status.strip():raise ValueError('Clean committed source required')
        data={'schema_version':1,'run_id':run.name+'-'+job,'status':'running','backend':'slurm','stage':a.stage,
            'git':{'repository':'TANG-LAB-WHU/HAP_CO2mineralization','branch':subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),'commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'dirty':False},
            'configuration':ref(a.input),'inputs':[ref(a.upstream)]+([ref(a.model)] if a.model else [])+([ref(a.restart_file)] if a.restart_file else []),
            'environment':{'container':ref(a.image),'software':{'host':platform.node(),'python':platform.python_version()},'scheduler_job_id':job,
                'allocation':{k:os.environ.get(k) for k in ['SLURM_CPUS_PER_TASK','SLURM_MEM_PER_NODE','SLURM_JOB_PARTITION','SLURM_JOB_GPUS']}},
            'outputs':[],'metrics':{},'validation':{'criteria':['Independent scientific and completeness checks are required after process termination'],'passed':False,'reviewer':None,'reviewed_at':None},
            'started_at':now(),'finished_at':None,'manuscript_eligible':False,'notes':'Unsubmitted template; normal execution does not imply production validity.'}
    else:
        data=json.loads(path.read_text());data.update(status='completed' if a.returncode==0 else 'failed',returncode=a.returncode,finished_at=now())
        data['outputs']=[ref(x) for x in a.output if Path(x).is_file()]
    path.write_text(json.dumps(data,indent=2)+'\n')
if __name__=='__main__':main()
