"""Host-side local Docker execution with immutable images, logs, hashes and resumption."""
import argparse,datetime,hashlib,json,os,shutil,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DP='deepmodeling/deepmd-kit@sha256:23d23b922459eb49b23f9d4977b6e453fd17feb2ee03a280c70fb7f4a6452e8c'
CP='cp2k/cp2k@sha256:e2912cc9500a42a08b3e057e383182b96f395b7477c8eb55e93db6e29cb1c2ae'
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def execute(args,log,timeout=1800):
    if log.exists():raise FileExistsError(f'Preserving execution log: {log}')
    log.parent.mkdir(parents=True,exist_ok=True);start=time.time()
    manifest={'command':args,'started_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'running','manuscript_eligible':False}
    mp=log.with_suffix('.execution.json');save(mp,manifest)
    try:
        with log.open('w',encoding='utf-8') as f:q=subprocess.run(args,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,timeout=timeout)
        manifest.update(returncode=q.returncode,status='completed' if q.returncode==0 else 'failed')
        if q.returncode:raise RuntimeError(f'Command failed ({q.returncode}); see {log}')
    except Exception as e:
        manifest['status']='failed';manifest['error']=repr(e)
        if isinstance(e,subprocess.TimeoutExpired) and '--name' in args:
            # Stop only this invocation's explicitly named container, never other workloads.
            subprocess.run(['docker','stop','-t','10',args[args.index('--name')+1]],capture_output=True,timeout=30)
        raise
    finally:manifest['elapsed_seconds']=time.time()-start;manifest['log_sha256']=sha(log);save(mp,manifest)
def common():
    if shutil.disk_usage(ROOT).free<40*2**30:raise RuntimeError('Less than 40 GiB free')
    return ['docker','run','--rm','--name',f'hap-refine-{os.getpid()}','--pull','never','--gpus','device=0','--cpus','4','--memory','10g','--mount',f'type=bind,source={ROOT},target=/work','-w','/work','-e','OMP_NUM_THREADS=4','-e','DP_INTRA_OP_PARALLELISM_THREADS=4','-e','DP_INTER_OP_PARALLELISM_THREADS=1']
def cp_jobs(r,jobs):
    from cp2k_validate import parse
    for job in jobs:
        d=r/job['path'];out=d/job['output'];mp=d/'job.json'
        if mp.exists():
            old=json.loads(mp.read_text());assert old['input_sha256']==sha(d/job['input']),'Input changed after recorded calculation'
            if old['status']=='completed':
                current=parse(out);assert current['output_sha256']==old['result']['output_sha256'],'Completed output changed'
                continue
            raise RuntimeError(f'Prior failed/incomplete job retained: {d}')
        save(mp,{'status':'running','input_sha256':sha(d/job['input'])})
        try:
            execute(common()+['-w','/work/'+d.relative_to(ROOT).as_posix(),CP,'cp2k','-i',job['input'],'-o',job['output']],d/'console.log',timeout=2400)
            result=parse(out);save(mp,{'status':'completed','input_sha256':sha(d/job['input']),'result':result})
            print('Completed',job['path'],flush=True)
        except Exception as e:save(mp,{'status':'failed','input_sha256':sha(d/job['input']),'error':repr(e)});raise
def main():
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['all','references','aimd','dft','fragments','train','lammps']);p.add_argument('--run',type=Path,required=True);a=p.parse_args();r=a.run.resolve()
    if a.stage=='all':
        r.mkdir(parents=True,exist_ok=True)
        def stage(name,gpu=False,allow_failed=False):
            marker=r/(name+'.manifest.json')
            if marker.exists():
                state=json.loads(marker.read_text())['status']
                if state=='completed' or (allow_failed and state=='failed'):return
                raise RuntimeError(f'Prior incomplete/failed stage: {name}')
            command=['workflow/refine_pipeline.py',name,'--run',r.relative_to(ROOT).as_posix()]
            try:
                if gpu:execute(common()+['--entrypoint','/work/.research-runs/venv/bin/python',DP]+command,r/(name+'.console.log'),timeout=7200)
                else:execute([sys.executable]+command,r/(name+'.console.log'),timeout=1800)
            except RuntimeError:
                if not allow_failed:raise
        def jobs(key):cp_jobs(r,json.loads((r/key).read_text()))
        stage('prepare');stage('mace',True,True);stage('mace-continue',True);stage('cp-prepare');jobs('reference_jobs.json');stage('references-audit')
        stage('fragment-mace',True);stage('fragments-prepare');jobs('fragment_jobs.json');stage('fragments-audit')
        stage('dft-prepare');jobs('dft_jobs.json');stage('dft-audit');stage('aimd-prepare');jobs('aimd_jobs.json');stage('export');stage('train-prepare')
        for label in ['overfit','compact','wider']:
            d=r/'models'/label
            record=d/'training_success.json'
            if record.exists():
                q=json.loads(record.read_text());assert q['input_sha256']==sha(d/'input.json') and q['model_sha256']==sha(d/'model.pth')
            else:
                execute(common()+['-w','/work/'+d.relative_to(ROOT).as_posix(),'--entrypoint','/bin/bash',DP,'-lc','set -e; dp --pt train input.json > train.log 2>&1; dp --pt freeze -o model.pth > freeze.log 2>&1'],d/'console.log',timeout=7200)
                save(record,{'input_sha256':sha(d/'input.json'),'model_sha256':sha(d/'model.pth'),'status':'completed'})
        stage('evaluate',True);stage('lammps-prepare')
        for label in ['mapping','continuous','first','restart','halfstep']:
            d=r/'lammps';log=d/(label+'.console.log')
            if log.exists():
                meta=json.loads(log.with_suffix('.execution.json').read_text());assert meta['status']=='completed' and meta['log_sha256']==sha(log)
            else:execute(common()+['-w','/work/'+d.relative_to(ROOT).as_posix(),'--entrypoint','lmp',DP,'-in',label+'.in','-log',label+'.log'],log,timeout=1800)
        stage('lammps-audit',True);stage('report')
        subprocess.run([sys.executable,'workflow/plot_refinement.py','--run',str(r),'--animate'],cwd=ROOT,check=True)
        return
    if a.stage=='train':
        for name in ['overfit','compact','wider']:
            d=r/'models'/name
            execute(common()+['-w','/work/'+d.relative_to(ROOT).as_posix(),'--entrypoint','/bin/bash',DP,'-lc','set -e; dp --pt train input.json > train.log 2>&1; dp --pt freeze -o model.pth > freeze.log 2>&1'],d/'console.log',timeout=7200)
    elif a.stage=='lammps':
        d=r/'lammps'
        for name in ['mapping','continuous','first','restart','halfstep']:
            execute(common()+['-w','/work/'+d.relative_to(ROOT).as_posix(),'--entrypoint','lmp',DP,'-in',name+'.in','-log',name+'.log'],d/(name+'.console.log'),timeout=1800)
    else:
        filename={'references':'reference_jobs.json','aimd':'aimd_jobs.json','dft':'dft_jobs.json','fragments':'fragment_jobs.json'}[a.stage]
        jobs=json.loads((r/filename).read_text());cp_jobs(r,jobs)
if __name__=='__main__':main()
