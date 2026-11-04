"""Read-only sampled resource monitor tied to one host orchestration process."""
import argparse,datetime,json,subprocess,time
from pathlib import Path
import psutil
def capture(args):
    try:
        p=subprocess.run(args,capture_output=True,text=True,encoding='utf-8',errors='replace',timeout=15)
        return {'returncode':p.returncode,'output':p.stdout.strip(),'error':p.stderr.strip()}
    except Exception as e:return {'error':repr(e)}
def main():
    p=argparse.ArgumentParser();p.add_argument('--pid',type=int,required=True);p.add_argument('--run',type=Path,required=True);p.add_argument('--minutes',type=int,default=180);a=p.parse_args()
    process=psutil.Process(a.pid);started=process.create_time()
    if 'run_refinement.py' not in ' '.join(process.cmdline()):raise ValueError('Expected the specific refinement orchestrator')
    path=a.run/'resource-monitor.jsonl';end=time.monotonic()+a.minutes*60
    while time.monotonic()<end:
        try:
            process=psutil.Process(a.pid)
            if process.create_time()!=started:break
        except psutil.NoSuchProcess:break
        q={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'host_orchestrator_pid':a.pid,
            'gpu_scope':'whole GPU including display and other workloads',
            'gpu':capture(['nvidia-smi','--query-gpu=memory.used,utilization.gpu,temperature.gpu,power.draw','--format=csv,noheader']),
            'containers':capture(['docker','stats','--no-stream','--format','{{.Name}} {{.MemUsage}} {{.CPUPerc}}'])}
        with path.open('a',encoding='utf-8') as f:f.write(json.dumps(q)+'\n')
        time.sleep(30)
    print('Resource monitor finished; retained',path)
if __name__=='__main__':main()
