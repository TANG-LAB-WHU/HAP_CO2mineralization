"""Read-only status of recorded local refinement jobs; never launches a calculation."""
import argparse,datetime,json
from pathlib import Path

def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();r=a.run
    replicas=[]
    for directory in sorted((r/'aimd').glob('replica-*')):
        job=directory/'job.json';item={'replica':directory.name,'recorded_status':load(job)['status'] if job.exists() else 'queued'}
        energy=directory/'md-1.ener'
        if energy.exists():
            rows=[x.split() for x in energy.read_text().splitlines() if x.strip() and not x.lstrip().startswith('#')]
            if rows:
                item.update(last_step=int(rows[-1][0]),time_fs=float(rows[-1][1]),last_step_seconds=float(rows[-1][-1]))
                item['energy_file_modified_UTC']=datetime.datetime.fromtimestamp(energy.stat().st_mtime,datetime.timezone.utc).isoformat()
        replicas.append(item)
    models=[]
    for name in ['overfit','compact','wider']:
        directory=r/'models'/name;item={'model':name,'recorded_status':'completed' if (directory/'training_success.json').exists() else 'not completed'}
        curve=directory/'lcurve.out'
        if curve.exists():
            rows=[x.split() for x in curve.read_text().splitlines() if x.strip() and not x.lstrip().startswith('#')]
            if rows:item['last_reported_update']=int(float(rows[-1][0]))
        models.append(item)
    report={'checked_at_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'run':str(r),'replicas':replicas,'models':models,
        'note':'Recorded file state only; running status does not prove a live process. Check the original runner and resource monitor before any recovery.'}
    if (r/'lammps_audit.json').exists():report['lammps_audit_complete']=True
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
