"""Capture final local pilot evidence without promoting publication eligibility."""
import argparse,csv,datetime,hashlib,importlib.metadata,json,platform,re,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
def cmd(args):
    q=subprocess.run(args,cwd=ROOT,capture_output=True,text=True,encoding='utf-8',errors='replace')
    return {'command':args,'returncode':q.returncode,'stdout':q.stdout,'stderr':q.stderr}
def append_row(path,row):
    with path.open(newline='',encoding='utf-8-sig') as stream:reader=csv.DictReader(stream);fields=reader.fieldnames;rows=list(reader)
    key=fields[0]
    if any(x[key]==row[key] for x in rows):return
    with path.open('a',newline='',encoding='utf-8') as stream:csv.DictWriter(stream,fieldnames=fields,lineterminator='\n').writerow(row)
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--verify',action='store_true');a=p.parse_args();r=a.run.resolve()
    if a.verify:
        record=json.loads((r/'reproducibility-manifest.json').read_text())
        for entry in record['artifacts']:
            file=(ROOT/entry['path']).resolve()
            if not file.is_relative_to(r):raise ValueError('Artifact path escapes the requested run')
            if not file.is_file() or file.stat().st_size!=entry['bytes'] or sha(file)!=entry['sha256']:
                raise ValueError(f'Captured artifact missing or changed: {file}')
        for entry in record.get('external_inputs',[]):
            file=(ROOT/entry['path']).resolve()
            if not file.is_relative_to(ROOT) or not file.is_file() or sha(file)!=entry['sha256']:
                raise ValueError(f'Captured source structure/model missing or changed: {file}')
        verification=json.loads((r/'verification.json').read_text())
        for name,item in verification['sources'].items():
            if sha(r/name)!=item['sha256'] or json.loads((r/name).read_text())!=item['data']:
                raise ValueError(f'Gate report source changed: {name}')
        print('Verified',len(record['artifacts']),'retained artifact hashes and',len(verification['sources']),'gate-report sources; scientific outcomes unchanged')
        return
    q=json.loads((r/'verification.json').read_text());assert q['functional_complete'],'Capture a terminal complete pilot, including failed scientific gates'
    if not re.fullmatch(r'[A-Za-z0-9_-]+',r.name):raise ValueError('Run name must be a simple identifier')
    tag=r.name.upper().replace('-','').replace('_','');identifier='LOCAL-HAP-'+tag
    target=ROOT/'workflow/validation'/f'refinement-{r.name}.json'
    record=ROOT/'workflow/validation/records'/f'{identifier}.json'
    manifest=r/'reproducibility-manifest.json'
    if any(p.exists() for p in [manifest,target,record]):raise FileExistsError('Preserving existing capture; review it instead of overwriting evidence')
    artifacts=[]
    for f in sorted(r.rglob('*')):
        if f.is_file() and f!=manifest:artifacts.append({'path':f.relative_to(ROOT).as_posix(),'bytes':f.stat().st_size,'sha256':sha(f),'manuscript_eligible':False})
    from run_refinement import DP,CP
    config_hash=json.loads((r/'mace.manifest.json').read_text())['config_sha256']
    snapshots=list((r/'source-snapshots').glob(config_hash+'-*.json'))
    if not snapshots or sha(snapshots[0])!=config_hash:raise ValueError('Missing original MACE configuration snapshot')
    config=json.loads(snapshots[0].read_text());model=ROOT/config['mace']['model']
    source=ROOT/'3.ML_AIMD/Step4_lammps_scaleup/HAP_ConventionalCell.cif'
    mace_record=json.loads((r/'mace.json').read_text());structure=json.loads((r/'structure_audit.json').read_text())
    if sha(model)!=mace_record['model_sha256']:raise ValueError('Pretrained MACE checkpoint changed since execution')
    if any(sha(source)!=x['source_sha256'] for x in structure['variants'].values()):raise ValueError('Starting CIF changed since preparation')
    external=[{'path':p.relative_to(ROOT).as_posix(),'sha256':sha(p),'role':role} for p,role in [(model,'pretrained MACE checkpoint'),(source,'starting crystallographic structure')]]
    result={'recorded_at_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'computer':platform.node(),'run_id':r.name,
        'git':cmd(['git','rev-parse','HEAD']),'branch':cmd(['git','branch','--show-current']),'working_changes':cmd(['git','status','--short']),
        'host_environment':{'python':sys.version,'executable':sys.executable,'platform':platform.platform(),
            'direct_packages':{name:importlib.metadata.version(name) for name in ['numpy','scipy','ase','pymatgen','matplotlib','Pillow']},
            'scope':'Actual host preparation/analysis environment; not a claim of clean-bootstrap reproducibility.'},
        'hardware':cmd(['nvidia-smi','--query-gpu=name,memory.total,driver_version','--format=csv']),
        'images':cmd(['docker','image','inspect',DP,CP,'--format','{{.Id}} {{json .RepoDigests}}']),
        'free_disk_bytes':shutil.disk_usage(ROOT).free,'run_bytes':sum(x['bytes'] for x in artifacts),'artifacts':artifacts,'external_inputs':external,
        'passed_all_pilot_gates':q['all_pilot_gates_pass'],'production_ready':False,'manuscript_eligible':False,
        'resource_note':'resource-samples.jsonl and resource-monitor.jsonl contain sampled Docker RAM and whole-GPU readings; the latter include display/other processes, not process-exclusive usage.'}
    save(manifest,result)
    shutil.copyfile(r/'verification.json',target)
    save(record,{'artifact_id':identifier,'local_raw_path':r.relative_to(ROOT).as_posix(),'status':'completed','manuscript_eligible':False,
        'notes':'Bounded structure/force-grid/small-model/engine validation; failed attempts and scientific gates retained. Not converged material predictions.',
        'reproducibility_manifest':manifest.relative_to(ROOT).as_posix(),'reproducibility_manifest_sha256':sha(manifest),'verification_report':target.relative_to(ROOT).as_posix()})
    append_row(ROOT/'literature/evidence.csv',{'evidence_id':'E-LOCAL-'+tag,'claim_id':'C-LOCAL-PILOT','status':'pending',
        'claim':'Local pilot diagnostics report structure checks, reference decomposition, short AIMD, compact model accuracy and engine transfer with explicit failed gates and limits.',
        'source_type':'computational_artifact','evidence_path':target.relative_to(ROOT).as_posix(),'citation_key':'','source_locator':'Actual gate measurements; diagnostic draft only',
        'run_id':r.name,'owner':'Shawn','verified_by':'','verified_at':'','notes':'Source path recorded for requested separate diagnostic draft. No material-property eligibility or human review inferred.'})
    append_row(ROOT/'research/artifact_inventory.csv',{'artifact_id':identifier,'stage':'analysis','system':'dry HAP001+CO2',
        'artifact_type':'local run record','path':record.relative_to(ROOT).as_posix(),'status':'completed','provenance_status':'path_only','manuscript_eligible':'no',
        'notes':'Actual bounded pilot and negative findings; see gate report. Thin dry slab, 4 fs replicas, small local model; not production evidence.'})
    for suffix,raw,status,note in [
        ('MACE-INITIAL','mace.json','partial','Initial 400-step adsorbed relaxation failed its force gate; later continuation is separate.'),
        ('D3-INPUT-FAIL','references/pbe-d3-500-attempt1','failed','Malformed first D3 section caused a parser failure; fixed input and rerun are separate.'),
        ('DFT-TIMEOUT','dft-relax-attempt1','failed','First DFT optimization hit 2400-second limit; original checkpoint, log and Hessian retained.'),
        ('DFT-STEP60','dft-relax-attempt2','partial','CP2K reached 60-step limit without geometry convergence; normal electronic termination is not an accepted minimum.'),
        ('DFT-STEP90','dft-relax-attempt3','partial','At 90 BFGS steps the force gate passed but displacement criteria failed; preserved separately from the successful limited-memory recovery.'),
        ('DFT-DIAGNOSIS','optimizer_diagnosis-step90.json','completed','Measured force/displacement/Hessian diagnosis preceding the bounded algorithm change; not proof of equilibrium.'),
        ('TRIAL-AUDIT-INITIAL-FIX','optimizer_trial_audit-attempt1.json','partial','First derived trial classification omitted the initial restart energy from its saved-frame comparison; corrected audit includes the previous saved state. Raw calculations were unchanged.'),
        ('DFT-LBFGS-RECOVERY','dft_lbfgs_recovery.json','completed','Fresh limited-memory history, unchanged electronic and convergence criteria; successful final audit is recorded separately.')]:
        if not (r/raw).exists():continue
        identifier='LOCAL-'+tag+'-'+suffix
        rec=ROOT/'workflow/validation/records'/f'{identifier}.json'
        save(rec,{'artifact_id':identifier,'local_raw_path':(r/raw).relative_to(ROOT).as_posix(),'status':status,'notes':note,
            'reproducibility_manifest':manifest.relative_to(ROOT).as_posix(),'manuscript_eligible':False})
        append_row(ROOT/'research/artifact_inventory.csv',{'artifact_id':identifier,'stage':'analysis','system':'dry HAP001+CO2',
            'artifact_type':'retained attempt or recovery diagnosis','path':rec.relative_to(ROOT).as_posix(),'status':status,'provenance_status':'path_only','manuscript_eligible':'no','notes':note})
    print('Captured',len(artifacts),'files and',result['run_bytes'],'bytes; all publication eligibility remains false')
if __name__=='__main__':main()
