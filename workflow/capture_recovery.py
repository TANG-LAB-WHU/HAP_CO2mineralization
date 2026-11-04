"""Freeze recovery diagnostics separately from the original failed pilot."""
import argparse,datetime,json,re,shutil,subprocess
from pathlib import Path
from capture_refinement import sha,save,append_row,ROOT
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--verify',action='store_true');a=p.parse_args();r=a.run.resolve()
    if not r.is_relative_to(ROOT/'.research-runs') or not re.fullmatch(r'[A-Za-z0-9_-]+',r.name):raise ValueError('Expected a named local recovery run')
    tag=r.name.upper().replace('-','').replace('_','');identifier='LOCAL-HAP-'+tag
    manifest=r/'reproducibility-manifest.json';anchor=ROOT/'workflow/validation/records'/f'{identifier}.json'
    target=ROOT/'workflow/validation'/f'recovery-{r.name}.json'
    if a.verify:
        registered=json.loads(anchor.read_text());assert sha(manifest)==registered['reproducibility_manifest_sha256'],'Recovery manifest changed'
        record=json.loads(manifest.read_text())
        for item in record['artifacts']:
            f=(ROOT/item['path']).resolve()
            if not f.is_relative_to(r) or not f.is_file() or f.stat().st_size!=item['bytes'] or sha(f)!=item['sha256']:raise ValueError(f'Changed recovery artifact: {f}')
        for item in record['external_inputs']:
            f=(ROOT/item['path']).resolve()
            if not f.is_relative_to(ROOT) or not f.is_file() or sha(f)!=item['sha256']:raise ValueError(f'Changed baseline input: {f}')
        summary=json.loads((r/'recovery-summary.json').read_text())
        for name,x in summary['sources'].items():
            assert sha(r/name)==x['sha256'] and json.loads((r/name).read_text())==x['data'],'Recovery summary source changed'
        print('Verified',len(record['artifacts']),'recovery artifact hashes; no original scientific gate changed');return
    if any(p.exists() for p in [manifest,anchor,target]):raise FileExistsError('Preserve prior recovery capture')
    summary=json.loads((r/'recovery-summary.json').read_text());assert summary['original_gates_changed'] is False and summary['production_ready'] is False
    snap=r/'source-snapshots';snap.mkdir(exist_ok=True)
    for name in ['recovery_fit_diagnostics.py','recovery_lbfgs.py','recovery_evaluate.py','recovery_report.py','recovery_dynamics.py','recovery_numerical_audit.py','plot_recovery_dynamics.py','instability_audit.py','animate_stable_demo.py','capture_recovery.py']:
        f=ROOT/'workflow'/name
        if f.exists():shutil.copyfile(f,snap/(sha(f)+'-'+name))
    artifacts=[{'path':f.relative_to(ROOT).as_posix(),'bytes':f.stat().st_size,'sha256':sha(f),'manuscript_eligible':False} for f in sorted(r.rglob('*')) if f.is_file() and f!=manifest]
    baseline=ROOT/'.research-runs/refine-001'
    inputs=[baseline/'reproducibility-manifest.json',baseline/'models/overfit/model.ckpt-1500.pt',baseline/'models/overfit/out.json',baseline/'models/wider/model.ckpt-6000.pt',baseline/'models/wider/out.json']
    for role in ['overfit','train','valid','test']:
        inputs.extend(sorted((baseline/'dataset'/role/'set.000').glob('*.npy')))
        inputs.append(baseline/'dataset'/role/'type.raw')
    from run_refinement import DP,CP
    record={'run_id':r.name,'recorded_at_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'images':{'deepmd':DP,'cp2k':CP},'artifacts':artifacts,'bytes':sum(x['bytes'] for x in artifacts),
        'external_inputs':[{'path':f.relative_to(ROOT).as_posix(),'sha256':sha(f)} for f in inputs],
        'original_gates_changed':False,'production_ready':False,'manuscript_eligible':False,
        'resource_scope':'Retained samples cover part of the DFT controls and optimizer continuation; whole GPU includes display/other workloads. All calculations used existing pinned runtimes on Shawn-PC.'}
    save(manifest,record);shutil.copyfile(r/'recovery-summary.json',target)
    save(anchor,{'artifact_id':identifier,'local_raw_path':r.relative_to(ROOT).as_posix(),'status':'completed','manuscript_eligible':False,
        'reproducibility_manifest':manifest.relative_to(ROOT).as_posix(),'reproducibility_manifest_sha256':sha(manifest),'verification_report':target.relative_to(ROOT).as_posix(),
        'notes':'Completed failure diagnosis and bounded recovery experiments. Original failures and all thresholds retained; no production-valid model claimed.'})
    append_row(ROOT/'literature/evidence.csv',{'evidence_id':'E-LOCAL-'+tag,'claim_id':'C-LOCAL-PILOT','status':'pending',
        'claim':'Separate local recovery diagnostics quantify force derivative consistency, DFT finite differences, optimizer response, energy/temperature units and repeat/restart sensitivity without altering original gates.',
        'source_type':'computational_artifact','evidence_path':target.relative_to(ROOT).as_posix(),'citation_key':'','source_locator':'Actual separate recovery measurements; diagnostic draft only',
        'run_id':r.name,'owner':'Shawn','verified_by':'','verified_at':'','notes':'No human verification, general material claim, or production eligibility inferred.'})
    append_row(ROOT/'research/artifact_inventory.csv',{'artifact_id':identifier,'stage':'analysis','system':'dry HAP001+CO2','artifact_type':'local recovery record',
        'path':anchor.relative_to(ROOT).as_posix(),'status':'completed','provenance_status':'path_only','manuscript_eligible':'no','notes':'Separate bounded optimizer and restart controls; retain original failures and retrospective-test caveat.'})
    failure=r/'fit-diagnosis-attempt1.log'
    if failure.exists():
        aid='LOCAL-'+tag+'-LOADER-DIAGNOSIS-FAIL';f=ROOT/'workflow/validation/records'/f'{aid}.json'
        save(f,{'artifact_id':aid,'local_raw_path':failure.relative_to(ROOT).as_posix(),'status':'failed','manuscript_eligible':False,
            'notes':'Diagnostic reconstruction initially used unresolved input defaults; strict state loading rejected it. Corrected helper uses original resolved out.json. Original training/checkpoints were unaffected.',
            'reproducibility_manifest':manifest.relative_to(ROOT).as_posix()})
        append_row(ROOT/'research/artifact_inventory.csv',{'artifact_id':aid,'stage':'analysis','system':'dry HAP001+CO2','artifact_type':'retained failed diagnostic helper',
            'path':f.relative_to(ROOT).as_posix(),'status':'failed','provenance_status':'path_only','manuscript_eligible':'no','notes':'Unresolved-default checkpoint reconstruction rejected; corrected resolved-config diagnostic retained separately.'})
    print('Captured',len(artifacts),'recovery files; scientific eligibility unchanged')
if __name__=='__main__':main()
