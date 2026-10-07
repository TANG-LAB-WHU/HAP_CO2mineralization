"""Package a complete, provenance-linked Overleaf source project after actual computation."""
import argparse,hashlib,json,re,shutil,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--recovery',type=Path,required=True);a=p.parse_args();r=a.run.resolve();d=ROOT/'drafts/local-validation';v=json.loads((r/'verification.json').read_text())
    if not v['functional_complete']:raise RuntimeError('Final package requires completed computational stages, even when scientific gates fail')
    if not (d/'main.pdf').exists():raise FileNotFoundError('Compile and inspect the real draft PDF first')
    text=(d/'main.tex').read_text(encoding='utf-8-sig');bib=(d/'references.bib').read_text(encoding='utf-8');keys=set(re.findall(r'@\w+\s*\{\s*([^,]+),',bib));cited=set()
    for group in re.findall(r'\\cite\w*\{([^}]+)\}',text):cited.update(group.split(','))
    if keys!=cited:raise ValueError('Missing or redundant bibliography entries')
    log=(d/'main.log').read_text(errors='replace')
    if re.search(r'Citation .* undefined|Reference .* undefined|There were undefined references|! LaTeX Error',log):raise ValueError('Unresolved LaTeX references/errors')
    qa=d/'visual_qa.json'
    if not qa.exists() or json.loads(qa.read_text()).get('pdf_sha256')!=sha(d/'main.pdf'):raise RuntimeError('Final PDF requires matching visual QA record')
    readme='''# Overleaf project: local HAP--CO2 validation draft

Import `HAP_local_validation_overleaf.zip` as a new Overleaf project. Select `main.tex`
as the main document and XeLaTeX as compiler. Use BibTeX (not Biber); Overleaf's normal
recompile resolves it automatically. All referenced figures and bibliography entries
are bundled. No login, upload, external image URL, custom font, or shell escape was used.

This is a scientific working draft about a bounded local workflow and small-model tests.
It does not establish a transferable reactive potential or CO2 mineralization. Read the
failed gates and limitations before using its results. Authorship, affiliations, funding,
and competing interests must be supplied by the research team. Raw data eligibility in
the original repository was not promoted; figure use is restricted to explicitly labeled
technical/negative-finding diagnostics requested for this draft.

`source_to_claim.json` maps numerical/figure claims to local source paths and checksums.
`verification.json` contains actual gate outcomes. `figure_sources.json` identifies the
plot inputs. `literature/` records primary-source access limits; `references.bib` contains
only entries cited by the manuscript. The separately delivered GIF shows actual saved
LAMMPS frames without interpolation. It is not a mineralization movie.

On the configured Shawn-PC checkout, verify the retained run and rebuild the document using:

```powershell
cd D:\\project\\HAP_CO2mineralization
.\\workflow\\run-refinement.ps1 -Stage all -Run refine-001
python workflow/build_validation_paper.py --run .research-runs/refine-001 --recovery .research-runs/recovery-001
```

A captured run is immutable: the all-stage command verifies hashes without rewriting it.
Use a new run ID for new calculations, or a separate review copy for figure edits.
These commands depend on the pinned local runtime and retained raw outputs,
documented in `workflow/REFINEMENT.md`; the Overleaf project itself is self-contained for
document compilation, not a redistribution of model weights or all trajectories.
'''
    (d/'README.md').write_text(readme,encoding='utf-8')
    selected=['main.tex','numbers.tex','results.tex','references.bib','README.md','source_to_claim.json','verification.json','figure_sources.json','visual_qa.json','recovery.tex','recovery_summary.json','recovery_sources.json','future_plan.tex','future_sampling_plan.md','future_sampling_profile.json']
    files=[d/x for x in selected]+sorted((d/'figures').glob('*.pdf'))+[d/'literature/evidence_notes.md',d/'literature/selection_notes.md']
    manifest={'files':[{'path':str(f.relative_to(d)).replace('\\','/'),'sha256':sha(f),'bytes':f.stat().st_size} for f in files],'compiled_pdf_sha256':sha(d/'main.pdf'),'cited_reference_count':len(cited),'raw_manuscript_eligible':False}
    (d/'package_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8');files.append(d/'package_manifest.json')
    out=ROOT/'deliverables';out.mkdir(exist_ok=True)
    with zipfile.ZipFile(out/'HAP_local_validation_overleaf.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for f in files:archive.write(f,str(f.relative_to(d)))
    shutil.copyfile(d/'main.pdf',out/'HAP_local_validation_draft.pdf')
    shutil.copyfile(r/'figures/08_lammps_actual.gif',out/'HAP_LAMMPS_actual_trajectory.gif')
    shutil.copyfile(r/'lammps/trajectory.dump',out/'HAP_LAMMPS_OVITO_trajectory.dump')
    shutil.copyfile(r/'lammps/trajectory.extxyz',out/'HAP_LAMMPS_OVITO_trajectory.extxyz')
    with zipfile.ZipFile(out/'HAP_validation_diagnostics.zip','w',zipfile.ZIP_DEFLATED) as archive:
        for f in sorted((r/'figures').glob('*')):
            if f.is_file():archive.write(f,'figures/'+f.name)
        archive.write(ROOT/'.research-runs/final-visual-review/visual_qa.json','visual_qa/animations_and_recovery_figures.json')
        archive.write(d/'visual_qa.json','visual_qa/final_pdf.json')
        for name in ['verification.json','structure_audit.json','reference_audit.json','force_grid_audit.json','dft_relaxation.json','cp2k_restart_audit.json','aimd_audit.json','dataset_audit.json','sampling_coverage.json','model_evaluation.json','lammps_audit.json']:
            archive.write(r/name,'records/'+name)
        for name in ['optimizer_diagnosis-step90.json','dft_lbfgs_recovery.json','optimizer_trial_audit.json','reference_interpretation_note.json','reproducibility-manifest.json','lammps_restart_diagnosis.json']:
            if (r/name).exists():archive.write(r/name,'records/'+name)
        for file in sorted((r/'source-snapshots').glob('*restart_boundary.py')):
            archive.write(file,'analysis_sources/'+file.name)
        for name in ['plot_refinement.py','refine_pipeline.py','cp2k_validate.py','aimd_to_deepmd.py','audit_cp2k_restarts.py']:
            archive.write(ROOT/'workflow'/name,'scripts/'+name)
    recovery=a.recovery.resolve();summary=json.loads((recovery/'recovery-summary.json').read_text())
    with zipfile.ZipFile(out/'HAP_validation_diagnostics.zip','a',zipfile.ZIP_DEFLATED) as archive:
        for name in ['recovery-summary.json','reproducibility-manifest.json','fit-diagnosis.json','dft-force-results.json','dynamics-units-audit.json','restart-control-results.json','restart-order-results.json','optimization-plan.json','optimization-continuation-plan.json','fullbatch-recovery-plan.json','fullbatch-evaluation.json','recovered-dynamics.json','09_recovery.pdf','09_recovery.png','technical_sources.json','numerical-followup-plan.json','numerical-followup-results.json','10_recovered_dynamics.pdf','10_recovered_dynamics.png','recovery-figure-manifest.json','final-instability-plan.json','final-instability-results.json','12_failure_dft.pdf','12_failure_dft.png','stable-demo-manifest.json']:
            if (recovery/name).is_file():archive.write(recovery/name,'recovery/'+name)
        for name in ['recovery_lbfgs.py','recovery_fit_diagnostics.py','recovery_evaluate.py','recovery_report.py','recovery_dynamics.py','recovery_numerical_audit.py','plot_recovery_dynamics.py','instability_audit.py','animate_stable_demo.py','capture_recovery.py']:
            archive.write(ROOT/'workflow'/name,'scripts/'+name)
        if (recovery/'lammps-recovered/trajectory.extxyz').exists():
            archive.write(recovery/'lammps-recovered/trajectory.extxyz','recovery/trajectory.extxyz')
    shutil.copyfile(recovery/'13_stable250fs_demo.gif',out/'HAP_stable_250fs_demonstration.gif')
    shutil.copyfile(recovery/'lammps-recovered/trajectory.extxyz',out/'HAP_stable_250fs_OVITO.extxyz')
    shutil.copyfile(recovery/'11_recovered_actual_1ps.gif',out/'HAP_recovered_actual_1ps.gif')
    shutil.copyfile(recovery/'lammps-recovered/trajectory-1ps.extxyz',out/'HAP_recovered_OVITO_1ps.extxyz')
    shutil.copyfile(recovery/'lammps-recovered/trajectory-1ps.dump',out/'HAP_recovered_OVITO_1ps.dump')
    overview=['# Local HAP validation deliverables','',
        '- HAP_local_validation_draft.pdf: compiled and visually checked working draft.',
        '- HAP_local_validation_overleaf.zip: self-contained XeLaTeX/BibTeX project.',
        '- HAP_LAMMPS_actual_trajectory.gif: original failed 250 fs run; actual saved frames every 5 fs; no interpolation.',
        '- HAP_LAMMPS_OVITO_trajectory.extxyz: cell and chemical elements included for OVITO.',
        '- HAP_LAMMPS_OVITO_trajectory.dump: original unwrapped trajectory; types 1=Ca, 2=P, 3=O, 4=H, 5=C.',
        '- HAP_validation_diagnostics.zip: plots, quantitative records and generating scripts.',
        '- HAP_stable_250fs_demonstration.gif: actual 51-frame short stable demonstration; same model fails near 395 fs. No long-time validity.',
        '- HAP_stable_250fs_OVITO.extxyz: all 101 actual frames from the separate stable short run.',
        '- HAP_recovered_actual_1ps.gif: separate 1 ps FAILURE animation, 101 actual frames every 10 fs; geometry and retrospective force gates fail.',
        '- HAP_recovered_OVITO_1ps.extxyz / .dump: all 401 saved recovery frames; unchanged atom IDs.','',
        'CO2 atom IDs are 45,46,47; its oxygen IDs are 46,47. Original and stable demonstrations cover 250 fs; the separate recovery-failure animation covers 1 ps. Both use dt=0.25 fs.',
        'This is a dry thin-slab and tiny-model pilot, not mineralization or production evidence.','',
        '| Pilot gate | Actual result |','|---|---|']
    overview.extend('| '+name+' | '+('PASS' if item['passed'] is True else 'FAIL')+' |' for name,item in v['gates'].items())
    sources={name:item['data'] for name,item in v['sources'].items()}
    model=sources['model_evaluation.json'];chosen=model['selected'];heldout=model['models'][chosen]['test'];md=sources['lammps_audit.json'];dataset=sources['dataset_audit.json']
    overview.extend(['','## Measured local outcomes','',
        f"- DFT maximum force: {sources['dft_relaxation.json']['fmax_eV_A']:.8f} eV/A; all optimizer criteria passed.",
        '- Dataset frames (training/validation/test): '+ '/'.join(str(dataset['frames_per_split'][x]) for x in ['train','valid','test'])+' from 4/2/2 separately seeded trajectories, each 4 fs. Statistical independence and equilibration are not established.',
        f"- Validation-selected model: {chosen}. Held-out force RMSE {heldout['force_rmse_eV_A']:.6f} eV/A; energy RMSE {heldout['energy_rmse_eV_atom']*1000:.6f} meV/atom.",
        f"- Python/LAMMPS same-state maximum force difference: {md['engine_force_max_difference_eV_A']:.6g} eV/A.",
        f"- NVE temperature: {md['temperature_K']['initial']:.1f} K initially, {md['temperature_K']['final']:.1f} K finally, maximum {md['temperature_K']['maximum']:.1f} K. Geometry failures: {md['geometry_failed_frame_count']}/101 saved frames. This does not establish physically valid heating or chemistry.",
        f"- Binary-restart maximum position difference: {md['restart_position_max_difference_A']:.6g} A; velocity difference: {md['restart_velocity_max_difference_A_per_ps']:.6g} A/ps.",
        f"- 250 fs NVE drift slopes: {md['energy_drift_meV_atom_ps']:.6g} (dt=0.25 fs) and {md['half_timestep_drift_meV_atom_ps']:.6g} (dt=0.125 fs), in meV/atom/ps. These are short-run diagnostics.",
        '', '## Verify on Shawn-PC','',
        'Raw inputs, logs, models, datasets and trajectories remain under `.research-runs/'+r.name+'` in `D:\\project\\HAP_CO2mineralization`.',
        '```powershell','cd D:\\project\\HAP_CO2mineralization',
        'python workflow/capture_refinement.py --run .research-runs/'+r.name+' --verify',
        '.\\workflow\\run-local.ps1 -Stage check','```',
        'The first command rechecks retained artifact hashes and gate-report sources; it does not grant scientific eligibility. The second runs static/evidence/integration contracts, not new GPU physics.'])
    overview.extend(['','## Separate bounded recovery','',
        f"The original one-frame force error {summary['overfit_original_eV_A']:.6f} eV/A was reduced to {summary['overfit_final_recovery_eV_A']:.6f} eV/A with the same architecture, cutoff and labels. Its unchanged recovery sanity gate {'passed' if summary['overfit_recovery_gate_pass'] else 'failed'}. Original-run gates above remain unchanged.",
        summary['limitation'],'',
        'Recovery raw data and controls: `.research-runs/recovery-001`. Verify using `python workflow/capture_recovery.py --run .research-runs/recovery-001 --verify`. See `workflow/RECOVERY.md` for exact training/restart limits and the retrospective-test caveat.',
        'Four DFT finite-difference calculations, model derivative checks, kinetic-energy reconstruction and repeated/order-controlled LAMMPS runs are retained. No new production-valid potential is claimed.'])
    if summary['full_data_recovery_executed']:
        e=summary['sources']['fullbatch-evaluation.json']['data'];overview.append(f"Full-data candidate validation force RMSE: {e['candidate']['valid']['force_rmse_eV_A']:.6f} eV/A; selected={e['recovery_selected']}.")
        if e['recovery_selected']:
            md=summary['sources']['recovered-dynamics.json']['data'];overview.append(f"Retrospective test force RMSE: {e['candidate']['test']['force_rmse_eV_A']:.6f} eV/A. Recovery trajectory: {md['geometry_failed_frame_count']}/{md['saved_frames']} geometry failures, maximum temperature {md['temperature_K']['maximum']:.1f} K. This does not override the failed original animation or establish physical validity.")
    n=summary['sources']['numerical-followup-results.json']['data'];long=n['one_ps']
    overview.append(f"Fixed-model numerical follow-up: strict restart gate={n['restart']['strict_original_gate_pass']}; 1 ps geometry failures {long['geometry_failed_frame_count']}/{long['saved_frames']}; T range {long['temperature_K']['minimum']:.1f}–{long['temperature_K']['maximum']:.1f} K; drift {long['energy_drift_meV_atom_ps']:.6g} meV/atom/ps. One ps is not equilibration or production convergence.")
    failure=summary['sources']['final-instability-results.json']['data'];halflong=failure['one_ps_halfstep'];probe=failure['first_failure_DFT_probe']
    overview.append(f"Full 1 ps half-timestep repeat: {halflong['geometry_failed_frame_count']}/401 geometry failures, starting at {halflong['first_geometry_failure_time_fs']} fs; drift {halflong['energy_drift_meV_atom_ps']:.6g} meV/atom/ps. A smaller timestep does not restore structural validity.")
    if probe['status']=='completed':overview.append(f"First-failure DFT check: force RMSE {probe['force_component_RMSE_eV_A']:.6f} eV/A; DFT energy change {probe['DFT_energy_change_from_start_eV']:.6f} eV versus learned-model change {probe['model_energy_change_from_start_eV']:.6f} eV. This point is diagnostic only and was not fitted.")
    shutil.copyfile(ROOT/'workflow/NEXT_SAMPLING_PLAN.md',out/'NEXT_SAMPLING_PLAN.md')
    shutil.copyfile(ROOT/'workflow/profiles/next-sampling-proposal.json',out/'next-sampling-proposal.json')
    overview.extend(['','## Future work only','', 'NEXT_SAMPLING_PLAN.md and next-sampling-proposal.json define a disabled, unexecuted campaign capped at 159 new DFT attempts, with fixed gates, targeted distortions, at most two acquisition rounds and a genuinely new sealed holdout. The 0.10788 eV/A reused-test result remains retrospective/exploratory.'])
    (out/'README_local_results.md').write_text('\n'.join(overview)+'\n',encoding='utf-8')
    print('Created Overleaf ZIP, checked PDF, actual GIF, OVITO dump and diagnostic bundle:',out)
if __name__=='__main__':main()
