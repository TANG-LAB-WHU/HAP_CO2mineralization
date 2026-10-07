"""Package a complete, provenance-linked Overleaf source project after actual computation."""
import argparse,hashlib,json,re,shutil,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();r=a.run.resolve();d=ROOT/'drafts/local-validation';v=json.loads((r/'verification.json').read_text())
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

On the configured Shawn-PC checkout, regenerate the calculations/figures using:

```powershell
cd D:\\project\\HAP_CO2mineralization
.\\workflow\\run-refinement.ps1 -Stage all -Run refine-001
.\\workflow\\run-refinement.ps1 -Stage plots -Run refine-001
python workflow/build_validation_paper.py --run .research-runs/refine-001
```

These simulation commands depend on the pinned local runtime and retained raw outputs,
documented in `workflow/REFINEMENT.md`; the Overleaf project itself is self-contained for
document compilation, not a redistribution of model weights or all trajectories.
'''
    (d/'README.md').write_text(readme,encoding='utf-8')
    selected=['main.tex','numbers.tex','results.tex','references.bib','README.md','source_to_claim.json','verification.json','figure_sources.json','visual_qa.json']
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
        for name in ['verification.json','structure_audit.json','reference_audit.json','force_grid_audit.json','dft_relaxation.json','cp2k_restart_audit.json','aimd_audit.json','dataset_audit.json','sampling_coverage.json','model_evaluation.json','lammps_audit.json']:
            archive.write(r/name,'records/'+name)
        for name in ['optimizer_diagnosis-step90.json','dft_lbfgs_recovery.json','reference_interpretation_note.json','reproducibility-manifest.json']:
            if (r/name).exists():archive.write(r/name,'records/'+name)
        for name in ['plot_refinement.py','refine_pipeline.py','cp2k_validate.py','aimd_to_deepmd.py','audit_cp2k_restarts.py']:
            archive.write(ROOT/'workflow'/name,'scripts/'+name)
    overview=['# Local HAP validation deliverables','',
        '- HAP_local_validation_draft.pdf: compiled and visually checked working draft.',
        '- HAP_local_validation_overleaf.zip: self-contained XeLaTeX/BibTeX project.',
        '- HAP_LAMMPS_actual_trajectory.gif: actual saved frames, every 5 fs; no interpolation.',
        '- HAP_LAMMPS_OVITO_trajectory.extxyz: cell and chemical elements included for OVITO.',
        '- HAP_LAMMPS_OVITO_trajectory.dump: original unwrapped trajectory; types 1=Ca, 2=P, 3=O, 4=H, 5=C.',
        '- HAP_validation_diagnostics.zip: plots, quantitative records and generating scripts.','',
        'CO2 atom IDs are 45,46,47; its oxygen IDs are 46,47. The duration is 250 fs at dt=0.25 fs.',
        'This is a dry thin-slab and tiny-model pilot, not mineralization or production evidence.','',
        '| Pilot gate | Actual result |','|---|---|']
    overview.extend('| '+name+' | '+('PASS' if item['passed'] is True else 'FAIL')+' |' for name,item in v['gates'].items())
    (out/'README_local_results.md').write_text('\n'.join(overview)+'\n',encoding='utf-8')
    print('Created Overleaf ZIP, checked PDF, actual GIF, OVITO dump and diagnostic bundle:',out)
if __name__=='__main__':main()
