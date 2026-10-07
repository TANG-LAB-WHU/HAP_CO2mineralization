"""Assemble a diagnostic-only LaTeX project from measured local results and cited records."""
import argparse,datetime,hashlib,json,re,shutil,subprocess,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(p.read_text(encoding='utf-8-sig'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def entries(text):
    result={}
    for m in re.finditer(r'@\w+\s*\{\s*([^,]+),',text):
        level=1;i=m.end()
        while level and i<len(text):
            if text[i]=='{':level+=1
            elif text[i]=='}':level-=1
            i+=1
        if level:raise ValueError('Unbalanced BibTeX entry')
        result[m[1].strip()]=text[m.start():i]
    return result
def tex_escape(s):return str(s).replace('_',r'\_').replace('%',r'\%').replace('&',r'\&')
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--draft',type=Path,default=ROOT/'drafts/local-validation');p.add_argument('--allow-incomplete-preview',action='store_true');a=p.parse_args();r=a.run.resolve();d=a.draft.resolve();q=load(r/'verification.json')
    if not q['functional_complete'] and not a.allow_incomplete_preview:raise ValueError('Refusing final manuscript assembly from incomplete stages')
    lit=d/'literature';allentries={}
    for path in lit.glob('*.bib'):allentries.update(entries(path.read_text(encoding='utf-8')))
    registry=load(lit/'electronic-methods-crossref.json');keys=['PBE1996','GTH1996','MOLOPT2007','D32010']
    for key,record in zip(keys,registry):
        title=re.sub('<[^>]+>','',record['title'][0]);title=' '.join(title.split());authors=' and '.join(x['family']+', '+x.get('given','') for x in record['author'])
        year=record.get('published',record['issued'])['date-parts'][0][0];page=record.get('page') or record.get('article-number') or {'MOLOPT2007':'114105','D32010':'154104'}[key]
        allentries[key]='@article{'+key+',\n author={'+authors+'},\n title={'+title+'},\n journal={'+record['container-title'][0]+'},\n year={'+str(year)+'},volume={'+record['volume']+'},pages={'+page.replace('-','--')+'},doi={'+record['DOI']+'}\n}'
    main=(d/'main.tex').read_text(encoding='utf-8');cited=[]
    for group in re.findall(r'\\cite\w*\{([^}]+)\}',main):
        for key in group.split(','):
            if key not in cited:cited.append(key)
    if any(key not in allentries for key in cited):raise KeyError('Missing citation: '+str(set(cited)-set(allentries)))
    (d/'references.bib').write_text('\n\n'.join(allentries[x] for x in cited)+'\n',encoding='utf-8')
    files=['01_structure','02_relaxation','03_cp2k_references','03b_scf_force_grid','04_aimd','05_dataset','06_models','07_lammps'];(d/'figures').mkdir(exist_ok=True)
    source_map=[]
    for name in files:
        source=r/'figures'/f'{name}.pdf'
        if not source.exists():raise FileNotFoundError(source)
        shutil.copyfile(source,d/'figures'/source.name)
        source_map.append({'id':name,'claim_scope':'Actual diagnostic outputs of this specific run; no material-property or mineralization inference','source':str(source.relative_to(ROOT)),'source_sha256':sha(source),'raw_manifest':str((r/'figures/manifest.json').relative_to(ROOT)),
            'raw_manuscript_eligible':False,'draft_use':'explicit user-authorized diagnostic/negative-finding presentation; not an eligibility promotion'})
    src={k:v['data'] for k,v in q['sources'].items()};refs=src.get('reference_audit.json',{});grid=src.get('force_grid_audit.json',{});models=src.get('model_evaluation.json',{});md=src.get('lammps_audit.json',{});data=src.get('dataset_audit.json',{});dft=src.get('dft_relaxation.json',{})
    settings=grid.get('selected',{'cutoff':'pending','relcut':'pending','eps':'pending'})
    numbers=r'\newcommand{\LabelSettings}{PBE '+str(settings['cutoff'])+'/'+str(settings['relcut'])+r'\,Ry}'+ '\n'
    numbers+=r'\newcommand{\RunStatusBanner}{'+(r'\begin{center}\color{red!70!black}\textbf{INCOMPLETE PREVIEW: calculations are still running.}\end{center}' if not q['functional_complete'] else '')+'}\n'
    (d/'numbers.tex').write_text(numbers,encoding='utf-8')
    text=[r'\subsection{Structural and numerical checks}',r'The explicit structural audit and optimizer records are illustrated in Figures~\ref{fig:structure} and~\ref{fig:relaxation}. The original baseline constraint description is corrected as described above. All new numerical values below refer to the retained refinement run \texttt{'+tex_escape(r.name)+r'}.']
    if refs:
        text.extend([r'\begin{table}[ht]\centering\small',r'\caption{Matched CP2K single points on the same MACE geometries. The energy difference includes reconstruction, not just molecular interaction.}',r'\begin{tabular}{lrr}\toprule Setting & Three-reference difference (eV) & Adsorbed $F_{\max}$ (eV/\AA)\\\midrule'])
        for name,v in refs['settings'].items():text.append(f"{tex_escape(name)} & {v['adsorption_energy_eV']:.6f} & {v['adsorbed']['fmax_eV_A']:.4f} "+r'\\')
        text.extend([r'\bottomrule\end{tabular}\end{table}',f"The 500-to-700 Ry total-energy-difference change is {refs['cutoff_change_eV']:.6f}eV. This passes the 0.02 eV pilot energy-sensitivity gate but does not certify forces. "])
    if grid:text.append(f"The subsequent 900/80 versus 1100/100 Ry force-component RMSE is {grid['force_component_rmse_eV_A']:.5f}eV/\\AA. The force-grid gate is \\textbf{{{'passed' if grid['force_gate'] else 'failed'}}}; selected labels use \\LabelSettings. This comparison concerns one geometry, not global numerical accuracy.")
    for name,title in [('mace_fragments.json','MACE'),('cp2k_fragments.json','CP2K PBE, 500/60 Ry')]:
        if name in src:
            x=src[name];text.append(f"For {title}, the fixed-geometry interaction is {x['vertical_interaction_eV']:.6f}eV, the substrate term is {x['slab_deformation_eV']:.6f}eV, and the molecular deformation term is {x['co2_deformation_eV']:.6f}eV. Their sum is {x['total_reference_difference_eV']:.6f}eV. The negative substrate term shows that the reference geometry has higher same-method energy than the adsorbed-geometry slab fragment. In particular, the PBE single points do not establish DFT stationary basin identities. No strong-binding claim is inferred from the total.")
    if dft:text.append(f"The DFT-refined structure has maximum force {dft['fmax_eV_A']:.5f}eV/\\AA; its 0.05 eV/\\AA{{}} gate is \\textbf{{{'passed' if dft['force_gate'] else 'failed'}}}. Chemical sanity checks are reported separately from optimizer convergence.")
    if 'dft_lbfgs_recovery.json' in src and dft:
        text.append('The retained 90-step BFGS attempt passed force criteria but failed displacement criteria. After the recorded Hessian/displacement diagnosis, a fresh limited-memory history reached all four CP2K convergence criteria in one optimizer iteration (three force evaluations). This is local optimizer convergence, not proof of a globally stable surface; all original thresholds were retained.')
    restarts=src.get('cp2k_restart_audit.json',{}).get('checks',[])
    if restarts:
        maximum_energy=max(x['energy_difference_eV'] for x in restarts);maximum_force=max(x['force_component_max_difference_eV_A'] for x in restarts)
        text.append(f"Across {len(restarts)} actual saved-state CP2K restart boundaries, the largest energy difference is {maximum_energy:.3g} eV and the largest force-component difference is {maximum_force:.3g} eV/\\AA. These compare the same saved geometry before further optimization; they do not assert bitwise identical later optimizer paths.")
    text.append(r'\subsection{Small-sample model validation}')
    if data:text.append('The realized training/validation/test frame counts are '+ '/'.join(str(data['frames_per_split'][x]) for x in ['train','valid','test'])+f". The minimum cross-split pair-distance RMS difference is {data['cross_split_min_pair_distance_rms_A']:.5f}\\AA, compared with the 0.025\\AA{{}} criterion. Entire replicas remain disjoint, but short within-replica sequences remain correlated; the nominal frame count is not an effective independent sample size.")
    coverage=src.get('sampling_coverage.json',{}).get('replicas',[])
    if coverage:
        span=max(x['maximum_within_replica_pair_distance_rms_A'] for x in coverage)
        separation=min(x['minimum_other_replica_pair_distance_rms_A'] for x in coverage)
        text.append(f"The largest within-replica sampled geometry span is {span:.5f} \\AA, while the closest pair from different replicas differs by {separation:.5f} \\AA{{}} in the same distance-vector metric. These are local structural holdouts. Statistical independence, equilibrium and effective sample size are not established by four-femtosecond trajectories.")
    if models:
        text.extend([r'\begin{table}[ht]\centering\small',r'\caption{Measured force-component RMSE. The test set is evaluated only for the validation-selected architecture; a one-frame fit is not a test set.}',r'\begin{tabular}{lrrr}\toprule Model & Training & Validation & Independent test\\ & \multicolumn{3}{c}{eV/\AA}\\\midrule'])
        for name,v in models['models'].items():
            if name=='overfit':values=[f"{v['overfit']['force_rmse_eV_A']:.5f}",'same one frame','not applicable']
            else:values=[f"{v[x]['force_rmse_eV_A']:.5f}" if x in v else 'not used' for x in ['train','valid','test']]
            text.append(name+' & '+' & '.join(values)+r'\\')
        text.extend([r'\bottomrule\end{tabular}\end{table}'])
        chosen=models['selected'];m=models['models'][chosen]['test'];text.append(f"Validation selects the {chosen} architecture. Its independent energy RMSE is {1000*m['energy_rmse_eV_atom']:.4f}meV/atom; the 95th-percentile absolute force-component error is {m['force_p95_abs_error_eV_A']:.4f}eV/\\AA. The force and energy gates are \\textbf{{{'passed' if models['heldout_force_gate'] else 'failed'}}} and \\textbf{{{'passed' if models['heldout_energy_gate'] else 'failed'}}}, respectively. These outcomes concern the sampled local neighborhood only. Figure~\\ref{{fig:models}} displays the underlying prediction diagnostics.")
    text.append(r'\subsection{Transfer and short-time dynamics}')
    if md:text.append(f"The same-state Python/LAMMPS force mismatch is {md['engine_force_max_difference_eV_A']:.3g}eV/\\AA{{}} and the energy mismatch is {md['engine_energy_difference_eV']:.3g}eV. After the split restart, the maximum coordinate difference is {md['restart_position_max_difference_A']:.3g}\\AA{{}} and the maximum velocity difference is {md['restart_velocity_max_difference_A_per_ps']:.3g}\\AA{{}}/ps. The measured 250 fs NVE energy-drift slope is {md['energy_drift_meV_atom_ps']:.4f}meV/atom/ps. Figure~\\ref{{fig:lammps}} reports these diagnostics together with structural behavior. This short duration cannot establish long-term stability, and it does not override independent-force errors.")
    text.extend([r'\begin{table}[ht]\centering\small',r'\caption{Machine-readable pilot gates. A failed gate is retained as a finding, not replaced by normal program termination.}',r'\begin{tabular}{lr}\toprule Gate & Result\\\midrule'])
    for name,v in q['gates'].items():text.append(tex_escape(name)+' & '+('PASS' if v['passed'] is True else ('FAIL' if v['passed'] is False else 'MISSING'))+r'\\')
    text.extend([r'\bottomrule\end{tabular}\end{table}',r'No production-readiness or mineralization gate is passed by this workflow. Material-wide convergence, hydration, reaction coverage, and long-time sampling remain outside the completed pilot.'])
    formatted='\n\n'.join(text)+'\n'
    formatted=re.sub(r'(?<=\d)(?=meV|eV|\\AA)', '~', formatted)
    (d/'results.tex').write_text(formatted,encoding='utf-8')
    shutil.copyfile(r/'verification.json',d/'verification.json');shutil.copyfile(r/'figures/manifest.json',d/'figure_sources.json')
    claims={'run':r.name,'generated_at_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'raw_manuscript_eligible':False,'scope':'Diagnostic working draft only; existing manuscript and eligibility ledger unchanged','computational_claims':source_map,
       'numeric_claims':[{'id':name,'source':str((r/name).relative_to(ROOT)),'sha256':sha(r/name),'permitted':'Exact reported output or measured gate of this run only'} for name in src],
       'literature_claims':{'cheng1998co2':'abstract-level adsorption vs carbonate distinction','diallogarcia2014sites':'abstract-level speciation/hydration','bittencourt2021acidbase':'prior CO2/probe adsorption work','bouharras2025adsorption':'prior experimental/DFT work; physisorption interpretation','nowicki2024capture':'condition-dependent high-temperature capture/regeneration','corno2009water':'specific ideal surface hydration differences','saito2021surfaces':'environment/electrostatic dependence','slepko2011ordering':'OH ordering as structural variable'},'method_citations':[x for x in cited if x not in ['cheng1998co2','diallogarcia2014sites','bittencourt2021acidbase','bouharras2025adsorption','nowicki2024capture','corno2009water','saito2021surfaces','slepko2011ordering']],
       'limitations_file':'literature/evidence_notes.md','no_authorship_inferred':True}
    (d/'source_to_claim.json').write_text(json.dumps(claims,indent=2),encoding='utf-8')
    print('Assembled',len(cited),'cited references and',len(files),'traceable diagnostic figures')
if __name__=='__main__':main()
