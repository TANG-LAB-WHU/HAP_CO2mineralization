"""Assemble a separate recovery record; never overwrite the original scientific gates."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args();r=a.run.resolve()
    if (r/'reproducibility-manifest.json').exists():raise RuntimeError('Captured recovery is immutable')
    required=['fit-diagnosis.json','dft-force-results.json','dynamics-units-audit.json','restart-control-results.json',
              'overfit-lbfgs/result.json','overfit-continuation/result.json','restart-order-results.json']
    records={name:load(r/name) for name in required}
    fitting=records['fit-diagnosis.json'];dft=records['dft-force-results.json'];units=records['dynamics-units-audit.json'];restart=records['restart-control-results.json']
    first=records['overfit-lbfgs/result.json'];last=records['overfit-continuation/result.json']
    full=None
    if last['original_force_gate_pass'] and dft['all_checks_pass']:
        full=load(r/'fullbatch-evaluation.json');required.append('fullbatch-evaluation.json');records[required[-1]]=full
        if full['recovery_selected']:
            records['recovered-dynamics.json']=load(r/'recovered-dynamics.json');required.append('recovered-dynamics.json')
            records['numerical-followup-results.json']=load(r/'numerical-followup-results.json');required.append('numerical-followup-results.json')
            records['final-instability-results.json']=load(r/'final-instability-results.json');required.append('final-instability-results.json')
    limitation=('One-frame recovery did not pass its unchanged gate within the stated budget; the conditional full-data fit was not executed. '
                'More DFT training examples cannot by themselves establish that this fitting setup reproduces even one existing label.') if not last['original_force_gate_pass'] else (
                'One-frame reproduction is an optimizer diagnosis, not a test of generalization. The original test has already been inspected; all recovery test results are retrospective diagnostics.')
    report={'run_id':r.name,'original_run':'refine-001','original_gates_changed':False,'production_ready':False,
        'overfit_original_eV_A':fitting['one_frame_force_rmse_eV_A'],'overfit_first_recovery_eV_A':first['final']['force_rmse_eV_A'],
        'overfit_final_recovery_eV_A':last['final']['force_rmse_eV_A'],'overfit_recovery_gate_pass':last['original_force_gate_pass'],
        'DFT_force_checks_pass':dft['all_checks_pass'],'full_data_recovery_executed':full is not None,'limitation':limitation,
        'sources':{name:{'sha256':sha(r/name),'data':records[name]} for name in required},'manuscript_eligible':False}
    (r/'recovery-summary.json').write_text(json.dumps(report,indent=2)+'\n')
    fmax=max(x['absolute_difference_eV_A'] for x in fitting['energy_force_finite_differences'] if x['delta_A']==1e-4)
    fdmax=max(abs(x['difference_from_original_label_eV_A']) for x in dft['components'])
    text=[r'\subsection{Bounded failure diagnosis and recovery}',
        'The original failed gates above are retained. A separate recovery campaign checks the force labels, model derivatives, optimization and restart behavior before considering further data collection. All atoms and Cartesian force components remain included; no constraints or force thresholds are changed.',
        f"The frozen and reconstructed trainable model agree to {fitting['frozen_vs_trainable_force_max_difference_eV_A']:.3g} eV/\\AA. The force-loss directional derivative agrees with a parameter finite difference to relative error {fitting['parameter_direction_derivative']['relative_difference']:.3g}; descriptor and fitting-network gradients are nonzero. Coordinate finite differences at $10^{{-4}}$ \\AA{{}} differ from analytic model forces by at most {fmax:.3g} eV/\\AA{{}} on the tested components. These checks disfavor a broken force-gradient implementation in the tested path; they are not an exhaustive code proof.",
        f"Four new PBE single points displace the original training frame by $\\pm0.001$ \\AA{{}} along P11-$y$ and H18-$z$. The cutoff and basis remain those of the labels; the SCF target is tightened to $10^{{-9}}$ for derivative noise. Central energy differences agree with the original labels within {fdmax:.3g} eV/\\AA, passing the predeclared 0.01 eV/\\AA{{}} check. This tests two force components, not global basis or force convergence.",
        f"Keeping the original resolved architecture, 5 \\AA{{}} cutoff, descriptor and training frame, full-batch L-BFGS first lowers force RMSE from {report['overfit_original_eV_A']:.5f} to {report['overfit_first_recovery_eV_A']:.5f} eV/\\AA{{}} in {first['closure_calls']} objective evaluations. A prospectively bounded continuation then reaches {report['overfit_final_recovery_eV_A']:.5f} eV/\\AA{{}} in {last['closure_calls']} additional evaluations. Its unchanged 0.02 eV/\\AA{{}} gate is \\textbf{{{'passed' if last['original_force_gate_pass'] else 'failed'}}}. The objective is force-component MSE plus per-atom energy MSE, in the original physical units. Only training labels enter optimization.",
        limitation,
        f"The saved velocities independently reconstruct LAMMPS kinetic energies within {units['kinetic_energy_max_reconstruction_error_eV']:.3g} eV and temperatures within {units['temperature_max_reconstruction_error_K']:.3g} K, using 138 degrees of freedom and the installed-version metal-unit constants. The run has no thermostat and uses a 0.00025 ps timestep. At peak temperature, model potential energy changes by {units['peak_temperature']['potential_change_eV']:.3f} eV while kinetic energy changes by {units['peak_temperature']['kinetic_change_eV']:.3f} eV. The heating is therefore associated with the learned potential-energy surface; it is not a hidden thermostat or temperature-unit conversion. A small fitted drift slope does not establish an adequate potential or a uniformly small instantaneous energy error."]
    for key,label in [('gpu_independent_repeat_vs_original','An independent GPU repeat without any restart'),('cpu_single_thread_continuous_vs_restart','The single-thread CPU continuous/restart control')]:
        x=restart[key];text.append(f"{label} gives endpoint differences of {x['position_max_difference_A']:.3g} \\AA{{}} in positions and {x['velocity_max_difference_A_ps']:.3g} \\AA/ps in velocity components. The original strict endpoint criterion is {'passed' if x['strict_original_gate_pass'] else 'failed'} in this control.")
    order=records['restart-order-results.json'];cpu_repeat=order['cpu_independent_uninterrupted_repeat'];ordered=order['cpu_no_sort_rebuild_every_step_restart'];growth=order['GPU_independent_repeat_difference_history']
    text.append(f"A second uninterrupted single-thread CPU run reproduces saved endpoint positions and velocities exactly (maximum differences {cpu_repeat['position_max_difference_A']:.3g} \\AA{{}} and {cpu_repeat['velocity_max_difference_A_ps']:.3g} \\AA/ps). Disabling atom sorting and rebuilding neighbors every step reduces the CPU restart discrepancies to {ordered['position_max_difference_A']:.3g} \\AA{{}} and {ordered['velocity_max_difference_A_ps']:.3g} \\AA/ps, but still fails the velocity criterion. The independent GPU pair starts with identical coordinates and velocities and a force difference of {growth[0]['force_max_difference_eV_A']:.3g} eV/\\AA; that small difference grows along the trajectory. This intervention and repeat evidence support amplification of roundoff/order differences, rather than loss of the saved state, as a contributor to the strict endpoint failure.")
    text.append('The repeat and CPU controls use byte-identical original LAMMPS physics inputs. Device and threading change jointly in the CPU control, so no single GPU kernel is identified as the cause. Numerical sensitivity is assessed without relaxing the original restart gate. The failed original trajectory remains visible in the animation and quantitative records.')
    if full:
        text.append(f"The conditional full-data continuation uses the original 16 training frames and wider architecture, with no additional test-label fitting. Its validation force RMSE is {full['candidate']['valid']['force_rmse_eV_A']:.5f} eV/\\AA; the recovery candidate is {'selected' if full['recovery_selected'] else 'not selected'} by the recorded validation-only comparison. No new prospectively untouched test set is claimed.")
        if full['recovery_selected']:
            md=records['recovered-dynamics.json'];test=full['candidate']['test']
            text.append(f"Retrospective test force RMSE is {test['force_rmse_eV_A']:.5f} eV/\\AA. The separately retained recovery trajectory has maximum temperature {md['temperature_K']['maximum']:.1f} K and {md['geometry_failed_frame_count']} failing geometry frames out of {md['saved_frames']}. These outcomes do not replace the original failures or establish a production potential.")
    if full and full['recovery_selected']:
        numerical=records['numerical-followup-results.json'];long=numerical['one_ps'];halved=numerical['half_timestep'];rs=numerical['restart']
        text.append(f"The full-data optimizer stops at its {full['training_fit']['closure_calls']}-evaluation budget, with training force RMSE {full['candidate']['train']['force_rmse_eV_A']:.5f} eV/\\AA. It has not reached its internal 0.02 eV/\\AA{{}} training target. The retrospective test force error remains above the original 0.1 eV/\\AA{{}} criterion; it is not relabeled as a pass.")
        text.append(f"With weights fixed, the recovered model's 500+500 restart gives maximum differences {rs['position_max_difference_A']:.3g} \\AA{{}} and {rs['velocity_max_difference_A_ps']:.3g} \\AA/ps; the unchanged strict gate is {'passed' if rs['strict_original_gate_pass'] else 'failed'}. The half-timestep run over 250 fs has drift {halved['energy_drift_meV_atom_ps']:.5g} meV/atom/ps and an endpoint minimum-image position RMS difference of {halved['endpoint_minimum_image_position_RMS_A']:.5g} \\AA. No post hoc coordinate-convergence threshold is assigned.")
        text.append(f"A separately predeclared 1 ps NVE run at 0.25 fs retains {long['saved_frames']} actual saved frames, of which {long['geometry_failed_frame_count']} fail the original geometry criteria. Its temperature range is {long['temperature_K']['minimum']:.1f}--{long['temperature_K']['maximum']:.1f} K, ending at {long['temperature_K']['final']:.1f} K; the total-energy drift is {long['energy_drift_meV_atom_ps']:.5g} meV/atom/ps. Figure~\\ref{{fig:recovered-dynamics}} compares predictions and dynamics. This bounded improvement does not establish equilibration, transferable forces or production convergence. Both original failed and recovery animations are retained.")
    if full and full['recovery_selected']:
        final=records['final-instability-results.json'];probe=final['first_failure_DFT_probe'];half_long=final['one_ps_halfstep']
        if probe['status']=='completed':
            text.append(f"A final post hoc CP2K point evaluates the first failing configuration at {probe['selected_frame']['time_fs']:.1f} fs. Relative to the same-method initial energy, DFT changes by {probe['DFT_energy_change_from_start_eV']:.4f} eV, whereas the learned model changes by {probe['model_energy_change_from_start_eV']:.4f} eV. Its force-component RMSE against this DFT point is {probe['force_component_RMSE_eV_A']:.4f} eV/\\AA. This directly measures a large local extrapolation error at the onset of failure; the point was selected for diagnosis and never used for fitting or presented as an unbiased test set.")
        else:text.append('The attempted DFT calculation at the first failed frame did not provide a converged diagnostic; its failure log is retained and no substitute labels are generated.')
        text.append(f"At half timestep over the same full 1 ps, {half_long['geometry_failed_frame_count']} of 401 saved frames fail geometry, beginning at {half_long['first_geometry_failure_time_fs']} fs. The maximum temperature is {half_long['temperature_K']['maximum']:.1f} K and drift is {half_long['energy_drift_meV_atom_ps']:.5g} meV/atom/ps. Thus timestep refinement does not restore geometric validity. The late energy-drift behavior and divergent paths are reported separately from the measured model error; energy conservation alone cannot validate the learned surface.")
        text.append(r'No production-ready learned potential was obtained. The implemented optimizer remedy improves local fitting but does not cure the narrow chemical coverage. Further fitting on these inspected test trajectories would not supply independent validation; a new, prospectively defined DFT sampling and holdout design is required. Figure~\ref{fig:failure-dft} reports the targeted DFT diagnostic.')
    text.append(r'Figure~\ref{fig:recovery} summarizes the measured controls. No hydration, reaction-path, broad-domain, or long-time validation was added by this bounded campaign.')
    (r/'recovery.tex').write_text('\n\n'.join(text)+'\n',encoding='utf-8')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(2,2,figsize=(10,7.8));fig.suptitle('recovery-001 | Separate diagnosis; original failures retained',x=.06,ha='left')
    h1=load(r/'overfit-lbfgs/history.json');h2=load(r/'overfit-continuation/history.json')
    x=[v['closure_calls'] for v in h1]+[h1[-1]['closure_calls']+v['closure_calls'] for v in h2[1:]]
    y=[v['force_rmse_eV_A'] for v in h1]+[v['force_rmse_eV_A'] for v in h2[1:]]
    ax[0,0].semilogy(x,y,label='Same training frame');ax[0,0].axhline(.02,color='#bd382b',ls='--',label='Unchanged sanity gate');ax[0,0].set(xlabel='Additional objective evaluations',ylabel='Force RMSE (eV/A)',title='Optimizer-only recovery');ax[0,0].legend(fontsize=8)
    labels=[f"{v['element']}{v['atom_id']}-{['x','y','z'][v['axis']]}" for v in dft['components']]
    ax[0,1].bar(labels,[abs(v['difference_from_original_label_eV_A']) for v in dft['components']],color=['#d68b26','#4a9c92']);ax[0,1].set_yscale('log');ax[0,1].axhline(.01,color='#bd382b',ls='--',label='Predeclared check');ax[0,1].set(ylabel='DFT force discrepancy (eV/A)',title='Actual central energy differences');ax[0,1].legend(fontsize=8)
    restart['cpu_ordered']=ordered
    keys=['gpu_original_continuous_vs_restart','gpu_independent_repeat_vs_original','cpu_single_thread_continuous_vs_restart','cpu_ordered'];short=['GPU restart','GPU repeat','CPU restart','CPU ordered']
    pos=[max(restart[k]['position_max_difference_A'],1e-16) for k in keys];vel=[max(restart[k]['velocity_max_difference_A_ps'],1e-16) for k in keys];xx=np.arange(4)
    ax[1,0].bar(xx-.17,pos,.34,label='Position (A)');ax[1,0].bar(xx+.17,vel,.34,label='Velocity (A/ps)');ax[1,0].set_xticks(xx,short,rotation=12);ax[1,0].set_yscale('log');ax[1,0].axhline(1e-8,color='#bd382b',ls='--');ax[1,0].set(ylabel='Max component difference',title='Strict 1e-8 gate retained');ax[1,0].legend(fontsize=8)
    baseline=ROOT/'.research-runs/refine-001/lammps/thermo.npy';t=np.load(baseline)
    ax[1,1].plot(t[:,1]*1000,t[:,3]-t[0,3],label='Potential');ax[1,1].plot(t[:,1]*1000,t[:,4]-t[0,4],label='Kinetic');ax[1,1].plot(t[:,1]*1000,t[:,5]-t[0,5],label='Total');ax[1,1].set(xlabel='Original trajectory time (fs)',ylabel='Energy change (eV)',title='Heating on the failed learned surface');ax[1,1].legend(fontsize=8)
    fig.text(.06,.02,'Actual diagnostic measurements. No test-set fitting, threshold relaxation, or production-validity claim.',fontsize=8)
    fig.tight_layout(rect=[0,.04,1,.94]);fig.savefig(r/'09_recovery.pdf');fig.savefig(r/'09_recovery.png',dpi=170);plt.close(fig)
    print('Created separate recovery summary, figure and draft section')
if __name__=='__main__':main()
