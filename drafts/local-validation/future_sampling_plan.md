# Proposed next sampling campaign — not executed

This plan defines future work. No calculations, active-learning loop, HPC connection,
or submission are started by this document. The current potential is unsuitable for
production. The 0.10788 eV/A reused-test result is retrospective/exploratory.

## Scope and bounded budget

Keep the initial campaign restricted to the same dry, ordered, neutral 47-atom
HAP(001)+CO2 cell. Do not infer hydrated, reactive, variable-composition or bulk
accuracy. Use the existing PBE/DZVP-MOLOPT-SR-GTH/GTH labels and initially 1100/100 Ry,
EPS_SCF=1e-8. First verify force-grid sensitivity on three deliberately distinct
configurations at 900/80 and 1100/100 Ry (six single points). Stop if their force
component RMSE exceeds 0.02 eV/A; revise the electronic-method study as a new version.
This is not a substitute for basis, surface, electrostatics or finite-size convergence.

The proposed maximum is **159 new energy/force calculations**:

- 6 numerical grid checks;
- 48 initial training configurations, 24 validation configurations and 24 sealed test
  configurations, all assigned before DFT labeling;
- at most two active-learning rounds of 24 training labels each;
- 9 preselected DFT spot checks on the final model's dynamical trajectories.

Failed or nonconverged calculations remain in the manifest and consume an attempted
calculation slot. No silent extra round or fabricated replacement labels are allowed.
If the budget is inadequate, stop and issue a new prospective plan. Existing old test
and failure-selected configurations may be identified as development data, but can
never become part of the new independent test set.

## Initial geometry design and split isolation

Construct the 48 training candidates as three balanced groups of 16:

1. OH restoring-force probes: vary each hydroxyl bond between 0.85 and 1.30 A and
   angular orientation within 15 degrees, including coupled and separate distortions.
2. Phosphate/Ca environment probes: P–O bond perturbations of 0.05–0.15 A in both
   directions, rotations up to 10 degrees, and Ca displacements up to 0.10 A.
3. CO2/interface probes: azimuth, lateral registry and height variations around the
   adsorption region; C–O bond distortions of 0.03–0.15 A and bending down to 160 degrees.

These are force-label probes, not asserted equilibrium configurations. Preserve atom
identity, cell, charge and all-mobile constraints. Reject numerical overlaps below
0.70 A, but do not discard deliberate OH elongations merely because they cross the
later dynamics geometry gate. Record the accepted design variables and every rejection.

Generate validation and sealed test geometries with separate design IDs and seeds
(training 20261021, validation 20261022, test 20261023), distinct parent perturbations
and CO2 registries. Assign whole parent families to one split. Require cross-split
pair-distance-vector RMS separation at least 0.025 A, as in the existing pilot. Resolve
duplicates before labeling using the predeclared generation sequence, not label values.
This geometric check does not establish thermal statistical independence.

Freeze coordinate hashes, split membership, seeds and electronic inputs before any
fit. Store the 24 new test labels separately from all optimizer and acquisition inputs;
do not inspect them until model selection and numerical-stability checks are complete.
Their selection is independent of the old inspected trajectory; no frame from either
old failing 1 ps trajectory belongs to the new test set.

## Fitting and acquisition rules

Use four fixed seeds (1101–1104) with the existing wider descriptor/fitting architecture
and explicit Ca/P/O/H/C mapping. Predeclare 20,000 Adam updates per member, start LR
1e-3, stop LR 1e-5, decay interval 1000; retain the existing Adam loss: energy preference 0.1 to 1, force preference 100 to 10, and virial preference zero. The L-BFGS objective is force-component MSE plus per-atom energy MSE with equal numerical weights in eV/A and eV/atom units.
Allow one training-only L-BFGS continuation per member, capped at 6000 objective
evaluations or 30 minutes, whichever occurs first. Do not add architectures after
seeing the test. Select the candidate with the lowest full-validation force RMSE;
report energy and per-element errors and the four-member spread separately.

Acquire only from development trajectories and designed perturbations. Flag a candidate
when the maximum atom force standard deviation across the committee is at least
0.15 eV/A, or when it approaches/exceeds a predeclared structural boundary. Committee
disagreement is a heuristic and may miss shared model error. Include diverse restoring
directions around OH, phosphate and CO2 distortions rather than only the single worst
frame. Select at most 24 per round with pair-distance RMS separation at least 0.025 A;
retain a manifest explaining each selection. Never acquire from the sealed test.

## Frozen acceptance gates and stop conditions

| Gate | Requirement, fixed before new labels |
|---|---|
| One-frame implementation sanity | Force-component RMSE <0.020 eV/A |
| New sealed-test force | RMSE <0.100 eV/A |
| New sealed-test energy | RMSE <0.010 eV/atom, no test-fitted offset |
| Python/LAMMPS same state | Max force difference <1e-5 eV/A; energy difference <1e-5 eV |
| Strict 250 fs restart | Max position <1e-8 A and velocity <1e-8 A/ps |
| NVE drift | Absolute fitted slope <5 meV/atom/ps at both 0.25 and 0.125 fs |
| Geometry | No saved frame with minimum contact <=0.70 A, P–O coordination other than four within 1.9 A, or nearest O–H distance not strictly between 0.8 and 1.2 A |
| DFT spot checks | Force RMSE <0.100 eV/A on each preselected dynamical snapshot; report all results |

Test three new dynamical initial conditions (velocity seeds 2101–2103, initially
200 K), all mobile and NVE. Run both timesteps over identical 1 ps durations first.
Only if all development gates pass, extend each to 5 ps at both timesteps. The nine
DFT spot checks are the 0, 0.5 and 1 ps frames for each of the three dt=0.125 fs runs,
fixed before inspecting results. Report minimum-image trajectory differences without
inventing a posterior tolerance or expecting chaotic long paths to be identical.

Open the sealed test once, after selection. Any failure ends this campaign without
changing a threshold. Further development requires a new plan and a genuinely new
holdout. Passing these bounded gates would establish only a larger local validation
domain, not equilibration, mineralization, transferability or production convergence.

## HPC preparation

Use the disabled templates in `workflow/hpc/` and the provenance/restart manifests.
Resolve the site's scheduler/account/partition/QoS, modules or Apptainer/image paths,
hardware/resource limits and shared/scratch paths before adaptation. The budget and
gates above are scientific planning constraints; no cluster runtime is assumed or
submitted here.
