# Recovery diagnostics on Shawn-PC

The local software chain runs, but **no production-ready HAP potential was obtained**.
Optimizer continuation repairs the one-frame fit and improves the first 250 fs. Both
1 ps trajectories fail structural validation; DFT at the first failed configuration
confirms a large model error. No threshold was relaxed and no inspected test frame
was added to training. The original `refine-001` failures remain immutable.

## Measured outcomes

| Check | Actual recovery result | Interpretation |
|---|---:|---|
| One-frame force RMSE | 0.016476 eV/A | Passes unchanged 0.020 sanity gate; same architecture/labels |
| Training force RMSE, 16 frames | 0.065398 eV/A | Budget ended before internal 0.020 training target |
| Validation force RMSE, 8 frames | 0.082894 eV/A | Improved from 0.137666; candidate selected before test read |
| Retrospective test force RMSE, 8 frames | 0.107883 eV/A | **FAIL** at unchanged 0.100 gate |
| Retrospective test energy RMSE | 0.787982 meV/atom | Pass at 10 meV/atom; does not override force failure |
| 250 fs geometry | 0/101 failed frames | Short-time pass only |
| 250 fs restart, max position / velocity | 1.39e-13 A / 1.04e-11 A/ps | Original strict 1e-8 criteria pass |
| 1 ps, dt=0.25 fs | 231/401 failed frames; first 395.0 fs | Geometry and drift fail |
| 1 ps, dt=0.125 fs | 233/401 failed frames; first 392.5 fs | Smaller timestep does not restore geometry |
| 1 ps half-step drift | -0.072749 meV/atom/ps | Energy conservation can pass for an invalid potential |
| DFT force discrepancy at first failure | 1.283081 eV/A | Post hoc diagnostic, never fitted |

At 395 fs, PBE raises energy by 5.965761 eV relative to
its initial label, while the model lowers its own energy by 2.289111 eV.
Each change uses its own method's initial energy for the same cell and composition.
This is not an adsorption-energy calculation. The diagnostic point was chosen because
the geometry failed; it is not an unbiased new test set.

## What was diagnosed

The exact resolved DeepMD `out.json` is required for checkpoint reconstruction: raw
`input.json` omits defaults such as descriptor `type_one_side` and fitting `resnet_dt`.
The initial diagnostic loader rejected that mismatch. Original training was correct;
the failed helper and corrected retry are retained.

Model energy/force finite differences and force-loss parameter directional derivatives
agree with analytic derivatives. Descriptor and fitting gradients are nonzero. Actual
DFT energy finite differences for P11-y and H18-z agree with original force labels within
0.000044 eV/A. All 47 atoms are mobile and fitted. The 5 A descriptor has adequate neighbor
capacity for the inspected structures, but no explicit long-range electrostatics.
Neither capacity accounting nor one-frame success establishes global representability.

The original Adam schedule ended after 1500 one-frame updates near learning rate 1.4e-5.
The same compact architecture subsequently passed the unchanged one-frame gate after
710 plus 5835 L-BFGS objective evaluations. This demonstrates an optimizer
remedy for that frame; it does not imply that the chosen network or training is converged.
The full-data candidate used 1800 additional objective evaluations, fixed before evaluation.
Its full validation RMSE determined selection. The old test had already been inspected,
so recovery results are explicitly retrospective. No offset or force correction was fitted.

Saved velocities independently reconstruct kinetic energy and temperature using 138
degrees of freedom and LAMMPS metal units. There is no thermostat or hidden constraint.
Original GPU repeat differences start near roundoff and grow during unstable dynamics.
A CPU uninterrupted repeat reproduces saved values exactly; ordering/rebuild controls
reduce but do not eliminate the original strict restart discrepancy. The recovered
model passes that same strict restart test over the stable first 250 fs. The binary
checkpoint itself is not shown to be corrupt.

## Fixed settings and limits

The model type order is Ca/P/O/H/C, double precision, `se_e2_a`, cutoff 5.0 A,
smoothing 4.5 A, `sel=[32,32,96,32,16]`, axis width 4. The full-data descriptor is
16/32/64 and fitting network 64/64/64. No virial loss is used. The diagnostic optimizer
minimizes force-component MSE plus per-atom energy MSE, with L-BFGS history 20, learning
rate 1, strong-Wolfe search, at most 20 inner iterations / 25 evaluations per outer step.
The first one-frame run used 30 outer steps; its final continuation used 250, capped at
6000 calls / 900 seconds. The full-data run allowed 90 outer steps, capped at 1800 calls /
900 seconds. Limits, resolved configuration, weights, predictions and exact commands
are retained under `.research-runs/recovery-001`.

New DFT diagnostics use the same plain PBE, DZVP-MOLOPT-SR-GTH and 1100/100 Ry label
settings. The two finite differences use EPS_SCF=1e-9 and +/-0.001 A. The first-failure
single point uses EPS_SCF=1e-8. They are diagnostic labels, not added training examples.

LAMMPS uses `units metal`, all-mobile NVE, 200 K initial velocities, seed 951, and
explicit CO2 IDs 45/46/47 (oxygen IDs 46/47). The full-duration controls use 4000 x
0.25 fs and 8000 x 0.125 fs. No thermostat, clipping or bonds were added to hide failure.

Compute containers were capped at 4 CPUs / 10 GiB RAM; optimizer Torch allocation at
35% VRAM. Only one GPU computation ran at a time. A single-thread CPU control briefly
coexisted with GPU fitting, under the Docker VM's memory capacity. Resource samples
include all GPU workloads, and are sampled maxima rather than continuous peaks.
The original `refine-001` resource profile remains separate. No system installation,
driver change, credential change, upload, remote connection or HPC submission was used
for this recovery.

## Verify preserved evidence

```powershell
cd D:\project\HAP_CO2mineralization
python workflow/capture_refinement.py --run .research-runs/refine-001 --verify
python workflow/capture_recovery.py --run .research-runs/recovery-001 --verify
.\workflow\run-refinement.ps1 -Stage all -Run refine-001
.\workflow\run-local.ps1 -Stage check
```

The captured-run entrypoint verifies and returns without rewriting calculations or
figures. Static tests are not GPU-physics validation. The new thermodynamic parser tests
reject truncated, malformed, duplicate-time and nonfinite output. Exact calculation
commands, exit status, elapsed time and log hashes are in each `*.execution.json`.

## Replay a known-failed diagnostic in a fresh folder

This repeats the actual selected model's 1 ps half-step case. It is useful for software
reproduction, **not production sampling**. Retain the existing raw runs and use a new name.

```powershell
cd D:\project\HAP_CO2mineralization
$replay='.research-runs/recovery-replay-001'
if (Test-Path -LiteralPath $replay) { throw 'Choose a new replay name' }
New-Item -ItemType Directory -Path $replay | Out-Null
foreach ($file in @('model.pth','system.data','long-half.in')) {
  Copy-Item -LiteralPath ".research-runs/recovery-001/lammps-recovered/$file" -Destination $replay
}
$image='deepmodeling/deepmd-kit@sha256:23d23b922459eb49b23f9d4977b6e453fd17feb2ee03a280c70fb7f4a6452e8c'
docker run --rm --pull never --gpus device=0 --cpus 4 --memory 10g --mount "type=bind,source=$PWD,target=/work" -w "/work/$replay" -e OMP_NUM_THREADS=4 -e DP_INTRA_OP_PARALLELISM_THREADS=4 -e DP_INTER_OP_PARALLELISM_THREADS=1 --entrypoint lmp $image -in long-half.in -log replay.log
if ($LASTEXITCODE -ne 0) { throw 'Replay failed; preserve its log' }
```

The deployed local runtime and existing pinned files were exercised. A pristine
reinstallation/bootstrap on an empty computer was not tested.

## Required next scientific work

Do not deploy these weights for production HAP/CO2 dynamics. The data cover only eight
four-femtosecond trajectories near one dry thin-slab structure; those frames do not
sample the excursions responsible for failure. A new prospective DFT sampling design
must cover relevant OH/phosphate/CO2 distortions and independent initial conditions,
with an untouched holdout and force/energy/stability gates fixed before fitting.
Longer optimizer budgets alone have not established that coverage. Surface thickness,
hydration, termination, basis/vacuum/cell convergence and long-range electrostatics also
remain unvalidated. Do not promote these artifacts to material-property evidence.

`PARAMETERS.md`, `profiles/hpc-candidate-production.json` and `HPC_HANDOFF.md` distinguish
actual smoke settings from unvalidated candidate production settings. HPC scripts remain
disabled and unsubmitted. Only site-specific scheduler/account/partition/QoS,
modules or Apptainer/image paths, hardware/resource limits and shared/scratch paths
remain necessary for site adaptation; scientific adequacy is a separate unresolved task.
