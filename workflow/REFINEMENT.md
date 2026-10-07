# Bounded local structure and model validation

This follow-up retains the earlier smoke baseline and tests the scientific meaning of
the interfaces. It is not a production mineralization calculation. Use the completed
run's verification JSON for measured outcomes; do not equate exit status with validity.

## One Windows entrypoint

```powershell
cd D:\project\HAP_CO2mineralization
# Resume/verify retained completed stages, refusing changed or failed evidence.
.\workflow\run-refinement.ps1 -Stage all -Run refine-001
# A fresh run uses a new directory; it can take substantial workstation time.
.\workflow\run-refinement.ps1 -Stage all -Run refine-002
# Read the recorded progress without starting another calculation.
python workflow/refinement_status.py --run .research-runs/refine-001
# Re-render only from actual completed outputs.
.\workflow\run-refinement.ps1 -Stage plots -Run refine-001
```

This uses the already tested project-local Linux venv and immutable local images.
It does not install a driver, change global Python, pull images, submit HPC jobs, or push Git.
The original baseline bootstrap guide remains `REPRODUCE_EXTENDED.md`. A pristine full
environment bootstrap is not claimed as tested. Do not run two calculation instances
at once. Completed job resumption checks input and output hashes and parses convergence.
A failed attempt must be reviewed and retained under a distinct attempt path before retry.

The pipeline stages are implemented in `refine_pipeline.py`; Docker orchestration is in
`run_refinement.py`. `plot_refinement.py` produces PDF/PNG diagnostics and a real-trajectory
GIF. Profiles are in `profiles/rtx4090-refine.json`. Inputs and per-stage source snapshots
inside the run record the exact executed settings, including the adaptive force-grid gate.

## What is different from the baseline

- The original intended FixAtoms constraint was lost on extxyz serialization. Its actual
  calculations were all-mobile. Documentation is corrected; old evidence is unaltered.
- Both parallel and antiparallel OH orderings have explicit occupancies, 44 atoms,
  Ca10(PO4)6(OH)2 composition, Ca/P=10/6, neutral formal charge and checked PO4/OH geometry.
  Antiparallel ordering cancels the initial formal ionic z dipole, not necessarily the
  electronic polarization. It is a test configuration, not an established ground state.
- MACE uses a tighter 0.03 eV/A force threshold. The original failed 400-step case and
  subsequent LBFGS continuation are retained separately.
- Matched PBE references compare 300/40, 500/60 and 700/60 Ry. A separate 500/60 Ry PBE-D3
  calculation changes the electronic method rather than tuning agreement with MACE.
- Fixed-geometry fragment energies separate molecular interaction from substrate
  reconstruction. Negative substrate contributions show a higher-energy isolated reference geometry
  and invalidate interpreting the total difference as pure adsorption strength. The
  PBE single points do not establish DFT stationary basin identities.
- The 500/700 Ry energy difference can pass while forces remain grid-sensitive. A separate
  900/80 versus 1100/100 Ry, SCF 1e-8 check gates acceptance of finer-grid PBE force labels.
- The selected DFT geometry must pass chemical sanity and maximum force <0.05 eV/A
  before sampling. All subsequent stages retain all-mobile boundaries and the same cell.

## Dataset and architecture scope

Eight distinct initial conditions vary seed, temperature (100--400 K), small atomic
displacements, CO2 azimuth and height. Four complete replicas are training, two validation,
and two test, assigned before labels. Each is bounded to 16 NVE steps at 0.25 fs (4 fs).
The fixed eight-start design rejects unsupported replica-count/timestep edits instead of silently ignoring them. Only steps 4/8/12/16 are selected. These trajectories are short local probes, not equilibrated
ensembles. Within-replica correlation and lack of broad coverage remain explicit.

Coordinates, forces, energy and cell records must match by step and atom identity. Energy
and force conversions are cross-checked against independently printed CP2K records.
Cross-split pair-distance-vector RMS differences must be >=0.025 A. Exact geometry hashes
alone are insufficient. No failed/hard test configurations are silently moved or dropped.
`sampling_coverage.json` compares within-replica geometry spans with the nearest other
replica and split. This is a local structural-separation diagnostic; neither an effective
sample size nor statistical independence is inferred from four-femtosecond sequences.

The one-frame 1,500-update overfit experiment is separate from 6,000-update compact/wider
models. The smooth `se_e2_a` descriptor uses rcut=5 A, smoothing at 4.5 A, explicit
Ca/P/O/H/C types, and neighbor-capacity verification with periodic images. Compact widths
are 8/16/32 plus 32/32 fitting; wider widths 16/32/64 plus 64/64/64 fitting. Both use double
precision, a 1e-3 to 1e-5 learning rate, explicit energy/force loss weights and fixed seeds.
The full validation force RMSE selects the architecture. Only then is its test set scored.
No test-set bias fitting is allowed. Rigid translations, rotations and same-species atom
permutations check model symmetries on real held-out coordinates.

Provisional gates: one-frame force RMSE <0.02 eV/A; independent force RMSE <0.1 eV/A;
independent energy RMSE <0.01 eV/atom. Thresholds are not relaxed to make a model pass.
They are pilot development criteria, not universal guarantees. This model has no explicit
long-range term and is not validated for changed compositions/cells, hydrated reactions,
isolated slab/CO2 reference energies, or broad material generalization. No virial labels
are trained. Printed pressure or zero virial error must not be interpreted as validated
stress. The earlier finite-strain stress conversion evidence remains in `stress-check`.

## LAMMPS and figures

The Python and LAMMPS engines must agree on identical 47-atom coordinates, types, forces
and energy before dynamics. CO2 IDs are 45/46/47; its oxygen IDs are 46/47, not all oxygen.
The main NVE run is 1,000 x 0.25 fs =250 fs. A 500+500 binary restart comparison and a
2,000 x 0.125 fs run test restart transfer and timestep sensitivity. Neither duration
establishes long-time stability. Failed model-accuracy gates remain failed even if MD runs.

Figures under `.research-runs/<run>/figures` show:

1. Initial OH arrangements and actual selected geometry.
2. MACE/DFT energy and maximum-force histories, including optimizer recovery.
3. Reference-energy sensitivity, SCF histories and force-grid differences.
4. AIMD temperature, conserved energy, force, displacement, contacts and autocorrelation.
5. Split energy/force distributions and cross-split geometry distances.
6. Learning curves, independent parity plots and error distributions.
7. LAMMPS drift, temperature, contacts, mapping and restart results.
8. Actual trajectory GIF, with physical time, element legend and unwrapped coordinates.

Missing stages are visibly marked by plotting code; no synthetic replacement data exist.
The GIF samples saved trajectory frames at 5 fs intervals without interpolating dynamics.
`lammps/trajectory.dump` is OVITO-readable; `trajectory.extxyz` preserves the cell and species.
All plots record source hashes in `figures/manifest.json`. They are diagnostic artifacts,
not publication-eligible material predictions. A separate user-requested working LaTeX
draft may show them only with this limited claim scope and factual captions.

## Resource and restart limits

Only one local GPU calculation container runs at a time, with 4 CPUs and 10 GiB RAM;
MACE and model-evaluation stages cap the Torch allocator at 35% VRAM; small DeepMD training uses the recorded architecture/batch settings and is monitored separately. The host requires at least 40 GiB free.
The run's intended disk budget is 8 GiB; logs and measured resource snapshots expose
actual usage. Do not overwrite completed calculations to retry. A named timed-out
container is stopped without touching other workloads. The existing open-webui is preserved.
HPC profiles remain disabled; no remote connection or submission is part of this work.

## Reviewed wall-time recovery

`resume-dft.ps1` accepts only a terminal wall-time failure, checks that no refinement
container remains active, and archives the entire first attempt before continuing.
The restart uses the saved CP2K coordinates, wavefunction and BFGS Hessian, retains the
60-step total bound and all convergence thresholds, and reduces the trust radius from
0.25 to 0.05 A after observed late-step force oscillations. Both attempts remain in the
plots. The exact checkpoint hashes are in `dft_continuation.json`.

```powershell
.\workflow\resume-dft.ps1 -Run refine-001
```

This is an explicit, narrowly guarded recovery command, not an instruction to rerun a
completed calculation. See the [CP2K BFGS input documentation](https://manual.cp2k.org/cp2k-2025_1-branch/CP2K_INPUT/MOTION/GEO_OPT/BFGS.html)
for Hessian restart and trust-radius units. A failed convergence gate remains failed.

The original Quarto manuscript can be rendered with the locally downloaded official
Quarto 1.10.19 portable archive (SHA256 in `.research-runs/tools/quarto-release.json`).
No system installation or shell PATH change is required. The separate LaTeX draft is
built with portable Tectonic 0.17.0; the Overleaf archive uses ordinary XeLaTeX/BibTeX.

If the 60-step attempt terminates without convergence while geometry remains valid,
`continue-dft.ps1` provides one additional reviewed 30-step allowance (90 total). It
archives the failed attempt and its audit and keeps all physical/convergence settings.
It refuses active jobs, existing archives, and already converged calculations. No
threshold is relaxed. `audit_cp2k_restarts.py` compares actual same-state energy and
force values across retained restart boundaries; final capture refreshes its hashes.

In `refine-001`, step 90 passed the force threshold but failed displacement convergence.
`optimizer_diagnosis-step90.json` records the positive stored BFGS Hessian approximation and its soft internal
modes, small center-of-mass drift, and remaining OH/CO2 motion. These observations did
not establish a singular or stalled optimizer. A separate bounded algorithm check,
`recover-dft-lbfgs.ps1`, archived that attempt and started fresh limited-memory history
(rank 5, at most 8 force evaluations per iteration, at most 20 additional iterations).
All force, displacement and electronic criteria remained unchanged and all atoms mobile.
It passed all four CP2K convergence tests in one iteration (three force evaluations),
with independently audited maximum force 0.032416 eV/A. This establishes local optimizer
convergence, not a globally optimal surface. The helper is guarded against active jobs
and repeated recovery; do not invoke it on the already completed relaxation.

Final evidence capture uses run-specific filenames and refuses to overwrite an existing
capture. After capture, `python workflow/capture_refinement.py --run .research-runs/refine-001 --verify`
rechecks every captured raw-file hash and the gate report's source hashes/data without
rerunning physics or changing gate outcomes. Failed/partial attempts, including the step-90 displacement failure, remain in
the inventory with publication eligibility disabled.

The [HPC handoff](HPC_HANDOFF.md) lists disabled scripts, actual pilot settings, restart
requirements, and the few site-specific details still needed. It makes no cluster-test
or production-convergence claim.
