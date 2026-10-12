# Closed HAP91 development diagnostic — 12 October 2026

**DFT data generation succeeded locally; the ML workflow executed; potential
generalization remains below the frozen requirements. This candidate is closed.**
All material here is DEVELOPMENT / manuscript_eligible=false. This package records completed calculations; it contains no new simulation or threshold change.

## What the evidence supports

The converged 91-atom parallel-OH reference and short direct-DFT trajectories
supplied real same-method labels. Six additional reference trajectories passed
recorded QA. A single coverage-v2 fit used 90 training frames and stopped normally
at its time admission boundary; training force RMSE 0.04334653 eV/A missed 0.03.
A separately authorized diagnostic on 16 R9/R10 frames gave 0.07058562 eV/A,
missing the original 0.05 gate. R9/R10 are development diagnostics, not a fresh
acceptance set. All prior optimizer, A1/R3 and 47-atom E2 failures remain retained.

The frozen v2 potential then completed 80 NVE steps / 20 fs / 81 frames at
0.25 fs, with seed 20261012, 91 mobile atoms and fixed triclinic XYZ cell.
Direct DeepMD/LAMMPS run0 agreed in energy exactly to reported precision and
in maximum force component to 1.32883288e-13 eV/A. All monitored geometry gates
passed. Maximum total-energy departure was 0.005717273 meV/atom; temperature
ranged 157.347–300 K, so 300 K initialization is not equilibration.

| 10-fs DFT comparison | Observed | Original gate | Exploratory screen |
|---|---:|---:|---:|
| Whole-system force-component RMSE (eV/A) | 0.08302348 | 0.05 FAIL | 0.10 pass |
| Maximum absolute component error (eV/A) | 0.34466946 | 0.20 FAIL | 0.50 pass |
| Largest element RMSE, H (eV/A) | 0.12032266 | 0.075 FAIL | 0.15 pass |
| CO2 force-component RMSE, nine components (eV/A) | 0.05185223 | no separate original CO2 gate | 0.05 FAIL |
| CO2 maximum component error (eV/A) | 0.09038251 | no separate original CO2 gate | 0.20 pass |
| Absolute energy error per atom (eV) | 0.00005176115 | 0.002 pass | 0.002 pass |

Reference and geometry QA passed at 10 fs. CO2 force RMSE alone failed the
purpose-specific exploratory screen, by 3.70%; the 20-fs DFT probe was therefore
skipped, not assigned a zero or inferred result. Ca/P/H also failed the original
element gate. Energy comparison uses potential energy, without kinetic energy
or a fitted energy offset; raw forces are uncorrected. Passing loose exploratory
criteria would not have reinstated model qualification. These warning criteria
were informed by already-seen development errors and are not independent acceptance.

Charged execution was 40.937 s for MD and 221.078 s for DFT (262.015 s overall),
within the 300/3600 s frozen budgets. Native run0/MD times were 6.380/15.801 s;
these are different timing scopes. Execution stopped for accuracy, not time.
No OOM was reported; R4 numeric files remain unopened and unhashed in this revision.

## One next scientific direction — proposed, not executed

Use a bounded, DFT-led comparison of distinct CO2 adsorption arrangements and
local surface distortions on the known HAP reference to test whether missing
configurational coverage explains the force failure. Predeclare physically
interpretable site/orientation/distance and local OH/PO4 perturbation families,
a finite label budget and unchanged accuracy criteria before any execution.
Keep entire configuration families out of fitting for a genuinely new assessment.
First assess reference energies and forces; only a separately approved decision
could commission another potential. Do not continue this candidate or adjust its
thresholds. Coverage deficiency is a hypothesis, not an exclusive causal diagnosis.

Independent initial velocity seeds isolate trajectory membership but start in
the same geometry/basin; correlated 20-fs samples do not establish independent
adsorption environments, equilibrium, barriers, reaction rates or mineralization.
There is only one learned trajectory and one DFT-checked endpoint here; no
confidence interval or broad transfer estimate is justified. The older 47-atom
PBE-only and current 91-atom PBE-D3(BJ)/dipole results are not a controlled size
or architecture comparison.

Scale-up requires converged DFT reference choices (basis/grid/slab/vacuum and
adsorption references), broader configuration/chemical coverage, a frozen model
passing predeclared independent family-level force/energy tests, reference checks
along longer trajectories, and finite-size/thermodynamic validation appropriate
to the claimed observable. Hydration, charge/proton transfer and reactive pathways
need explicit coverage before mineralization claims. GPU throughput and stable
integration alone satisfy none of these scientific requirements.

## Package contents

- `summary.json` and `data/` preserve measured results, frozen gates and plotting data.
- `figures/` contains the two PNG diagnostics and their captions.
- `inputs/` contains geometry, the actual 10-fs CP2K input, and the frozen MD/probe specification.
- `data/endpoint-10fs-forces.csv` retains all 273 DFT and predicted force components without correction.
- `history/prior-hap91-summary.json` retains the earlier A1/R3 force failure.
- `cuda-cpu-singlepoint/` records the completed same-input GPU/CPU numerical comparison.
- `runtime.json`, `recorded-commands.json` and `provenance.json` retain versions, identities and origin hashes.
- `verify_bundle.py` checks the package and recomputes endpoint errors from the exported components; `make_figures.py` replays the charts.

Scientific validity is unchanged: original force qualification failed, the 20-fs DFT probe was skipped, and R4 numerical files remain sealed and excluded. Short-time integration and GPU throughput do not establish transferability or reaction chemistry. No weights, checkpoints, WFN, raw trajectories, training arrays, machine telemetry, chat records or document artifacts are included. `omitted-artifacts.json` identifies the main retained local assets without publishing them.

See `REPRODUCE.md` for the portable report-replay commands and exact scope.
