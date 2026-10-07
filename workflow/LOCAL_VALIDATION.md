# Shawn-PC local validation

Checkout: `D:\project\HAP_CO2mineralization`, branch `integration/shawn-4090-reproducible`.
Fresh fetch confirmed `shawn_dev` 028ca5ac5f6bdfa5b5f72715110edbb94f5997e1 contains
`feat-mattergen-embedded` c2b961103dfa7b6821e8e156f60d1be13f075280 (4 ahead, 0 behind).
No merge was needed; no remote branch/default was changed or pushed.

## Verified environment

SHAWN-PC: RTX 4090, 24,564 MiB VRAM; 68,447,625,216 bytes physical RAM.
D: initially had 356,942,217,216 bytes free. Docker Desktop 4.40.0 / Engine 28.0.4
was already running, with a 32 GiB Linux VM. Only docker-desktop WSL was installed.
Host Anaconda Torch 2.5.1 is CPU-only; GPU calculations use existing Linux images.
Existing open-webui was preserved. No Windows driver, security, global Python, WSL
distribution, or existing image was replaced.

Pinned images are in `profiles/rtx4090-smoke.json`; MACE package overlay is locked in
`environment/requirements-mace.lock`. The optional Dockerfile recreates this overlay
inside an image; the tested runtime is the project-local Linux venv at
`.research-runs/venv` over the pinned DeepMD image. The Dockerfile itself is not build-tested.
MatterGen has a separate `.research-runs/mattergen-venv` and official source commit
92423660a8bd70e83679086e88f88596d484dc16. Large weights, datasets, caches and raw logs are ignored.
The installed environments and Windows-to-Docker entrypoint were exercised. A pristine
full bootstrap from the final locks was not repeated; do not claim that as tested.

## Actual verification matrix

| Stage | Verified execution | Scientific interpretation |
|---|---|---|
| HAP preparation | 44-atom ordered bulk, (001) slab, isolated CO2, 47-atom adsorbed configuration; source checksum, composition and minimum-distance checks | Explicit OH ordering is an assumption, not reviewed truth; one thin dry termination only |
| MACE CUDA | Matched references, float64, fmax 0.05 eV/A at 78/1/77 steps; saved-geometry energy repeatability; clean final-code rerun | Model/surface/site/coverage convergence unvalidated |
| CP2K | All three matched PBE single points SCF-converged; actual installed CUDA image works on 4090 | Final energy difference +0.108263 eV vs MACE -0.176288 eV; **ranking accuracy not validated** |
| Stress | +/-0.1% isotropic strain: numerical derivative 0.2881297 eV vs -trace(virial) 0.2885053 eV | Confirms conversion sign/volume/unit implementation on this test, not pressure convergence |
| AIMD | Three CO2 and three HAP+CO2 runs, 5 steps each, dt 0.25 fs; real CP2K energies/forces/cells | 1.25 fs per run is only an execution test |
| Dataset / DeepMD | Exact coordinate/force/cell step alignment; shared initial geometry removed; 15 unique HAP+CO2 configurations across three splits; 20 GPU updates, freeze, held-out test | Held-out force RMSE 0.375753 eV/A and energy RMSE 0.030304 eV/atom fail provisional production gates; close replicas remain correlated |
| LAMMPS | Actual HAP+CO2 GPU potential execution; 5+restart+5 vs uninterrupted10, high-precision unwrapped state difference ~1.1e-14 | 2.5 fs is not dynamics validation; trained model is smoke-only |
| MatterGen | One actual 1000-step sample, seed42, CH2CaO4; MatterSim relaxation and reference evaluation completed; chemistry filter rejects 1/1 and saves zero candidates | Generated structure lacks P and is not an apatite candidate; n=1 cannot establish generator quality |
| Repository checks | 154 tests pass in local Linux runtime; manuscript/evidence checker passes | Static/evidence checks supplement but do not replace calculations |
| HPC | Disabled parameter profiles, job/dependency templates, syntax checks | No HPC connection, submission, or cluster validation |

Raw results: `.research-runs/smoke-final`, `stress-check`, `aimd-adsorbed-{0,1,2}`,
`adsorbed-dataset-v2`, and `mattergen-smoke`. The earlier `smoke-001`, `relax-200`,
`co2-dataset`, `adsorbed-dataset`, and restart-check-attempt1 files are preserved.
`validation/verification.json` is the compact machine-checkable result; the local
`.research-runs/reproducibility-manifest.json` records artifact hashes and host evidence.

## Run on this installed machine

From PowerShell:

```powershell
Set-Location D:\project\HAP_CO2mineralization
# Existing artifacts: inspect/verify, without redoing calculations.
python workflow/verify_local_run.py .research-runs
.\workflow\run-local.ps1 -Stage check

# Fresh matched-reference test; choose a new run name to preserve earlier evidence.
.\workflow\run-local.ps1 -Stage prepare -Run my-smoke
.\workflow\run-local.ps1 -Stage mace -Run my-smoke -Steps 200
.\workflow\run-local.ps1 -Stage cp2k -Run my-smoke
```

All Docker operations use `--pull never`, 4 CPUs and 8 GiB RAM for these stages.
MACE limits its allocator to 35% of VRAM. The helper refuses existing run manifests/logs.
To reproduce the model download after intentional environment setup:
`python workflow/download_mace_smoke.py` verifies the exact 32,581,838-byte release hash.

The additional executable stage recipes are in `REPRODUCE_EXTENDED.md`. The candidate
production parameter/convergence plan is in `PARAMETERS.md`; it is deliberately disabled.

## Baseline provenance correction (follow-up audit)

The intended bottom-half FixAtoms constraint did not survive the original extxyz
save/reload. Actual MACE calculations were all-mobile, as independently inspected
from saved files and trajectory behavior. The old raw records are preserved; this
correction supersedes any earlier fixed-bottom description.

## Failures and recovery retained

- First model-fetch attempt preceded completed installation and failed import; retry succeeded.
- First MatterGen dependency resolution selected incompatible emmet-core; pinned 0.84.2,
  mendeleev 0.19.0 and protobuf 5.28.3 in its isolated venv, then generation succeeded.
- First MatterGen evaluation hit LMDB double-open behavior; the wrapper reuses its existing
  read-only transaction and streams identical decompressed bytes to avoid a RAM spike.
  Subsequent retained failures required a pymatgen serialization alias and SMACT 2.5.5
  instead of its incompatible 4.x API; the fourth evaluation attempt completed.
- The first MLIP split included identical initial geometries across replicas. These datasets
  are retained and superseded by the corrected 15-frame dataset; no validity is retroactively claimed.
- First HAP restart comparison used wrapped six-digit coordinates. Preserved files show the
  false cell-vector difference; high-precision unwrapped comparison passes.
- Native Windows aggregate checks fail on POSIX directory-open safeguards and GBK Unicode.
  The supported check command uses Linux. A namespace collision with an installed `scripts`
  package was fixed by making the repository utilities an explicit Python package.
- Eight legacy CP2K/XYZ conflict-marker files remain quarantined, without guessing their intended inputs.

No scientific result or literature claim was promoted. The manuscript remains unsupported
by new production evidence. No enhanced-sampling calculation was executed.
