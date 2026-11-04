# Extended executable recipes

Run these from the checkout on Shawn-PC. They use only the local Docker engine.
Existing output directories are intentionally preserved; use the verification command for
already-completed runs. The examples below describe rebuilding on a fresh checkout/run area.

## Common PowerShell setup

```powershell
Set-Location D:\project\HAP_CO2mineralization
New-Item -ItemType Directory -Force .research-runs/environment | Out-Null
$image='deepmodeling/deepmd-kit@sha256:23d23b922459eb49b23f9d4977b6e453fd17feb2ee03a280c70fb7f4a6452e8c'
$cp2k='cp2k/cp2k@sha256:e2912cc9500a42a08b3e057e383182b96f395b7477c8eb55e93db6e29cb1c2ae'
$common=@('run','--rm','--pull','never','--gpus','device=0','--cpus','4','--memory','8g','--mount',"type=bind,source=$PWD,target=/work",'-w','/work','-e','OMP_NUM_THREADS=4')
```

Stop on any nonzero Docker/Python exit code. Do not treat a final file's existence as success.
The installed venv paths are Linux paths and must be executed inside these containers.

## MACE environment recreation

The tested installation used a system-site-packages venv over the immutable image. Exact
overlay versions are in `workflow/environment/requirements-mace.lock`:

```powershell
docker run --rm --pull never --cpus 4 --memory 8g --mount "type=bind,source=$PWD,target=/work" -w /work --entrypoint /bin/bash $image -lc 'python -m venv --system-site-packages .research-runs/venv && .research-runs/venv/bin/python -m pip install --no-cache-dir -r workflow/environment/requirements-mace.lock'
python workflow/download_mace_smoke.py
```

Then `run-local.ps1` performs prepare/MACE/CP2K. A fresh `smoke-001` with these three stages
provides inputs and matching wavefunctions for the following AIMD smoke.

## Actual HAP+CO2 AIMD → labels → DeepMD → LAMMPS

```powershell
python workflow/prepare_aimd.py .research-runs --system adsorbed --steps 5
docker @common $cp2k /bin/bash -lc 'set -e; for i in 0 1 2; do cd /work/.research-runs/aimd-adsorbed-$i; /usr/local/bin/entrypoint.sh cp2k -i md.inp -o md.log; done'
python workflow/aimd_to_deepmd.py .research-runs --system adsorbed --output .research-runs/adsorbed-dataset-v2
python workflow/prepare_lammps_smoke.py .research-runs --system adsorbed --dataset .research-runs/adsorbed-dataset-v2
docker @common -w /work/.research-runs/adsorbed-dataset-v2 -e DP_INTRA_OP_PARALLELISM_THREADS=4 -e DP_INTER_OP_PARALLELISM_THREADS=1 --entrypoint /bin/bash $image -lc 'set -e; dp --pt train input.json > train.log 2>&1; dp --pt freeze -o adsorbed-smoke.pth > freeze.log 2>&1; dp --pt test -m adsorbed-smoke.pth -s test -n 5 > test.log 2>&1; lmp -in lammps-first.in -log lammps-first.log > lammps-first.console.log 2>&1; lmp -in lammps-restart.in -log lammps-restart.log > lammps-restart.console.log 2>&1; lmp -in lammps-continuous.in -log lammps-continuous.log > lammps-continuous.console.log 2>&1'
```

The final exporter removes the shared step-zero geometry and checks exact coordinate,
force and cell step alignment plus duplicate-coordinate hashes. This avoids exact leakage;
nearby short replicas still do not provide a statistically meaningful held-out benchmark.

## MatterGen

Official source is pinned to `92423660a8bd70e83679086e88f88596d484dc16`. Clone it to
`.research-runs/mattergen` with LFS smudging disabled, check out that exact commit, then fetch
only `checkpoints/chemical_system/**` and `data-release/alex-mp/reference_TRI2024correction.gz`.
The wrapper verifies the published LFS SHA256 hashes before loading either.
The exact final overlay is `workflow/environment/requirements-mattergen.lock`. Its local
MatterGen source reference assumes the same `/work/.research-runs/mattergen` mount location.
Install it into a separate system-site-packages venv, using the official graph-wheel index:

```powershell
docker run --rm --pull never --cpus 4 --memory 8g --mount "type=bind,source=$PWD,target=/work" -w /work --entrypoint /bin/bash $image -lc 'python -m venv --system-site-packages .research-runs/mattergen-venv && .research-runs/mattergen-venv/bin/python -m pip install --no-cache-dir -f https://data.pyg.org/whl/torch-2.2.0+cu121.html -r workflow/environment/requirements-mattergen.lock'
python workflow/download_mace_smoke.py --mattersim
docker @common --memory 10g -e XDG_CACHE_HOME=/work/.research-runs/mattergen-cache --entrypoint /bin/bash $image -lc 'timeout 1200 .research-runs/mattergen-venv/bin/python workflow/mattergen_smoke.py generate > .research-runs/environment/mattergen-generate.log 2>&1'
docker @common --memory 12g -e XDG_CACHE_HOME=/work/.research-runs/mattergen-cache --entrypoint /bin/bash $image -lc 'timeout 1200 .research-runs/mattergen-venv/bin/python workflow/mattergen_smoke.py evaluate > .research-runs/environment/mattergen-evaluate.log 2>&1'
python 0.Generative_Discovery/post_filter.py --results-dir .research-runs/mattergen-smoke
```

The evaluation wrapper documents three compatibility adaptations: streaming lossless gzip,
reuse of an open read-only LMDB transaction, and a serialization module-path alias for
ComputedEntry classes. No reference energies or structures are fabricated or replaced.
`is_novel` labels refer to the evaluation reference, not to novelty of an apatite material.
The generated smoke sample CH2CaO4 contains no phosphorus and must not become a HAP slab.

## Evidence and offline HPC handoff

```powershell
python workflow/verify_local_run.py .research-runs
.\workflow\run-local.ps1 -Stage check
python workflow/capture_evidence.py .research-runs
```

The HPC profiles are disabled and require the site's actual scheduler/allocation/runtime/
hardware limits. `workflow/hpc/dependencies.sh` only prints a dependency plan; it never
submits. The Slurm template must be adapted and cluster-tested by an authorized operator.
