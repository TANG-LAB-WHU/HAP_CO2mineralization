# Reproducing this snapshot

There are three levels, and only the first works from a plain clone.

## 1. Check the package and recompute its tables (standard library only)

From the repository root:

```sh
python results/pipeline-v2-20261104/scripts/summarize.py --check
```

It checks every file against `manifest.json` (size and SHA-256), checks that no host path is left in the records,
recomputes the four tables in `tables/` from the bundled JSON/CSV files, and compares them byte for byte. It loads
no model and runs no calculation. Without `--check` it rewrites `tables/` and prints the headline numbers quoted in
the README.

## 2. Regenerate the snapshot from the run directory (Shawn-PC, or a copy of it)

```sh
python workflow/finite_size_compare.py --labels .research-runs/pipeline-v2/finite-size/labels --out finite_size_pairs.csv
python workflow/export_pipeline_v2_results.py --run .research-runs/pipeline-v2 --finite-size finite_size_pairs.csv --out results/pipeline-v2-<date>
python results/pipeline-v2-<date>/scripts/summarize.py
```

`finite_size_compare.py` needs ASE and NumPy. `export_pipeline_v2_results.py` uses the standard library, plus the
Docker CLI to record image identities and package versions. Copy `scripts/`, `README.md` and `REPRODUCE.md` into a
new package directory before running the exporter. The exporter copies only stage records that the orchestrator
has marked done in `pipeline_state.json`, counts label steps still running, and records excluded assets by
SHA-256. It rewrites absolute paths under the repository, and the container mount `/work`, to relative paths; each
rewritten file is flagged with `host_paths_rewritten` and keeps its source SHA-256 in `manifest.json`. Against an
unchanged run directory and the same images, a rerun reproduces every record except the snapshot time and the
progress counts of label steps that are still running (`state/progress.json`, and the manifest entries for them).

## 3. Rerun the computation

This needs the run's inputs and environment, and it is not bitwise reproducible.

- **Code:** the commit that adds this package. The run itself used commit `de20ed8` of
  `integration/shawn-4090-reproducible`, with the uncommitted changes listed in `manifest.json`
  (`run_checkout_uncommitted_code`). Those changes are part of this commit: the CUDA allocator setting in
  `workflow/pipeline_v2.py`, which is active from the orchestrator restart on 2026-11-03 at 18:21 local time, so
  from train-3 on, and documentation.
- **Images** (`environment.json`): `cp2k/cp2k:2026.2_mpich_generic_psmp` pulled by the recorded digest, and
  `hap-mlip:2026.10` built locally from `workflow/environment/Dockerfile.mlip` on the
  `ghcr.io/deepmodeling/deepmd-kit:3.1.3_cuda129` base. That Dockerfile does not pin its pip packages, so a rebuild
  can differ. The versions in the image that ran are listed in `environment/hap-mlip-2026.10-packages.txt`: DeePMD-kit
  3.1.3, PyTorch 2.10.0, TensorFlow 2.19.1, LAMMPS 29 Aug 2024, PLUMED 2.9.2, mace-torch 0.3.17, MatterSim 1.2.5 and
  ASE 3.29.0. Pin them with `pip install -c` when rebuilding. `hap-cp2k:2026.2-cuda-sm89` was used only for
  benchmarks.
- **Inputs:** the exploration start structures are in `exploration/` (w8 and s3, plus the seed-777 sealed-test
  sets). The MatterGen survivors come from `0.Generative_Discovery/runs/Ca-P-O-H-C_gfscan_RTX5090_20260329/`.
- **Command** (orchestrator; it resumes from `pipeline_state.json`):

  ```sh
  python workflow/pipeline_v2.py --run .research-runs/pipeline-v2 --cutoff 800 --rel-cutoff 60 --iters 6 \
      --gpu-jobs 2 --explore-jobs 3 --warm-steps 100000 --extra-systems s3
  ```

  The initial MACE exploration, frame selection, finite-size labels and S3 extras ran as separate steps, with the
  commands and seeds in `workflow/PIPELINE_V2.md`. Seeds: 20261029 for exploration, frame selection and LAMMPS
  velocities (plus the temperature), 777 for the sealed-test trajectories, and DeePMD seeds 11, 23, 37 and 53 for
  models 0–3 (`iter-*/train/model-*/input.json`).
- **Not bitwise reproducible:** GPU training and GPU MD are not deterministic, and MD trajectories diverge from
  any small difference, so a rerun gives different frames, candidates and models. Compare a rerun statistically,
  for example against the spread of the four ensemble members, not number for number.
- **Hardware:** i9-14900K, 64 GB, RTX 4090 24 GB (driver in `environment.json`), Docker Desktop on WSL2. Timings
  in `workflow/PIPELINE_V2.md` hold only for this machine and its concurrent load.
