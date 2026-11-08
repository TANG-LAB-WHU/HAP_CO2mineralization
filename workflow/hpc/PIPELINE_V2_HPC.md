# Pipeline v2 on a Slurm cluster

Handover of the running pipeline v2 campaign (`workflow/pipeline_v2.py`, protocol and results in
`workflow/PIPELINE_V2.md`) from the local RTX 4090 workstation to a cluster, prepared 2026-11-08. Nothing has
been transferred or submitted, and none of the cluster paths below has been run yet: the Apptainer, Slurm and
MPICH settings are checked locally only as far as stated.

## State at handover

| Item | State on 2026-11-08 12:20 (local time) |
|---|---|
| `pipeline_state.json` | done: label-0 … explore-4, `label-extra-s3`, sealed test. The next step is label-5 |
| label-5 (`iter-005/labels`, 661 frames from explore-4) | local run stopped at 144/661 done; 2 frames recovered by the CG fallback (`candidates_00072`, `00110`); `candidates_00144` was running and restarts from the atomic guess. The 144 finished labels use the same inputs and are kept: `run_labels.sh` skips them |
| Biased-exploration labels (`extra/metad/labels`) | 165/180 done, `..._00165` was running. Labelled, screened and collected by the cpu phase before train-5 (`label-extra-metad`) |
| `extra/metad-test/labels` | 24/24 done, not yet collected (same step) |
| Iterations | `--iters` raised 6 → 8 (`pipeline_v2_args.sh`): iterations 5, 6 and 7 remain, then the final ensemble and stage 7 |

The local stage-7 preview of the iteration-4 models (`preview-iter4/`) keeps running on the workstation; it is a
diagnostic and not part of the pipeline state.

## Files in this directory

| File | Purpose |
|---|---|
| `pipeline_v2_args.sh` | Science arguments of every orchestrator job (cutoffs, iterations, steps, NPT, extra sets) |
| `pipeline-v2.sbatch` | Example job: one phase of the orchestrator under Apptainer, with a run manifest |
| `cp2k-labels-array.sbatch` | Optional example array: CP2K labels of one prepared label directory, sharded over tasks |
| `pipeline_v2_record.py` | Writes `<run>/manifests/job-<id>.manifest.json` (schema `workflow/schemas/run_manifest.example.json`) |
| `pipeline_v2_transfer.rsync` | rsync filter for the part of the run directory the cluster needs |

## Images

| Image | Source | Apptainer |
|---|---|---|
| CP2K 2026.2 (MPICH, generic, psmp) | `cp2k/cp2k:2026.2_mpich_generic_psmp`, digest `sha256:28e35f63a5253b95cc0f50e786370669185ad9f58eacacf33170cb1257216f36` | `apptainer pull cp2k_2026.2.sif docker://cp2k/cp2k@sha256:28e35f63a5253b95cc0f50e786370669185ad9f58eacacf33170cb1257216f36` |
| `hap-mlip:2026.10` | built 2026-10-29 from `workflow/environment/Dockerfile.mlip` on `ghcr.io/deepmodeling/deepmd-kit:3.1.3_cuda129`; local image id `sha256:9c1be5f68117760f0e4d08d91a511ee720ab8573aba211256eaa11f69b6376a8`, no registry digest | see below |

The Dockerfile leaves pip versions open, so a rebuild can differ from the image that trained iterations 0–4. To
keep the same software, export the local image (about 20 GB) and convert it on the cluster:

```bash
docker save hap-mlip:2026.10 -o hap-mlip_2026.10.tar
apptainer build hap-mlip_2026.10.sif docker-archive://hap-mlip_2026.10.tar
```

Pipeline v2 uses DeePMD-kit 3.1.3 (PyTorch for training, TensorFlow `.pb` models for the energy screen),
LAMMPS with the DeePMD plugin, ASE and pymatgen from this image. The CUDA 12.9 build needs a matching driver on
the GPU nodes. `pipeline_v2_record.py` stores each image's SHA256 beside it (`<image>.sha256`) after the first job.

## Host requirements

- Apptainer with `--nv` on GPU nodes, `git`, and Python ≥ 3.10 with numpy and ase outside the images. The
  orchestrator runs `cp2k_label.py` (prepare, collect, relabel) and `dp_ensemble.py write` with the host Python.
  Stage 7 also needs pymatgen and scipy for `cp2k_hull.py cells`.
- A clean checkout of the pushed branch; the sbatch examples refuse a dirty tree.

## Transfer

From the workstation (the operator runs this; no transfer was made):

```bash
rsync -a --filter="merge workflow/hpc/pipeline_v2_transfer.rsync" .research-runs/pipeline-v2/ <cluster>:<run-dir>/
```

The run directory may be outside the repository: the orchestrator then mounts it at `/hap-run` in the MLIP
image (`/work` is the repository). The paths stored in `pipeline_state.json` `info` are local records only; no step
reads them.

## Phases

`--phase cpu` runs CP2K labels with their energy screens, relabels and the recheck before each training step;
`--phase gpu` runs DeePMD training and the sealed test, LAMMPS exploration and stage-7 validation. Each run stops
before the first step of the other phase and writes that phase to `<run>/next_phase` (`done` at the end). From the
current state the sequence is:

| Job | Steps |
|---|---|
| cpu | rest of label-5; `label-extra-metad`; recheck-5 |
| gpu | train-5, test-5, explore-5 (120k steps) |
| cpu | label-6, recheck-6 |
| gpu | train-6, test-6, explore-6 (140k steps) |
| cpu | label-7, recheck-7 |
| gpu | train-7, test-7, explore-7 (160k steps; it measures the train-7 ensemble, its candidates are not labelled) |
| cpu | recheck-8 (every training set, with the iteration-7 models) |
| gpu | train-final (4 × 300k steps from scratch), test-final, validate-w8 (500 ps), validate-s3 (200 ps) |

The loop ends earlier if an exploration reaches the 0.99 accurate fraction with every extra set explored. A job
whose phase has nothing to do exits at once, so an alternating chain can be submitted ahead:

```bash
export HAP_REPO=<repo> HAP_RUN=<run-dir> HAP_CP2K_SIF=<cp2k.sif> HAP_MLIP_SIF=<hap-mlip.sif>
c=$(sbatch --parsable --export=ALL,HAP_PHASE=cpu --cpus-per-task=<cores> workflow/hpc/pipeline-v2.sbatch)
g=$(sbatch --parsable --dependency=afterok:$c --export=ALL,HAP_PHASE=gpu --gpus=<n> workflow/hpc/pipeline-v2.sbatch)
```

`afterok` only checks the exit code; the gates in `PIPELINE_V2.md` are reviewed from the recorded results.

### What the Apptainer runtime changes

- CP2K: `apptainer exec <cp2k.sif> /opt/cp2k/install/bin/launch bash run_labels.sh …`. Docker runs the image's
  entrypoint `launch`, which puts `mpirun` and `cp2k.psmp` on PATH; `apptainer exec` does not, so it is called
  explicitly.
- MPICH 5.0.1 in the image starts its ranks through `srun` when it sees a Slurm job, and the image has no `srun`.
  With fake Slurm variables in the local container, `mpirun -np 2` failed with "execvp error on file srun"; with
  `HYDRA_LAUNCHER=fork HYDRA_RMK=user` it ran. The orchestrator and the array example set both. One job therefore
  uses one node; `--ranks` is the MPI rank count on it (20 on the workstation).
- Training: each of the `--gpu-jobs` models at a time gets one of the job's GPUs (`CUDA_VISIBLE_DEVICES`) in turn,
  so `--gpus=4 --gpu-jobs 4` trains the ensemble at once. The LAMMPS runs of one exploration share a GPU
  (`--explore-jobs`).
- Seed containers (`seed-*`) and the pause scheduler are Docker-only; under Apptainer the extra sets are
  labelled by the cpu phase.
- The cpu phase screens labels with the TF models in the MLIP image without a GPU (`--nv` only warns on CPU
  nodes). On the workstation the screen ran on the GPU, so its CPU time on the cluster is not measured.

### Sharded labels (optional)

`run_labels.sh <dir> <ranks> <threads> <shard> <shards>` labels every `<shards>`-th frame directory starting at
`<shard>`; done and failed frames are skipped. To spread label-5 over nodes:

```bash
python3 workflow/cp2k_label.py scripts --label-dir $HAP_RUN/iter-005/labels   # current run_labels.sh
sbatch --array=0-7 --export=ALL,HAP_LABEL_DIR=$HAP_RUN/iter-005/labels workflow/hpc/cp2k-labels-array.sbatch
```

then a cpu-phase job screens, relabels and collects. Label directories of later iterations exist only after the
gpu phase has explored (`cp2k_label.py prepare` runs at the start of the cpu phase), so arrays need a prepare
step first. Never run an array beside a cpu-phase job on the same directory.

## Measured local costs

From the workstation runs (i9-14900K, 20 MPI ranks; RTX 4090), as a starting point for walltimes, not cluster
sizing:

| Step | Measured |
|---|---|
| CP2K label, median wall time per frame (label-4) | 80 atoms 148 s, 120 atoms 191 s, 160 atoms 253 s, 169 atoms 289 s, 224 atoms 334 s |
| label-5 remainder | 517 frames, about 30 h of CP2K wall time at those medians |
| train-4 (100k warm-start steps, 2 models at a time on one GPU) | 4.4 h per pair, 8.8 h for the ensemble |
| explore-4 (100k steps, 3 LAMMPS runs at once) | 14.2 h |
| Stage-7 preview, S3 systems (200 ps, 3 at once, preview-iter4) | 5.0–5.5 h per system (first four) |

## Not settled by this handover

- Site details: account, partition, QoS, GPU type and count, walltime limits, Apptainer module.
- `AGENTS.md` still says "HPC integration is disabled", and `workflow/profiles/whu-hpc.example.yaml` has
  `submit_enabled: false`. Enabling submission is the operator's decision; these files were not changed.
- Labels computed on the cluster use the same image and inputs as the local ones; mixed datasets are expected.
  A cluster rerun of a few local frames would show whether the two machines agree to SCF precision.
