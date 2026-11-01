# Pipeline v2: GPU-first MatterGen → CP2K → DeePMD → LAMMPS

Status: **in progress (started 2026-10-29)**. Nothing produced by this pipeline is
manuscript-eligible until its gate report and run manifest are reviewed.

Run directory: `.research-runs/pipeline-v2/` (not in Git). Orchestrator: `python workflow/pipeline_v2.py
--run .research-runs/pipeline-v2 --cutoff 800 --rel-cutoff 60 --iters 6 --gpu-jobs 2 --explore-jobs 3
--warm-steps 100000 --extra-systems s3`; it resumes from `pipeline_state.json`. Log: `logs/pipeline_v2.log`. The DFT hull of the
MatterGen survivors runs beside it (`0.Generative_Discovery/screening/cp2k_hull.py`, log `logs/dft_hull.log`).

| Stage | State on 2026-10-29 |
| --- | --- |
| Images | `hap-cp2k:2026.2-cuda-sm89` and `hap-mlip:2026.10` built; CP2K labels use `cp2k/cp2k:2026.2_mpich_generic_psmp` |
| CP2K protocol | done: 800/60 Ry; ≥ 8 Å Γ cells with the recorded bias |
| MatterGen S1/S2 | done: 363 → 149 (see the run README) |
| MatterGen S3 | done: 149 → 10 (see the run README) |
| Exploration | done: 8 systems × 3 temperatures (MACE-MPA-0, float32); independent seed-777 run for the sealed test |
| Iteration-0 labels | done: 216/216 (20 MPI ranks after the 2026-10-29 Docker crash). The 28 CaO frames failed the first pass and were relabelled after the cell fix below; they enter training from iteration 1 |
| Sealed test labels | done: 83/84 (one portlandite frame with SCF not converged is excluded) |
| Iteration-0 ensemble | trained (4 × 300k steps, one model at a time, about 6.8 h each); sealed test passes, see below |
| Iteration-0 LAMMPS exploration | first run invalid (see "LAMMPS neighbour skin" below), kept in `iter-000/explore_INVALID_lammps_skin2`; rerun with skin 0 started 2026-10-31 |
| Finite-size check | 5/20 paired labels done, paused to free the CPU; resume with `run_labels.sh` in `finite-size/labels` |
| Iteration 1 | 153/153 labels; 4 × 100k warm-start steps (2 models at a time); sealed test passes; LAMMPS exploration (40k steps) started 2026-11-01 00:50, restarted with 3 jobs at once after its 6th job |
| MatterGen S3 survivors in training | MACE exploration done; 90 + 40 frames labelling on the CPU from 2026-11-01 05:00, expected to join at train-3 (see "MatterGen survivors in the pipeline") |
| MatterGen S4 DFT hull | 44 structures queued 2026-11-01 01:29; by 05:50, 7 of 10 survivors relaxed (cell optimisation converged, final single point done), references not started (see "MatterGen survivors in the pipeline") |

## Why v2

The v1 local pilot failed its production gates: retrospective force RMSE 0.108 eV/Å on
32 frames from four-femtosecond trajectories, and 1 ps NVE runs that broke apart (see
`deliverables/README_local_results.md`). The cause was too few, too correlated labels and
a tiny model, not the software stack. v2 changes the data strategy and uses the hardware
for what it is good at.

## Hardware split (Shawn-PC: i9-14900K 24C/32T, 64 GB, RTX 4090 24 GB)

| Work | Where | Why |
| --- | --- | --- |
| CP2K labelling (PBE-D3, FP64) | CPU, MPI ranks on all cores; GPU build benchmarked before use | RTX 4090 FP64 is 1/64 of FP32 (≈1.3 TFLOPS), so CP2K GPU offload gains are not guaranteed |
| Foundation-MLIP screening and exploration MD (MACE, MatterSim, SevenNet) | GPU | FP32/FP64 small cells, GPU-bound |
| DeePMD training (4-model ensemble) and model deviation | GPU | FP32 |
| LAMMPS + DeePMD production MD | GPU | FP32 |

CP2K and the GPU stages run at the same time so neither side idles.

### Measured CP2K GPU vs CPU (2026-10-29, HAP bulk single point, PBE-D3(BJ), 600/60 Ry, Γ)

Same input; other GPU jobs paused during the runs. Energies agree to 1e-12 Ha.

| Cell | `hap-cp2k:2026.2-cuda-sm89`, 4 ranks × 8 threads + RTX 4090 | same, 8 × 4 | `cp2k/cp2k:2026.2_mpich_generic_psmp`, 24 × 1 (CPU only) | 16 × 2 | 32 × 1 |
| --- | ---: | ---: | ---: | ---: | ---: |
| 88 atoms | 251 s | 267 s | **97 s** | 116 s | 154 s |
| 352 atoms | 902 s | — | **476 s** | — | — |

Docker Hub `cp2k/cp2k` also has CUDA tags, but only `generic_cuda_P100_psmp` builds, the newest
being 2025.1. They contain sm_60 cubins plus sm_60 PTX, so they run on the 4090 through driver JIT.
Measured on the same 88-atom input: `2025.1_mpich_generic_cuda_P100_psmp` took 349 s (4 × 8) and
368 s (8 × 4), slower than both the sm_89 build and the CPU image.

**Decision:** CP2K labelling runs CPU-only (24 MPI ranks × 1 thread) on the official 2026.2 image.
The sm_89 GPU build works but is 1.9–2.6× slower because the 4090's FP64 throughput is low.
It is kept for re-testing on FP64-capable GPUs (A100/H100) or larger cells. The GPU's full
capacity goes to the MLIP stages.

## Images

| Image | Source | Notes |
| --- | --- | --- |
| `hap-cp2k:2026.2-cuda-sm89` | `environment/Dockerfile.cp2k-cuda-sm89` | CP2K v2026.2 psmp, CUDA sm_89, toolchain `-march=native`. Upstream ships CUDA images only to 2025.1 and only for P100 (sm_60); GitHub releases ship source and a CPU-only static `.ssmp` |
| `hap-mlip:2026.10` | `environment/Dockerfile.mlip` | DeePMD-kit 3.1.3 CUDA 12.9 base (includes LAMMPS) + MACE, MatterSim, SevenNet venvs, phonopy, pymatgen |

## CP2K labelling protocol (one protocol for every label)

PBE + D3(BJ), DZVP-MOLOPT-SR-GTH, GTH-PBE pseudopotentials (Ca q10, P q5, O q6, C q4, H q1),
Γ point only on cells whose perpendicular widths all exceed the converged minimum,
OT/DIIS with FULL_SINGLE_INVERSE, EPS_SCF 1e-7. Cutoff and minimum width are fixed by
`cp2k_protocol.py` before labelling starts:

- cutoff: smallest value whose energy is within 1 meV/atom and force RMS within 0.01 eV/Å of 1000 Ry;
- width: smallest Γ supercell whose energy is within 1 meV/atom of the 16 Å cell.

Labels from v1 (PBE without dispersion, 300 Ry) are not mixed into v2 training data.

### Measured protocol tests (2026-10-29, HAP bulk; `cp2k-protocol/protocol_report.json`)

| Cutoff / rel. cutoff (Ry), 88 atoms | ΔE vs 1000 Ry (meV/atom) | ΔF RMS (eV/Å) | Gate |
| --- | ---: | ---: | --- |
| 400 / 60 | −0.01 | 0.055 | fail |
| 500 / 60 | −1.11 | 0.030 | fail |
| 600 / 40, 600 / 60, 600 / 80 | −0.06, −0.05, −0.05 | 0.011 | fail (just above 0.01) |
| **800 / 60** | −0.02 | 0.0025 | **pass → chosen** |

| Γ cell (600/60 Ry) | ΔE vs 528-atom cell (meV/atom) | Gate |
| --- | ---: | --- |
| 88 atoms (min width 8.2 Å) | +2.67 | fail |
| 352 atoms (12 Å) | +0.04 | pass |

Between the 88- and 352-atom cells, force RMS differs by 0.017 eV/Å (max 0.040).

**Decision (deviation from the width gate, recorded):** active-learning labels use cells with every
perpendicular width ≥ 8 Å (88–170 atoms). The 12 Å cells would make the interface systems
355–676 atoms and multiply single-workstation labelling time by roughly 3–4×. The cost is a
systematic Γ-sampling energy bias of about 2.7 meV/atom; the force effect (0.017 eV/Å) is
far below the 0.10 eV/Å target. The bias is measured directly on 10 HAP frames labelled at
both widths (`finite-size/`). Energy-RMSE gates are evaluated against 8 Å labels, and the
bias is reported next to them. Labels move to ≥ 12 Å when HPC becomes available.

## Accuracy gates (from published Deep Potential work on the same chemistry)

Reference: the SCAN-ML calcium carbonate Deep Potential ("Ab initio machine-learning simulation of
calcium carbonate from aqueous solutions to the solid state", PNAS,
[doi:10.1073/pnas.2415663122](https://www.pnas.org/doi/10.1073/pnas.2415663122); metadata not yet
entered in `literature/references.csv`) reports energy
RMSE 0.4 meV/atom and force RMSE 101 meV/Å. It validates with 1 ns NPT lattice constants and
10 ns interface and solvation runs. Other mineral MLPs typically report 0.5–6.5 meV/atom and
60–200 meV/Å ([C-S-H MLP, Adv. Sci.](https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.77942)).
This project targets the same order. It does not chase higher accuracy:

| Gate | Threshold |
| --- | --- |
| Sealed test set (independent trajectories, never used for selection) | force RMSE ≤ 0.10 eV/Å overall and ≤ 0.12 eV/Å on interface frames; energy RMSE ≤ 2 meV/atom |
| Bulk HAP lattice (NPT 300 K, ≥ 200 ps) | a, c within 1 % of the CP2K cell optimisation |
| Stability | 300 K NVT/NPT ≥ 500 ps with no geometry-gate failures; NVE drift ≤ 1 meV/atom/ps |
| Active-learning convergence | ≥ 99 % of explored frames below trust_lo |

### Iteration-0 sealed test (2026-10-31, `iter-000/train/sealed_test.json`)

83 frames from the independent seed-777 trajectories. Overall force RMSE is atom-weighted; energy RMSE is frame-weighted.

| Model | Force RMSE, all (eV/Å) | Force RMSE, HAP(001) interfaces (eV/Å) | Energy RMSE (meV/atom) |
| --- | ---: | ---: | ---: |
| 0 | 0.071 | 0.101 | 1.64 |
| 1 | 0.070 | 0.099 | 1.74 |
| 2 | 0.070 | 0.099 | 1.72 |
| 3 | 0.070 | 0.101 | 1.65 |

All four models pass the sealed-test gate (≤ 0.10, ≤ 0.12 and ≤ 2 meV/atom). Weakest systems:
HAP(001) + water + CO2 has a force RMSE of 0.106–0.109 eV/Å, and portlandite has an energy RMSE of 3.8–4.0 meV/atom.
CaO was absent from iteration-0 training but still scores 0.020 eV/Å and 1.2–1.7 meV/atom.
This passes one gate only. The lattice, stability and convergence gates are still open.

Cell fix: CP2K 2026.2 aborts with "No meaningful cell information" on ASE's rock-salt primitive
vectors, so all CaO frames failed at first. `cp2k_label.py prepare` now rotates every new frame to
the lower-triangular standard form, and it leaves finished labels untouched.

### DeePMD training throughput and warm start (2026-10-31)

Short benchmark: iteration-0 data, 4000 steps per model, CP2K running on 20 cores at the same time.

| Models trained at once | Aggregate steps/s | Peak GPU utilisation |
| ---: | ---: | ---: |
| 1 | 4.8 | 60 % |
| 2 | 8.1 | 95 % |
| 4 | 11.3 | 98 % |

The same single model ran at 12.3 steps/s when the CPU was idle during iteration 0, so DeePMD-PyTorch training
is partly CPU-bound. Running the 4 models at once recovers most of that loss. No nvlddmkm events were logged
during the 55 min at 95–98 % GPU load.

Iteration-1 training with `--gpu-jobs 4` (2026-10-31 17:09 local) used 23.7 of 24.5 GB VRAM. Within 6 min
it logged about 200 "Graphics FECS Exception" nvlddmkm events, about 40 per minute and not decreasing. The same
error preceded the Docker VM crash on 2026-10-29. The 4 containers were removed and the run restarted with
`--gpu-jobs 2`. That run uses 19.2 GB, has logged no nvlddmkm events, and trains each model at 6.5 steps/s
(1000 steps in 153 s, CPU idle). The orchestrator now runs with `--gpu-jobs 2`.

### LAMMPS exploration throughput (2026-11-01)

One LAMMPS-DeePMD run (4 compressed models, skin 0, 1 MPI rank × 2 OpenMP threads) left the GPU at about
40 % utilisation, and exploration ran its jobs one after another. In iteration 1 that meant about 1 h per job
while MACE also used the GPU. `dp_explore.py --parallel N` now runs N jobs at once and analyses them in input
order afterwards. Benchmark: iteration-1 models, 2000 steps, with the iteration-1 exploration and 24 CP2K ranks
running at the same time.

| Jobs at once | Systems (atoms) | Loop time per job | Job-steps/s |
| ---: | --- | --- | ---: |
| 1 | monetite (224) | 228 s | 8.8 |
| 3 | monetite (224), HAP(001) + water + CO2 (169), HAP bulk (88) | 352 s, 368 s, 350 s | 16.8 |

The run time hardly depends on cell size, so a job is bound by per-step latency (the neighbour list is rebuilt
every step) rather than GPU arithmetic. The 3-job run reached 95 % GPU utilisation and 14.1 GB VRAM with the 4th
process, and no nvlddmkm events. Monetite gave the same result in both runs (accurate 0.902, 4 candidates;
Python cross-check 1.7e-5 eV/Å). The orchestrator runs with `--explore-jobs 3` from the iteration-1 exploration
on. That exploration was restarted after its 6th job (2026-11-01). The 5 completed jobs and the 6th
(hap001_co2_T500) are reused; the 7th had not written `model_devi.out` when it was stopped and reran from the start. On the same hardware the 3-job exploration runs at about 6 steps/s per job (2026-11-01 05:44–05:50: 2000 steps in about 5.5 min, 83–93 % GPU, 9.0 GB).

Iterations ≥ 1 warm-start each model from the previous iteration's checkpoint (`--init-model`), with sel fixed
to the iteration-0 values (`--sel-from`). They train for 100k steps (`--warm-steps`) instead of 300k.
Iteration-1 frames reach at most [39, 15, 53, 64, 16] neighbours, within sel [44, 16, 60, 72, 20].
When the loop ends, a final ensemble is retrained from scratch for 300k steps on all labels (`final/train`) and
sealed-tested again. Only that ensemble is meant for production MD.

### LAMMPS neighbour skin (2026-10-31)

The first iteration-0 exploration flagged every frame of every system as "failed": max_devi_f was 0.7–3.6 eV/Å
already at step 0, and HAP bulk lost atoms. In Python (`deepmd.infer.DeepPot`), the same models on the same
structure agree to 0.019 eV/Å. Measured on HAP bulk (88 atoms):

| LAMMPS `neighbor` skin | uncompressed model, RMS(F_LAMMPS − F_Python) | compressed model |
| --- | ---: | ---: |
| 0.0 Å | 0.0000 eV/Å | 0.0000 eV/Å |
| 0.3 Å | — | 2.74 eV/Å |
| 0.5 Å | 0.0000 eV/Å | — |
| 1.0 Å | 0.35 eV/Å | — |
| 2.0 Å (metal default) | 1.13 eV/Å | 2.78 eV/Å |

An isolated 30-atom cluster agrees exactly at any skin. The likely explanation is that the DeePMD-kit 3.1.3
PyTorch LAMMPS plugin in `hap-mlip:2026.10` does not drop periodic neighbours between rcut and rcut + skin.
This is an inference and has not been checked upstream. Fix: `neighbor 0.0 bin`, which rebuilds the list
every step. `dp_explore.py` now re-evaluates 5 frames per run in Python and stops if max_devi_f differs by
more than 0.005 eV/Å. Production LAMMPS-DeePMD inputs must use the same setting.

Active learning follows DP-GEN: a 4-model ensemble, with frames whose maximum force deviation
falls between trust_lo = 0.10 and trust_hi = 0.30 eV/Å sent to CP2K. The DP-GEN water example
uses 0.05/0.15 eV/Å; the looser window fits the larger target error for mineral–water interfaces.

## MatterGen survivors in the pipeline (2026-11-01)

Stage 3 below calls for the S3 survivors in the exploration set. The w8 set has none of them: its `CaCO3_20` is the
Alex-MP phase mp-556235. The 10 survivors now join as the extra system set `s3`.

- **Exploration:** `explore_systems.py --only-extra --extra exploration/s3_survivors.extxyz`. Settings are the
  same as w8: MACE-MPA-0 (float32), Langevin NVT at 300/500/700 K, 3000 steps of 0.5 fs, a frame every 40 steps,
  cells with every width ≥ 8 Å (80–160 atoms), seed 20261029. The sealed-test run uses seed 777 at 300/500 K.
- **Frames:** `select_frames.py` with 5 frames dropped at the start of each block. Training: 2 per temperature
  plus 3 strained copies (±4 %) = 9 per system, 90 in all (`--per-block 2`, seed 20261029). Sealed test, from the
  seed-777 run: 1 per temperature plus 2 strained = 4 per system, 40 in all. In every training trajectory (every
  15th frame checked) the shortest interatomic distance stays at 0.87–0.95 Å for the H-bearing cells and
  1.21–1.24 Å for the two CaCO3 cells, so no MD run collapsed. s3_055 still has 22 C neighbours in the selected
  frames (sel 20); every other count is within sel.
- **Labels:** the same CP2K protocol (DFT section identical to iteration 1). The `seed-s3` container labels both
  sets at 16 ranks (`extra/s3/labels`, then `extra/s3-test/labels`). A label costs about 4 min, so the 130 frames
  take about 9–10 h. Blocking a training step on them would leave the GPU idle that long, so they are labelled on
  the CPU while train-2 runs on the GPU. With `--extra-systems s3` the orchestrator:
  - defers s3 to the next training step while `seed-s3` is running or paused;
  - pauses every `seed-*` container while it labels its own frames;
  - labels any remaining s3 frames and collects them once `seed-s3` has exited;
  - makes the final ensemble wait for `seed-s3`;
  - does not stop on convergence until every extra system has been in an exploration.

  From then on the S3 frames are in training, their test frames are in the sealed test (reported per system), and
  their 8 Å cells join the LAMMPS exploration. The S3 systems are new, so active learning starts converging again;
  `--iters` goes from 4 to 6.
- **sel:** warm starts keep the iteration-0 sel [44, 16, 60, 72, 20]. One survivor, s3_055 (CaCO3, 10-atom cell),
  has 22 C neighbours within 6 Å at its MACE geometry. The 2 beyond sel are at 5.97 Å, where the se_e2_a switch
  function is effectively zero, and DeePMD drops the farthest neighbours. The final from-scratch ensemble
  recomputes sel from all data.

**S4 DFT hull** (`0.Generative_Discovery/screening/cp2k_hull.py`, run directory `dft-hull/`):

- **Method:** the labelling protocol (PBE-D3(BJ), DZVP-MOLOPT-SR-GTH, 800/60 Ry, Γ point, widths ≥ 8 Å). Each
  structure gets a cell optimisation (BFGS, 1 atm, pressure tolerance 100 bar, CP2K default force and displacement
  criteria), starting from the S2 MACE geometry. A final single point follows at the relaxed geometry. Survivors and
  tier-1 references with ≤ 450 atoms also get a ≥ 12 Å single point at the start geometry, which gives the Γ
  finite-size correction. CP2K 2026.2 has removed `CELL_OPT/TYPE`; cell optimisations are always direct.
- **References (34):**
  - Tier 1 (10) is every S2 decomposition product of a survivor in any of the three MLIP hulls: dry-ice CO2, two
    ice H2O polymorphs, calcite, portlandite, Ca2P2O7, monetite, HAP (mp-721624), CaH2P2O7, and
    Ca(H2PO4)2·H2O.
  - Tier 2 (24, including 3 added by hand) is one polymorph of every other formula in the CaO–P2O5–CO2–H2O space
    with Alex-MP E_hull ≤ 30 meV/atom: within 10 meV/atom of the lowest polymorph, the smallest cell. The 3 added
    by hand (`--include 187 188 281`) are the two Pnma CaCO3 entries (aragonite candidates) and hexagonal HAP
    (mp-41472).
- **Hull:** E_hull is the distance to the lower convex envelope of the computed references, found by linear
  programming. Every phase is formally neutral in the CaO–P2O5–CO2–H2O space, so no elemental references are needed.
  The gate is E_hull ≤ 50 meV/atom, as in S2, with the decomposition recorded.
- **Scheduling:** 3 jobs × 8 MPI ranks while the pipeline is neither labelling nor training, and 1 × 8 otherwise
  (`seed-s3` counts as labelling; the switch at 05:00 on 2026-11-01 preempted 2 jobs as intended).
  A stopped cell optimisation resumes from its restart file. Tested on CaO: it continued at step 2 with the
  wavefunction restart.

## Stages

1. **MatterGen screening (GPU):** S1 chemistry gate → S2 three-model hull committee
   (E_hull ≤ 50 meV/atom in ≥ 2 models, references from Alex-MP relaxed with the same model)
   → S3 phonons (min frequency ≥ −0.3 THz) and 300 K NPT (5 ps in the 2026-10-29 run, as `screening_summary.json` records; the `mlip_dynamics.py` default is 10 ps). Scripts: `0.Generative_Discovery/screening/`.
2. **CP2K protocol tests (CPU, then GPU timing):** `cp2k_protocol.py`.
3. **Exploration (GPU):** foundation-MLIP MD of HAP bulk, HAP(001) + CO2, HAP(001) + water + CO2,
   calcite, and S3 survivors at 300–600 K; farthest-point selection for the first labels.
4. **Labelling (CPU):** CP2K single points with the fixed protocol.
5. **DeePMD (GPU):** se_e2_a, rcut 6 Å, as in the CaCO3 reference; 4 seeds; compressed for MD.
6. **Active learning (GPU MD + CPU labels):** LAMMPS-DeePMD exploration with model deviation until converged.
7. **Validation (GPU):** the gates above; a failure is recorded, not hidden.
