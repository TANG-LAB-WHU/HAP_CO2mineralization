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
| Finite-size check | done 2026-11-03 09:42: 20/20 paired labels. 8 Å labels lie 2.46 ± 0.24 meV/atom above 12 Å labels; forces differ by 0.016 eV/Å RMS (see "CP2K labelling protocol"). The last 14 ran as container `seed-fs` (24 ranks), paused during train-2 because they halved its speed and run only during explore-2, gated by a watcher logged as `fs-gate` in the pipeline log |
| Iteration 1 | 153/153 labels; 4 × 100k warm-start steps (2 models at a time, no virial term, see below); sealed test passes; exploration (40k steps, restarted with 3 jobs at once after its 6th job): accurate 0.723, 253 candidates |
| Iteration 2 | labels 251/253, done 2026-11-02 07:40 (about 5 min each at 20 ranks while the hull held 8). Two interface frames failed with SCF not converged after 1050 steps (`candidates_00056` hap001_co2, `candidates_00121` hap001_water_co2); kept and excluded. train-2 (virial term on) started 07:40. Models 0–1 trained at 1.0 s/batch while `seed-s3` and the hull held 28 CPU ranks, 0.50–0.56 with the hull's 8 ranks, and 0.32 with the CPU idle (see "DeePMD training and CPU load"); models 2–3 ran 23:58–08:37 at 0.30–0.31 with the CPU idle. Sealed test passes (see "Iterations 1–2 sealed test"). Exploration (60k steps, first run on the TF backend, 3 jobs at once) started 2026-11-03 08:41. The CaCO3_20 300 K LAMMPS process crashed at step 50 600 (09:18) with a segmentation fault (kernel log: signal 11 in `lmp`, null-pointer read). Its largest force deviation just before the crash was 0.055 eV/Å, so the trajectory had not blown up. Its 50 600 steps are analysed and recorded as an incomplete run (`lammps_complete: false`). Exploration finished at 12:10 (3.5 h, the other 15 runs complete) with accurate 0.838 and 208 candidates; label-3 started at 12:10 and was resumed at 104/208 after the orchestrator restart at 18:21 (see "GPU memory fragmentation"). Runs were 60k steps against 40k in iteration 1, so the accurate fractions are not strictly comparable. Per run, the two blow-ups of iteration 1 are gone: hap001_co2 at 500 K went from 0.057 to 0.883 accurate, and monetite at 500 K from 0.069 to 0.506, with a maximum deviation of 0.59 eV/Å (13.6 before). HAP(001) + water + CO2 is still the gap: 0.008 accurate at 300 K and 0.106 at 500 K. At 300 K, 19 % of frames are above trust_hi, mostly after step 24 000. Every 3000th frame of that run still has both CO2 intact, all 25 H2O, the same 4 singly protonated O and no geometry-gate failure, so the deviation comes from unsampled configurations, not from broken chemistry |
| Iteration 3 | labels running: 187/208 converged at the 2026-11-04 snapshot (07:27 UTC, `results/pipeline-v2-20261104`). `candidates_00127` (HAP(001) + water + CO2, 500 K, step 35 500, max_devi_f 0.104 eV/Å) failed with SCF not converged after 1050 steps; its geometry is normal (169 atoms, shortest distance 0.92 Å, both CO2 intact), and it is kept and excluded like the two iteration-2 failures. train-3 follows with the S3 systems (about 41 training systems, 2 models at once with expandable segments) |
| MatterGen S3 survivors in training | MACE exploration done; training frames 90/90 and test frames 40/40 labelled (`seed-s3` exited 2026-11-02 about 11:00); deferred by train-2, joins at train-3 (see "MatterGen survivors in the pipeline") |
| MatterGen S4 DFT hull | done 2026-11-02 22:51: 44/44 cell optimisations converged, every final single point converged. 8 of 10 survivors pass E_hull ≤ 50 meV/atom; s3_262 (60.9) and s3_180, the only carbonate-phosphate (68.0), fail. None is on the hull. Table and comparison with S2 in the MatterGen run README |
| Validation (stage 7) | `dp_validate.py` written and smoke-tested on iteration-1 models; runs automatically after the final ensemble (see "Validation, stage 7") |

## Why v2

The v1 local pilot failed its production gates: retrospective test force RMSE 0.108 eV/Å on
8 frames (gate 0.100; 32 frames in all, from four-femtosecond trajectories), and both 1 ps trajectories
failed structural validation (see `workflow/RECOVERY.md`). The cause was too few, too correlated labels and
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

Measured 2026-11-03 at 800/60 Ry (`finite-size/labels`, all 20 SCF converged). Each 352-atom frame is its
88-atom frame repeated 2 × 2 × 1, so the two labels differ only in Γ sampling:

| 8 Å minus 12 Å | Mean over 10 frames | Range |
| --- | ---: | ---: |
| Energy (meV/atom) | 2.46 (SD 0.24) | 2.06 – 2.84 |
| Force RMS difference (eV/Å) | 0.016 | 0.014 – 0.016 |
| Largest force-component difference (eV/Å) | | up to 0.047 |
| Stress RMS difference (GPa) | 0.012 | 0.007 – 0.018 |

The offset is close to the 2.67 meV/atom found at 600/60 Ry. Only its spread across frames, 0.24 meV/atom, adds
to an energy RMSE measured against 8 Å labels. The constant part shifts every HAP bulk energy per atom by the same
amount. It cancels in energy differences between HAP bulk frames, but the other systems' offsets have not been
measured, so it does not cancel between systems.

## Accuracy gates (from published Deep Potential work on the same chemistry)

Reference: the SCAN-ML calcium carbonate Deep Potential ("Ab initio machine-learning simulation of
calcium carbonate from aqueous solutions to the solid state", PNAS,
[doi:10.1073/pnas.2415663122](https://www.pnas.org/doi/10.1073/pnas.2415663122); registered as the
unverified candidate `piaggi2025caco3` in `literature/references.csv`; the numbers below are not yet checked against the full text) reports energy
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

### Iterations 1–2 sealed test (`iter-001/train/sealed_test.json`, `iter-002/train/sealed_test.json`)

Same 83 frames and weighting as iteration 0. Iteration 1 trained without the virial term, and iteration 2 with it
(see "Virial loss was off in iterations 0–1").

| Model | Force, all: it 1 → it 2 (eV/Å) | Force, interfaces: it 1 → it 2 (eV/Å) | Energy: it 1 → it 2 (meV/atom) |
| --- | ---: | ---: | ---: |
| 0 | 0.072 → 0.067 | 0.107 → 0.098 | 1.67 → 1.59 |
| 1 | 0.070 → 0.066 | 0.103 → 0.096 | 1.58 → 1.64 |
| 2 | 0.071 → 0.066 | 0.105 → 0.097 | 1.80 → 1.33 |
| 3 | 0.071 → 0.066 | 0.106 → 0.098 | 1.82 → 1.67 |

Both iterations pass the gate. In iteration 1, interface force RMSE rose by 0.004–0.006 eV/Å over iteration 0.
In iteration 2 every model is below iteration 0 on both force measures. HAP(001) + water + CO2 is still the
weakest system for force (0.104–0.107 eV/Å). It is now also the weakest for energy (2.8–3.9 meV/atom), in place of
portlandite.

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

### DeePMD inference speed: TensorFlow backend for LAMMPS (2026-11-01)

The PyTorch backend of DeePMD-kit 3.1.3 needs 80–100 ms per energy-and-force evaluation of HAP bulk (88 atoms,
one compressed model). The time is the same in LAMMPS and in Python, at float32 or float64, and with the CPU
nearly idle. A `torch.profiler` trace attributes about 75 ms per evaluation to the TorchScript interpreter
itself (`forward` self time). The trace also shows about 1800 kernel launches and 27 ms of GPU kernels. On the
CPU alone the same evaluation takes 80 ms. The model is therefore overhead-bound and gains little from the GPU.

`dp convert-backend graph.pth graph.pb` turns the same model into a TensorFlow graph. Measured with LAMMPS on
HAP bulk, model 0 of iteration 1, while CP2K used 28 cores:

| Backend | 1 model, ms/step | 4 models (model deviation), ms/step |
| --- | ---: | ---: |
| PyTorch, compressed (used until now) | 93 | 134 |
| TensorFlow, uncompressed | 21–25 | 39 |
| TensorFlow, compressed | 22 | 56 |
| TensorFlow, uncompressed, CPU only | 55 | — |

The converted model reproduces the PyTorch model. On 12 trajectory frames (HAP(001) + water + CO2, monetite,
HAP(001) + CO2) the largest differences are 1.3e-10 meV/atom in energy and 2.6e-13 eV/Å in force. TF
compression keeps the forces (2.6e-7 eV/Å), but the energy is off by about 510 eV/atom: the per-type energy bias
is lost. The compressed TF graph is not faster either, so it is not used. Three uncompressed TF 4-model jobs at
once run at about 57 ms/step each (52 job-steps/s). The PyTorch exploration managed 18 job-steps/s.

`dp_explore.py` now defaults to `--backend tf`. It converts `graph.pth` (next to the `graph-compress.pth` it is
given) for each of the 4 models into `<out>/models/model-k.pb` and evaluates both on every starting structure.
It stops if they differ by more than 1e-6 eV/atom or 1e-6 eV/Å, and records the check in `explore_summary.json`.
In a 1000-step side-by-side test (3 systems, same seed), the TF and PyTorch explorations gave identical
summaries. Along the trajectories, max_devi_f differed by at most 3e-6 eV/Å, and the TF runs took 57–61 s
against 152–161 s. TensorFlow runs with 2 intra-op threads and 1 inter-op thread per job (`DP_*_PARALLELISM_THREADS`),
which leaves the CPU to CP2K at the same speed.

Training stays on PyTorch. A TensorFlow run of the iteration-0 model-0 input (same data, seed, sel and schedule) trained at 9.6
steps/s against 3.9 for PyTorch under the same CPU load. It diverged, though: the energy RMSE rose from 0.7 to
13 eV/atom between steps 1000 and 5000, and it stopped with NaN at step 6000
(`benchmarks/tf-train-ab/model-0`, kept as a negative result). The two backends also apply this loss
differently (next section), so the run is not a like-for-like comparison.

### Virial loss was off in iterations 0–1 (found 2026-11-01)

`dp_ensemble.py` asked for `start_pref_v 0.0, limit_pref_v 0.02`, a virial weight that rises from 0 to 0.02 as
the learning rate decays. DeePMD-kit 3.1.3 PyTorch keeps a loss term only when both of its prefactors are nonzero
(`deepmd/pt/loss/ener.py`, `has_v = start_pref_v != 0.0 and limit_pref_v != 0.0`). TensorFlow keeps it when
either is nonzero. The iteration-0 and iteration-1 ensembles were therefore trained on energies and forces only:
their `lcurve.out` files have no `rmse_v` column. That holds even though the bulk labels include stress, and the
strained copies were added for NPT. The sealed-test numbers above are force and energy RMSEs and still stand.
How well those models predict stress has not been measured. In a 400-step check, iteration-1 model 0 has a
virial RMSE of about 0.015 eV/atom on its training data. For HAP bulk that is roughly 0.2 GPa.

Fix, effective from train-2: `start_pref_v 2e-4`, which keeps the intended ramp to within 1 % and switches the
term on. A warm-started PyTorch run now reports `rmse_v`. `tests/test_pipeline_v2.py` checks that every term
has two nonzero prefactors. The final ensemble trains from scratch with the virial term from the start.

### DeePMD training and CPU load (2026-11-02)

Each DeePMD-PyTorch training process keeps one CPU thread busy (`docker stats` about 99 %, GPU about 43 %), so
any CP2K work on the CPU slows it. Per model, train-2 with two models at once:

| Other CPU work | Batches | s/batch |
| --- | --- | ---: |
| `seed-s3` (20 ranks) + hull (8 ranks) | first ~10k | 1.0 |
| hull (8 ranks) | 29k–90k | 0.50–0.56 |
| finite-size labels (24 ranks) | 90k–91k | 1.06–1.14 |
| none | 91k–92k | 0.32–0.33 |

So the plan of labelling on the CPU while the GPU trains cost train-2 about a factor of 2. From now on CPU jobs
outside the pipeline's own labels run only during LAMMPS exploration. The TF-backend timings above were measured
while CP2K used 28 cores, though how much that load slows exploration has not been measured. The hull is finished and `seed-s3` has exited, which leaves only the
finite-size labels, gated to exploration.

With the CPU idle, train-2 still takes 0.32 s/batch against 0.154 in train-1. The difference is not isolated.
train-2 adds the virial term and has 23 training systems instead of 15, with a mean of 133 atoms per batch
instead of 124. DeePMD computes the virial in the forward pass whether or not the loss uses it
(`deepmd/pt/model/model/transform_output.py`), so the loss term alone should add little. The benchmark below
isolates the cause.

### GPU memory fragmentation, not virial or data (2026-11-03, `bench-virial/`)

The benchmark ran during label-3, with CP2K holding 20 ranks throughout. Every run warm-starts from iteration-2
model 0 for 2000–4000 steps, and the timing excludes the first 1000 batches. GPU memory was sampled every 30 s
from the Windows `GPU Process Memory` counters (`bench-virial/gpu_memory.log`).

| Run | Data | Virial | Allocator | Run with | s/batch | GPU memory |
| --- | --- | --- | --- | --- | ---: | ---: |
| A / B | iteration 2 (23 systems) | on / off | default | each other | 0.771 / 0.792 | 24 013 MiB in use (not split) |
| B2 alone | iteration 2 | off | default | nothing | 0.212 | not sampled |
| C alone | iteration 1 (15 systems) | off | default | nothing | 0.213 | 9.0 GiB |
| C + B | 1 and 2 | off | default | each other | stalled: neither reached batch 1000 in 40 min, GPU 100 % | 13.8 + 9.0 GiB, 24 010 of 24 564 MiB |
| C + B | 1 and 2 | off | `expandable_segments` | each other | **0.253 / 0.255** | **5.8 + 4.1 GiB** |

- The virial term costs nothing measurable (A against B).
- Training data does not set the speed: 15 and 23 systems run at the same speed alone.
- The cost is PyTorch's default CUDA allocator. Two trainings each cache far more memory than they need, the 24 GB
  GPU fills, and both slow down or stall. An earlier pair (C with B2) also broke down: C ran at 3.2–4.2 s/batch, and B2
  reached batch 1000 only after 4 h, then ran at 0.21 once C had exited.
- With `PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True`, the same pair uses 9.8 GiB instead of 22.8 GiB and runs at
  0.25 s/batch each, 7.9 batches/s in all against 4.7 for one model alone.

The allocator setting changes only how memory is reserved, not the arithmetic. The orchestrator now sets it in every
`hap-mlip` container (`docker_mlip` in `pipeline_v2.py`). It was restarted at 18:21 during label-3, at 104/208.
The frame running then (`candidates_00104`) restarts, and finished labels are kept. Iteration 3 has about 41
training systems, so with the default allocator its 2-model pairs would most likely have stalled.

The CPU-load table above was measured with the default allocator. With expandable segments and 20 CP2K ranks, a pair
runs at 0.25 s/batch, faster than the 0.32 that the default allocator gave with the CPU idle. So part of the slowdown
there blamed on CPU load may have come from the allocator. Its CPU sensitivity with expandable segments has not been
measured with the CPU idle, so the rule of keeping other CPU jobs to exploration windows stays.

### Validation, stage 7 (`dp_validate.py`, 2026-11-01)

The orchestrator runs the gates after `test-final`, with model 0 of the final ensemble on the TF backend
(converted and checked as above):

- **Bulk phases:** 0 K `box/relax tri` at 1 atm. Then NPT at 300 K and 1 atm (`tri`, `flip no`), followed by NVE.
  Both use the exploration cell: the ≥ 8 Å label cell, the same cell as the CP2K reference.
- **Slabs** (`hap001*`): NVT at 300 K, followed by NVE.
- **Run lengths:** 500 ps for the w8 systems, which carry the formal gates, and 200 ps for the extra (S3) systems.
  Every NVE segment is 50 ps.
- **Lattice references:** `cp2k_hull.py cells`. It picks the DFT-hull structure with the same composition, atom
  count and lattice-vector order (lengths within 5 %, angles within 3°) that also passes a StructureMatcher fit.
  It then exports that structure's CP2K-relaxed cell (`opt-1.restart`) with its convergence flag. The pipeline
  `hap_bulk` matches ref-281 (mp-41472, P6_3), calcite ref-186, portlandite ref-218, monetite ref-272 and
  CaO ref-118; each S3 survivor matches its own candidate. CaCO3_20 (mp-556235) and the two slabs have no
  reference.
- **Gates** (thresholds as in the gate table above):
  - lattice: lengths, and also volume, of the NPT average after 100 ps, per input cell, within 1 % of the CP2K
    cell (the formal one is HAP a and c);
  - geometry on every dumped frame (1 ps): the v1 gate of `refine_pipeline.py`, which needs a minimum distance
    > 0.7 Å, 4 O within 1.9 Å of every P, and the nearest O of every H at 0.8–1.2 Å. Every C must also have 2 or
    3 O within 1.6 Å, so both CO2 and CO3 pass;
  - NVE drift ≤ 1 meV/atom/ps.
- **Reported, not gated:** the 4-model deviation on 100 dumped frames.

A smoke test with the iteration-1 models (2 ps NPT/NVT, 1 ps NVE; `validation-smoke/iter-001`) ran every stage
for calcite, portlandite, HAP bulk, HAP(001) + CO2 and s3_055. It found two problems, both fixed:

- NPT flipped the rhombohedral calcite cell into another set of basis vectors. `flip no` keeps the basis, and
  the volume comparison does not depend on it.
- The default 8-digit thermo output rounded the NVE total energy to 0.01 eV per cell. The output now has
  12 digits.

The smoke-test numbers do not test the gates.

Cost: with ≥ 14 Å supercells (640–720 atoms), one model on the TF backend took 50–100 ms per step with 3 jobs at
once, about 22 h per 500 ps run. On 88-atom HAP bulk, one model ran at 24 ms/step alone. With 3 jobs at once it ran
at 47 ms/step each (64 job-steps/s, 71 % GPU), and with 6 jobs at 82 ms/step (74 job-steps/s, 91 % GPU); CP2K
used 28 cores throughout. The validation therefore runs on the exploration cells, 3 jobs at once, and the
lattice comparison is cell for cell. The w8 set (8 × 550 ps) and the S3 set (10 × 250 ps) together take an
estimated 2.5 days.

The orchestrator was restarted in label-2 (2026-11-01 13:45) to load this code. One frame in progress was
recomputed. The restart rules are now:

- a training container still running when the orchestrator restarts is waited for, not restarted;
- an interrupted LAMMPS exploration run (no "Loop time") is moved to `*.interrupted-<k>` and rerun;
- a `seed-*` container left paused by a stopped label step is resumed when the next label step ends.

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
  wavefunction restart. The 8 ranks during training were not free: they slowed train-2 from 0.32 to 0.50–0.56
  s/batch (see "DeePMD training and CPU load").
- **Result:** see the MatterGen run README (S4 row and table) and `screening/s4_dft_hull.csv` there.

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
