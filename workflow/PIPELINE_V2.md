# Pipeline v2: GPU-first MatterGen → CP2K → DeePMD → LAMMPS

Status: **in progress (started 2026-10-29)**. Nothing produced by this pipeline is
manuscript-eligible until its gate report and run manifest are reviewed.

Run directory: `.research-runs/pipeline-v2/` (not in Git). Orchestrator: `python workflow/pipeline_v2.py
--run .research-runs/pipeline-v2 --cutoff 800 --rel-cutoff 60 --iters 6 --gpu-jobs 2 --explore-jobs 3
--warm-steps 100000 --extra-systems s3`, with `--npt-temperatures 300` and `--extra-systems s3 metad` added later
(see below); it resumes from `pipeline_state.json`. From label-5 the campaign moves to a Slurm cluster with `--iters 8`
(arguments in `workflow/hpc/pipeline_v2_args.sh`, handover in `workflow/hpc/PIPELINE_V2_HPC.md`). A snapshot
from 2026-11-04 07:27 UTC (iterations 0–2 done, label-3 at 187/208) is published in `results/pipeline-v2-20261104`. Log: `logs/pipeline_v2.log`. The DFT hull of the
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
| Iteration 3 | labels 207/208, done 2026-11-04 02:39 (resumed at 104/208 after the 18:21 restart). `candidates_00127` (HAP(001) + water + CO2, 500 K, step 35 500, max_devi_f 0.104 eV/Å) failed with SCF not converged after 1050 steps; its geometry is normal (169 atoms, shortest distance 0.92 Å, both CO2 intact), and it is kept and excluded like the two iteration-2 failures. train-3 started 02:39 with 40 training systems (the 10 S3 survivors included) and 21 validation systems, 2 models at once with expandable segments: 0.157 s/batch for each model, 13.0 GB of GPU memory and 95 % utilisation for the pair. Models 0–1 ran 02:39–07:06 and models 2–3 07:06–11:27, both pairs with the CPU idle. Sealed test at 11:36: on the 83 original frames all 4 models pass; on all 123 frames (S3 test frames included) model 0 fails the energy gate at 2.10 meV/atom (see "Iteration-3 sealed test"). Exploration (18 systems × 2 temperatures = 36 runs, 80k steps, 3 jobs at once) ran 11:34–21:42 (10.1 h), every run complete: mean accurate 0.496 and 746 candidates. The fall from 0.838 is the S3 systems in their first exploration: their 20 runs average 0.193 and give 571 candidates, while the 16 original runs average 0.875 (0.838 in iteration 2) and give 175. Accurate fractions in the S3 runs run from 0.001 to 0.997 (s3_070 at 300 K). Five S3 runs have most frames above trust_hi (s3_180 at both temperatures 95–98 %, s3_294 at 300 K 57 % and at 500 K 98 %, s3_273 at 500 K 59 %), and their candidates come from the early steps only (s3_180 at 500 K: steps 50–3500). The largest deviations are 6.96 eV/Å (s3_070, 500 K), 5.17 (s3_273, 500 K) and 4.82 (s3_294, 500 K). All 746 candidates pass a geometry check (shortest distance 0.87–1.27 Å, every C with 2 or 3 O within 1.6 Å, every P with 4 O within 1.9 Å). 15 frames, from monetite at 500 K and from s3_039, s3_180 and s3_250, have one H more than 1.25 Å from any O, which looks like a proton between two O. HAP(001) + water + CO2 at 300 K rises from 0.008 to 0.468 accurate, from 0.29 in its first 10k steps to 0.61 in its last 10k. At 500 K it falls from 0.106 to 0.021 (largest deviation 2.34 eV/Å, 1.7 % above trust_hi). In every 500th frame of both runs, both CO2 keep 2 O and no H, and all 25 H2O are intact, so again no conversion |
| Iteration 4 | label-4 started 2026-11-04 21:42: 746 frames, 54–224 atoms, 20 ranks. The first two (160 atoms) took 213–218 s. From iteration-3 timings per atom count, the whole set should take about 45 h, so roughly until the evening of 2026-11-06. On 2026-11-05 the orchestrator was restarted twice (01:13 and 01:16, at 65/746; `candidates_00065` restarts): `--extra-systems s3 metad`, the SCF fallback in `run_labels.sh` and the label energy screen (see "Wrong SCF states"). The GPU runs the first biased exploration meanwhile. Outer SCF capped at 6 loops from `candidates_00068` (01:45); by 2026-11-05 21:54, 359/746 done and 3 HAP(001) + water + CO2 frames (`00085`, `00091`, `00110`) recovered by the CG fallback, each 0.33–0.36 Ha below where DIIS had stalled. Remaining 387 frames are all S3 bulk phases. 2026-11-06 06:40: 487/746, none failed; at the 00:29–06:40 rate (about 16 per hour) the rest takes about 16 h. 17:49: 688/746, none failed. A fourth CG recovery, the first outside the interface: `candidates_00388` (s3_070 CaCO3 bulk, 160 atoms), OT DIIS stalled at −2902.880 Ha at the 6-loop cap, CG reached −2903.218 Ha (0.338 Ha lower). Done 2026-11-07 00:10: 746/746, none failed. From 2026-11-06 17:35 the step shared the CPU with `seed-metad`, which the seed scheduler resumed because of the NPT smoke test run by hand (see "NPT exploration of bulk cells"): 7–9 frames per hour instead of 16–23, so the last 58 frames took 6.4 h (17:49–00:10) instead of about 3 h. Screen: 14 flagged. `candidates_00065` was a wrong state (see "Wrong SCF states"); the other 13 (s3_273 ×2, s3_294 ×11, 10–21 meV/atom from their system medians) reproduced to within 1e-8 meV/atom when relabelled and stay in. Relabels 00:13–02:12 at 20 ranks; second screen 0 flagged; 746 collected. train-4 started 02:14 (models 0 and 1 first). The orchestrator was restarted at 02:24 with `--npt-temperatures 300` and waited for `dp-train-4-0` and `-1` instead of restarting them. Models 0–1 done 06:38, 2–3 at 11:02 (4 × 100k warm-start steps, 59 training systems). Sealed test 11:08: on the 83 frames of iterations 0–2 all four models pass and improve; on all 123 frames only model 2 is at the energy gate, because S3 energies got worse while S3 forces improved (see "Iteration-4 sealed test"). Recheck with the iteration-4 models: `candidates_00304` was a wrong state the label screen had missed; relabelled and `iter-004/data` collected again (see "Wrong SCF states"). explore-4 ran 11:08 to 01:22 on 2026-11-08 (14.2 h; 36 runs of 100k steps, bulk systems in NPT at 300 K for the first time, 3 at once), every run complete: mean accurate 0.721, 661 candidates. The 16 original runs average 0.884 (0.875 in explore-3) and give 156 candidates; the 20 S3 runs average 0.591 (0.192) and give 505. HAP(001) + water + CO2 is still the gap: 0.349 accurate at 300 K (0.468 in explore-3) and 0.026 at 500 K (11 % of frames above trust_hi). Largest deviations: 3.59 eV/Å (s3_039, 500 K) and 2.65 (HAP(001) + water + CO2, 500 K). NPT cells and the s3_055 transformation: see "NPT exploration of bulk cells" |
| Iteration 5 | label-5 started 2026-11-08 01:22: 661 frames, 80–224 atoms (240 of them 120 atoms), 20 ranks. From the label-4 timings per atom count without CPU sharing it takes about 40 h, so until about 17:30 on 2026-11-09. The orchestrator was restarted at 01:24 with the same arguments to load the recheck before train-5 (see "Wrong SCF states"); no frame was done, and the label container the stopped orchestrator left running was stopped first so that no frame is labelled twice (`candidates_00000` restarts from the atomic guess). Frames 00000–00002 took 3.6–3.9 min. The GPU runs the stage-7 preview on the iteration-4 ensemble meanwhile (`preview-iter4/`, see "Stage-7 preview on the iteration-4 ensemble"). Stopped on the workstation 2026-11-08 12:20 at 144/661, none failed, 2 recovered by the CG fallback (`candidates_00072` and `00110`, both HAP(001) + water + CO2); `candidates_00144` was running and restarts from the atomic guess. The campaign continues on a Slurm cluster from here, and `--iters` goes from 6 to 8 (see "Iterations 6 and 7"). `pipeline_state.json` is unchanged (next step label-5); the 144 labels are kept and the cluster's label-5 skips them (`workflow/hpc/PIPELINE_V2_HPC.md`) |
| Biased exploration, round 1 | option A approved 2026-11-05. 9 runs on the iteration-3 models, 00:52–05:37, all complete: 180 training + 24 sealed-test frames selected, 38 + 6 with a C bonded to a third O. CP2K labels in `seed-metad`, paused except during exploration (see "First round" under "Biased exploration of CO2 conversion"). Beside label-4 (2026-11-06 17:35 to 2026-11-07 02:14, see Iteration 4) it labelled all 24 sealed-test frames (17:47–22:09) and then 24/180 training frames (22:47–02:24, the last in the minute the orchestrator restart resumed it); not in train-4. Ran again during explore-4 (11:08 to 01:22 on 2026-11-08) and reached 165/180 training frames; paused by label-5 since. `seed-metad` was stopped with label-5 on 2026-11-08 12:20 (`metad_selected_00165` was running and restarts). On the cluster the cpu phase labels the other 15 before train-5, then screens and collects the set (the wrong state `metad_selected_00000` is relabelled then) |
| Label fixes, iterations 0–3 | 6 converged labels in a wrong SCF state and 3 SCF failures found 2026-11-05; `candidates_00127` recovered by the CG fallback (ground state, 2.35 meV/atom from the ensemble); the other 8 are being rerun at 4 ranks, then iterations 0–3 are collected again before train-4 |
| MatterGen S3 survivors in training | MACE exploration done; training frames 90/90 and test frames 40/40 labelled (`seed-s3` exited 2026-11-02 about 11:00); deferred by train-2, joined train-3 on 2026-11-04 (see "MatterGen survivors in the pipeline") |
| MatterGen S4 DFT hull | done 2026-11-02 22:51: 44/44 cell optimisations converged, every final single point converged. 8 of 10 survivors pass E_hull ≤ 50 meV/atom; s3_262 (60.9) and s3_180, the only carbonate-phosphate (68.0), fail. None is on the hull. Table and comparison with S2 in the MatterGen run README. In NPT at 300 K the iteration-4 model turns s3_055 into the aragonite-type CaCO3 of ref-187 within 25 ps (see "NPT exploration of bulk cells") |
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

### Iteration-3 sealed test (2026-11-04, `iter-003/train/sealed_test.json`)

From iteration 3 the sealed test also holds the 40 S3 test frames (4 per survivor, from the seed-777 MACE runs; see
"MatterGen survivors in the pipeline"). Same weighting as above. The first three columns are the 83 frames used in
iterations 0–2, so they compare across iterations.

| Model | Force, all 83: it 2 → it 3 (eV/Å) | Force, interfaces: it 2 → it 3 (eV/Å) | Energy, 83: it 2 → it 3 (meV/atom) | S3, 40 frames: force / energy | All 123 frames: force / energy |
| --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.067 → 0.064 | 0.098 → 0.093 | 1.59 → 1.81 | 0.084 / 2.60 | 0.072 / 2.10 |
| 1 | 0.066 → 0.064 | 0.096 → 0.092 | 1.64 → 1.59 | 0.082 / 1.73 | 0.071 / 1.64 |
| 2 | 0.066 → 0.064 | 0.097 → 0.092 | 1.33 → 1.47 | 0.084 / 1.53 | 0.072 / 1.49 |
| 3 | 0.066 → 0.064 | 0.098 → 0.093 | 1.67 → 1.50 | 0.082 / 1.78 | 0.071 / 1.60 |

On the 83 frames all four models pass the gate, and both force measures drop again. Energy moves both ways (two
models up, two down). On all 123 frames, model 0 fails the energy gate (2.10 meV/atom against 2.00), so 3 of 4
models pass the sealed test as it now stands. It is an intermediate ensemble, and the gate applies formally to the
final one. HAP(001) + water + CO2 improves in force (0.104–0.107 → 0.100–0.101 eV/Å) but not in energy
(2.8–3.9 → 3.2–4.3 meV/atom).

Per S3 system (range over the 4 models), with only 4 frames each and 9 training frames per system so far:

| System | Force (eV/Å) | Energy (meV/atom) |
| --- | ---: | ---: |
| s3_014 Ca2CO5H2 | 0.056–0.061 | 0.92–1.67 |
| s3_039 Ca2CH6O7 | 0.094–0.098 | 0.73–1.91 |
| s3_055 CaCO3 | 0.036–0.040 | 0.51–0.83 |
| s3_070 CaCO3 | 0.070–0.074 | 0.82–1.78 |
| s3_180 Ca2CH4O11P2 (carbonate-phosphate) | 0.101–0.107 | 0.74–1.50 |
| s3_250 Ca2H5O7P | 0.084–0.088 | 0.77–2.14 |
| s3_262 Ca3H2O9P2 | 0.063–0.069 | 0.42–1.25 |
| s3_273 CaH6O9P2 | 0.097–0.101 | 1.47–4.18 |
| s3_294 CaH4O8P2 | 0.093–0.102 | 1.59–5.71 |
| s3_302 CaH3O5P | 0.069–0.074 | 0.55–3.37 |

The hydrogen-rich phosphates (s3_273, s3_294, s3_302) carry most of the energy error, and the 4 models disagree widely on
them, so these numbers are mostly a measure of how few S3 frames training has seen. s3_180, the only
carbonate-phosphate, has the largest force error of any system in the test.

Cell fix: CP2K 2026.2 aborts with "No meaningful cell information" on ASE's rock-salt primitive
vectors, so all CaO frames failed at first. `cp2k_label.py prepare` now rotates every new frame to
the lower-triangular standard form, and it leaves finished labels untouched.

### Iteration-4 sealed test (2026-11-07, `iter-004/train/sealed_test.json`)

Same 123 frames and weighting as iteration 3. Train-4 added the 746 iteration-4 labels, 571 of them from S3 phases
(DP-explored with the iteration-3 models).

| Model | Force, 83: it 3 → it 4 (eV/Å) | Force, interfaces: it 3 → it 4 (eV/Å) | Energy, 83: it 3 → it 4 (meV/atom) | S3, 40: force it 3 → it 4 | S3, 40: energy it 3 → it 4 | All 123: force / energy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.064 → 0.063 | 0.093 → 0.089 | 1.81 → 1.14 | 0.084 → 0.067 | 2.60 → 3.47 | 0.064 / 2.19 |
| 1 | 0.064 → 0.063 | 0.092 → 0.089 | 1.59 → 1.29 | 0.082 → 0.066 | 1.73 → 3.70 | 0.064 / 2.36 |
| 2 | 0.064 → 0.063 | 0.092 → 0.090 | 1.47 → 1.19 | 0.084 → 0.067 | 1.53 → 3.05 | 0.064 / 1.995 |
| 3 | 0.064 → 0.062 | 0.093 → 0.089 | 1.50 → 1.14 | 0.082 → 0.068 | 1.78 → 3.40 | 0.064 / 2.15 |

On the 83 frames of iterations 0–2 every measure improves and all four models pass. HAP(001) + water + CO2 goes
from 0.100–0.101 to 0.096–0.098 eV/Å and from 3.2–4.3 to 2.2–2.8 meV/atom; portlandite from 1.2–1.6 to 0.4–0.5
meV/atom. S3 forces improve in every system (s3_180 0.101–0.107 → 0.079–0.081 eV/Å). S3 energies get worse, so on
all 123 frames only model 2 passes the energy gate (1.995 meV/atom; iteration 3: 3 of 4).

**Why S3 energies got worse** (`diagnostics/s3-energy-it4/bias.py`, CPU). Per frame, the iteration-4 error on the S3
test frames is a constant offset per system, with little scatter: mean +0.0 to +7.3 meV/atom (model above DFT), standard
deviation 0.2–1.2 within a system. The largest offsets are s3_273 (+6.2 to +6.9), s3_294 (+5.1 to +7.3) and s3_180 (+2.9
to +4.9). Model 0 shows the same offsets on the 9 MACE training frames of each S3 phase (s3_273 +7.3, s3_294 +7.6,
s3_180 +4.7), but not on the DP-explored training frames of the same phases (means −1.6 to −0.1). Those 35–60
DP-explored frames per phase outnumber the 9 MACE frames, and train-4 places the MACE group up to 7 meV/atom too
high relative to them. The sealed-test frames come from MACE trajectories (seed 777) and share that offset. A likely
reason, not checked: the iteration-3 models, extrapolating on S3, drove the DP-explored frames into geometries far
from the MACE trajectories. s3_055 (+1.5 to +3.4; iteration 3 −0.2 to +0.6) also had a wrong-state label in
train-4, `candidates_00304`, 0.33 Ha too high (see "Wrong SCF states"), which may account for part of its offset.
Explore-4 now samples S3 cells in NPT at 300 K, closer to the MACE region; train-5 shows whether that removes the
offset. Nothing is changed for it.

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

The two HAP(001) slab systems have no stress labels on purpose: their cells have vacuum along c, so `cp2k_label.py
collect` writes no `virial.npy` for them. DeePMD drops the virial term for a batch whose system has no virial, and
prints `rmse_v = nan` for that batch in the log and in `lcurve.out`. That is expected, not a diverged loss. The
slabs hold 223 of 562 training frames in train-2 (40 %) and 330 of 840 in train-3 (39 %). In the train-2 model-0
`lcurve.out`, `rmse_v_trn` is `nan` in 31 of 102 rows and `rmse_v_val` in 89 of 102, so the validation virial curve
is mostly empty and says little.

### Wrong SCF states in converged labels (found 2026-11-05)

While preparing the SCF retry that the biased-exploration plan asks for, the energies of all finished labels were
compared with the iteration-3 ensemble (`workflow/screen_labels.py`; residual (E_model − E_DFT)/N against the
median of its system). The other labels of a system scatter by a median absolute deviation of 0.1–1.0 meV/atom.
Six labels sit 50–86 meV/atom too high, while the 4 models agree on them to about 1 meV/atom:

| Label | System | Deviation (meV/atom) | Outer SCF loops |
| --- | --- | ---: | ---: |
| `iter-000/hap001_water_co2_selected_00020` | HAP(001) + water + CO2 | −52.1 | 6 |
| `iter-000/portlandite_CaOH2_selected_00021` | portlandite | −85.8 | 16 |
| `iter-001/candidates_00074` | HAP(001) + water + CO2 | −53.9 | 6 |
| `iter-002/candidates_00086` | HAP(001) + water + CO2 | −53.8 | 15 |
| `iter-003/candidates_00094` | HAP(001) + water + CO2 | −58.8 | 7 |
| `iter-003/candidates_00136` | HAP(001) + water + CO2 | −56.5 | 7 |

Each passed every check in `cp2k_label.py collect`: SCF converged to EPS_SCF 1e-7, forces parsed, no force above
15 eV/Å. For a 169-atom frame the offset is about 9.5 eV. These are converged SCF solutions in a wrong electronic
state, and their forces are off too: `candidates_00136` differs from the ensemble by 0.160 eV/Å RMS, against 0.073
for its neighbour `candidates_00135`. No sealed-test, S3 or finite-size label is flagged. Each was collected into
the data of its iteration (training or validation split), so the data of ensembles 0–3 held 2, 3, 4 and 6 of them.
They are 5 of the 179 HAP(001) + water + CO2 labels outside the sealed test (2.8 %), which may account for part of that system's weakness
in the sealed test.

The SCF failure `iter-003/candidates_00127` is the same problem. Its OT DIIS run stalled at −2161.065 Ha, close to
the wrong-state labels (for example −2161.053 Ha for `candidates_00136`). Restarted from its last wavefunction with
OT CG + FULL_ALL (`.research-runs/scf-retry/`, 4 ranks), it converged in 59 CG steps to −2161.413 Ha, 9.4 eV lower.
That is 2.35 meV/atom from the ensemble, against a system median of 1.41, so the ground state.

**Fix, from 2026-11-05 (2 orchestrator restarts during label-4, at 65/746).**
- `run_labels.sh` (written by `cp2k_label.py prepare`, also installed in `iter-004/labels`): if OT DIIS fails, the frame
  runs again from its last wavefunction with OT CG + FULL_ALL. Functional, basis, cutoffs and EPS_SCF are unchanged.
  The first output is kept as `output.diis.out`. A frame that fails both gets a `failed` marker.
- After each label step from iteration 4 on (and for extra sets), the orchestrator screens the new labels with the
  newest explored ensemble. It flags |deviation| > 10 meV/atom where the 4 models agree to within 5 meV/atom, in
  systems with ≥ 5 frames. Flagged frames are relabelled once from the atomic guess with OT CG + FULL_ALL; the first
  output is kept as `output.suspect.out`. Then they are screened again.
- A relabel within 1 meV/atom of the first result confirms that result, and the frame stays. Anything still flagged is
  left out by `collect --exclude-screen` with status `energy_outlier` in `label_audit.csv`.
- Iterations 0–3 are fixed in place before train-4. The six frames are relabelled from the atomic guess, and the two
  iteration-2 failures (`candidates_00056`, `candidates_00121`) are rerun with the CG fallback (container
  `relabel-batch`, 4 ranks). Each `iter-00k/data` is then collected again, and the old directory kept as
  `data.pre-20261105`. The new collection changes the train/valid split of those iterations. Results: below.
- Label inputs cap the outer SCF at 6 loops instead of 20 (`OUTER_MAX_SCF` in `cp2k_label.py`; convergence criteria
  unchanged). The 457 converged labels of iterations 2–3 needed 2 outer loops (455) or 3 (2). A stuck DIIS run
  therefore reaches the CG fallback after 300 steps rather than 1000 (at 5–9 s per step under full CPU load, 1–2 h
  saved per stuck frame): `iter-004/candidates_00065` sat at −2161.2326 Ha through 8 outer loops. The CG reruns so far
  converged within 2 outer loops (`00127`, `00136`, `iter-002/00056`), so the cap carries over to `input.cg.inp`.
  Applied in place on 2026-11-05 at 01:45 to the 678 iteration-4 inputs from `candidates_00068` on; frames 00000–00067
  keep 20. `candidates_00065` itself converged by DIIS after 14 outer loops, at −2161.2326 Ha. It was a wrong state: the
  iteration-4 screen flagged it at −53.5 meV/atom from the HAP(001) + water + CO2 median, and the relabel from the
  atomic guess reached −2161.5679 Ha (0.335 Ha lower, like the CG recoveries), 0.4 meV/atom from the median. So DIIS
  can also converge, not only stall, in the wrong state; the screen is what catches that case.

**Results, iterations 0–3 (2026-11-05, 01:00–04:04).** All 9 frames converged. Each relabel from the atomic guess
lowered the energy by about the deviation that flagged it, and every repaired frame now lies within 1.5 meV/atom of
its system median (`iter-00k/labels/energy_screen.json`, iteration-3 ensemble):

| Label | Fix | Energy change (meV/atom) | Deviation after (meV/atom) |
| --- | --- | ---: | ---: |
| `iter-000/hap001_water_co2_selected_00020` | relabel from atomic guess | −51.3 | −0.76 |
| `iter-000/portlandite_CaOH2_selected_00021` | relabel from atomic guess | −84.4 | −1.41 |
| `iter-001/candidates_00074` | relabel from atomic guess | −54.1 | +0.18 |
| `iter-002/candidates_00086` | relabel from atomic guess | −53.2 | −0.48 |
| `iter-003/candidates_00094` | relabel from atomic guess | −57.6 | −1.17 |
| `iter-003/candidates_00136` | relabel from atomic guess | −56.1 | −0.37 |
| `iter-002/candidates_00056` (failed) | CG from last wavefunction | — | +0.01 |
| `iter-002/candidates_00121` (failed) | CG from last wavefunction | — | −1.15 |
| `iter-003/candidates_00127` (failed) | CG from last wavefunction | — | +0.87 |

The second screen flags none of the 830 labels (216, 153, 253, 208). Collected again with `--exclude-screen`, the data
of iterations 0–3 holds 216, 153, 253 and 208 frames, all `ok`; before, 216, 153, 251 and 207, with the 3 failures
listed as `missing` and the 6 wrong-state labels inside. Iterations 0 and 1 keep their frame lists and so their split;
iterations 2 and 3 gain one or two frames, which changes their train/valid split. Train-4 is the first training on
the repaired data; ensembles 0–3 were trained with the wrong-state labels in.

`relabel.sh` first ran its list with `while read`, and `mpirun` read the rest of the list from stdin, so only the
first frame of each list ran; the portlandite frame was skipped without a `failed` marker and rerun separately. The
script now loops over the list with `for` and gives `mpirun` `/dev/null` as stdin.

**Iteration 4 and the recheck before each training step (2026-11-07).** In label-4, 6 of 746 frames were in a wrong
state: four stalled OT DIIS runs recovered by the CG fallback (`00085`, `00091`, `00110`, `00388`) and two DIIS runs that
converged normally (`00065`, `00304`). The label screen caught `00065`. It missed `candidates_00304` (s3_055 CaCO3, 80
atoms, 300 K, converged by DIIS in 4 outer loops at −1451.4175 Ha): its deviation was −110.2 meV/atom, but the
iteration-3 models, which had seen S3 phases only in 90 MACE frames, disagreed on it by 5.95 meV/atom, just over the
5 meV/atom the screen requires, so it was not flagged and train-4 trained on it. The iteration-4 models flag it at
−109.5 meV/atom with a spread of 0.5 (`diagnostics/s3-energy-it4/screen_iter004_labels_it4models.json`). Relabelled from
the atomic guess with OT CG (20 ranks, 11:35–11:39) it reached −1451.7430 Ha, 0.3255 Ha (110.7 meV/atom) lower, in 2
outer loops; the first output is kept as `output.suspect.out`. `iter-004/data` was collected again (old directory
kept as `data.pre-20261107`): same frames and split, and only that frame changed (energy by 8.86 eV, force components by
up to 0.79 eV/Å). Train-5 is the first training without it.

The same iteration-4 models screened every other label set (`diagnostics/s3-energy-it4/screen_*_it4models.json`):
nothing in iterations 0–3, the sealed test, the S3 training set or the metad sealed test. S3 test sets have 4 frames
per system, below the screen's minimum of 5; per frame they scatter by 0.2–1.2 meV/atom around their system offset
(`bias.json`), so none holds a wrong state. In the biased-exploration labels, `hap001_water_co2_metad_selected_00000`
is flagged at −54.2 meV/atom (spread 1.5); that set has not been collected yet, and its own label step screens and
relabels it.

From the next training step on, the orchestrator repeats this check (`recheck` in `pipeline_v2.py`, step
`recheck-<k>` in `pipeline_state.json`): before train-k (and before the final training), every training set collected
so far (iterations and extra sets, not the test sets) is screened with the models that screened label-k, into
`energy_screen.recheck-<k>.json`. If they flag a frame that the screens at labelling time did not already leave out,
the newly flagged frames that were never relabelled are relabelled once from the atomic guess (`relabel.recheck-<k>.txt`)
and screened again. The set is then collected again, leaving out what either screen still flags; the previous data is
kept as `data.pre-recheck-<k>`. Test sets are not rechecked automatically, because leaving out frames on the models'
word would bias the test; they are checked by hand as above. Effective from the orchestrator restart during label-5 (2026-11-08 01:24).

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
The frame running then (`candidates_00104`) restarts, and finished labels are kept. Iteration 3 has 40 training
systems, so with the default allocator its 2-model pairs would most likely have stalled. With expandable segments its
first pair holds 13.0 GB in all.

The CPU-load table above was measured with the default allocator. With expandable segments and 20 CP2K ranks, a pair
runs at 0.25 s/batch, faster than the 0.32 that the default allocator gave with the CPU idle. So part of the slowdown
there blamed on CPU load may have come from the allocator. With the CPU idle, the first train-3 pair runs at
0.157 s/batch each (iteration-3 data, 40 systems, so not the same training set as the benchmark). That suggests CP2K
load still costs about a third of the speed with expandable segments, so the rule of keeping other CPU jobs to
exploration windows stays.

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
    3 O within 1.6 Å, so both CO2 and CO3 pass. The H limit is 1.45 Å since 2026-11-06 (Shawn's decision after the
    iteration-3 preview, before the formal run; see "Stage-7 preview on the iteration-3 ensemble");
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

### Stage-7 preview on the iteration-3 ensemble (2026-11-05, `preview-iter3/`)

An early look while the GPU waits for label-4, not the stage-7 result: `dp_validate.py` unchanged, iteration-3
models (conversion check ≤ 2e-13), 200 ps for every system instead of the formal 500 ps for w8, 3 jobs at once
(`preview-iter3/run_preview.sh`). w8 ran 05:39–19:00; the S3 set follows.

| System | MD | Geometry gate | NVE drift (meV/atom/ps) | Max devi (eV/Å) | Frames ≥ 0.10 | ΔV 300 K | ΔL 300 K (a, b, c) | ΔV 0 K | ΔL 0 K (a, b, c) |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- | ---: | --- |
| CaCO3_20 | 200 ps | pass | 2e-6 | 0.043 | 0.00 | — | no reference | — | — |
| calcite | 200 ps | pass | 4e-6 | 0.058 | 0.00 | +0.69 % | +0.70, +0.78, +0.75 % | −0.01 % | −0.05, −0.02, −0.01 % |
| HAP(001) + CO2 | 200 ps | pass | 3e-6 | 0.102 | 0.01 | — | slab | — | — |
| HAP(001) + water + CO2 | stopped at 154 ps | not finished | — | 0.188 | 0.53 | — | slab | — | — |
| HAP(001) + water + CO2, rerun alone | 200 ps | pass | −2.5e-5 | 0.197 | 0.56 | — | slab | — | — |
| HAP bulk | 200 ps | pass | −6e-5 | 0.076 | 0.00 | +0.84 % | +0.17, +0.17, +0.51 % | −0.31 % | −0.22, −0.16, +0.10 % |
| lime | 200 ps | pass | 2e-6 | 0.021 | 0.00 | +1.35 % | +0.41, +0.45, +0.46 % | −0.03 % | −0.03, −0.01, −0.01 % |
| monetite | 200 ps | **fail**, 105 of 201 frames | −4e-7 | 0.133 | 0.03 | +1.21 % | +0.25, +0.44, +0.40 % | 0.00 % | −0.04, +0.04, +0.25 % |
| portlandite | 200 ps | pass | 1.3e-4 | 0.059 | 0.00 | +1.50 % | +0.29, +0.30, +0.92 % | +0.76 % | +0.01, −0.01, +0.77 % |

- **Lattice.** The formal gate (HAP a and c within 1 % at 300 K) passes: +0.17 % and +0.51 %. The code's lattice flag
  tests the three lengths; the volume is reported only, and lime, monetite and portlandite are 1.2–1.5 % larger in
  volume at 300 K than the 0 K CP2K cell. The 0 K relaxations with the same model agree with CP2K to ≤ 0.31 % in
  volume except portlandite (+0.76 %, all of it in c, the interlayer axis), so most of the 300 K offset is thermal
  expansion, which the 0 K reference does not contain.
- **HAP(001) + water + CO2 stopped** at step 309200 (154.6 ps, 14:29). Its log ends mid-line, with no LAMMPS error and
  no failed frame before it; `dp_validate.py` sent LAMMPS output to /dev/null. Even before it stopped, 53 % of its
  frames had a deviation ≥ 0.10 eV/Å, so the iteration-3 ensemble is not yet converged on this system, which fits
  its sealed-test errors.
  - **Cause: the WSL GPU driver, not the model.** The docker-desktop kernel log has, at 14:28:57 (uptime 618405.6 s,
    mapped with the current uptime), a vmalloc warning inside `dxgvmb_send_submit_command_hwqueue`, then
    `dxgkio_submit_command_to_hwqueue: Ioctl failed: -12` (out of memory in the GPU paravirtualisation driver) and a
    fatal signal 6 in `lmp`. The CaCO3_20 crash of iteration 2 (signal 11, 2026-11-03 09:17) is in the same log,
    with no driver warning in the 12 s before it.
  - **Rerun alone** (same `in.lmp` and `conf.lmp`, GPU, 03:55–12:33 UTC = 21:55–06:33, output kept,
    `preview-iter3/w8-rerun/`): exit code 0, 200 ps NVT and 50 ps NVE. `analyse_rerun.py` applies the gates of
    `dp_validate.py` with the same models (`rerun_summary.json`): 201 frames, no geometry failure at 1.2 or
    1.45 Å, NVE drift −2.5e-5 meV/atom/ps, maximum deviation 0.197 eV/Å on 100 frames, 56 % ≥ 0.10 (58 % up to
    154 ps, 48 % after). GPU runs are not bitwise reproducible, so this is a second trajectory from the same start,
    not a replay. The preview result for this system is the rerun; the stopped run stays in the table.
  - `dp_validate.py` now writes LAMMPS output to `lmp.out` and the exit code to `exit_code.txt` in each run
    directory, and reports `lmp_exit_code`. A crashed run still fails the gate ("the run must also finish"); a
    rerun is a separate, recorded decision.
- **Monetite fails the geometry gate on the H rule, not on structure.** Every failure is an H whose nearest O lies
  at 1.21–1.27 Å, above the v1 limit of 1.2 Å. In the CP2K-optimised cell (ref-272) the H sit 1.04–1.06 Å from their
  donor O, in hydrogen bonds with O···O of 2.47 and 2.58 Å. Along the MD (200 frames × 32 H) the nearest O of an H
  differs from its frame-0 donor in 1993 H-frames: the protons move between donor and acceptor. At each failure the
  second O is 1.20–1.45 Å away (median 1.25 Å) and O···O is 2.32–2.66 Å (median 2.45 Å). The H is shared in a short
  hydrogen bond. P–O, Ca–O and every other check pass. The ensemble agrees there (max 0.133 eV/Å, 3 % of frames
  ≥ 0.10). The CP2K-labelled monetite frames, taken from MLIP exploration at 300–600 K, have 3.4 % of their H beyond 1.2 Å
  (maximum 1.31 Å). Whether this proton motion matches DFT dynamics is not tested here: classical nuclei, no AIMD.
  Changing the H rule for the formal run is Shawn's decision (thresholds are fixed before the final run, never
  relaxed inside `dp_validate.py`).
- **Decision (Shawn, 2026-11-06): the stage-7 H rule is 0.8–1.45 Å** for every system, taken after this preview and
  before the formal run. 1.45 Å is the limit of the metad selection gate; an H with no O within 1.45 Å still fails.
  `dp_validate.py` has it as `H_O_MAX` (test: an H at 1.25 Å passes, one at 1.5 Å fails). The preview results above
  stay as they were run, with the 1.2 Å rule. Re-checking the same w8 frames with `geometry_gate` (no new MD): at
  1.2 Å monetite again fails 105 of 201 frames, at 1.45 Å it fails none, and the other 7 systems pass both. The
  largest nearest-O distance of any H is 1.270 Å in monetite and 1.094–1.139 Å in the four other H-bearing systems.
  The S3 preview was started before the change and uses 1.2 Å; its frames are re-checked the same way when it ends.

**S3 set** (10 MatterGen survivors, 200 ps NPT + 50 ps NVE each, ran 19:35 on 2026-11-05 to 17:29 on 2026-11-06,
`preview-iter3/s3/validation_summary.json`). Geometry is given as run (H limit 1.2 Å) and re-checked at 1.45 Å on
the same frames. Δ are relative to each survivor's CP2K-relaxed cell; angles in degrees.

| System | Geometry fails, 1.2 / 1.45 Å | NVE drift | Max devi | Frames ≥ 0.10 | ΔV 300 K | ΔL 300 K (a, b, c) | Δangles 300 K | ΔV 0 K | ΔL 0 K (a, b, c) |
| --- | --- | ---: | ---: | ---: | ---: | --- | --- | ---: | --- |
| s3_014 C1Ca2H2O5 | 0 / 0 | −1.9e-5 | 0.178 | 0.99 | +4.21 % | +1.42, +3.68, +2.77 % | +4.4, +0.7, −1.0 | +0.93 % | +1.45, −0.24, +0.01 % |
| s3_039 C1Ca2H6O7 | 100 / 0 | −3.1e-5 | 0.235 | 0.99 | +2.08 % | +1.55, +2.24, −1.89 % | +0.6, −0.1, −0.8 | +0.27 % | +1.56, −0.07, −1.46 % |
| s3_055 CaCO3 | 0 / 0 | 1.7e-5 | 0.193 | 0.84 | −7.99 % | −0.65, −5.05, −0.87 % | 0.0, 0.0, −1.6 | −0.59 % | −0.39, −0.36, +0.15 % |
| s3_070 CaCO3 | 0 / 0 | 1.0e-5 | 0.141 | 0.10 | +1.93 % | −0.27, +0.76, +0.66 % | 0.0, +1.7, 0.0 | +0.21 % | −0.36, −0.07, −0.02 % |
| s3_180 C1Ca2H4O11P2 | 15 / 0 | −6.1e-5 | 1.283 | 0.99 | +0.69 % | −1.75, +11.49, −5.18 % | −3.3, +3.0, −1.9 | −0.85 % | +1.76, +0.22, −2.57 % |
| s3_250 Ca2H5O7P1 | 93 / 0 | 6.3e-5 | 0.236 | 0.96 | +2.07 % | +0.73, +0.17, +1.18 % | −0.3, +1.0, −0.6 | +0.36 % | 0.00, +0.58, −0.24 % |
| s3_262 Ca3H2O9P2 | 0 / 0 | 2.9e-5 | 0.297 | 0.99 | −2.98 % | −0.41, +0.04, +4.79 % | +0.2, +7.9, −1.0 | −0.10 % | +0.08, +0.01, +0.92 % |
| s3_273 Ca1H6O9P2 | 5 / 0 | 2.7e-5 | 0.529 | 0.99 | −3.91 % | +0.80, −3.72, −1.36 % | −2.3, +0.5, +4.6 | −5.92 % | −4.42, −2.23, −1.11 % |
| s3_294 Ca1H4O8P2 | 37 / 0 | −6.9e-7 | 0.737 | 0.99 | −1.82 % | −1.68, +18.38, −4.27 % | +7.8, −0.8, +8.0 | −2.27 % | −1.55, −1.96, +0.13 % |
| s3_302 Ca1H3O5P1 | 24 / 0 | −1.4e-5 | 0.282 | 0.97 | −3.07 % | +2.33, +0.51, +2.78 % | −12.9, +11.0, +8.3 | +1.24 % | −1.52, +0.88, +2.05 % |

- **Gates as run.** Every run finished and every NVE drift passes. All geometry failures are the H rule (nearest O
  1.20–1.28 Å); at 1.45 Å none fails. The shortest distance in any frame is 0.89 Å; no P, C or contact failure.
  Only s3_070 passes the lattice flag (lengths within 1 %); for the S3 set it is reported, not a formal gate.
- **The iteration-3 ensemble extrapolates on these phases.** Apart from s3_070 (10 %), 84–99 % of the frames have a
  deviation ≥ 0.10 eV/Å, up to 1.28 eV/Å (s3_180). The S3 phases entered exploration only in iteration 3; their
  571 candidates are in label-4. So these lattice and stability numbers do not yet test the phases themselves.
- **Cell shape over time** (thermo averages over 0–2, 2–20, 20–100, 100–150 and 150–200 ps). s3_014, s3_055 and
  s3_262 move to a new cell within 20 ps and stay there (s3_055 b −5 %; s3_262 c +4.8 % and β +7.9°). s3_180,
  s3_294 and s3_302 change over 20–100 ps (b +11.5 % and +18.7 %, angles up to 12.9°), then stay. s3_039 and
  s3_273 drift slowly through the run. s3_070 and s3_250 stay within 1.3 % and 1.7°. Whether any of these is a
  transformation of the hypothetical phase or an artefact of extrapolation cannot be told from these runs.
- **0 K.** The model's own cell relaxation is within 1.3 % in volume of CP2K for 8 of 10; s3_273 is −5.9 % (a −4.4 %)
  and s3_294 −2.3 %.
- **Exploration gap.** `dp_explore.py` ran NVT only, in the fixed starting cell. Training has the starting cells and
  ±4 % strained copies, but no frame from a cell that moved like these. From explore-4 on, bulk systems explore in
  NPT at 300 K (see "NPT exploration of bulk cells").

### Stage-7 preview on the iteration-4 ensemble (2026-11-08, `preview-iter4/`)

Same early look as `preview-iter3/` on the iteration-4 models (`graph-compress.pth`): `dp_validate.py` unchanged,
200 ps for every system, 3 jobs at once (`preview-iter4/run_preview.sh <groups>`). Started 2026-11-08 01:33 with
the S3 group only, about 22 h at the iteration-3 pace; w8 (about 13 h) follows if label-5 still leaves the GPU free,
so that the preview does not share the GPU with train-5. Not the stage-7 result. Results go here when it ends.

### NPT exploration of bulk cells (from explore-4, decided 2026-11-06)

- **Change.** `dp_explore.py --npt-temperatures 300` (orchestrator `--npt-temperatures 300`): at 300 K every bulk
  system (not `hap001*`) runs NPT at 1 atm after `change_box all triclinic`, with every length and angle free (`tri`,
  `flip no`, Tdamp 0.05 ps, Pdamp 0.5 ps, the barostat of stage 7). At 500 K bulk systems stay NVT; slabs stay NVT at
  both. Run directories are `<system>_T300_npt`, and each candidate carries `ensemble` and its own cell.
- **Why.** The S3 preview above: in NPT the S3 cells moved far from the fixed cells that NVT exploration samples.
- **Cost.** The number of runs (36) and the cap of 30 candidates per run are unchanged, so the upper bound on labels
  per iteration is too. The mean accurate fraction now mixes NPT and NVT runs and is not strictly comparable with
  iterations 0–3.
- **Smoke test** (`npt-smoke/`, iteration-3 models, 2 ps, hap_bulk, lime, s3_294 and the slab hap001_co2, not used as
  data). Python re-evaluation of the NPT frames matches LAMMPS to ≤ 1.2e-5 eV/Å (tolerance 0.005), so the moving
  cell is read correctly. hap_bulk and lime: 100 % accurate in NPT and in NVT. s3_294: in NPT 1.2 % accurate, 86 %
  candidate, 12 % failed; in NVT at 500 K 44 % candidate, 54 % failed. The slab ran NVT at both temperatures.
  LAMMPS warns that the triclinic skew is large for the hexagonal and rhombohedral cells (speed only).
- **Start.** Restarted 2026-11-07 02:24 with `--npt-temperatures 300`, ten minutes into train-4: the orchestrator
  waited for `dp-train-4-0` and `-1` (nothing recomputed). Its start-up resumed `seed-metad`; the seed scheduler
  paused it again at 02:24:58.
- **Side effect of the smoke test.** The seed scheduler counted any `dp_explore.py` container as an exploration
  step, so the smoke test resumed `seed-metad` at 17:35 during label-4, which then ran at half speed until train-4
  (see Iteration 4 in the status table). Since 2026-11-07 the scheduler (`metad/seed_scheduler.sh`) also keeps
  `seed-metad` paused while any other CP2K container (a label or relabel step) runs.
- **Result (explore-4, 2026-11-07 11:08 to 2026-11-08 01:22).** Every NPT run finished. In the original bulk
  systems the cells stayed close to their start: second-half mean volume −1.6 to +2.1 % of the start, and within
  1.3 % of where a linear fit of the DFT potential pressure of their labels against volume crosses −P_kin
  (hap_bulk: DFT +2.55 GPa at the start cell, NPT +2.0 %, fit +2.0 %, which also checks the sign of the virial).
  S3 accurate fractions at 300 K rose from 0.001–0.997 in NVT (explore-3) to 0.004–0.999 in NPT (s3_014 0.388 →
  0.978, s3_250 0.038 → 0.964); s3_180 (0.004) and s3_302 (0.097) are still mostly candidates.
- **The S3 start cells are not the DFT cells.** Exploration starts from the MACE-MPA-0 geometry of each survivor,
  1.7–6.8 % larger in volume than its CP2K-relaxed cell in the DFT hull (`preview-iter4/s3/reference_cells.extxyz`).
  Against the start cell the S3 cells shrank by up to 9 %; against the CP2K cell, six of ten are 0.5–2.8 % larger
  at 300 K (thermal expansion of the size seen in the old bulk systems in `preview-iter3/`). The other four are
  s3_055 −7.5 %, s3_302 −3.0 %, s3_180 −2.1 % and s3_014 +5.9 %; all four also moved in the iteration-3 preview.
  The negative DFT pressure of the S3 labels at the start cells (−0.4 to −2.6 GPa) is therefore expected and is no
  evidence of a transformation by itself.
- **Model pressure against DFT** (`diagnostics/npt-volume-it4/pressure.py`, mean of the four iteration-4 models,
  every labelled bulk frame): mean difference −0.26 to +0.16 GPa and RMSE 0.04–0.28 GPa per system (largest for
  s3_055). The same script also scans one S3 frame affinely from 0.88 to 1.06 of its volume; without relaxing the
  internal coordinates that scan compresses the C–O and P–O bonds too and overstates the stiffness, so it says
  nothing about how the cells moved. What the trajectories do show is that the S3 cells changed shape as well as
  volume (s3_302: α 103° → 90°, γ 104° → 114°; s3_180: c −10 %), mostly within the first 5 ps, and then held.
- **s3_055 becomes aragonite-type CaCO3** (`diagnostics/s3_055-transform/`). In the NPT run its Ca go from 6- to
  8.75–8.94-coordinated (Ca–O < 3.0 Å) within the first 25 ps, at 58.3–59.0 Å³ per CaCO3; calcite in the same run
  type stays 6-coordinated. Relaxed at 0 K with model 0 (cell free), the last frame is Pnma with 9-coordinated Ca
  at 57.65 Å³ per CaCO3, and pymatgen's StructureMatcher matches it to the DFT-hull reference ref-187 (Pnma,
  9-coordinated, 57.40 Å³; rms 0.012 Å), which is not s3_055 (P6₃22, 6-coordinated, 63.48 Å³). The same
  relaxation of the start frame returns to s3_055 (rms 0.002 Å), so s3_055 is a local minimum of the model and the
  change at 300 K crosses a barrier the model puts within reach of a few ps. In the DFT hull (our protocol) ref-187
  is 6.4 meV/atom above the hull phase ref-186 (R-3c, 6-coordinated) and s3_055 is 15.8 above it, so the change
  lowers the DFT energy by about 9.4 meV/atom; the model puts the relaxed structures 10.1 meV/atom apart and within
  0.4 meV/atom of the DFT energies of s3_055 and ref-187. The energies are DFT; the rate and the path are the
  model's and are not checked with DFT. Label-5 holds 11 candidates from this run (and 9 from the 500 K NVT run),
  whose DFT forces and stress test the model on the transformed structure.

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
- **Iterations 6 and 7 (decided 2026-11-08).** `--iters` goes from 6 to 8. With 6, the last labelled candidates
  were those of explore-4: explore-5 only measured the train-5 ensemble, while explore-4 left the interface at 0.349
  (300 K) and 0.026 (500 K) accurate and s3_180 and s3_302 far from 0.99. With 8, explore-5 and explore-6 are
  labelled (label-6, label-7); explore-7 measures the train-7 ensemble. The interface did not improve from
  iteration 3 to 4, so two more iterations are not expected to close that gap on their own.
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

## Biased exploration of CO2 conversion (drafted 2026-11-04, option A approved 2026-11-05, first round explored, labels queued)

Status: Shawn chose option A on 2026-11-05 and asked for the plan to be adjusted as needed. The first round runs on
the iteration-3 models. What was built and how it departs from the draft is under "First round" at the end of this
section. The first round's runs finished 2026-11-05 05:37; its frames wait for CP2K labels. No biased frame has
entered training yet.

**Aim.** Collect CP2K labels along CO2 + H2O/OH⁻ → HCO3⁻/CO3²⁻ at the hydrated HAP(001) surface, so that the
potential has seen the reactive configurations before any mechanism study. Unbiased exploration cannot provide
them. In iterations 0–2, HAP(001) + water + CO2 never converted a CO2, and the iteration-2 runs (30 ps each) show
why. In both, each C stays at least 2.0 Å from any O other than its own two.

- **300 K:** one CO2 stays at the water–vacuum surface, 8.7–11 Å from the nearest slab Ca. The other moves into the
  film and stops about 6 Å from the slab.
- **500 K:** in the sampled frames both CO2 are above the film by step 7 500 (3.75 ps). From step 22 500 (11 ps)
  on, they sit mostly 2.2–4.1 Å from O atoms of the bare bottom face of the periodic image; one frame has a C at
  6.0 Å.

The bare bottom face was found while drafting this plan (`traj.dump` every 150th frame;
scratch scripts, not in Git). The 500 K run therefore samples CO2 on a dry HAP face, which is the chemistry of
`hap001_co2`, not of the hydrated interface. That the CO2 never reached the hydrated surface also bears on
the "no CO2 conversion" statement: it was not tested at the surface.

**Literature basis** (all `candidate` in `literature/references.csv`, not verified, abstracts only):
- metadynamics-driven active learning of a reactive NN potential in water (`yang2022metadynamics`);
- WT-metadynamics of CO2 mineralization on portlandite under a nanometric water film (`zhu2025portlandite`);
- a DFT-level CaCO3 potential covering bicarbonate chemistry (`piaggi2025caco3`).

See `literature/novelty_check_2026-11.md`.

**Software (checked 2026-11-04).**
- LAMMPS in `hap-mlip:2026.10` lists both `PLUMED` and `COLVARS` among its installed packages (`lmp -h`).
- The PLUMED binary is `/opt/deepmd-kit/bin/plumed`, version 2.9.2. The `bias`, `colvar` and `function` modules
  are on; `opes` and `ves` are off. Well-tempered `METAD` (`BIASFACTOR`) is therefore available; OPES would need a
  PLUMED rebuild.
- A biased run is the `dp_explore.py` input plus `fix pl all plumed plumedfile plumed.dat outfile plumed.log`.
  That keeps `neighbor 0.0 bin` (see "LAMMPS neighbour skin"), the 4 TF models with `out_freq`, NVT (Tdamp
  0.05 ps) and the 0.5 fs step.

**Atom groups** come from `exploration/w8/hap001_water_co2_initial.extxyz`, with C2H54Ca20O81P12 in 169 atoms.
PLUMED indices are 1-based:

| Group | Atoms |
| --- | --- |
| C | 146, 152 |
| CO2 O (C–O < 1.4 Å) | 147, 148, 153, 154 |
| O other than CO2 | 77: 48 phosphate, 25 water, and the 4 channel OH⁻ (O 43, 44, 87, 88) |
| H | 54 |
| Height reference | top slab Ca 46 (z 13.88 Å) |

In the initial structure the slab spans z 0.33–13.88 Å (Ca) and the water film z 15.7–24.6 Å (O). The cell is
33.74 Å along z, so about 10 Å of vacuum separates the top of the film from the bare bottom face of the next image.
The 4 OH⁻ lie in the HAP channels at z 1.9–11.5 Å, not on the top surface. The film holds no surface OH⁻, so the
fast OH⁻ + CO2 route is open only after a proton transfer, or with a hydroxylated termination, which is a separate
structure to build.

**Collective variables (candidates).**
1. `cn`: `COORDINATION GROUPA=C GROUPB=O_other R_0=1.6 NN=6 MM=12`, the number of new C–O bonds. It is near 0 for
   CO2 and near 1 per C for HCO3⁻ or CO3²⁻. The main biased CV.
   - In the iteration-2 runs the per-C value stays at or below 0.15 at 300 K and 0.23 at 500 K. The sum over the two C
     starts at 0.26 because the 6/12 switching function has a long tail.
   - To test before production: `D_MAX`, or a steeper function (`NN=8 MM=16`), so that the baseline drops and a
     formed bond counts near 1.
2. `ch`: `COORDINATION GROUPA=O_CO2 GROUPB=H R_0=1.2`, the protonation of the former CO2 oxygens (HCO3⁻ or H2CO3).
   It is a second biased CV in a 2D run, or printed only.
3. Height of each C above Ca 46: `DISTANCE ATOMS=46,146 COMPONENTS`, then `.z`. It is used in walls, not biased.

A single-C form of `cn` (one METAD per C) may separate the two CO2. With the summed form, the bias cannot tell
which C reacted.

**Walls.**
- `UPPER_WALLS` on both heights at about 10.5 Å, which is below the film top at about 10.7 Å above Ca 46. This keeps
  the CO2 in the water instead of letting it desorb onto the dry image face, as at 500 K.
- `LOWER_WALLS` at about 2 Å stops a C from being pushed into the slab.
- An optional second set of runs uses a moving restraint (`MOVINGRESTRAINT`) on the height. It brings one CO2 into
  contact with the surface water layer first, so that conversion is sampled at the surface and not only in bulk-like
  water.

**WT-METAD settings to test** (trial values, not taken from a paper):

| Setting | Trial value | Note |
| --- | --- | --- |
| Temperature | 300 K and 400 K | 500 K desorbed the CO2 in the unbiased run; the walls remove that, but 400 K is the safer upper value |
| `HEIGHT` | 0.02 eV | about 0.8 kT at 300 K |
| `SIGMA` | 0.05 (CN units) | to be reset from the unbiased CV fluctuation after the `D_MAX` change |
| `PACE` | 50 steps (25 fs) | the same as the dump and `model_devi` stride, so each dumped frame carries its current bias |
| `BIASFACTOR` | 20 | ΔT about 5700 K at 300 K; to be revised once the first barrier crossing is seen |
| Length | 100 000 steps (50 ps) per run | stop early if the geometry gate fails |

**Frame selection and what the bias means.**
- Selection is unchanged: max_devi_f between trust_lo 0.10 and trust_hi 0.30 eV/Å, at most 30 frames per run, evenly
  spread in time.
- Every candidate also records `cn`, `ch`, both heights and `mtd.bias` (from `COLVAR` at the dump step). The labelled
  set can then be checked for coverage along the CV.
- Frames above trust_hi are counted, not labelled, as now. A biased run may cross trust_hi early. The run is not
  discarded: its candidates before the crossing still count.
- No free energies or rates come from these runs. The potential changes between iterations, and the bias is
  deposited on a model that is not yet converged, so WT reweighting is not valid. Free energies come only from
  production runs with the final ensemble.

**When and how many.**
- **Option A, recommended:** start after train-3, as an extra exploration beside explore-3. That gives 4 runs:
  300/400 K × 2 seeds, with the summed `cn` and walls.
- **Option B:** wait until the unbiased hydrated interface reaches about 0.5 accurate at 300 K. This risks never
  starting, because the unbiased runs do not reach the surface.

Either way it starts as a separate script (`dp_explore.py` with an optional `--plumed` template, or a sibling
script). It writes its own `explore_summary.json` and `candidates.extxyz`, merged into the next label step.

**Cost**, from measured numbers:
- **GPU:** 3 TF 4-model jobs at once run at about 57 ms/step each (see "DeePMD inference speed"), and the PLUMED
  overhead is negligible (smoke test). One 100k-step run is about 1.6 h, so 4 runs at 3 at once take about
  3.2 h of GPU.
- **CPU:** labels are up to 30 per run, 120 for 4 runs. At about 5 min each with 20 ranks, that is about 10 h of
  CP2K.

**Risks.**
- **SCF failures.** Reactive and high-energy frames are the ones most likely to fail SCF. Three interface frames
  already have: `candidates_00056` (hap001_co2), and `candidates_00121` and `candidates_00127` (hap001_water_co2),
  each after 1050 steps. Excluding them biases the training set against exactly these frames.
  - Before biased labels start, one of the failed frames is retried with a more robust SCF setting (for example
    another OT minimizer or preconditioner, or more steps from the saved wavefunction). The original failure is
    kept in any case.
- **Unphysical frames.** A frame can be over-biased into an unphysical state, such as C–O < 1.1 Å or a fragmented
  phosphate. The trust_hi cap and the stage-7 geometry gate (every C with 2–3 O within 1.6 Å, every H with its
  nearest O at 0.8–1.2 Å; 1.45 Å since 2026-11-06) are checked on each candidate before labelling. A failure is
  recorded.
- **Box size.** The 2 CO2 and 25 H2O in a 9.4 Å × 9.4 Å cell give a high local CO2 concentration, and the ions that
  form carry a charge. Conversion in this box need not represent the dilute case. This limits any mechanism claim,
  not label coverage.
- **Wrong CV.** If the CV is incomplete (for example, the proton leaves through a different path), the bias may
  produce hysteresis without crossings. The check is to watch `ch` and the O–H distances along the run.

**Acceptance.** The interface counts as "covered" for reactive configurations only when all of the following hold:
1. At least one C reaches `cn` ≥ 0.8 with a geometry-gate pass in at least two independent biased runs.
2. Labelled frames span the CV range from CO2 to HCO3⁻/CO3²⁻, with at least 5 labels in each CV bin of width 0.2.
3. A separate seed-777-style biased run, never used for selection, supplies about 20 sealed-test frames. The sealed
   test reports them as their own group, against the same gate (≤ 0.12 eV/Å interface force RMSE).
4. A later biased exploration with the next ensemble shows max_devi_f below trust_hi for at least 90 % of frames.

Until then, nothing about the mechanism is claimed.

**Smoke test (2026-11-04, `.research-runs/metad-smoke/`, not in Git).**
- **Setup:** the iteration-2 TF models with `make_inputs.py` (in that directory) on hap001_water_co2 at 300 K, run
  while the GPU was idle during label-3.
- **What works:**
  - `fix plumed` loads in `hap-mlip:2026.10`, and a 400-step run completes (`run-400`).
  - At step 0, `cn` is 0.265 (sum over both C), `ch` 0.185, and the heights 9.55 and 9.34 Å.
  - After 400 steps (10 hills) the bias is 0.124 eV, as expected for 10 hills of 0.02 eV at an almost unchanged CV.
  - Pair time is 29.5 ms/step alone, and PLUMED (`Modify`) takes 0.6 % of the run.
  - The largest max_devi_f over the 9 output steps is 0.188 eV/Å, already in the candidate window.
- **Bug found and fixed in the smoke input:** `DISTANCE ATOMS=a,b` returns b − a, so the first input gave negative
  heights, and the upper wall at +10 Å would never have acted. With the reference atom first (`run-200-sign`),
  the heights are +9.55 and +9.34 Å.
- **Warnings:**
  - LAMMPS warns that the triclinic skew is large (the hexagonal cell). This is already the case in the unbiased runs.
  - PLUMED warns that the `DISTANCE` components are not periodic. That is harmless while the walls keep the C within
    about 10 Å of Ca 46 (half the cell is 16.9 Å).
- **What it does not show:** a 0.2 ps run says nothing about crossings or settings.

**First round (2026-11-05, `3.ML_AIMD/Step3_mlip_deepmd/dp_explore_metad.py`, `.research-runs/pipeline-v2/metad/iter-003`).**
- **CV, chosen from data instead of the draft's trial value.** The per-C coordination of the iteration-3 unbiased
  runs of this system was computed offline, with the same groups.
  - The draft's `NN=6 MM=12` sits at 0.038 ± 0.022 per C at 300 K.
  - A `D_0`-shifted form is 0 beyond 2 Å, so it would give the bias no gradient to pull with.
  - `R_0=1.6 NN=8 MM=16` sits at 0.008 ± 0.006 (300 K) and 0.010 ± 0.015 (500 K), gives 0.69–0.86 for an O at
    1.27–1.45 Å, and still has a slope at 2.4–2.8 Å. It is used per C (`cn1`, `cn2`), so the bias knows which C
    reacted.
  - PLUMED's value at step 0 (0.0386, 0.0326) equals the offline one.
- **Bias:** one well-tempered METAD per C: SIGMA 0.04, HEIGHT 0.02 eV, PACE 50, BIASFACTOR 20, grid −0.5 to 3.0.
  Walls at 2.0 and 10.5 Å above Ca 46 (KAPPA 2 eV/Å²). `ch` and both heights are printed, not biased.
- **Surface runs** (the draft's optional moving restraint): C1 is pulled to 3.5 Å above Ca 46 over 10 ps
  (KAPPA 0.5 eV/Å²) and held there. Its METAD keeps acting.
- **Runs:** 100k steps each on the iteration-3 TF models, 3 at once, at 57.8 ms/step each (2000-step test,
  `metad/smoke-iter-003`).
  - Training: free film at 300 and 400 K with seeds 1 and 2, and surface at 300 and 400 K with seed 1.
  - Sealed test: seed 777, free film at 300 and 400 K and surface at 300 K.
- **Selection.**
  - Training frames use the dp_explore window (0.10–0.30 eV/Å, ≤ 30 per run, evenly spread in time).
  - Every selected frame passes a geometry check: shortest distance > 0.75 Å, 4 O on each P, 2–3 O on each C with
    none closer than 1.1 Å, and each H within 1.45 Å of an O. The H limit was looser than stage 7's 1.2 Å, so a
    proton in transfer is kept; the stage-7 result is stored on the frame. Since 2026-11-06 stage 7 uses 1.45 Å too.
  - Test frames: 8 per run, evenly spread in time among frames that pass the check, whatever their deviation.
  - Frames go to `exploration/metad/selected` and `exploration/metad-test/selected` with their CV values. System
    name `hap001_water_co2_metad` (a slab, so no stress); the sealed test reports it as its own system.
- **Into the pipeline:** the orchestrator runs with `--extra-systems s3 metad` since 2026-11-05 01:13.
  - A `seed-metad` container labels the frames when the CPU is free. Iteration 4 defers the set while it is running,
    and the last iteration waits for it (changed at the 01:16 restart). So the frames join train-4 or train-5, and the
    iteration-5 exploration tests an ensemble that has seen them.
  - Label energy screen and relabelling as for every other set (see "Wrong SCF states").
- **First 20 ps:**
  - In the free runs `cn` per C reaches 0.7–0.8, and `ch` 0.85–1.4: a C bonded to a water O, with the CO2 oxygens
    protonated. The C stay 8.9–10.6 Å above Ca 46.
  - max_devi_f reaches 2.7 eV/Å, so these configurations are far outside the training data and are not selected.
    The selected frames come from the approach to them.
  - None of this is a mechanism result: the bias drives an unconverged model.
- **Result, 100 ps per run (finished 2026-11-05 05:37, `metad/iter-003/metad_summary.json`).** All 9 LAMMPS runs
  completed. The Python check of model 0 against LAMMPS agrees to ≤ 7.5e-4 eV/Å. Fractions are of the 2001 dumped
  frames (accurate < 0.10 eV/Å ≤ candidate ≤ 0.30 eV/Å < failed):

  | Run | Accurate | Candidate | Failed | Max devi (eV/Å) | Geometry fail | `cn` max | `ch` max | C height (Å) | Picked | Picked with `cn` > 0.8 |
  | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: |
  | free 300 K s1 | 0.006 | 0.056 | 0.938 | 2.88 | 0.000 | 0.94 | 1.58 | 5.5–10.8 | 30 | 0 |
  | free 300 K s2 | 0.034 | 0.181 | 0.784 | 2.73 | 0.003 | 1.83 | 3.45 | 6.0–10.8 | 30 | 18 |
  | free 400 K s1 | 0.011 | 0.056 | 0.933 | 6.62 | 0.018 | 2.70 | 5.86 | 2.7–10.8 | 30 | 1 |
  | free 400 K s2 | 0.005 | 0.142 | 0.853 | 3.94 | 0.009 | 2.34 | 6.18 | 5.5–10.8 | 30 | 19 |
  | surface 300 K s1 | 0.009 | 0.071 | 0.919 | 2.88 | 0.000 | 0.94 | 1.96 | 2.4–10.8 | 30 | 0 |
  | surface 400 K s1 | 0.008 | 0.059 | 0.933 | 2.92 | 0.003 | 1.83 | 3.65 | 1.7–10.7 | 30 | 0 |
  | test free 300 K | 0.005 | 0.052 | 0.943 | 2.98 | 0.000 | 0.95 | 2.60 | 7.7–10.7 | 8 | 1 |
  | test free 400 K | 0.010 | 0.129 | 0.861 | 2.87 | 0.000 | 1.09 | 2.09 | 1.7–10.0 | 8 | 3 |
  | test surface 300 K | 0.014 | 0.062 | 0.923 | 2.74 | 0.000 | 0.94 | 1.83 | 3.0–10.8 | 8 | 2 |

  - 78–94 % of every run lies beyond the trust window, so the bias moves the system mostly where the iteration-3
    ensemble does not hold. The 180 training frames come from the 5–18 % in the window; 38 of them have a C bonded to a
    third O (`cn` > 0.8), all from the two seed-2 free runs and one from free 400 K s1. In the surface runs and free
    300 K s1, every picked frame has `cn` ≤ 0.09: their reacted configurations were all beyond 0.30 eV/Å.
  - The 24 sealed-test frames span max_devi 0.09–2.45 eV/Å (median 1.07); 6 have `cn` > 0.8. Every selected frame,
    training and test, also passes the stage-7 geometry gate with its H limit of the time (1.2 Å).
  - The `cn` maxima of 1.8–2.7 and `ch` up to 6.2 come from frames beyond the trust window; they are not read as
    chemistry. A second round on an ensemble trained with these frames is what can move the window outward.
- **Labels (from 2026-11-05 05:39):** `cp2k_label.py prepare` wrote 180 + 24 inputs in `extra/metad/labels` and
  `extra/metad-test/labels`, identical to the iteration-4 HAP(001) + water + CO2 inputs apart from the coordinates
  (6 outer SCF loops). Container `seed-metad` (20 ranks from `extra/RANKS`, test set first) was started and paused at
  once. `metad/seed_scheduler.sh` resumes it only while LAMMPS exploration runs and no training does, so it does not
  slow label-4 or the trainings. When label-4 ends, the orchestrator resumes it and defers the set (iteration 4 < the
  last); the scheduler pauses it again when train-4 starts. Expected cost: 204 frames × about 4–5 min at 20 ranks.

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
