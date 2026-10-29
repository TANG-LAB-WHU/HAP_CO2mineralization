# Pipeline v2: GPU-first MatterGen → CP2K → DeePMD → LAMMPS

Status: **in progress (started 2026-10-29)**. Nothing produced by this pipeline is
manuscript-eligible until its gate report and run manifest are reviewed.

Run directory: `.research-runs/pipeline-v2/` (not in Git). Orchestrator: `python workflow/pipeline_v2.py
--run .research-runs/pipeline-v2 --cutoff 800 --rel-cutoff 60 --iters 4`; it resumes from
`pipeline_state.json`. Log: `logs/pipeline_v2.log`.

| Stage | State on 2026-10-29 |
| --- | --- |
| Images | `hap-cp2k:2026.2-cuda-sm89` and `hap-mlip:2026.10` built; CP2K labels use `cp2k/cp2k:2026.2_mpich_generic_psmp` |
| CP2K protocol | done: 800/60 Ry; ≥ 8 Å Γ cells with the recorded bias |
| MatterGen S1/S2 | done: 363 → 149 (see the run README) |
| MatterGen S3 | running on the GPU (2 shards) |
| Exploration | done: 8 systems × 3 temperatures (MACE-MPA-0, float32); independent seed-777 run for the sealed test |
| Iteration-0 labels | running: 216 frames, 24 MPI ranks |
| Sealed test labels | queued: 84 frames, labelled while iteration-0 trains |
| Finite-size check | 5/20 paired labels done, paused to free the CPU; resume with `run_labels.sh` in `finite-size/labels` |

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

Active learning follows DP-GEN: a 4-model ensemble, with frames whose maximum force deviation
falls between trust_lo = 0.10 and trust_hi = 0.30 eV/Å sent to CP2K. The DP-GEN water example
uses 0.05/0.15 eV/Å; the looser window fits the larger target error for mineral–water interfaces.

## Stages

1. **MatterGen screening (GPU):** S1 chemistry gate → S2 three-model hull committee
   (E_hull ≤ 50 meV/atom in ≥ 2 models, references from Alex-MP relaxed with the same model)
   → S3 phonons (min frequency ≥ −0.3 THz) and 10 ps 300 K NPT. Scripts: `0.Generative_Discovery/screening/`.
2. **CP2K protocol tests (CPU, then GPU timing):** `cp2k_protocol.py`.
3. **Exploration (GPU):** foundation-MLIP MD of HAP bulk, HAP(001) + CO2, HAP(001) + water + CO2,
   calcite, and S3 survivors at 300–600 K; farthest-point selection for the first labels.
4. **Labelling (CPU):** CP2K single points with the fixed protocol.
5. **DeePMD (GPU):** se_e2_a, rcut 6 Å, as in the CaCO3 reference; 4 seeds; compressed for MD.
6. **Active learning (GPU MD + CPU labels):** LAMMPS-DeePMD exploration with model deviation until converged.
7. **Validation (GPU):** the gates above; a failure is recorded, not hidden.
