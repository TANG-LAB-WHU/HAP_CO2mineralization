# AI-Driven Materials Discovery for Biomimetic CO2 Mineralization

**Current work (4 November 2026):** the GPU-first [pipeline v2](workflow/PIPELINE_V2.md) is running active
learning on Shawn-PC. Iterations 0–2 are complete, and every ensemble passes the sealed-test gate. The
lattice, stability and convergence gates have not run, the final ensemble does not exist yet, and unbiased
dynamics has sampled no CO2 conversion. Measured state, with hash-checked records:
[pipeline-v2 snapshot](results/pipeline-v2-20261104/README.md).

Earlier development results: the [17 October 2026 release](results/RELEASE-20261017.md) (conditional DFT
references and failed dynamics/model-validation diagnostics) and the [MatterGen generation record](results/mattergen-generation-20261011/README.md)
(an 8-atom smoke and a 44-atom out-of-distribution pilot; all five samples rejected). No production-ready
potential or mineralization mechanism is established, and nothing is manuscript-eligible.

![Project Status](https://img.shields.io/badge/Status-Active-brightgreen)
![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Machine Learning](https://img.shields.io/badge/AI4Science-MatterGen%20%7C%20MACE%20%7C%20DeepMD%20%7C%20CP2K%20%7C%20LAMMPS-orange)

## 📖 Introduction & Scientific Background

This repository organizes a developing, multi-scale computational study of CO2 interactions with hydroxyapatite (HAP) and related apatite-family surfaces. The intended workflow spans candidate generation, surface screening, first-principles calculations, machine-learning interatomic potentials, molecular dynamics, analysis, and manuscript preparation. Its stages have different implementation and evidence maturity, summarized below.

The project investigates mineral carbonation as a possible route to durable carbon storage. Hydroxyapatite ($Ca_{10}(PO_4)_6(OH)_2$) is used as a calcium-rich, compositionally tunable model material. HAP-specific surface reactivity and CO2-mineralization performance are research questions here, not established findings of this repository.

The planned study examines composition and facet effects, including possible F, Cl, Sr, and Mg substitutions. It is designed to combine generative models for proposing bulk candidates with machine-learning interatomic potentials for preliminary screening and first-principles calculations for validation. No composition or facet is currently qualified as superior.

---

## 🔁 From the original workflow to pipeline v2

The original plan was CP2K geometry optimization of each HAP facet, then short CP2K AIMD, a DeePMD potential
trained on those trajectories, LAMMPS scale-up with that potential, and post-processing. Its scripts are still in the
numbered directories (`3.ML_AIMD/workflow_documentation.md`, `3.ML_AIMD/run_workflow.bat`). The local v1 run of that
route failed its production gates: retrospective test force RMSE 0.108 eV/Å on 8 frames against a 0.100 gate, and both
1 ps trajectories failed structural validation ([v1 recovery record](workflow/RECOVERY.md)). The 32 frames in all came from
short four-femtosecond AIMD replicas, too few and too correlated to cover the configurations that MD later visits.
Pipeline v2 keeps the same stages but replaces the training data strategy with DFT-labelled active learning:

| Original step | Pipeline v2 implementation | State (4 Nov 2026) |
| --- | --- | --- |
| Initial structures (`data/`, 22 facets) | HAP bulk with ordered OH (`3.ML_AIMD/Step4_lammps_scaleup/HAP_ConventionalCell.cif`), HAP(001) dry + CO2 and hydrated + 2 CO2 slabs, calcite, CaCO3, portlandite, monetite and CaO cells ≥ 8 Å (`workflow/explore_systems.py`), plus 10 MatterGen survivors | done |
| Candidate generation (Phase 0) | MatterGen Ca-P-O-H-C scan → S1 chemistry → S2 three-MLIP hull committee → S3 phonons + NPT → S4 CP2K hull (`0.Generative_Discovery/screening/`) | done: 2,324 → 363 → 149 → 10 → 8 within 50 meV/atom, none on the hull |
| CP2K geometry optimization | one fixed labelling protocol, tested before use (`workflow/cp2k_protocol.py`): PBE-D3(BJ), DZVP-MOLOPT-SR-GTH, 800/60 Ry, Γ point, widths ≥ 8 Å | done; 8 Å bias measured at 2.46 ± 0.24 meV/atom |
| Short CP2K AIMD | MACE-MPA-0 exploration MD, then CP2K single points on selected frames (`workflow/cp2k_label.py`); sealed test from independent seed-777 trajectories | iterations 0–2 labelled (620 frames), iteration 3 running |
| DeepMD training | 4-model se_e2_a ensemble, rcut 6 Å, warm-started between iterations (`3.ML_AIMD/Step3_mlip_deepmd/dp_ensemble.py`) | iterations 0–2 trained; final from-scratch ensemble pending |
| LAMMPS with DeepMD | model-deviation exploration, trust 0.10–0.30 eV/Å (`dp_explore.py`), then the stage-7 lattice, stability and NVE-drift gates (`dp_validate.py`) | exploration iterations 0–2 done; stage 7 not run |
| Enhanced sampling | biased exploration of CO2 → HCO3⁻/CO3²⁻ | PLANNED, not implemented |
| Post-processing and figures | `4.PostProcess_Analysis/`, `5.Figure_Publication/` | not applied to v2 |

The orchestrator `workflow/pipeline_v2.py` runs labelling (CPU), training and exploration (GPU) in a resumable loop.
Method, every deviation, and the benchmarks behind the hardware choices are in [`workflow/PIPELINE_V2.md`](workflow/PIPELINE_V2.md).

## 📊 Current results (all development-only)

| Result | Value | Record |
| --- | --- | --- |
| MatterGen screen | 2,324 unique S.U.N. → S1 363 → S2 149 → S3 10 → S4 8 within 50 meV/atom of the CP2K hull; none on the hull; the only carbonate-phosphate fails (68.0 meV/atom) | [run README](0.Generative_Discovery/runs/Ca-P-O-H-C_gfscan_RTX5090_20260329/README.md) |
| MLIP vs DFT hull | three MLIPs agree within 2.5–3.3 meV/atom (median pairwise) but differ from CP2K by −12.8 to +20.7 meV/atom | same |
| Γ finite-size bias | 8 Å minus 12 Å: 2.46 meV/atom (SD 0.24, 10 HAP frames); forces 0.016 eV/Å RMS | [snapshot](results/pipeline-v2-20261104/README.md) |
| Sealed test, iteration 2 | force RMSE 0.066–0.067 eV/Å overall, 0.096–0.098 eV/Å on interfaces; energy 1.33–1.67 meV/atom; gate passed by 4/4 models | same |
| Weakest system | HAP(001) + water + CO2: force 0.104–0.107 eV/Å, energy 2.8–3.9 meV/atom (fails the 2 meV/atom gate on its own) | same |
| Active learning | mean accurate fraction 0.673 → 0.723 → 0.838; hydrated interface 0.008 (300 K) and 0.106 (500 K) in iteration 2; target ≥ 0.99 | same |
| Reactive events | no CO2 conversion observed in unbiased MD | same |

Open gates: active-learning convergence, bulk HAP lattice (NPT ≥ 200 ps), 500 ps stability and NVE drift, all on the final
ensemble. A literature novelty check of the planned contributions is in [`literature/novelty_check_2026-11.md`](literature/novelty_check_2026-11.md);
it found partial precedents for an HAP machine-learning potential and for generative-screen attrition.

## 🔍 Reproducibility

| Check | Command | Needs |
| --- | --- | --- |
| Pipeline-v2 snapshot: file hashes and recomputed tables | `python results/pipeline-v2-20261104/scripts/summarize.py --check` | Python standard library |
| 17 October release inventory | `python -B results/verify_release_20261017.py` | Python standard library |
| Registry, evidence and inventory rules | `make check` | Python |
| Workflow unit tests | `python -m pytest tests` | pytest, NumPy, ASE |
| Regenerate a snapshot from the run directory | see [`results/pipeline-v2-20261104/REPRODUCE.md`](results/pipeline-v2-20261104/REPRODUCE.md) | the run directory on Shawn-PC |
| Rerun pipeline v2 | `python workflow/pipeline_v2.py --run .research-runs/pipeline-v2 --cutoff 800 --rel-cutoff 60 --iters 6 --gpu-jobs 2 --explore-jobs 3 --warm-steps 100000 --extra-systems s3` | Docker images and inputs listed in the snapshot's `REPRODUCE.md`; NVIDIA GPU |

Run directories, labels, models and trajectories stay outside Git (`.research-runs/`); result packages record their SHA-256.
GPU training and MD are not bitwise reproducible, so a rerun reproduces the statistics, not identical frames.

---

## 🏗️ Architecture & Project Structure

The repository is organized as a nominal pipeline from Phase 0 (Generative Discovery) to Phase 5 (Publication). The directory order describes intended data flow; it does not imply that every stage is implemented or that existing artifacts passed scientific validation.

| Stage | Directory | Intended role | Current status |
| :--- | :--- | :--- | :--- |
| **Data** | `data/` | Initial PDB/XYZ structures for 22 HAP facets. | **Present, provenance pending** |
| **Phase 0** | `0.Generative_Discovery/` | Propose and filter apatite-like bulk candidates with MatterGen. | **Screened, development-only:** imported Ca-P-O-H-C scan through S1–S4 (10 survivors labelled; they join training at iteration 3); the earlier bounded smoke and 44-atom pilot rejected all five samples |
| **Phase 1** | `1.HT_Screening/` | Enumerate compositions, build slabs, place CO2, screen, and rank candidates. | **Placeholder:** several core scripts and READMEs are empty; no ranking exists |
| **Pre-Relax** | `2.GeoOpt_PreRelax/` | CP2K geometry optimization and planned MACE pre-relaxation. | **Legacy artifacts, unqualified:** completion states are mixed and some inputs contain conflict markers |
| **Phase 2-4** | `3.ML_AIMD/`, `workflow/` | AIMD data generation, MLIP training, LAMMPS scale-up, and enhanced sampling. | **Pipeline v2 in progress:** DeePMD ensemble, model-deviation exploration and validation code in `3.ML_AIMD/Step3_mlip_deepmd/`; sealed test passed in iterations 0–2; no validated model, scale-up or enhanced-sampling result. Legacy AIMD outputs remain ineligible |
| **Phase 4** | `4.PostProcess_Analysis/` | Analyze qualified LAMMPS outputs and compare materials. | **Legacy analysis, rejected for manuscript use:** source data and adequate quality evidence are missing |
| **Phase 5** | `5.Figure_Publication/` | Produce provenance-complete publication figures. | **Placeholder:** no qualifying publication figure exists |
| **Results** | `results/` | Compact, hash-checked result packages. | **Development records:** pipeline-v2 snapshot, 17 October release, MatterGen records |
| **Writing** | `paper/`, `literature/`, `research/`, `workflow/` | Maintain the manuscript, reference registry, bibliography, evidence relationships, state, and provenance contracts. | **Phase 2A validated:** eight seeds human-verified and cleanly BBT-exported; 20 further references are unverified candidates; all claim-specific full-text locators remain pending |
| **Shared** | `utils/` | Shared structure and format utilities. | **Partial:** only selected utilities contain implementations |

---

## 🚀 Detailed Workflow Guide

### Phase 0: Generative Discovery (`0.Generative_Discovery/`)
The generation stage uses **MatterGen** to propose structures in selected chemical systems.
- **Imported scan**: a Ca-P-O-H-C guidance-factor scan (gf 0.001–2, 16,000 structures, MatterSim-5M relaxation). Its tracked summary, manifest and shortlist are in `runs/Ca-P-O-H-C_gfscan_RTX5090_20260329/`; raw outputs stay untracked in `results/`. The MatterGen commit and seeds were not recorded.
- **Screening** (`screening/`): S1 chemistry gate, S2 hull committee of MACE-MPA-0, MatterSim-5M and SevenNet-MF-OMPA against 288 Alex-MP references, S3 MACE phonons and 5 ps NPT, S4 CP2K PBE-D3(BJ) hull with 34 computed references (`cp2k_hull.py`).
- **Bounded generation records**: an 8-atom smoke and four 44-atom samples, all rejected before HAP qualification. The 44-atom pilot extends the atom-count distribution locally beyond the recorded native domain; it does not establish exact-composition generation, fine-tuning, or a validated candidate campaign.
- **Scan summary**: `summarize_gf_scan.py` classifies S.U.N. structures by CO3/PO4 motifs and charge balance, deduplicates across runs, and writes checksums.
- **Original job scripts**: `legacy_siqi_jobscripts/` (unrelated Ca-V-O-Ru and La-C-O scripts are in `unrelated_systems/`).
- **Batch scanning interface**: `batch_scan.sh` enumerates configured chemical systems and calls the generation, evaluation, statistics, and filtering scripts.
- **Filtering interface**: `post_filter.py` reads precomputed detailed metrics and filters by an energy-above-hull threshold and novelty flag. It does not independently establish thermodynamic validity.

### Phase 1: High-Throughput Screening (`1.HT_Screening/`)
This directory reserves a pipeline for ranking generated materials. The end-to-end workflow is not implemented in the current repository state.
1. **`01_composition_enumeration/`**: Placeholder for apatite-family composition rules; the implementation and known-composition dataset are empty.
2. **`02_slab_generation/`**: Contains command-line placeholders for slab generation, CO2 placement, and LAMMPS conversion; no structures are currently produced by these scripts.
3. **`03_mace_screening/`**: Contains a MACE-to-MLIAP conversion utility and a LAMMPS minimization template; run orchestration and adsorption-energy calculation are empty.
4. **`04_ranking_and_selection/`**: Ranking and Top-N export files are empty.

### Phase 2-4: High-Fidelity Validation (`3.ML_AIMD/` and `workflow/`)
**Pipeline v2 (current).** `workflow/pipeline_v2.py` drives the loop below. It resumes from `pipeline_state.json` in its run directory.
1. **Protocol** (`workflow/cp2k_protocol.py`): cutoff and Γ-width tests on HAP bulk fix the labelling protocol before any label.
2. **Exploration and selection** (`workflow/explore_systems.py`, `workflow/select_frames.py`): MACE-MPA-0 Langevin MD at 300–700 K, even-in-time frame selection, strained bulk copies; an independent seed-777 run gives the sealed test.
3. **Labels** (`workflow/cp2k_label.py`): CP2K single points in `cp2k/cp2k:2026.2_mpich_generic_psmp` on the CPU; SCF failures are kept and excluded.
4. **Ensemble** (`3.ML_AIMD/Step3_mlip_deepmd/dp_ensemble.py`): four DeePMD-kit 3.1.3 PyTorch models with energy and force terms, plus the virial term from iteration 2 on; sealed test per system.
5. **Active learning** (`dp_explore.py`): LAMMPS + DeePMD (TensorFlow backend, checked against PyTorch) with model deviation; candidates go back to step 3.
6. **Validation** (`dp_validate.py`): lattice, geometry and NVE-drift gates on the final ensemble.

**Legacy route.** The command below is the original Windows-oriented orchestration interface. It must not be treated as evidence that every stage has run successfully.
```bat
cd 3.ML_AIMD
run_workflow.bat [HKL_PARAM_TRAIN] [HKL_PARAM_LAMMPS]
```
1. **Step 1 - AIMD (CP2K)**: Inputs exist for 22 facets. Only three facets have legacy MD outputs, and none is currently eligible as training evidence.
2. **Step 2 - Dataset synthesis**: A conversion script exists, but the tracked DeepMD pool is empty.
3. **Step 3 & 5 - MLIP training**: A DeepMD template and preparation script exist; no validated model exists. The MACE fine-tuning files are placeholders.
4. **Step 4 & 6 - Scale-up and enhanced sampling**: A LAMMPS input template exists, but no qualifying scale-up output is present. The enhanced-sampling implementation is empty.

---

## ⚡ Accelerated Pre-relaxation (`2.GeoOpt_PreRelax/`)

The directory reserves a MACE-MH-1 pre-relaxation route for comparison with CP2K geometry optimization. Its scientific validity has not yet been established for this system.
- Navigate to `2.GeoOpt_PreRelax/mace_prerelax/`
- `run_mace_prerelax.bat` and a GPU-oriented Docker configuration reserve a batch interface. Do not use generated structures as evidence until the model domain, convergence, and CP2K comparison are reviewed.

---

## 💻 Runtime Dependencies & Setup

- **Local writing baseline**: Git, Python 3.10+, Make, and Quarto. Run `make check` and `make paper` from the repository root.
- **Pipeline v2 compute** (tested on Shawn-PC: i9-14900K, 64 GB, RTX 4090 24 GB, Docker Desktop on WSL2):
  - `hap-mlip:2026.10`, built from `workflow/environment/Dockerfile.mlip` (DeePMD-kit 3.1.3 base with LAMMPS and PLUMED, plus MACE, MatterSim and SevenNet environments). The Dockerfile does not pin its pip packages; the versions that ran are recorded in `results/pipeline-v2-20261104/environment/`.
  - `cp2k/cp2k:2026.2_mpich_generic_psmp` from Docker Hub, pinned by the digest recorded in the snapshot.
  - `hap-cp2k:2026.2-cuda-sm89` (`workflow/environment/Dockerfile.cp2k-cuda-sm89`), used only for the GPU benchmark.
- **Python 3.10+**: Core dependencies include:
  - `pymatgen`
  - `ase`
  - `mattergen`
  - `mace-torch`
  - `lammps` (Python wrapper, optional if using Docker)

*Container images for the legacy routes are configured in the `docker-compose.yml` files within the phase directories.*

---

## 🧪 Earlier Shawn-PC local validation (v1)

See [local validation and commands](workflow/LOCAL_VALIDATION.md), [extended reproduction](workflow/REPRODUCE_EXTENDED.md), and [smoke versus proposed production parameters](workflow/PARAMETERS.md). All local results remain manuscript-ineligible; MACE/CP2K adsorption disagreement and inadequate MLIP coverage are explicitly retained. HPC profiles are disabled and no jobs were submitted.

The separate [bounded structure and model refinement](workflow/REFINEMENT.md) adds explicit OH ordering, force-grid checks, reference decomposition, disjoint small-sample learning tests, and real trajectory diagnostics. Its [working LaTeX draft](drafts/local-validation/main.tex) concerns local technical validation only; it does not replace the living manuscript or promote raw artifact eligibility.

## 📄 License

This project is licensed under the MIT License.
See [LICENSE](LICENSE) for details.
