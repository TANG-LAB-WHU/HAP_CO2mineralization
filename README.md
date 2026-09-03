# AI-Driven Materials Discovery for Biomimetic CO2 Mineralization

![Project Status](https://img.shields.io/badge/Status-Active-brightgreen)
![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Machine Learning](https://img.shields.io/badge/AI4Science-MatterGen%20%7C%20MACE%20%7C%20DeepMD%20%7C%20CP2K%20%7C%20LAMMPS-orange)

## 📖 Introduction & Scientific Background

This repository organizes a developing, multi-scale computational study of CO2 interactions with hydroxyapatite (HAP) and related apatite-family surfaces. The intended workflow spans candidate generation, surface screening, first-principles calculations, machine-learning interatomic potentials, molecular dynamics, analysis, and manuscript preparation. Its stages have different implementation and evidence maturity, summarized below.

The project investigates mineral carbonation as a possible route to durable carbon storage. Hydroxyapatite ($Ca_{10}(PO_4)_6(OH)_2$) is used as a calcium-rich, compositionally tunable model material. HAP-specific surface reactivity and CO2-mineralization performance are research questions here, not established findings of this repository.

The planned study examines composition and facet effects, including possible F, Cl, Sr, and Mg substitutions. It is designed to combine generative models for proposing bulk candidates with machine-learning interatomic potentials for preliminary screening and first-principles calculations for validation. No composition or facet is currently qualified as superior.

---

## 🏗️ Architecture & Project Structure

The repository is organized as a nominal pipeline from Phase 0 (Generative Discovery) to Phase 5 (Publication). The directory order describes intended data flow; it does not imply that every stage is implemented or that existing artifacts passed scientific validation.

| Stage | Directory | Intended role | Current status |
| :--- | :--- | :--- | :--- |
| **Data** | `data/` | Initial PDB/XYZ structures for 22 HAP facets. | **Present, provenance pending** |
| **Phase 0** | `0.Generative_Discovery/` | Propose and filter apatite-like bulk candidates with MatterGen. | **Partial:** scripts exist; no tracked result summary or run manifest |
| **Phase 1** | `1.HT_Screening/` | Enumerate compositions, build slabs, place CO2, screen, and rank candidates. | **Placeholder:** several core scripts and READMEs are empty; no ranking exists |
| **Pre-Relax** | `2.GeoOpt_PreRelax/` | CP2K geometry optimization and planned MACE pre-relaxation. | **Legacy artifacts, unqualified:** completion states are mixed and some inputs contain conflict markers |
| **Phase 2-4** | `3.ML_AIMD/` | AIMD data generation, MLIP training, LAMMPS scale-up, and enhanced sampling. | **Partial/planned:** three short legacy AIMD outputs exist; no validated dataset, model, scale-up, or enhanced-sampling result |
| **Phase 4** | `4.PostProcess_Analysis/` | Analyze qualified LAMMPS outputs and compare materials. | **Legacy analysis, rejected for manuscript use:** source data and adequate quality evidence are missing |
| **Phase 5** | `5.Figure_Publication/` | Produce provenance-complete publication figures. | **Placeholder:** no qualifying publication figure exists |
| **Writing** | `paper/`, `literature/`, `research/`, `workflow/` | Maintain the manuscript, bibliography, evidence ledger, state, and provenance contracts. | **Implemented baseline; validation is being strengthened** |
| **Shared** | `utils/` | Shared structure and format utilities. | **Partial:** only selected utilities contain implementations |

---

## 🚀 Detailed Workflow Guide

### Phase 0: Generative Discovery (`0.Generative_Discovery/`)
The planned generation stage uses **MatterGen** to propose structures in selected chemical systems. Current scripts do not constitute a completed candidate campaign.
- **Batch scanning interface**: `batch_scan.sh` enumerates configured chemical systems and calls the generation, evaluation, statistics, and filtering scripts.
- **Filtering interface**: `post_filter.py` reads precomputed detailed metrics and filters by an energy-above-hull threshold and novelty flag. It does not independently establish thermodynamic validity.

### Phase 1: High-Throughput Screening (`1.HT_Screening/`)
This directory reserves a pipeline for ranking generated materials. The end-to-end workflow is not implemented in the current repository state.
1. **`01_composition_enumeration/`**: Placeholder for apatite-family composition rules; the implementation and known-composition dataset are empty.
2. **`02_slab_generation/`**: Contains command-line placeholders for slab generation, CO2 placement, and LAMMPS conversion; no structures are currently produced by these scripts.
3. **`03_mace_screening/`**: Contains a MACE-to-MLIAP conversion utility and a LAMMPS minimization template; run orchestration and adsorption-energy calculation are empty.
4. **`04_ranking_and_selection/`**: Ranking and Top-N export files are empty.

### Phase 2-4: High-Fidelity Validation (`3.ML_AIMD/`)
This directory contains inputs and partial tooling for a future validation pipeline. The command below is a Windows-oriented orchestration interface and must not be treated as evidence that every stage has run successfully.
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
- **Docker + Docker Compose**: Used only by selected compute workflows. Existing containers request NVIDIA GPUs and are not part of the macOS writing baseline.
- **NVIDIA GPU runtime**: Required for reasonable inference and training times for MatterGen, MACE, and DeepMD.
- **Python 3.10+**: Core dependencies include:
  - `pymatgen`
  - `ase`
  - `mattergen`
  - `mace-torch`
  - `lammps` (Python wrapper, optional if using Docker)

*Container images for CP2K, DeepMD-kit, and MACE are configured in respective `docker-compose.yml` files located within the phase directories.*

---

## 📄 License

This project is licensed under the MIT License.
See [LICENSE](LICENSE) for details.
