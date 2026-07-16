# AI-Driven Materials Discovery for Biomimetic CO2 Mineralization

![Project Status](https://img.shields.io/badge/Status-Active-brightgreen)
![Python](https://img.shields.io/badge/Python-3.10%2B-blue)
![Machine Learning](https://img.shields.io/badge/AI4Science-MatterGen%20%7C%20MACE%20%7C%20DeepMD%20%7C%20CP2K%20%7C%20LAMMPS-orange)

## 📖 Introduction & Scientific Background

This repository provides an integrated, multi-scale computational workflow to discover and study **accelerated CO2 mineralization (AAM-CO2)** on the surfaces of hydroxyapatite (HAP) and related Apatite-family materials. 

**Biomimetic CO2 Mineralization** is a promising negative-emission technology that aims to permanently sequester CO2 by converting it into stable solid carbonates, mimicking natural biological and geological processes. Hydroxyapatite ($Ca_{10}(PO_4)_6(OH)_2$), the primary inorganic component of human bones and teeth, presents a highly active surface for such catalytic conversions. 

However, exploring the vast chemical space of potential dopants (e.g., F, Cl, Sr, Mg) to optimize the CO2 capture efficiency has traditionally been prohibitively slow using Density Functional Theory (DFT) alone. To overcome this bottleneck, this project bridges the gap between **Generative AI** for theoretical bulk discovery and **Machine Learning Interatomic Potentials (MLIPs)** for high-throughput surface catalytic/adsorption validation.

---

## 🏗️ Architecture & Project Structure

The workflow is organized into a strict, linear pipeline mapping from Phase 0 (Generative Discovery) to Phase 5 (Publication).

| Stage | Directory | Description |
| :--- | :--- | :--- |
| **Data** | `data/` | Initial structural data (PDB/XYZ/CIF) representing perfect HAP facets (e.g., 002, 004, 100-513). Serves as the static ground truth. |
| **Phase 0** | `0.Generative_Discovery/` | **Generative AI phase space search.** Uses the diffusion model *MatterGen* to hallucinate novel, stable Apatite-like compositions (e.g., Ca-P-O-H-C, Ca-Sr-P-O-H). |
| **Phase 1** | `1.HT_Screening/` | **High-throughput screening pipeline.** Cleaves bulk candidates into slabs via `pymatgen`, places CO2 adsorbates, and evaluates adsorption energies using the universal foundation model *MACE-MH-1* in LAMMPS. |
| **Pre-Relax** | `2.GeoOpt_PreRelax/` | Geometry optimization for generated candidates and perfect facets. Supports classical CP2K relaxations and ultra-fast *MACE-MH-1* pre-relaxations. |
| **Phase 2-4** | `3.ML_AIMD/` | **Main high-fidelity pipeline.** For the Top-N candidates: runs *CP2K AIMD* ➔ Trains custom *DeepMD/MACE* potentials ➔ Executes *LAMMPS* scale-up ➔ Performs Enhanced Sampling (Metadynamics). |
| **Phase 4** | `4.PostProcess_Analysis/` | Deep post-processing of LAMMPS outputs: extracts thermodynamics, computes radial distribution functions (RDFs), evaluates structural dynamics, and performs cross-material comparisons. |
| **Phase 5** | `5.Figure_Publication/` | Scripts and assets for rendering publication-quality scientific figures. |
| **Shared** | `utils/` | Shared utility scripts (e.g., structure IO, LAMMPS data converters, and element mappings). |

---

## 🚀 Detailed Workflow Guide

### Phase 0: Generative Discovery (`0.Generative_Discovery/`)
Instead of randomly guessing dopant ratios, we use the generative diffusion model **MatterGen** to intelligently explore the compositional phase space of Apatites.
- **Batch Scanning**: Use `./batch_scan.sh` to run bulk generation across multiple chemical systems simultaneously (Ca-P-O-H, Ca-Sr-P-O-H, etc.).
- **Thermodynamic Filtering**: The script `post_filter.py` automatically computes the Energy above Hull ($E_{hull}$) against the `TRI2024` dataset to extract only the most thermodynamically stable and novel candidates for the next stage.

### Phase 1: High-Throughput Screening (`1.HT_Screening/`)
A fast, automated pipeline designed to rank the generated materials based on their CO2 adsorption capabilities.
1. **`01_composition_enumeration/`**: Enumerates and applies Apatite super-family constraints to the candidates.
2. **`02_slab_generation/`**: Programmatically cleaves optimal Miller index surfaces and systematically places CO2 adsorbates via `pymatgen`.
3. **`03_mace_screening/`**: Utilizes the **MACE-MH-1** universal foundation model to rapidly perform molecular mechanics minimization of the slab+CO2 systems in LAMMPS, bypassing expensive DFT optimizations.
4. **`04_ranking_and_selection/`**: Extracts final energies, calculates adsorption energies ($E_{ads}$), and exports the Top-N structures directly into the AIMD pipeline.

### Phase 2-4: High-Fidelity Validation (`3.ML_AIMD/`)
The traditional rigorous validation pipeline exclusively reserved for the Top-N screened candidates.
```bat
cd 3.ML_AIMD
run_workflow.bat [HKL_PARAM_TRAIN] [HKL_PARAM_LAMMPS]
```
1. **Step 1 - AIMD (CP2K)**: Short, highly accurate *ab-initio* molecular dynamics sampling to generate real reactive trajectories.
2. **Step 2 - Dataset synthesis**: Converts CP2K trajectories into ML-ready training datasets.
3. **Step 3 & 5 - MLIP Training**: Train a custom potential from scratch using DeepMD (`Step3_mlip_deepmd/`) or fine-tune the existing MACE foundation model (`Step5_mace_finetune/`) on the specific Apatite-CO2 interactions.
4. **Step 4 & 6 - Scale-up Simulation**: Perform large-scale, nanosecond-long LAMMPS simulations (`Step4_lammps_scaleup/`) and use Plumed for metadynamics reaction pathway analysis (`Step6_enhanced_sampling/`).

---

## ⚡ Accelerated Pre-relaxation (`2.GeoOpt_PreRelax/`)

For standard HAP facets or bulk candidates requiring quick geometric convergence before deep DFT simulations, we integrate a specialized MACE-MH-1 pre-relaxation script.
- Navigate to `2.GeoOpt_PreRelax/mace_prerelax/`
- Run `run_mace_prerelax.bat` to instantly optimize all 22 HAP perfect facets using a GPU-accelerated LAMMPS Docker container. This saves thousands of CPU core-hours compared to CP2K geometry optimization.

---

## 💻 Runtime Dependencies & Setup

- **Docker + Docker Compose**: Highly recommended. Container orchestrations request `gpus: all` for seamless GPU acceleration without dealing with complex CUDA libraries.
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
