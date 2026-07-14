# HAP_CO2mineralization

This repository provides an integrated, multi-scale computational workflow to discover and study accelerated CO2 mineralization (AAM-CO2) on the surfaces of hydroxyapatite (HAP) and related Apatite-family materials.

The workflow spans from **generative AI** for material discovery, to **high-throughput screening** with universal machine learning potentials, and finally to **high-fidelity AIMD** and **custom MLIPs** (Machine-Learning Interatomic Potentials) for precise scale-up simulations.

---

## Project Structure

The repository is organized into distinct phases (0 to 5) that represent the full material discovery and validation lifecycle:

| Directory | Description |
|-----------|-------------|
| **mattergen/** | **[Phase 0]** Generative AI phase space search using MatterGen. Generates and evaluates novel Apatite-like compositions (e.g., Ca-P-O-H-C, Ca-Sr-P-O-H). |
| **0.InitialStructureConfig/** | Initial HAP structures (PDB/XYZ) for multiple facets (e.g. 002, 004, 100-513). |
| **1.GeoOpt/** | Geometry optimization for each facet. Supports classical CP2K relaxations and ultra-fast **MACE-MH-1 pre-relaxations** (`mace_prerelax/`). |
| **5.HT_Screening/** | **[Phase 1]** High-throughput screening pipeline. Evaluates MatterGen candidates for CO2 adsorption capability using MACE-MH-1 and LAMMPS. |
| **2.ML_AIMD/** | **[Phase 2-4]** Main high-fidelity pipeline (CP2K AIMD ➔ DeepMD/MACE Fine-tuning ➔ LAMMPS scale-up ➔ Enhanced Sampling). |
| **3.Data_PostProcess/** | Post-processing of LAMMPS outputs: thermodynamics, structure/dynamics metrics, and cross-material comparisons. |
| **4.Figure_Making/** | Scripts/assets for publication-quality figures. |
| **utils/** | Shared utility scripts (e.g., `xyz_to_lammps.py`, element mappings). |

---

## Workflow Guide

### Phase 0: Generative Discovery (`mattergen/`)
Use the generative diffusion model (MatterGen) to explore the compositional phase space of Apatites.
- **Batch Scanning**: Use `./batch_scan.sh` to run bulk generation across multiple chemical systems.
- **Filtering**: `post_filter.py` automatically extracts stable (low energy above hull) and novel candidates for the next stage.

### Phase 1: High-Throughput Screening (`5.HT_Screening/`)
A fast, automated pipeline to rank the generated materials based on CO2 adsorption energy.
1. **`01_composition_enumeration/`**: Enumerates Apatite super-family constraints.
2. **`02_slab_generation/`**: Cleaves surfaces and places CO2 adsorbates via `pymatgen`.
3. **`03_mace_screening/`**: Uses the MACE-MH-1 universal foundation model to rapidly minimize the slab+CO2 systems in LAMMPS.
4. **`04_ranking_and_selection/`**: Ranks candidates by adsorption energy and exports the Top-N structures to the AIMD pipeline.

### Phase 2-4: High-Fidelity Validation (`2.ML_AIMD/`)
The traditional rigorous pipeline for the Top-N candidates.
Run from `2.ML_AIMD`:
```bat
run_workflow.bat [HKL_PARAM_TRAIN] [HKL_PARAM_LAMMPS]
```
1. **Step 1 - AIMD (CP2K)**: Short ab-initio molecular dynamics sampling.
2. **Step 2 - Dataset synthesis**: Converts CP2K trajectories into ML training sets.
3. **Step 3 & 5 - MLIP Training**: Train a custom DeepMD model (`Step3_mlip_deepmd/`) or fine-tune the MACE foundation model (`Step5_mace_finetune/`).
4. **Step 4 & 6 - Scale-up**: Large-scale LAMMPS simulations (`Step4_lammps_scaleup/`) and metadynamics for reaction pathways (`Step6_enhanced_sampling/`).

---

## Accelerated Pre-relaxation (`1.GeoOpt/mace_prerelax/`)
To bypass expensive DFT geometry optimizations for standard surfaces, we integrate MACE-MH-1 pre-relaxation:
- Navigate to `1.GeoOpt/mace_prerelax/`
- Run `run_mace_prerelax.bat` to instantly optimize all 22 HAP perfect facets using a GPU-accelerated LAMMPS container.

---

## Runtime Dependencies

- **Docker + Docker Compose**: Required for CP2K, DeepMD, and LAMMPS (MACE/MLIAP) services.
- **NVIDIA GPU runtime**: Container orchestrations request `gpus: all`.
- **Python 3.10+**: With `pymatgen`, `ase`, and `mattergen` dependencies installed.

*Container images are configured in respective `docker-compose.yml` files throughout the repository.*

---

## License

This project is licensed under the MIT License.
See [LICENSE](LICENSE) for details.
