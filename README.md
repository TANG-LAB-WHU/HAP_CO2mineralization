# HAP_CO2mineralization

This repository combines **AIMD** (ab initio molecular dynamics) and **machine-learning interatomic potentials** to study accelerated mineralization (AAM-CO2) on hydroxyapatite (HAP) surfaces across multiple Miller indices.

The workflow links:
- CP2K short AIMD sampling
- DeepMD dataset synthesis and model training/freezing
- LAMMPS large-scale simulation with the trained model

---

## Project structure

| Directory | Description |
|-----------|-------------|
| **0.InitialStructureConfig** | Initial HAP structures (PDB/XYZ) for multiple facets (e.g. 002, 004, 100-513). |
| **1.GeoOpt** | CP2K geometry optimization for each facet (`HAP_xxx_Perfect`) and optimized structures. |
| **2.ML_AIMD** | Main four-stage pipeline (CP2K AIMD -> DeepMD data -> DeepMD train/freeze -> LAMMPS scale-up). |
| **3.Data_PostProcess** | Post-processing of LAMMPS outputs: thermodynamics, structure/dynamics metrics, plotting. |
| **4.Figure_Making** | Scripts/assets for publication-quality figures. |

Inside `2.ML_AIMD`:
- `Step1_aimd_cp2k_runs`: per-facet CP2K MD inputs/outputs
- `Step2_dataset_synthesis`: CP2K -> DeepMD conversion scripts and pooled datasets
- `Step3_mlip_deepmd`: DeepMD training/freeze inputs and model outputs
- `Step4_lammps_scaleup`: slab generation + LAMMPS run inputs
- `run_workflow.bat`: unified workflow entry point
- `docker-compose.yml`: container orchestration for CP2K/DeepMD/LAMMPS

---

## Unified workflow in `2.ML_AIMD`

Run from `2.ML_AIMD`:

```bat
run_workflow.bat [HKL_PARAM_TRAIN] [HKL_PARAM_LAMMPS]
```

### Parameter behavior

- `HKL_PARAM_TRAIN` (default `112`): facet used for Stage 1-3 (CP2K + DeepMD training data/model)
- `HKL_PARAM_LAMMPS` (default = `HKL_PARAM_TRAIN`): facet used for Stage 4 slab generation and LAMMPS run

Examples:

```bat
run_workflow.bat
run_workflow.bat 213
run_workflow.bat 112 431
```

### Stage-by-stage summary

1. **Stage 1 - AIMD (CP2K)**
   - Working folder: `Step1_aimd_cp2k_runs/HAP_<HKL_PARAM_TRAIN>_Perfect`
   - If `<prefix>_opt-pos-1.xyz` exists but `<prefix>_final.xyz` is missing, `extract_last_frame.py` builds final coordinates.
   - Runs service `cp2k_run` only if `<prefix>_md-pos.xyz` is absent.

2. **Stage 2 - Dataset synthesis**
   - `generate_deepmd_data.py` writes per-system data to `Step2_dataset_synthesis/deepmd_pool/HAP_<HKL_PARAM_TRAIN>`.
   - `prepare_multisystem_input.py` generates `Step3_mlip_deepmd/input.multisystem.json`.

3. **Stage 3 - DeepMD train/freeze**
   - Runs `deepmd_train` then `deepmd_freeze` from `docker-compose.yml`.
   - Expected model artifact: `Step3_mlip_deepmd/model/hap_model.pth`.

4. **Stage 4 - LAMMPS scale-up**
   - `create_hap_slabs.py` builds slab/data for `HKL_PARAM_LAMMPS`.
   - Runs `lammps_run` with `Step4_lammps_scaleup/in.deepmd.lammps`.

---

## Runtime dependencies

- **Docker + Docker Compose** (required for CP2K / DeepMD / LAMMPS services)
- **NVIDIA GPU runtime** (the compose file requests `gpus: all`)
- **Python 3.10+** available in active conda env or system `PATH`
- Python packages for helper scripts (see subfolder requirements if provided)

Container images are configured in `2.ML_AIMD/docker-compose.yml` (CP2K and DeepMD-kit based images).

---

## Inputs and outputs (quick reference)

- **Step1 input**: `Step1_aimd_cp2k_runs/HAP_xxx_Perfect/HAP_xxx_md.inp`
- **Step1 key outputs**: `*_md-pos.xyz`, `*_md-1.ener`, `*_forces.dat`, `*_virial.dat`
- **Step2 pooled data**: `Step2_dataset_synthesis/deepmd_pool/HAP_xxx`
- **Step3 model**: `Step3_mlip_deepmd/model/hap_model.pth`
- **Step4 run files**: `Step4_lammps_scaleup/hap_<HKL>.data`, LAMMPS logs/trajectory outputs

---

## Notes and troubleshooting

- `run_workflow.bat` performs strict pre-checks and exits on missing scripts/templates.
- If CP2K trajectory already exists (`*_md-pos.xyz`), Stage 1 is skipped intentionally.
- In `docker-compose.yml`, shell variables inside `deepmd_test` use `$$` escaping to avoid Docker Compose interpolation warnings.
- If containers start but jobs fail, inspect workflow output plus `docker logs <container_id>` for the failed stage.

---

## License

This project is licensed under the MIT License.
See [LICENSE](LICENSE) for details.
