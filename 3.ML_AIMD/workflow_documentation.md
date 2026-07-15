# 2.ML_AIMD Workflow

## Four-stage unified workflow

This workflow is organized as a fixed four-stage pipeline:

1. AIMD generation (CP2K)
2. Training dataset synthesis (multi-structure pool)
3. DeepMD training and model export
4. LAMMPS large-scale simulation

## Stage 1: AIMD (CP2K)

- `run_workflow.bat` accepts HKL as parameter (default: `112`)
- Stage 1 working directory:
  - `Step1_aimd_cp2k_runs/HAP_<hkl>_Perfect`
- CP2K outputs trajectory, energy, force, and virial files in this folder

## Stage 2: Training dataset synthesis

- `Step2_dataset_synthesis/generate_deepmd_data.py` converts one Step1 system into DeepMD format
- Output is written directly into pooled system path:
  - `Step2_dataset_synthesis/deepmd_pool/HAP_<hkl>/`

## Stage 3: DeepMD training and export

- `Step3_mlip_deepmd/prepare_multisystem_input.py` scans `Step2_dataset_synthesis/deepmd_pool/*` and generates:
  - `Step3_mlip_deepmd/input.multisystem.json`
- `deepmd_train` reads `Step3_mlip_deepmd/input.multisystem.json`
- `training.training_data.systems` uses pooled multi-structure systems
- sampling policy defaults to:
  - `auto_prob: "prob_sys_size"`
- `deepmd_freeze` exports:
  - `Step3_mlip_deepmd/model/hap_model.pth`

## Stage 4: LAMMPS large-scale simulation

- `lammps_run` loads `Step3_mlip_deepmd/model/hap_model.pth`
- executes `Step4_lammps_scaleup/in.deepmd.lammps`
- writes trajectories and analysis outputs in `Step4_lammps_scaleup/`

## Data pool convention

Each pooled system directory should contain:

- `set.000/`
- `type.raw`
- `type_map.raw`

Example:

- `Step2_dataset_synthesis/deepmd_pool/HAP_100/`
- `Step2_dataset_synthesis/deepmd_pool/HAP_112/`
- `Step2_dataset_synthesis/deepmd_pool/HAP_300/`

## DeepMD multi-system probability options

In `Step3_mlip_deepmd/input.multisystem.json`, use one of:

- `auto_prob: "prob_sys_size"` (default)
- `auto_prob: "prob_uniform"`
- `auto_prob: "prob_sys_size;0:4:0.5;4:8:0.5"`
- `sys_probs: [...]` for manual weights
