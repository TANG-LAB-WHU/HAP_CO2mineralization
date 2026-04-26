# HAP_hkl_template Workflow

## Four-stage unified workflow

This template is now organized as a fixed four-stage pipeline:

1. AIMD generation (CP2K)
2. Training dataset synthesis (multi-structure pool)
3. DeepMD training and model export
4. LAMMPS large-scale simulation

## Stage 1: AIMD (CP2K)

- `run_workflow.bat` detects HKL from folder name (`HAP_XXX_template`)
- `generate_config.py` produces dynamic files:
  - `cp2k_run/HAP_<hkl>_md.inp`
  - `lammps_run/in.deepmd.lammps`
  - `docker-compose.yml`
- CP2K outputs trajectory, energy, force, and virial files in `cp2k_run/`

## Stage 2: Training dataset synthesis

- `generate_deepmd_data.py` converts CP2K outputs into `deepmd_data/`
- Workflow syncs `deepmd_data/` into pooled system path:
  - `deepmd_pool/HAP_<hkl>/`
- `prepare_multisystem_input.py` scans `deepmd_pool/*` and generates:
  - `deepmd_train/input.multisystem.json`

## Stage 3: DeepMD training and export

- `deepmd_train` reads `input.multisystem.json`
- `training.training_data.systems` uses pooled multi-structure systems
- sampling policy defaults to:
  - `auto_prob: "prob_sys_size"`
- `deepmd_freeze` exports:
  - `deepmd_model/hap_model.pth`

## Stage 4: LAMMPS large-scale simulation

- `lammps_run` loads `deepmd_model/hap_model.pth`
- executes `lammps_run/in.deepmd.lammps`
- writes trajectories and analysis outputs in `lammps_run/`

## Data pool convention

Each pooled system directory should contain:

- `set.000/`
- `type.raw`
- `type_map.raw`

Example:

- `deepmd_pool/HAP_100/`
- `deepmd_pool/HAP_112/`
- `deepmd_pool/HAP_300/`

## DeepMD multi-system probability options

In `deepmd_train/input.multisystem.json`, use one of:

- `auto_prob: "prob_sys_size"` (default)
- `auto_prob: "prob_uniform"`
- `auto_prob: "prob_sys_size;0:4:0.5;4:8:0.5"`
- `sys_probs: [...]` for manual weights
