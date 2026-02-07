# HAP_hkl_template Workflow

- **Legend**
  - **Paths**: relative, wrapped in backticks (e.g., `cp2k_run/*`)
  - **Containers**: from `docker-compose.yml` (e.g., `cp2k_run`)
  - **Entry point**: `run_workflow.bat`

### 1) User Config (`run_workflow.bat`)
- **Surface**: `MILLER=(1,1,2)`, `LAYERS=4`, `SUPERCELL=(3,3)`, `VACUUM=15 Å`
- **CO2**: `ADD_CO2=true`, `N_CO2=1`, `CO2_HEIGHT=10 Å`
- **Method**: `SLAB_METHOD=ase`
- **Logs**: `logs/*`
- Proceeds to surface and initial LAMMPS structure generation

### 2) Surface & LAMMPS Data Prep (container: `python:3.10-slim`)
- **Inputs**: `lammps_run/HAP_ConventionalCell.cif`, `lammps_run/create_hap_slabs.py`
- **Run**: install `ase/pymatgen` → `create_hap_slabs.py` (optional `--add-co2`)
- **Output**: `lammps_run/hap_hkl.data`; ensure `lammps_run/model/`

### 3) CP2K stage (service: `cp2k_run`)
- **Input**: `cp2k_run/HAP_hkl_md.inp`
- **Run**: `cp2k.psmp` (GPU, MPI)
- **Outputs**:
  - `cp2k_run/HAP_hkl_md-pos.xyz`
  - `cp2k_run/HAP_hkl_md-1.ener`
  - `cp2k_run/HAP_hkl_md-1.cell`
  - `cp2k_run/HAP_hkl_forces.dat`
  - `cp2k_run/HAP_hkl_virial.dat`
  - Final frame: `cp2k_run/HAP_hkl_final.xyz`

### 4) DeepMD dataset generation (container: `python:3.10-slim`)
- **Script**: `generate_deepmd_data.py` (installs `requirements.txt`)
- **Parses**: XYZ/Cell/Energy/Forces/Virial (with robust fallbacks)
- **Outputs**: `deepmd_data/`
  - `set.000/coord.npy(.raw)`, `box.npy(.raw)`, `force.npy(.raw)`, `energy.npy(.raw)`, `virial.npy(.raw)`
  - `type.raw`, `type_map.raw`
- Helper: `deepmd_prep` copies `type.raw` if needed

### 5) DeepMD training (service: `deepmd_train`)
- **Image**: `deepmodeling/deepmd-kit:3.1.0_cuda129`
- **Inputs**: `deepmd_data/*`, `deepmd_train/input.json`
- **Run**: `dp --pt train input.json -o model`
- **Outputs**: `deepmd_train/checkpoint`, `model.ckpt-*.pt`, `lcurve.out`

### 6) Freeze to PyTorch model (service: `deepmd_freeze`)
- **Run**: `dp --pt freeze -o /mnt/shared/model/hap_model.pth`
- **Output**: `deepmd_model/hap_model.pth`
- Optional test: `deepmd_test` (`dp --pt test`)

### 7) LAMMPS run with DeepMD model (service: `lammps_run`)
- **Prereqs**: `lammps_run/hap_hkl.data`, `deepmd_model/hap_model.pth`
- **Actions**: copy model to `lammps_run/model/` → run `lmp -in lammps_run/in.deepmd.lammps`
- **Outputs**: `lammps_run/*` (trajectories/logs)

### Global logs (from `run_workflow.bat`)
- `logs/workflow_log_...md`
- `logs/workflow_errors_...txt`
- `logs/terminal_output_...txt`
- `logs/combined_output_...txt`

### Artifact map
- **CP2K outputs**: `cp2k_run/*`
- **DeepMD dataset**: `deepmd_data/set.000/*`, `deepmd_data/type*.raw`
- **Training checkpoints**: `deepmd_train/*`
- **Frozen model**: `deepmd_model/hap_model.pth`
- **LAMMPS I/O**: `lammps_run/*`

## For example, if hkl was changed to 112, then update the names in the corresponding scripts in line with the instructions listed below to maintain consistency.

### 1. **File Renames Required**
```
cp2k_run/HAP_hkl_md.inp → cp2k_run/HAP_112_md.inp
cp2k_run/HAP_hkl_opt-pos-1.xyz → cp2k_run/HAP_112_opt-pos-1.xyz
HAP_hkl_specification.txt → HAP_112_specification.txt
```

### 2. **Content Modifications Required**

#### A. **workflow_documentation.md** (7 locations)
- Line 1: `### HAP_hkl_Workflow` → `### HAP_112_Workflow`
- Line 18: `lammps_run/hap_hkl.data` → `lammps_run/hap_112.data`
- Line 21: `cp2k_run/HAP_hkl_md.inp` → `cp2k_run/HAP_112_md.inp`
- Lines 24-29: All `HAP_hkl_*` references → `HAP_112_*`
- Line 51: `lammps_run/hap_hkl.data` → `lammps_run/hap_112.data`

#### B. **run_workflow.bat** (11 locations)
- Line 89: `--output hap_hkl` → `--output hap_112`
- Lines 93-97: All `hap_hkl.data` references → `hap_112.data`
- Lines 113-116: All `hap_hkl.data` references → `hap_112.data`
- Line 139: `hap_hkl.data` → `hap_112.data`
- Lines 211-212: `HAP_hkl_opt-pos-1.xyz` → `HAP_112_opt-pos-1.xyz`
- Line 225: `HAP_hkl_opt-pos-1.xyz` and `HAP_hkl_final.xyz` → `HAP_112_*`
- Lines 233, 272-273, 286-287: All `HAP_hkl_md-pos.xyz` → `HAP_112_md-pos.xyz`
- Line 458: `hap_hkl.data` → `hap_112.data`

#### C. **lammps_run/in.deepmd.lammps** (1 location)
- Line 14: `read_data hap_hkl.data` → `read_data hap_112.data`

#### D. **cp2k_run/HAP_hkl_md.inp** (2 locations)
- Line 1: `@SET PROJECT HAP_hkl` → `@SET PROJECT HAP_112`
- Line 122: `@include 'HAP_hkl_final.xyz'` → `@include 'HAP_112_final.xyz'`

#### E. **lammps_run/create_hap_slabs.py** (1 location)
- Line 640: `workflow_file = "hap_hkl.data"` → `workflow_file = "hap_112.data"`

#### F. **generate_deepmd_data.py** (1 location)
- Line 738: `project_name = "HAP_hkl"` → `project_name = "HAP_112"`

#### G. **docker-compose.yml** (1 location)
- Line 25: `HAP_hkl_md.inp -o HAP_hkl_md.log` → `HAP_112_md.inp -o HAP_112_md.log`

#### H. **lammps_run/VMD visuaization.md** (1 location)
- Line 297: `hap_hkl.data` → `hap_112.data`

## Summary

**Total files requiring modification: 8 files**
- **3 files** need to be renamed
- **8 files** need content modifications
- **Total locations to change: 25+ references**

The changes follow a consistent pattern:
- `HAP_hkl` → `HAP_112` (for project identifiers and file prefixes)
- `hap_hkl` → `hap_112` (for LAMMPS data files)
- All file references and paths need to be updated accordingly

This systematic approach ensures complete consistency across the entire workflow when changing from "hkl" to "112" naming convention.