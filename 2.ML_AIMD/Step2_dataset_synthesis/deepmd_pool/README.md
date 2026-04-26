# DeepMD Multi-Structure Data Pool

This directory stores pooled DeepMD systems used by the unified training stage.

Expected structure:

- `deepmd_pool/HAP_100/`
- `deepmd_pool/HAP_112/`
- `deepmd_pool/HAP_300/`

Each system directory should include:

- `set.000/`
- `type.raw`
- `type_map.raw`

Workflow behavior:

1. `run_workflow.bat` calls `Step2_dataset_synthesis/generate_deepmd_data.py`.
2. CP2K outputs from `Step1_aimd_cp2k_runs/HAP_<hkl>_Perfect` are converted into one DeepMD system.
3. The system is written to `Step2_dataset_synthesis/deepmd_pool/HAP_<hkl>/`.
4. `Step3_mlip_deepmd/prepare_multisystem_input.py` scans this pool and writes `Step3_mlip_deepmd/input.multisystem.json`.
5. DeepMD training reads `input.multisystem.json` from `Step3_mlip_deepmd`.
