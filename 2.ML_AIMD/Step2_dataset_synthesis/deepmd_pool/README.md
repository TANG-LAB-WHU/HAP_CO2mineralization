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

1. `run_workflow.bat` generates `deepmd_data/` for current structure.
2. The dataset is synced to `deepmd_pool/<HAP_PREFIX>/`.
3. `prepare_multisystem_input.py` scans `deepmd_pool/` and builds `deepmd_train/input.multisystem.json`.
4. DeepMD training reads the generated multi-system input.
