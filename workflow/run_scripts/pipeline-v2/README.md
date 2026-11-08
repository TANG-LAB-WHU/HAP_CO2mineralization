# Hand-written scripts of the local pipeline v2 run

Copies, as run, of the scripts written by hand inside the run directory `.research-runs/pipeline-v2/` (not in Git),
kept here so that the results they produced in `workflow/PIPELINE_V2.md` have their source. Each sits at its path
relative to the run directory and uses the paths of the workstation (Docker mounts `/work` = repository,
`/data` = the directory given to the CP2K image). Copied 2026-11-08; the orchestrator (`workflow/pipeline_v2.py`)
does not call any of them.

| Script | Used for (section of PIPELINE_V2.md) |
|---|---|
| `cp2k-protocol/run_all.sh` | CP2K cutoff and cell-width tests that set the labelling protocol |
| `benchmarks/tf-train-ab/run.sh` | TensorFlow vs PyTorch training A/B on the iteration-0 model-0 configuration |
| `dft-hull/run_job.sh` | CP2K cell optimisations of the S4 DFT hull (resumes a stopped optimisation) |
| `metad/seed_scheduler.sh` | Workstation-only: pauses the `seed-metad` label container during training and labelling |
| `preview-iter3/run_preview.sh`, `preview-iter3/w8-rerun/analyse_rerun.py` | Stage-7 preview on the iteration-3 ensemble and the hap001_water_co2 rerun |
| `preview-iter4/run_preview.sh` | Stage-7 preview on the iteration-4 ensemble |
| `diagnostics/s3-energy-it4/bias.py` | S3 energy offset per system in the iteration-4 sealed test |
| `diagnostics/npt-volume-it4/pressure.py` | Model against DFT pressure on the labelled bulk frames (NPT exploration) |
| `diagnostics/s3_055-transform/relax.py`, `compare.py` | s3_055 → aragonite-type CaCO3: DP relaxation and structure match against the DFT hull |
