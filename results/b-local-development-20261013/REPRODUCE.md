# Reproduce the archived diagnostic plots

## Requirements

Use Python 3 with Matplotlib installed. The plotting script uses only the Python standard library and Matplotlib; it does not require DeepMD, PyTorch, CP2K, LAMMPS, CUDA or a GPU. It uses Matplotlib's noninteractive `Agg` backend.

## Run from the repository root

```sh
python results/b-local-development-20261013/plot_results.py
```

Alternatively, run from this package directory:

```sh
python plot_results.py
```

The script resolves its input relative to its own location, so both commands read the same `training-curve-snapshots.csv`. They write or replace `training-snapshots.png` and `training-snapshots.svg` in this directory. Font and rendering differences between Matplotlib versions or platforms may change image bytes; the numerical input is the archived CSV.

## What this reproduces

The CSV contains **50 existing score-snapshot records**, extracted from completed training outputs. The script groups B16 continuations by cumulative accepted-update count and plots the old development40 fit separately. Repeated records at the same accepted-update count are represented by one plotted point. It plots recorded pooled all-atom and CO2 component RMSE values and fixed threshold reference lines.

This is a redraw only: no model loading, new prediction, objective evaluation, rescoring, optimization, trajectory generation or label processing occurs. Pooled curves do not replace the per-family qualification checks. The B-H independent test failed, the old development40 fit failed training qualification, and the new joint40 CPU proposal was not executed. These diagnostic figures are not automatically eligible for manuscript use.

## Provenance and limits of computational replay

`source-provenance.csv` records the original byte SHA256 and location of each source summary used for the local-development tables and curves. `basin-potential-summary.json` separately records the sources for the preceding 142-frame fit and optimization diagnostic. Paths beginning `.research-runs/` refer to external artifacts relative to the original project root; `workspace/` denotes external preparation metadata. These files need not exist in a fresh Git clone.

Model/checkpoint hashes in `model-lineage.json` are values reported by the archived terminal records. The model and checkpoint binaries were neither copied nor rehashed for this publication package. Raw labels, trajectories, wavefunctions, runtime releases and the complete training/evaluation worker source bundles remain external.

Consequently, this package supports numerical-summary inspection and plot replay, **not a standalone rerun of training or inference**. Such a rerun would require the matching external data and model assets, original recipe/protocol, pinned software environment, and exact worker/helper source bytes verified against the original immutable runtime-release SHA256 manifest. A summary-file hash is not a worker-source hash; installing a current worker or library version is not a verified substitute. No new training, evaluation or sealed-label access is part of these reproduction commands.
