# Reproduce this report

This bundle is a static export of completed development records. Reading it requires no CP2K, GPU, model checkpoint, or sealed data.

1. Use `manifest.json` to verify the published byte hashes with a SHA-256 utility. The manifest excludes itself.
2. Use `source-records.json` to locate each retained original under `${RESEARCH_RUNS}` and verify its original SHA-256 if that archive is available. JSON copies normalize local paths and omit operational owner/container identifiers; their original and publication hashes therefore differ. Inputs and endpoint EXTXYZ files are byte-for-byte copies.
3. Inspect the existing summaries in `summaries/`. `data/summary-values.csv` gives the exact JSON pointer for each tabulated value. These values were copied, not recalculated from energies, forces, or structures. The README rounds them for readability.
4. Check `history/` and `data/execution-summary.json` alongside the scalar point records. In particular, preserve the fragment driver's postprocessing failure and the separate zero-SCF repair summary.

The original A91, clean88 and gas3 endpoints and the earlier B partial stages remain in the sibling `adsorption-20261012` package. This package adds the completed follow-ups; it does not replace that history. `inputs/` contains the exact frozen text inputs and coordinates used by the new completed jobs, including each bounded B/clean optimization segment. `dependencies.json` records available original seed/Hessian hashes. CP `centers.extxyz` contains all basis centers; the matching input KIND/GHOST definitions identify real nuclei and ghost centers.

This is not a self-contained scientific rerun bundle. Raw native logs, trajectories, WFN/Hessian binaries, runtime containers, pseudopotential/basis assets, telemetry, private dependency maps, and third-party materials are intentionally excluded. Their retained provenance appears in `source-records.json`, `omitted-artifacts.json`, the input dependency records and existing summary references. A separate authorized rerun requires the exact original dependencies and runtime/resource contract; do not treat omitted files as permission to substitute an atomic guess, Hessian, cell, method, or ghost mask. Final energy-only records must never become force labels. No H/R4 sealed numerical files were read or hashed to produce this package.
