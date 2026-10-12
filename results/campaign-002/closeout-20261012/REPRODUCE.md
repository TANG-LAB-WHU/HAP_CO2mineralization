# Reproduce saved results and figures

This is a portable report package for a completed, unsuccessful force-qualification campaign. It is not a self-contained training or simulation rerun.

From this directory:

```text
python -B verify_bundle.py
python -B make_figures.py
```

Verification uses Python's standard library and checks the listed file hashes, all 273 endpoint force components, energy error, and the 81-frame thermo table against the recorded summary. Chart replay needs NumPy and Matplotlib; the preparation environment used Python 3.12.4, NumPy 1.26.4 and Matplotlib 3.8.4. No installation, model evaluation or scientific execution is performed by either script.

`recorded-commands.json` preserves the actual historical Docker command structure with `${CAMPAIGN_ROOT}`, `${GPU_UUID}` and `${RUN_ID}` placeholders replacing local identities. These commands are provenance, not a ready-to-run procedure. They depend on omitted local guards, weights and restart assets. Scientific settings, pinned image digests, time/memory limits and thread counts remain in that record. `inputs/protocol.json` uses the same root placeholder and preserves the frozen scientific values; its original source hash is retained in `provenance.json`.

The runtime record reports DeepMD 3.1.0, PyTorch 2.6.0 and native executable hashes. CP2K reports version 2025.1, revision git:9635df4. The CP2K comparison uses pinned official CPU/CUDA images recorded in the historical commands and source references. `sources.json` identifies original local result records by campaign-relative paths and SHA256; `provenance.json` distinguishes copied files from publication transcriptions. Source hash identities do not claim that the original large assets are included.

The final potential did not pass original force criteria. R9/R10 are development diagnostics, R4 stays sealed, and the missing 20-fs DFT result must remain absent. No extra fitting, inference or dynamics is authorized by this package. A scientific rerun requires all matching original assets, complete reviewed execution code and a separately bounded run plan; this package only reproduces reporting and endpoint arithmetic.
