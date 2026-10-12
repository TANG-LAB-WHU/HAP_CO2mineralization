# Local computational result packages

These packages preserve measured results, data, figure-generation inputs and limitations. They do not release a qualified interatomic potential or establish mineralization mechanisms.

- [Campaign 002: HAP91 closeout](campaign-002/closeout-20261012/README.md): completed DFT/data/software workflow, failed original force qualification, one bounded 20-fs learned trajectory with a failed 10-fs CO2 force screen, and a same-input GPU/CPU single-point comparison. Includes diagnostic PNG figures, captions, numerical tables and report replay.
- [Campaign 003: adsorption relaxation](campaign-003/adsorption-20261012/README.md): converged A/clean88/gas3 and a conditional electronic adsorption energy of -0.3461139532 eV per CO2. B remains partial, so no A/B ranking is reported. Inputs, bounded execution history and raw-evidence retention paths accompany the result.

Each package provides its own instructions, limitations and file-hash manifest. Large raw trajectories, wavefunctions, Hessians and model/checkpoint binaries remain outside ordinary Git. R4 numerical files stay sealed. Missing results and failed criteria remain explicit.

Campaign 003 also includes [four figures with captions](campaign-003/adsorption-20261012/analysis-figures/figure-captions.md) and [portable redraw instructions](campaign-003/adsorption-20261012/analysis-figures/REPRODUCE.md).
