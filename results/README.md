# Local computational result packages

These packages preserve measured results, data, figure-generation inputs and limitations. They do not release a qualified interatomic potential or establish mineralization mechanisms.

- [MatterGen generation records and correction](mattergen-generation-20261011/README.md): the early 8-atom smoke and the later four-candidate 44-atom pilot are separate completed generation records. All five samples were rejected. The 44-atom run is an out-of-distribution atom-count adaptation, not fine-tuning or fixed-composition CSP training.
- [17 October release and replay scope](RELEASE-20261017.md): later Campaign 003 optimization/reference/CP results; B-local development and the earlier failed 40-frame fit; 182-atom static, direct-DFT and failed learned-trajectory diagnostics; and the completed T2 training-fit diagnostic. The **new joint 40-frame fit was not executed**.
- [Campaign 003 final follow-up](campaign-003/final-followup-20261013/README.md): later completion states, conditional reference comparisons, and inconclusive cases, preserved separately from the original snapshot.
- [B-local development](b-local-development-20261013/README.md): candidate training and failed validation/coverage gates.
- [182-atom diagnostics](b182-diagnostics-20261013/README.md): geometry failure and DFT force mismatch are reported separately.
- [T2 training-fit diagnostic](t2-fit-diagnostic-20261014/README.md): eight original training frames fitted to diagnostic thresholds; no independent model qualification.
- [Campaign 002: HAP91 closeout](campaign-002/closeout-20261012/README.md): completed DFT/data/software workflow, failed original force qualification, one bounded 20-fs learned trajectory with a failed 10-fs CO2 force screen, and a same-input GPU/CPU single-point comparison. Includes diagnostic PNG figures, captions, numerical tables and report replay.
- [Campaign 003: original adsorption-relaxation snapshot](campaign-003/adsorption-20261012/README.md): converged A/clean88/gas3 and a conditional electronic adsorption energy of -0.3461139532 eV per CO2. B was partial at that snapshot; the final-follow-up package records later work. Inputs, bounded execution history and raw-evidence retention paths accompany the result.

Each package provides its own instructions, limitations and file-hash manifest. Large raw trajectories, wavefunctions, Hessians and model/checkpoint binaries remain outside ordinary Git. R4 numerical files stay sealed. Missing results and failed criteria remain explicit.

From the repository root, `python -B results/verify_release_20261017.py` checks the new release's file inventory without starting a scientific calculation. Follow each package's `REPRODUCE.md` to redraw its figures from bundled data.

Campaign 003 also includes [four figures with captions](campaign-003/adsorption-20261012/analysis-figures/figure-captions.md) and [portable redraw instructions](campaign-003/adsorption-20261012/analysis-figures/REPRODUCE.md).
