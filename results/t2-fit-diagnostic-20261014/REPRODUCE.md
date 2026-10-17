# Reproducing the evidence presentation

This bundle reproduces inspection and plotting of saved training metrics. It
does not rerun training, load a model, or recreate an independent validation set.

## Portable commands

From the repository root:

```sh
python results/t2-fit-diagnostic-20261014/scripts/plot_and_verify.py --verify-only
python results/t2-fit-diagnostic-20261014/scripts/plot_and_verify.py --output t2-training-curves.png
```

Python's standard library is sufficient for verification. Plotting additionally
requires Matplotlib. The script verifies bundled content hashes, record counts,
accepted-step order, finite scalar metrics, final reported metrics and the first
simultaneous crossing of the two force thresholds. It reads only these CSV/JSON
files and performs no model inference. Its optional four-panel PNG is a new
layout of the saved metrics, not a byte-identical recreation of the archived
eight-panel figure. No smoothing, interpolation or synthetic points are added.

## Included files and transformations

| File | Purpose and transformation |
|---|---|
| `data/curve.csv` | Original canonical scalar history, byte-identical; steps 0–10,978. |
| `data/discarded-*.csv` | Three original discarded scalar histories, byte-identical and separate. |
| `figures/t2-fit-diagnostic.png` | Archived figure, byte-identical; explicitly labels training, not validation. |
| `summary.json` | Selected result/provenance/accounting fields; local paths removed. |
| `run-parameters.json` | Actual T2 objective, optimizer, precision, frame selection, resources and stopping rules. |
| `model-architecture.json` | The model object extracted from the original recipe, without legacy data/system paths. |
| `scripts/plot_and_verify.py` | New portable evidence-only presentation code. |

The original recipe also contains an earlier Adam schedule, system list and
20,000-step setting. Those fields were not the active T2 LBFGS procedure and are
deliberately not presented as reproduction instructions. Architecture seeds
describe inherited initialization, not a new random restart on T2.

`summary.json` records original result, score, completion, worker, helper, plot,
exported-model and checkpoint hashes. The original run-specific worker/helper
and plotting script depend on operational mounts and predecessor records and
are not advertised as portable executables here. No source is silently replaced:
the supplied presentation script is separately identified as new code.

## What this package cannot reproduce

Exact training continuation requires the original checkpoint, model, fixed
eight-frame energy/force labels, original environment and provenance records.
They are not included. The absence of these artifacts is explicit; publishing
the scalar curve does not create a fully reproducible training package. No
additional data or model evaluation is required to inspect the reported fit.

Discarded branches preserve the negative engineering history. Retained accepted
steps overlap with those branch step numbers after recovery; never concatenate
all CSV files as one optimization trajectory. Likewise, do not sum native
LBFGS function-evaluation counters with physical complete-objective counts.

These files are diagnostic repository evidence, not manuscript figures or
validated model outputs. They do not change manuscript or independent-validation
eligibility, and the 40-frame joint fit remains unexecuted.
