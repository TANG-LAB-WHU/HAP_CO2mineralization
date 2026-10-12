# Reproduce the campaign-003 development figures

This folder is the self-contained reporting snapshot for repository publication.
It has measured compact data, provenance, scripts and four PNG figures. It
does not launch CP2K, train a potential, run LAMMPS, install dependencies or access
sealed R4 labels. Original runtime logs and restart/checkpoint files remain under
the campaign-relative source paths recorded in the data and summary files.

## Render from included compact data

Use an existing Python environment with NumPy and Matplotlib. From this folder:

```text
python -B plot_structures.py
python -B plot_optimization.py
python -B plot_reference_and_workflow.py
```

The original snapshot used Python 3.12.4, NumPy 1.26.4 and Matplotlib 3.8.4.
The published generators emit PNG only. The default rendering path reads the
included data and starts no scientific engine. Structure rendering also rewrites
its derived coordinate/measurement tables. Plot-library differences may change
image bytes. `${CAMPAIGN_003_ROOT}` denotes the original raw-run directory;
it is provenance metadata, not needed for default rendering.

## Data and processing conventions

- results-summary.json and compact-data.csv preserve the A/clean/gas final-energy
  convention and the separate complete-force evaluations used for convergence.
- data/structures-* contains actual initial/endpoint coordinates, cell and atom
  order. Periodic images are derived from the recorded triclinic cell. Plot labels
  use one-based atom IDs; dashed contacts are geometric distances, not bond orders.
- data/optimization-* preserves the initial B segment, engineering recovery and
  final continuation. Duplicate restart evaluations remain in the data but are
  not counted as new geometry updates. The interrupted SCF-only step and final
  energy-only records do not become invented force labels or zero-valued points.
- data/reference-workflow.json records the actual A subtraction and the earlier
  campaign-002 v2 metrics. The workflow diagram is an execution summary, not a
  causal explanation of model error or a claim that failed gates were passed.
- Forces are eV/angstrom; displacement criteria are bohr; total reference energies
  are hartree and adsorption differences eV per cell/CO2. No interpolation,
  fitted energy offset, force correction or synthetic scientific values is used.

## Scope of the four figures

B remains a partial optimization; no two-site ranking or distinct-minimum claim
is supported. A's value is conditional on the specified dry periodic slab and
coverage, with no isolated-gas, BSSE or vibrational-minimum qualification. The v2
diagram keeps integration/reference QA separate from failed force accuracy and
shows the separately authorized diagnostic/exploration transitions explicitly.

The source commit identifies the original checkout baseline. Source snapshot
identity and the publication-only path/output adaptations are recorded in
publication-provenance.json. FILE_MANIFEST.sha256 identifies this portable copy.
No scientific execution follows from replaying the plots.
