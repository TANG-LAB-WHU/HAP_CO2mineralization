# Campaign 003: conditional A adsorption energy; B remains partial

These are measured development results from the dry parallel-OH, two-repeat
91-atom HAP-CO2 model. This report adds no scientific calculation or authorization.
All source artifacts retain manuscript_eligible=false; evidence reconciliation,
method convergence and scientific publication qualification remain separate.

A91, clean88 and gas3 each passed the same four optimization criteria: maximum
force component 0.010 eV/A, RMS force 0.005 eV/A, maximum displacement 0.003 bohr
and RMS displacement 0.0015 bohr. The matched periodic-cell PBE-D3(BJ)/Z-dipole
reference difference is

`E_ads(A) = E(A91, relaxed) - E(clean88, relaxed) - E(gas3, relaxed)`

**E_ads(A) = -0.34611395321154986 eV per cell and per CO2.** The electronic result
is conditional on the specified slab, OH arrangement and coverage. The gas
reference remains periodic and its isolated limit has not been verified. There
is no BSSE correction, vibrational minimum proof, finite-temperature free energy
or inference of mineralization. The earlier 47-atom fixed-geometry interaction
energies use a different model and cannot be combined with this value.

B91's final continuation completed 20 geometry updates and 21 complete force
evaluations, ending at 2026-10-12 22:43:26 UTC after 1864.563 s for that segment.
It stopped partial_budget. Its final maximum/RMS forces are approximately
0.0482595/0.00620484 eV/A; maximum/RMS displacements are
0.0944863066/0.0091575896 bohr. All four criteria remain unmet. Reference and
geometry QA passed, which does not establish an optimized endpoint. B is
preserved as a partial optimization, excluded from the adsorption result and
from any A/B ranking. No B convergence or ranking is inferred from its energy.

## Compact data and provenance

- results-summary.json retains exact scalar results, gates, method limitations,
  campaign-relative source paths and SHA256 values.
- compact-data.csv has one row per system/recorded segment. It uses energies in
  hartree, force components/RMS in eV/A and displacements in bohr.
  B's counters describe its final continuation, not all historical attempts.
- A/clean/gas energies are the final energy-only evaluations, each with converged
  SCF and zero geometry difference from the optimized endpoint. Forces and gates
  come from the last complete force evaluations; the energy-only frames remain
  excluded from force labels. These roles match the original conditional record.
- Source roots are .research-runs/campaign-003/first-batch-run/jobs/A91/native,
  references-run/jobs/clean88/native, references-run/jobs/gas3/native and
  B-continuation1-run/jobs/B91/native. Original logs, trajectories and endpoint
  files remain there. The conditional record is
  references-run/A-conditional-adsorption-energy.json. The older initial B
  result.json is a mid-run snapshot and is not used as the final B record.

The claim ledger uses C-C003-A-ADSORPTION and C-C003-B-PARTIAL. These remain
pending development claims, with no change to human verification or eligibility.
Four source-backed development figures now accompany this record: initial/final
local structures, the segmented optimization history, A's reference subtraction,
and the earlier v2 validation workflow. They retain all partial/failed scientific
criteria. None is an A/B energy-ranking plot or a publication-eligibility decision.
Existing campaign-002 diagnostic files remain unchanged. See REPRODUCE.md and
figure-captions.md for rendering instructions and the intended reading of each
figure; figure_manifest.csv and FILE_MANIFEST.sha256 identify this snapshot.

## Next scientific priority (proposed only)

Resolve the B endpoint under the same scientific method and four criteria before
attempting a two-configuration comparison. First use the complete saved tail to
decide whether another bounded continuation is justified; do not automatically
extend time, relax criteria or declare a distinct minimum. If B converges to A's
basin, report that outcome; if it remains partial, retain that limitation without
ranking. Method/periodic-reference/basis and BSSE controls follow only for the
claims that require them. Larger cells, ML retraining, MatterGen and pathways
remain downstream questions, not actions authorized by this reporting update.
