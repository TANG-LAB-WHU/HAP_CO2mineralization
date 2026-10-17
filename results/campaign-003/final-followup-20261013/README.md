# Campaign 003: completed follow-up record

B91 and the B-derived clean88 reference completed the original four optimization criteria. The A tilt test remains inconclusive. Four frozen own-basis fragments and four counterpoise (CP) points completed; their existing conditional summaries are preserved here.

This is a reporting-only continuation of [adsorption-20261012](../adsorption-20261012/README.md). It **supersedes the old B-partial completion status**, while preserving that package as the historical record. It does not establish a stable-minimum ranking or manuscript eligibility. No simulation, model evaluation, new energy/force calculation, or new scientific interpretation was performed to create this package.

| Completed follow-up | Existing outcome | Record |
|---|---|---|
| A rigid -5/+5 degree pair | Both singlepoints completed; directional test **inconclusive** | [A tilt](summaries/A-tilt.json) |
| B direction checkpoint | Five updates completed; still partial at this stage | [Checkpoint](summaries/B-direction-checkpoint.json) |
| B adaptive continuation | Converged after 11 further updates across two segments | [Common-reference energies](summaries/B-final-common-reference.json) |
| Four own-basis fragments | Four accepted singlepoints; subsequent summary parser repaired without new SCF/force calls | [Decomposition](summaries/fragment-decomposition-postprocessed.json) |
| B-derived clean88 | Converged after 31 updates across four segments; original clean reference retained | [Reference comparison](summaries/new-clean-common-reference.json) |
| Four ghost-fragment CP points | Four completed conditional finite-basis controls | [CP correction](summaries/cp-correction.json) |

The following values are copied from those existing summaries, in eV per periodic cell (one CO2 in each complex). They are conditional electronic-energy bookkeeping, not accuracy bounds or new calculations.

| Named reference/control | A | B | B minus A |
|---|---:|---:|---:|
| Original clean88, raw adsorption energy | -0.3461139532 | -0.9559156773 | -0.6098017241 |
| B-derived clean88, raw adsorption energy | -0.0426654984 | -0.6524672225 | -0.6098017241 |
| Original clean88, CP-corrected adsorption energy | -0.3276443481 | -0.9129776768 | -0.5853333286 |
| B-derived clean88, CP-corrected adsorption energy | -0.0241958933 | -0.6095292219 | -0.5853333286 |

Full existing scalar precision and field pointers are in [summary-values.csv](data/summary-values.csv) and the source summaries. [Point results](data/point-results.json) separate final energy-only records from complete force evaluations. [Execution records](data/execution-summary.json) retain exit/OOM and resource evidence; completion alone is not scientific validation.

Conditions are unchanged: the specified dry parallel-OH slab, fixed cell and coverage, PBE-D3(BJ), and Z dipole correction. The periodic gas reference is not a demonstrated isolated-gas limit. A's tilt remains inconclusive; neither optimization nor these finite-basis CP controls proves a vibrational/global minimum, basis-set limit, or functional accuracy. D3 is already included in total energies and is not added a second time. Ghost-center derivatives are not ML force labels. Negative reference-relative slab contributions are retained, not clipped into positive deformation costs. All records remain development evidence with `manuscript_eligible=false`.

See [PROCESS.md](PROCESS.md) for failures and repairs, [REPRODUCE.md](REPRODUCE.md) for report replay and omitted dependencies, and [source-records.json](source-records.json) for original byte hashes.
