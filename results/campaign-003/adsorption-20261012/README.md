# Campaign 003: bounded adsorption-relaxation record

**A, clean88 and gas3 converged. B remains partial and cannot be ranked against A.**

The recorded electronic adsorption energy for A is **-0.34611395321154986 eV per periodic cell and per CO2**, from the three converged energies in `data/A-conditional-adsorption-energy.json`. This is conditional on the specified dry parallel-OH slab, coverage, common periodic cell and PBE-D3(BJ)/Z-dipole method. The periodic gas reference is not a verified isolated-gas limit. No BSSE or vibrational correction is included; optimization convergence is not proof of a vibrational minimum.

This package preserves successful results, B's interruption and bounded continuations, frozen thresholds, inputs, small endpoint structures, scalar optimization snapshots and last complete force evaluations. It does not silently convert `partial_budget`, `reference_QA_passed` or exit zero into scientific convergence. The B energy is retained as partial evidence, not used for an A/B preference or adsorption-energy claim.

`METHODS.md` describes the model and criteria. `PROCESS.md` explains the stages, monitoring failure and budget accounting. `INSTRUCTIONS.md` and `REPRODUCE.md` define safe report replay and the requirements for a separate scientific rerun. `omitted-artifacts.json` gives portable original retention paths and hashes for raw evidence excluded from ordinary Git. All results remain development records with `manuscript_eligible=false` and no production qualification.
Four measured [figures and captions](analysis-figures/figure-captions.md),
[compact plot data and generators](analysis-figures/REPRODUCE.md) accompany the
record. The plot data reconstruct the complete original B tail from the saved
CP2K log and position/force files: 72 + 8 + 21 complete force evaluations,
98 completed geometry updates after excluding two restart re-evaluations.
`data/optimization-snapshots.csv` separately preserves the earlier parser
snapshot (53 complete evaluations in the initial B segment); it is not the
complete initial native trajectory. See the plot provenance for the distinction.
