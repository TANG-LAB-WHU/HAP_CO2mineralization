# B182 diagnostics and short trajectories

This bundle preserves three distinct evidence classes: static CP2K representation checks, an independent 4 fs CP2K trajectory, and an 85 fs exploratory trajectory made with a potential that failed qualification, followed by a full-system CP2K force diagnostic. It contains compact numerical evidence, selected structures, input records, and reproducible plots. **All results remain manuscript-ineligible and do not qualify the potential for production dynamics.**

The starting B182 structure is an exact 2 × 1 × 1 repetition of B91: **176 slab atoms + 2 CO₂ molecules = 182 atoms**, with composition Ca₄₀P₂₄O₁₀₈H₈C₂. This preserves CO₂ coverage while increasing lateral area. It does not increase slab thickness or establish thickness convergence. It is distinct from the earlier 179-atom construction containing one CO₂ molecule. Some original job folders used the operational alias `B91`; the records here identify the actual atom count.

The common periodic cell, in Å, is:

```text
18.848   0.000000000000   0.0000000
-4.712   8.161423405265   0.0000000
 0.000   0.000000000000  37.7401146
```

All three directions are periodic. CO₂ atom IDs are 89–91 and 180–182, using one-based indexing.

## What each calculation establishes

| Evidence | Recorded result | Interpretation and limits |
| --- | --- | --- |
| Static B182 and clean176 | One energy/force evaluation per structure; zero geometry updates. Both completed the recorded reference and physical checks. | Checks the expanded representation and forces at supplied coordinates. Neither structure was optimized in these jobs. The pair alone does not establish adsorption energy, thermal stability, or cell-size convergence. |
| Direct CP2K B182 dynamics | NVE, initial 300 K, 0.25 fs step, 16 steps, 17 frames, **4 fs** total. Recorded physical and energy-drift checks passed. | A short direct-DFT response. It is independent of the failed ML route and does not demonstrate equilibration or ML accuracy. |
| Exploratory ML B182 dynamics | NVE, initial 300 K, 0.25 fs step, 341 frames through step 340, **85 fs** total. An O–H distance alarm first occurred at 82.75 fs. | Diagnostic exploration with an unqualified potential. It stopped at 85 fs; it was not a completed 200 fs trajectory. |
| CP2K at the ML endpoint | One full-182-atom energy/force evaluation at the exact 85 fs coordinates; no optimization. | Measures force disagreement at one problematic configuration. It is not a new validation set or evidence for a reaction pathway. |

The static maximum atomic force norms are 0.0417311 eV/Å for B182 and 0.0412289 eV/Å for clean176. Their respective maximum component differences between the two replicated copies are approximately 5.21 × 10⁻⁶ and 4.88 × 10⁻⁶ eV/Å. These comparisons support numerical consistency of the replicated structures, not exact energy scaling or a hardware performance claim.

The direct 4 fs trajectory has a recorded maximum absolute total-energy departure of 0.00742467 meV/atom and temperatures of 256.718–300 K. Its short duration is a material limitation even though the recorded checks passed.

## The 85 fs failure and endpoint diagnostic

The H109–O135 distance starts at 0.977993 Å, first exceeds the 1.4 Å diagnostic threshold at step 331 (82.75 fs; 1.405479 Å), and reaches 1.498646 Å at step 340 (85 fs). The recorded CO₂ and phosphate checks remained intact through these frames. The small recorded energy drift, at most 0.00738466 meV/atom, does not establish the accuracy of the potential.

The stop procedure also had an engineering failure: an included control file used `jump SELF stop_exploration`, while the target label belonged to the main input. LAMMPS exited with code 1 and the enclosing gate with code 2. The completed 341 frames are retained, and the failure is not relabeled as a successful production run. The archived control input preserves this defect.

At the same 85 fs coordinates, the ML-minus-DFT force-component RMSE is **1.257604 eV/Å across all 182 atoms** and **7.603828 eV/Å for H109 and O135**. The largest component difference is 17.925892 eV/Å. The local neighborhood was frozen before the comparison as the union of atoms within 3.5 Å of either endpoint under minimum-image distances: 20 atoms including the O–H pair, or 18 excluding it.

Both methods give a restoring relative force projection along the O135 → H109 axis: −1.296922 eV/Å for ML and −20.200598 eV/Å for DFT. These are projections of the atoms' total forces in their environment, not isolated O–H bond forces. They do not establish opposite force directions, a reaction mechanism, or the cause of the trajectory's departure. No mean-force subtraction or force correction was applied.

## Files and figures

| Location | Contents |
| --- | --- |
| `summary.json` | Compact findings and explicit validity limits. |
| `provenance.json` | Original source hashes, portable source identifiers, and extraction or path-replacement records. |
| `inputs/run-parameters.json` | Run settings and external asset requirements. |
| `tables/static-singlepoints.csv` | B182 and clean176 static results. |
| `tables/dft4fs-curve.csv` | All 17 direct-DFT trajectory records. |
| `tables/ml85fs-curve.csv`, `tables/ml85fs-OH-distances.csv` | All 341 exploratory trajectory records and eight tracked O–H distances. |
| `tables/step340-force-metrics.csv` | Endpoint force-error summaries by group. |
| `tables/step340-all-atom-forces.csv`, `tables/step340-local-forces.csv` | Raw ML and DFT force components, their differences, and neighborhood membership. |
| `snapshots/` | Selected structures serialized from recorded coordinates. |
| `inputs/cp2k/`, `inputs/lammps/` | Portable input records, with qualifications in [REPRODUCE.md](REPRODUCE.md). |
| `scripts/plot_results.py` | Plotting and input-verification entry point. |

Figures: [direct DFT, 4 fs](figures/dft4fs.png), [exploratory ML, 85 fs](figures/ml85fs.png), [O–H diagnostic](figures/OH-diagnostic.png), and [endpoint force errors](figures/step340-force-errors.png).

Large trajectories, complete native logs, model files, and wavefunction files are excluded. The exploratory model is identified by SHA-256 `57c4fc2d7b6890dccdb1346c9d5641f18047feeeb3800df94028c0bbd993fe54`; its omission does not change its failed qualification status. A repository publication of these diagnostic artifacts is not approval to use them as manuscript evidence.
