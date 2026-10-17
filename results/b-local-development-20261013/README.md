# B-local potential development and coverage expansion

This package preserves completed development results and negative outcomes from 13 October 2026. It is a result archive, not a release of a validated potential. Source records retain `manuscript_eligible=false`. No model, checkpoint, wavefunction, raw trajectory or sealed validation-label arrays are included.

**The old 40-frame coverage fit was executed and exported an unqualified model. The new joint40 CPU comparison/refit was prepared but has not executed.** They are separate experiments.

## Results

Force errors below are component RMSE in eV/Å. Numerical values are copied from existing completed score summaries; no inference or rescoring was performed for this package.

| Experiment | Completed work | Force RMSE | Outcome |
|---|---|---:|---|
| Broad 142-frame fit | 20,000 retained Adam updates, with engineering recoveries | 0.141052068 | Original training gates failed; no holdout release. |
| Optimization diagnostic | Additional 4000 updates; 24,000 retained total | 0.128271480 | Pooled improvement 9.06%, below the required 20%; original training gates still failed. |
| B16 initial fit | 300 accepted LBFGS updates; 630 complete objectives | 0.040489081 | Training gates failed; retained as a negative result. |
| B16 continuation 1 | 600 additional accepted updates; cumulative 900 / 1901 objectives | 0.028458995 | Family training gates still failed despite pooled RMSE <0.03; model and LBFGS history preserved. |
| B16 continuation 2 | 197 additional accepted updates; cumulative 1097 / 2317 objectives | 0.027019878 | Training gates passed at the first qualifying accepted state. This establishes fit quality on the 16 development frames only. |
| Independent B-H test of that B16 model | One trajectory; four selected correlated frames; one consumed score attempt | 0.079154306 | Test failed. Maximum force-component error 0.428435005 eV/Å; CO2 force RMSE 0.121064967 eV/Å. No learned MD release. |
| Old development40 coverage fit | 569 accepted updates; 1208 complete objectives; exact optimizer recovery after an engineering interruption | 0.058128689 | Export completed, but training gates failed at the fixed optimization time limit. CO2 force RMSE 0.076647458 eV/Å. |
| New joint40 CPU comparison/refit | Prepared only; no candidate comparison or training started | — | **PLANNED / NOT EXECUTED**. No result can be attributed to this run. |

The broader 142-frame fit and the subsequent 4000-update optimization diagnostic are documented in [basin-potential-summary.json](basin-potential-summary.json), including the exact source-summary hashes. They precede B16 and remain part of its initialization history; a successful local B16 fit does not erase their broader-domain failure. The 20,000-update lineage includes four non-bitwise sampling resets and 768 lost updates. The optimization diagnostic improved the pooled RMSE and each macrofamily RMSE, but the B-family maximum component error increased from 0.351523 to 0.563608 eV/Å; improvement must not be generalized to every metric.

## Data and method

B91 contains Ca20 P12 O54 H4 C1: an 88-atom slab plus one CO2. The periodic cell is recorded in [parameters.json](parameters.json). Labels use PBE-D3(BJ), GPW 1200/80 Ry, inner/outer SCF tolerance 1e-8, DZVP-MOLOPT-SR-GTH / GTH-PBE, neutral singlet and the frozen campaign-003 method assets.

B16 uses 12 B-optimization frames and four B-T frames. Its deterministic objective weights these families 0.5 / 0.5. LBFGS uses learning rate 1, history 20, strong-Wolfe line search, microbatch two, and a complete objective accumulated over every selected frame before each optimizer update. Authorized continuations preserve the earlier failures and exact accepted model/optimizer state.

The old 40-frame dataset comprises B-optimization (12), B-T (4), **B-H-development (4)**, T1 (8), T2 (8), and paired perturbations (4). Group weights are 0.25 for optimization and 0.15 for each other group, uniform within a group. The pooled RMSE is an unweighted component SSE/count aggregate, not the weighted optimization objective. Complete group metrics, thresholds and dataset identifiers are retained in the CSV/JSON files.

B-H was independent when its single test failed. Its four selected frames were subsequently exposed and explicitly promoted to development data for the 40-frame fit. They must not be described as an independent test of that later model. H2/H3 were planned as new tests but were not generated after the training failure; their labels and R4 remain outside this package.

T1 and T2 each use 32 NVE steps at 0.25 fs (8 fs), initialized at 300 K with independent seeds 2026101401 and 2026101402. Each retains 33 raw frames and selects steps 4, 8, …, 32. P+/P− perturb the adsorbate and S+/S− perturb the substrate; the exact deterministic definitions are in `parameters.json`. All six canonical development jobs passed their recorded reference/physical QA. The interrupted first T1 attempt was retained separately and was not concatenated with its same-seed rerun.

Training qualification requires pooled and each-group force RMSE ≤0.03 eV/Å, each-group CO2 force RMSE ≤0.05 eV/Å, and CO2 maximum component error ≤0.2 eV/Å. In the old 40-frame terminal model, only B-optimization passes all its group checks. Low energy error alone does not qualify either failed model.

The 40-frame run started with fresh LBFGS history from the accepted B16 model, keeping its statistics and energy bias fixed. After 246 accepted updates, a telemetry artifact access failure caused an engineering interruption. The continuation restored its model and optimizer exactly, added 323 accepted updates, and stopped at 569 total. Recorded effective fit use was 3455.221851 s, including the prior segment. The exported model SHA256 is `57c4fc2d7b6890dccdb1346c9d5641f18047feeeb3800df94028c0bbd993fe54`; this identifier is copied from its terminal record, not a fresh binary hash.

The subsequent direct B182 DFT fallback is a separate calculation, not a successful MLMD deployment; its trajectory is documented by the dedicated B182 result package. This package does not duplicate that trajectory or the later T2-only diagnostics.

## Files and provenance

- [results-summary.csv](results-summary.csv): terminal states and accepted/objective counters. B16 continuation counters are cumulative; development40 counters refer to that one logical fit.
- [aggregate-metrics.csv](aggregate-metrics.csv): original pooled and family force/energy/CO2 metrics, with exact score-source hashes.
- [training-curve-snapshots.csv](training-curve-snapshots.csv): existing checkpoint-score snapshots, not newly sampled or recomputed scores.
- [training-snapshots.png](training-snapshots.png) and [SVG](training-snapshots.svg): diagnostic pooled curves from those archived snapshots. The CO2 group threshold is shown only as a reference; each family must pass separately, and pooled curves alone do not establish qualification. These are not manuscript-qualified figures or holdout learning curves.
- [holdout-test-summary.json](holdout-test-summary.json): once-only B-H test, failed checks and provenance linkage.
- [parameters.json](parameters.json): frozen method, data membership, seeds, optimizer and thresholds.
- [model-lineage.json](model-lineage.json): terminal-reported model/checkpoint hashes and export checks; binaries remain external.
- [engineering-and-negative-results.json](engineering-and-negative-results.json): interruptions, recoveries and scientific failures retained separately.
- [joint40-planned-status.json](joint40-planned-status.json): prepared-only CPU proposal, 8 CPU / decimal 32 GB / four-hour effective budget. `run1` did not exist at the publication check.
- [source-provenance.csv](source-provenance.csv): byte hashes and sizes of the existing source summaries. `.research-runs/…` paths are relative to the project root; `workspace/…` denotes preparation metadata retained outside Git. These hashes identify source records, not the curated derivatives.
- [plot_results.py](plot_results.py): portable plotting of the included aggregate snapshots only. Requires Python and Matplotlib; run `python plot_results.py` from this directory.

All trajectories are very short local responses. They do not establish thermal equilibration, adsorption free energies, reaction rates, mineralization, long-time stability or transferability to a larger system. No independent H2/H3 success, successful learned B182 trajectory or joint40 result is implied.
