# MatterGen generation: smoke and HAP44 pilot

**MatterGen was run beyond the earlier one-sample smoke: four 44-atom candidates were generated on 2026-10-11. All four failed the exact-composition gate; none advanced to relaxation, hull evaluation or DFT.** This package corrects a smoke-only account of progress. Generation completed; accepted HAP discovery did not.

## Two separate records

The earlier [smoke record](smoke/context.json) used a pretrained chemical-system checkpoint, seed 42 and Ca-P-O-H-C conditioning. It generated one 8-atom CH2CaO4 structure in 176.47064661979675 seconds. The evaluation completed in 341.66393518447876 seconds using **MatterSim-v1.0.0-1M with TRI2024**, not CP2K or DFT. Its necessary HAP chemistry filter rejected the sample because phosphorus was absent, saving zero candidates. Both the generated structure and recorded MatterSim-relaxed structure are preserved in `smoke/`; their stored outputs are not DFT labels. The archive timestamp in `smoke/context.json` is not a newly inferred generation timestamp.

The later HAP44 pilot used seeds 4201–4204, Ca-P-O-H conditioning, guidance 2.0, 1000 sampling steps and at most two concurrent workers. Each generated structure has 44 atoms. The original [execution summary](hap44/execution-summary.json) records native start **2026-10-11T00:59:13.827641065Z**, finish **2026-10-11T01:05:39.622580135Z**, exit code 0 and no OOM. Host elapsed time was 387.68700000003446 seconds; cumulative worker charging was 769.8762434759992 seconds. These measure different scopes and must not be added or substituted for one another.

| Candidate | Seed | Recorded element counts | Recorded outcome |
|---|---:|---|---|
| MG44-01 | 4201 | Ca4 P3 O19 H18 | exact_composition_mismatch |
| MG44-02 | 4202 | Ca6 P4 O23 H11 | exact_composition_mismatch |
| MG44-03 | 4203 | Ca7 P4 O23 H10 | exact_composition_mismatch |
| MG44-04 | 4204 | Ca4 P4 O20 H16 | exact_composition_mismatch |

The frozen required composition was **Ca10 P6 O26 H2**. [filter_summary.json](hap44/filter_summary.json) records `qualified_for_DFT=0` and `automatic_relaxations=0`. Local-geometry and topology statuses remain `pending_review`; their metrics were not evaluated after rejection. Relaxation, hull and novelty are `not_evaluated`. These missing evaluations are not passes, failures of stability, or zero-valued measurements.

## Scope of the adaptation

[hap44-screening.json](hap44/hap44-screening.json) is the original byte-preserved contract. Its `approved_scope_execution_not_verified` status describes the frozen pre-execution contract, not the final run. Use the independent terminal and filter records for the outcome.

The contract explicitly discloses that 44 atoms is outside the model's native 1–20 atom sampling/training-size domain. The worker introduced a **process-local atom-count distribution**, retained the pretrained checkpoint and restricted allowed elements. It did not fine-tune weights or train a fixed-composition crystal-structure-prediction (CSP) model. Chemical-system conditioning does not enforce stoichiometry or apatite topology. The HAP CIF/ordered bulk named in the contract is a control/reference, not evidence that a generated candidate entered the later CP2K/DeepMD/LAMMPS chain.

Subsequent CSP discussion was an investigation; these records establish no verified CSP weights, training dataset or completed training result. This bounded negative pilot does not prove that MatterGen adaptation or CSP is impossible. No automatic extra seeds, threshold relaxation or candidate repair is implied.

All records remain development evidence with `manuscript_eligible=false`. This reporting export performed no generation, model loading, relaxation, hull calculation, DFT, training or molecular dynamics. No sealed numerical evidence was read.

See [REPRODUCE.md](REPRODUCE.md), [source-provenance.json](source-provenance.json) and [manifest.json](manifest.json) for byte-level provenance and the deliberately limited verifier.
