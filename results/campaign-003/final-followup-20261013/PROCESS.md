# Preserved completion and failure history

The sequence continues after the earlier B-continuation1 partial endpoint. B-direction1 completed five bounded updates without convergence; the subsequent B-adaptive1 run converged in its second segment. The final B endpoint is in `endpoints/B91-final.extxyz`. The B-derived clean optimization converged in cleanB1 segment04; its endpoint is `endpoints/clean88-B-derived.extxyz`. Intermediate segments remain explicitly partial in `history/` and `data/point-results.json`.

The first A-tilt attempt failed before obtaining an accepted point. The repair retained the original tilt deadline and completed the two singlepoints; its scientific classification remained inconclusive. Both initial and repair driver records are retained.

The first fragment attempt stopped on a missing/stale dedicated CPU thermal observation. The repair completed all four native singlepoints. Its driver then stopped with `Native final energy does not match selected energy` during postprocessing. This terminal state is intentionally retained in `history/fragments1-repair1-run.json`. The later existing `fragment-decomposition-postprocessed.json` identifies the energy-block parser repair, preserves the old analysis hashes, and records `new_SCF_or_force_calls=0`. It does not rewrite the failed driver state or claim a fresh successful calculation.

The CP run completed four ghost-fragment points after the new clean reference became available. Its summary names the original and B-derived clean references separately. Ghost-center derivatives do not use the physical real-atom net-force acceptance gate and are not training labels.

Charged driver walltime and native command time remain separate, as recorded by the original jobs. Phase totals can contain prior umbrella costs and must not simply be summed. In particular, a missing duration in the failed initial fragment attempt is unknown, not zero. No new resource measurements were taken for this publication.
