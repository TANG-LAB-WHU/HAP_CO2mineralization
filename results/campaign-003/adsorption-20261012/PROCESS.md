# Recorded process and negative evidence

1. The first batch refines A under the common criteria; A converges after one geometry update and two complete force evaluations.
2. During initial B optimization, a CPU-temperature telemetry process-start timeout stops the supervisor. A container-name guard also rejects the stop request; the resulting unmonitored interval and later explicit container stop are preserved. The early result snapshot is not a complete record of all subsequent native work. Later recovery accounting records 73 started / 72 complete evaluations and checkpoint step71; its orphan interval and repair time remain charged. No temperature value is invented for the gap.
3. One engineering recovery resumes B from its native checkpoint/WFN/Hessian, allowing seven updates/eight evaluations. It ends partial, stopping the original queue.
4. A separately released reference stage relaxes clean88 and gas3 under the unchanged scientific criteria. They converge after 33 and two geometry updates respectively. The three converged endpoints support only A's conditional electronic adsorption energy.
5. The explicitly bounded B continuation permits at most 20 additional updates/21 evaluations and the remaining cumulative time. It finishes these updates without convergence and stops permanently as partial. Recorded cumulative B charge is 10214.589597999933 seconds against 10800; this stop is not claimed to have consumed the entire clock allowance. No automatic continuation or extra configuration follows.

The original global allowance is 28800 seconds including startup, failure, export and cleanup. Initial per-system limits and subsequent accounting are recorded in process.json. Per-segment restarted step indices and repeated checkpoint evaluations are not distinct trajectory configurations and must not be naively concatenated or summed. `optimization-snapshots.csv` preserves absent fields as blank, never zero. The original full logs, results and native checkpoint assets remain local at the paths listed in omitted-artifacts.json.

The package preserves the earlier campaign-002 potential's failed force qualification in its separate result bundle. No new training, inference, MD, HPC work or scientific rerun is performed to produce this report.

The finalized plot package independently reconstructed the saved native tail:
initial B has 72 complete force evaluations and 71 updates, followed by 8/7
in recovery and 21/20 in continuation. This gives 101 complete evaluations,
99 unique evaluated geometries, and 98 updates. Its original step 72 stopped
during SCF and contributes no invented force point. These complete native-log
counts coexist with, and do not rewrite, the preserved earlier parser snapshot.
See `analysis-figures/data/optimization-provenance.json` for source hashes and
restart-coordinate checks. B remains unconverged after all three segments.
