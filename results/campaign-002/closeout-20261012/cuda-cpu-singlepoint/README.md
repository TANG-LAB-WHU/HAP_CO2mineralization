# Same-input HAP91 GPU/CPU single-point comparison

One sequential matched pair completed on 2026-10-11, each with one MPI rank, eight OpenMP threads, 24 GiB RAM limit, no swap and the same frozen geometry, CP2K input and initial WFN. Both converged in 73 SCF steps and exited zero without OOM. GPU grid collocation/integration and FFT execution were present in the native output.

GPU and CPU command times were 411.582818 and 569.599221 seconds: 1.383924x observed speed ratio for this pair only. Total-energy difference was 7.275958e-12 eV; maximum force-component difference was 3.162443e-12 eV/angstrom, with no correction or fitted offset. `comparison.json` records all predeclared limits, force RMSE and timing scopes.

The GPU device-wide sampled memory peak was 4222 MiB including a 1900 MiB desktop baseline. It is a sampled lower bound, not exact per-process allocation; WDDM did not expose the browser's individual allocation. This is numerical-equivalence and performance evidence, not a general speedup estimate, model qualification or an adsorption reference comparison.
