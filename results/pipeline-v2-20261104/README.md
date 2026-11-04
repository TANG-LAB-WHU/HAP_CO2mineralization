# Pipeline v2 active-learning snapshot — 4 November 2026 (run in progress)

This package records the measured state of the GPU-first pipeline v2 (MatterGen → CP2K → DeePMD → LAMMPS) on
Shawn-PC at the time in `state/progress.json` (2026-11-04 07:26:58 UTC). The run is **not finished**:
iterations 0–2 of active learning are complete, iteration-3 labelling is running, and the final ensemble and the
stage-7 validation gates have not run. Nothing here is a qualified potential or a mineralization result.
`manuscript_eligible=false` and `production_qualified=false` apply to every file.

Method, decisions and measured benchmarks: [`workflow/PIPELINE_V2.md`](../../workflow/PIPELINE_V2.md).
MatterGen screening (S1–S4) and its DFT hull: [`0.Generative_Discovery/runs/Ca-P-O-H-C_gfscan_RTX5090_20260329/README.md`](../../0.Generative_Discovery/runs/Ca-P-O-H-C_gfscan_RTX5090_20260329/README.md).
Every number below is recomputed from the bundled records by `scripts/summarize.py` (see [REPRODUCE.md](REPRODUCE.md)).

## Protocol

CP2K 2026.2, PBE + D3(BJ), DZVP-MOLOPT-SR-GTH, GTH-PBE pseudopotentials, 800/60 Ry, Γ point on cells with every
perpendicular width ≥ 8 Å (`protocol/protocol_report.json`). The 8 Å cells deviate from the protocol's width gate,
and the deviation is measured on 10 HAP bulk frames labelled at both widths (`finite-size/finite_size_pairs.csv`):

| 8 Å (88 atoms) minus 12 Å (352 atoms) | Value |
| --- | --- |
| Energy | 2.46 meV/atom (SD 0.24, range 2.06–2.84) |
| Force RMS difference | 0.016 eV/Å (largest component 0.047) |
| Stress RMS difference | 0.012 GPa |

Only the spread (0.24 meV/atom) adds to energy errors within HAP bulk. The offsets of the other systems are not
measured, so the bias does not cancel between systems.

## Sealed test (83 frames, 8 systems, independent seed-777 trajectories)

Atom-weighted force RMSE and frame-weighted energy RMSE (`tables/sealed_test_overall.csv`). Gate: force ≤ 0.10 eV/Å
overall and ≤ 0.12 eV/Å on the HAP(001) interfaces, energy ≤ 2 meV/atom. One of 84 sealed frames (portlandite)
failed SCF and is excluded.

| Iteration | Training | Force, all (eV/Å) | Force, interfaces (eV/Å) | Energy (meV/atom) | Gate |
| --- | --- | --- | --- | --- | --- |
| 0 | 300k steps from scratch, energy + force | 0.070–0.071 | 0.099–0.101 | 1.64–1.74 | pass, 4/4 models |
| 1 | 100k warm-start steps, energy + force | 0.070–0.072 | 0.103–0.107 | 1.58–1.82 | pass, 4/4 |
| 2 | 100k warm-start steps, virial term on | 0.066–0.067 | 0.096–0.098 | 1.33–1.67 | pass, 4/4 |

Per system (`tables/sealed_test_by_system.csv`), HAP(001) + water + CO2 is the weakest: in iteration 2 its force
RMSE is 0.104–0.107 eV/Å and its energy RMSE 2.8–3.9 meV/atom, so it fails the energy gate on its own. In iterations
0–1 the virial term was silently off (`workflow/PIPELINE_V2.md`, "Virial loss was off"), so those models were never
fitted to stress.

## Active learning

DP-GEN model deviation, trust window 0.10–0.30 eV/Å, at most 30 candidates per run, LAMMPS NVT at 300 and 500 K with
0.5 fs steps (`tables/exploration_runs.csv`, `tables/labels_by_system.csv`):

| Iteration | CP2K labels | Exploration per run | Mean accurate fraction | Candidates |
| --- | --- | --- | --- | --- |
| 0 | 216/216 (28 CaO frames relabelled after a cell-format fix) | 20k steps (10 ps) | 0.673 | 153 |
| 1 | 153/153 | 40k steps (20 ps) | 0.723 | 253 |
| 2 | 251/253 (2 SCF failures) | 60k steps (30 ps) | 0.838 | 208 |
| 3 | 187/208 at the snapshot, 1 SCF failure so far | not started | — | — |

The exploration length grows with the iteration, so the accurate fractions are not strictly comparable. The
convergence gate (≥ 99 % of explored frames below 0.10 eV/Å) is far from met. HAP(001) + water + CO2 is the gap:

| Run | Accurate fraction, iteration 0 → 1 → 2 | Largest deviation in iteration 2 (eV/Å) |
| --- | --- | --- |
| HAP(001) + water + CO2, 300 K | 0.000 → 0.013 → 0.008 | 0.98 |
| HAP(001) + water + CO2, 500 K | 0.000 → 0.004 → 0.106 | 1.85 |
| HAP(001) + CO2, 500 K | 0.042 → 0.057 → 0.883 | 0.49 |
| monetite, 500 K | 0.005 → 0.069 → 0.506 | 0.59 |

No CO2 conversion to bicarbonate or carbonate was observed in any unbiased exploration. In the 300 K
water + CO2 run of iteration 2, every 3000th frame still has both CO2 molecules and all 25 H2O intact
(`workflow/PIPELINE_V2.md`).

The 10 MatterGen survivors of S3 are labelled (90 training and 40 sealed-test frames, all converged;
`state/progress.json`) and join training at iteration 3.

## Negative and partial results kept in this package

- `iter-000/explore_INVALID_lammps_skin2_summary.json`: the first iteration-0 exploration, invalid because a LAMMPS
  neighbour skin of 2 Å gave wrong DeePMD forces (every frame above trust_hi). It is rerun with skin 0, and the rerun
  is the valid `iter-000/explore_summary.json`.
- Iteration 2, CaCO3_20 at 300 K: LAMMPS segmentation fault at step 50 600 of 60 000, recorded with
  `lammps_complete: false`; the frames before the crash are analysed.
- SCF failures after 1050 steps, kept and excluded: `candidates_00056` (HAP(001) + CO2) and `candidates_00121`
  (HAP(001) + water + CO2) in iteration 2, `candidates_00127` (HAP(001) + water + CO2, 500 K) in iteration 3, and
  one sealed-test portlandite frame. Three of the four are interface frames, so the hardest interface configurations
  are slightly under-sampled.

## Not included

Training and sealed-test labels, DeePMD data sets, models, checkpoints, LAMMPS trajectories and CP2K outputs stay in
`.research-runs/pipeline-v2/` on Shawn-PC. `manifest.json` records their sizes and SHA-256 (tree hashes for data
directories) under `excluded_assets`, so a later copy can be checked against this snapshot. Benchmarks described
in `workflow/PIPELINE_V2.md` (GPU/CPU CP2K timing, TensorFlow versus PyTorch inference, allocator benchmark) are
documented there and not packaged.

## Interpretation boundaries

1. The sealed-test gate is one of four gates. Lattice, stability (≥ 500 ps) and active-learning convergence are
   open, and they apply only to the final ensemble, which has not been trained.
2. Sealed-test RMSEs are measured against 8 Å Γ labels; the Γ bias above is not included in them.
3. Unbiased exploration has sampled no reactive event. A mechanism claim needs biased sampling, which is planned
   but not implemented.
