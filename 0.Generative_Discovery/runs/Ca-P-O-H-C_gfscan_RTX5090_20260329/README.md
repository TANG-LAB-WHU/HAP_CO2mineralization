# Ca-P-O-H-C guidance-factor scan (RTX 5090, siqi, 2026-03 → 2026-10)

Imported on 2026-10-28 from the group Google Drive folder `results20260329_RTX5090`.
Raw outputs (18.1 GB, 45,174 files) are kept out of Git at
`0.Generative_Discovery/results/Ca-P-O-H-C_gfscan_RTX5090_20260329/`; this directory
holds only the tracked summary. Regenerate it with:

```bash
python -I 0.Generative_Discovery/summarize_gf_scan.py \
  --scan-dir 0.Generative_Discovery/results/Ca-P-O-H-C_gfscan_RTX5090_20260329 \
  --out-dir  0.Generative_Discovery/runs/Ca-P-O-H-C_gfscan_RTX5090_20260329
```

## Configuration (from `legacy_siqi_jobscripts/CommandLine.txt`)

| Item | Value |
| --- | --- |
| Script | `legacy_siqi_jobscripts/predict_ca_p_o_h_c.py` → `mattergen-generate` + `mattergen-evaluate` |
| Model | pretrained `chemical_system`, condition `{'chemical_system': 'Ca-P-O-H-C'}` |
| Sampling | batch 64 × 50 batches = 3,200 per run; guidance factor ∈ {0.001, 0.5, 1, 1.5, 2} |
| Relaxation / energy | MatterSim-v1.0.0-5M |
| Hull reference | `reference_MP2020correction.gz` (Alex-MP), `structure_matcher=disordered` |
| Classification | `legacy_siqi_jobscripts/extracting_structure.py` (identical in all five runs) |
| Not recorded | MatterGen commit, seeds, MatterSim package version, wall time |

## Per-run results

| gf | stable (E_hull<0.1) | S.U.N. | avg E_hull (eV/atom) | all 5 elements | Ca+CO3+PO4 in S.U.N. |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.001 | 1,593 (49.8 %) | 727 | 0.160 | 4 | 2 |
| 0.5 | 1,036 (32.4 %) | 704 | 0.212 | 29 | 10 |
| 1.0 | 661 (20.7 %) | 521 | 0.275 | 180 | 28 |
| 1.5 | 429 (13.4 %) | 368 | 0.334 | 544 | 54 |
| 2.0 | 277 (8.7 %) | 240 | 0.372 | 1,066 | 56 |

All structures have 1–20 atoms (median 9–10), as expected for the Alex-MP-20 training domain.

**Trade-off:** raising the guidance factor makes MatterGen respect the full five-element
system (4 → 1,066 of 3,200) but lowers stability (49.8 % → 8.7 %). Low gf mostly yields
Ca-O-H or binary/ternary subsystems. The yield of DAC-relevant Ca–CO3–PO4 phases peaks at gf 1.5–2.

## Cross-run S.U.N. set

- 2,560 S.U.N. structures → **2,324 unique** after StructureMatcher dedup within each reduced formula.
- Classes: Ca-O-H only 480, other P/C 931, Ca carbonate 354, Ca phosphate 412, **Ca carbonate-phosphate 147**.
- **52 priority candidates** (`priority_Ca_CO3_PO4_shortlist.csv`): Ca + isolated CO3 + isolated PO4,
  every C in CO3 and every P in PO4, no shared O, formal-charge neutral (Ca²⁺ P⁵⁺ C⁴⁺ H⁺ O²⁻).
  Most frequent: Ca2(HPO4)(CO3) ×16 polymorphs, CaPH5CO8 ×9, Ca2P2H4CO11 ×7.
  Chemically notable: Ca4(PO4)2CO3 (carbonate-apatite-like stoichiometry) and Ca3PHCO8
  (consistent with Ca3(PO4)(CO3)(OH), an apatite-like anion set).

## Data-quality flags

1. `gf1/metrics.json` reports `avg_rmsd_from_relaxation = 0.00056 Å`, but recomputing
   generated → relaxed gives a median of 0.108 Å, in line with the other runs (0.06–0.13 Å).
   The gf1 metrics file may come from a different evaluation pass. gf1 files are dated 2026-04-21,
   its log 2026-03-28, and its trajectory archive is named `generated_trajectories-002.zip`.
2. `gf2/generated_trajectories.zip` is **0 bytes**, so its denoising trajectories are missing.
   The gf1.5 run was regenerated on 2026-10-17.
3. Extraction logs exist only as `.txt` for gf0.5/1.5/2. gf0.001 has a Drive-exported PDF and gf1 a `.docx`.
4. Up to 16 "novel" polymorphs per formula within 0.1 eV/atom points to a flat MatterSim landscape
   or a loose threshold, not to real polymorphism. Expect heavy attrition under DFT.

## Screening progress (pipeline v2, 2026-10-29)

Scripts: `0.Generative_Discovery/screening/`. Tracked results: `screening/screening_summary.json`, `screening/s2_hull.csv`.

| Gate | Rule | In → out |
| --- | --- | --- |
| S1 chemistry | Ca carbonate/phosphate classes, formal-charge neutral, every C in an isolated CO3 and every P in an isolated PO4, no shared O, every H on O | 2,324 → **363** |
| S2 hull committee | MACE-MPA-0, MatterSim-5M and SevenNet-MF-OMPA each relax the candidates and 288 Alex-MP reference phases (DFT E_hull ≤ 0.1 eV/atom); E_hull ≤ 50 meV/atom in ≥ 2 models | 363 → **149** (131 unanimous; 70 carbonates, 59 phosphates, 20 carbonate-phosphates) |
| S3 dynamics | MACE phonons (min frequency ≥ −0.3 THz on an 8³ mesh) and 5 ps 300 K NPT with CO3/PO4 intact | running |

The three models agree closely: the median pairwise |ΔE_hull| is 2.5–3.3 meV/atom. That agreement
does not replace DFT, because all three are trained on similar PBE data. The DFT hull with the v2 CP2K
protocol (`workflow/PIPELINE_V2.md`) remains the next gate for S3 survivors.

## Status and next gates

`legacy_unverified` (external import): no result here is manuscript-eligible (see `research/decision_log.md`, 2026-10-28).

1. Re-relax the 52 priority structures with a second MLIP (MACE-MP-0 or SevenNet) and keep only
   those with E_hull < 50 meV/atom under both.
2. Rebuild the hull with **the same DFT setup used downstream** (CP2K PBE-D3 or VASP r2SCAN) for the
   survivors plus known phases (calcite, aragonite, HAP, brushite, monetite, portlandite).
3. For survivors, compute ΔG of HAP + CO2 (+H2O) → phase at 298 K / 420 ppm CO2 / RH 50 %.
4. Only then feed selected phases into slab generation (`1.HT_Screening/02_slab_generation`).
