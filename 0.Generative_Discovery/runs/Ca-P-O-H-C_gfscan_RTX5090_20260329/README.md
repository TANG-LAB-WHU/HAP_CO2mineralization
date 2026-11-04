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
| S3 dynamics | MACE phonons (min frequency ≥ −0.3 THz on an 8³ mesh) and 5 ps 300 K NPT with CO3/PO4 intact; calcite, portlandite, monetite and CaO controls pass | 149 → **10** (`screening/s3_dynamics.csv`) |
| S4 DFT hull | CP2K PBE-D3(BJ) cell optimisation (v2 labelling protocol) of the survivors and 34 reference phases; E_hull ≤ 50 meV/atom against the computed references | 10 → **8** (`screening/s4_dft_hull.csv`) |

The three models agree closely: the median pairwise |ΔE_hull| is 2.5–3.3 meV/atom. That agreement
does not replace DFT, because all three are trained on similar PBE data.

S3 survivors: CaCO3 ×2 (10 and 20 atoms), Ca2H2CO5, Ca2H6CO7, Ca2PH5O7, CaP2(H2O3)3, Ca3P2H2O9,
CaP2(HO2)4, CaPH3O5, and one carbonate-phosphate, Ca2P2H4CO11 (formally Ca2(H2PO4)2CO3).
Six of the ten have a lowest phonon frequency between −0.30 and −0.05 THz, close to the threshold.
The median over all 149 tested structures is −0.96 THz.

S4 (finished 2026-11-02; method and references in `workflow/PIPELINE_V2.md`): all 44 cell optimisations
converged and every final single point converged. DFT E_hull in meV/atom (with the Γ finite-size correction
in brackets):

| Candidate | Formula | DFT E_hull | S2 MLIP E_hull (MACE / MatterSim / SevenNet) | Lowest decomposition |
| --- | --- | ---: | --- | --- |
| s3_055 | CaCO3 (10 atoms) | 15.8 (14.9) | 14.1 / 9.9 / 21.4 | calcite |
| s3_273 | CaP2(H2O3)3 | 26.3 (26.1) | 21.9 / 20.6 / 21.9 | Ca(H2PO4)2·H2O |
| s3_014 | Ca2H2CO5 | 32.8 (33.8) | 44.3 / 45.6 / 46.4 | calcite + portlandite |
| s3_070 | CaCO3 (20 atoms) | 36.1 (37.4) | 40.3 / 43.3 / 47.3 | calcite |
| s3_039 | Ca2H6CO7 | 39.7 (39.8) | 42.1 / 44.3 / 44.9 | calcite + portlandite + ice |
| s3_294 | CaP2(HO2)4 | 41.7 (41.3) | 41.9 / 38.7 / 40.8 | Ca(H2PO4)2 |
| s3_302 | CaPH3O5 | 42.5 (42.2) | 35.8 / 31.3 / 34.7 | monetite + ice |
| s3_250 | Ca2PH5O7 | 47.0 (46.9) | 35.6 / 35.0 / 39.2 | HAP + portlandite + ice |
| s3_262 | Ca3P2H2O9 | **60.9** (60.9) | 42.1 / 38.1 / 42.2 | HAP + monetite + ice |
| s3_180 | Ca2P2H4CO11 | **68.0** (67.2) | 47.3 / 46.9 / 48.7 | monetite + dry ice + ice |

No survivor is on the DFT hull. The only carbonate-phosphate (s3_180) fails, so no Ca carbonate-phosphate
from this scan survives to S4. DFT minus the median MLIP E_hull is +4.1 meV/atom on average (−12.8 to
+20.7). The MLIP hulls used Alex-MP references relaxed with each model, so this compares two different
reference sets, not only two energy models. The two CaCO3 survivors lie 15.8 and 36.1 meV/atom above calcite,
against 6.4 and 20.3 for the two Pnma CaCO3 references (aragonite candidates) in the same hull.
The finite-size correction moves no survivor's E_hull by more than 1.3 meV/atom. Hexagonal HAP (mp-41472, the
pipeline's HAP bulk) is 2.6 meV/atom above mp-721624.

## Status and next gates

`legacy_unverified` (external import): no result here is manuscript-eligible (see `research/decision_log.md`, 2026-10-28).

1. Re-relax the 52 priority structures with a second MLIP (MACE-MP-0 or SevenNet) and keep only
   those with E_hull < 50 meV/atom under both.
2. Rebuild the hull with **the same DFT setup used downstream** (CP2K PBE-D3 or VASP r2SCAN) for the
   survivors plus known phases (calcite, aragonite, HAP, brushite, monetite, portlandite).
   Done with CP2K PBE-D3 for the S3 survivors (S4 above; the CaHPO4·2H2O composition is ref-275 in
   `screening/s4_dft_hull.csv`).
3. For survivors, compute ΔG of HAP + CO2 (+H2O) → phase at 298 K / 420 ppm CO2 / RH 50 %.
4. Only then feed selected phases into slab generation (`1.HT_Screening/02_slab_generation`).
