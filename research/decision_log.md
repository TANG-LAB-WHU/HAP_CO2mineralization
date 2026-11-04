# Research decision log

Record decisions that change scope, assumptions, data, models, or interpretation. Do not rewrite old entries; append a superseding decision.

## 2026-09-03 — local-first research workflow

- **Decision:** Begin with a macOS-local literature and manuscript workflow on `shawn_dev`.
- **Reason:** The immediate bottleneck is structuring claims, evidence, and required results rather than launching more computation.
- **Boundary:** HPC is represented only by a disabled profile template and run-manifest contract.
- **Revisit when:** The team approves the research questions, evidence gaps, and validation gates for the next compute campaign.

## 2026-09-03 — legacy artifacts require eligibility review

- **Decision:** Treat every existing computational artifact as manuscript-ineligible until its provenance, completeness, and scientific quality are explicitly reviewed in `research/artifact_inventory.csv`.
- **Reason:** The repository contains incomplete and failed optimizations, short or unstable AIMD trajectories, unresolved conflict markers, missing raw postprocessing inputs, and legacy code that could generate synthetic RDF-like fallback data.
- **Consequence:** Existing files are preserved as historical evidence, including negative results, but none may support a manuscript finding at this stage.
- **Boundary:** This decision does not repair, delete, move, rerun, or reinterpret any historical calculation.
- **Revisit when:** A result has a run manifest, checksums, explicit units and acceptance criteria, and an independent eligibility review.

## 2026-09-03 — separate reference verification from claim evidence

- **Decision:** Use `literature/references.csv` for one-row-per-publication bibliographic status and `literature/evidence.csv` for one-row-per-claim-source relationship.
- **Reason:** A verified publication identity does not establish that the publication supports a particular manuscript claim, and one claim may require several sources.
- **Consequence:** `evidence_id` is unique while `claim_id` may repeat; active manuscript citations require both a human-verified reference and a supported or partial relationship with a full-text locator.
- **Automation boundary:** Machine metadata checks stop at `metadata_partial` and never populate a human reviewer identity.

## 2026-09-03 — Phase 2A bibliography and Zotero gate

- **Decision:** Keep `paper/references.bib` free of candidate and metadata-partial entries. Configure Better BibTeX only through Zotero's GUI and do not imitate its output before a real automatic export is observed.
- **Local data boundary:** Zotero databases, profiles, managed attachments, PDFs, and XPI installers remain outside Git.
- **Current state:** Eight seeds are machine-checked but await human review; eight additional candidates fill the authorized first-round cap. Zotero 10.0.1 is installed, and the official Better BibTeX 9.0.63 XPI was downloaded with SHA-256 `3a4d080ec94153a8c24bf8275689c5f946d99e314b6ae0fab5060d6970f56388`; GUI installation, collections, key pinning, and auto-export remain `blocked_by_manual_zotero_setup`.

## 2026-09-05 — Phase 2A human reconciliation

- **Human confirmation:** Shawn Zhang explicitly approved all eight seed citation keys with verification date 2026-09-04; the registry records only those eight as `verified`.
- **Evidence boundary:** Bibliographic verification does not change claim-source status. All literature relationships remain `pending` because no legal full-text locator was supplied.
- **Zotero/BBT state:** The Zotero 10.0.1 and Better BibTeX 9.0.63 GUI gate, collections, native `Citation Key` fields, and automatic export were confirmed. Zotero 10 plus BBT 9 uses the native field, so the legacy pin action is no longer required.
- **Export finding:** The first automatic export contained six local Zotero attachment paths. They were removed from the repository copy; the BBT field-omission setting must be corrected before the next export. Better BibTeX Git integration remains disabled.
- **Final export validation:** The user set `Fields to omit from export` to `file,attachment,note,annote,abstract,keywords` and reran the verified-collection auto-export. The regenerated file was checked on 2026-09-05 and contains no forbidden fields, local paths, or Zotero attachment identifiers; the earlier export blocker is resolved.

## 2026-09-05 — Phase 2B-1 machine extraction remains below the human gate

- **Decision:** Rank all 16 registry references by project relevance and use at most four sources to validate the full-text evidence-card workflow.
- **First wave:** Qomi covers general mineralization/interfacial-water background, Astala is a direct HAP surface/hydration candidate, Nowicki is only a 500 °C scope boundary, and Kühne supports an actually used CP2K method capability.
- **Acquisition:** Three official PDFs passed fail-closed validation; the RSC PDF request for Nowicki returned HTTP 403 and remains `blocked_fulltext`. No third-party mirror or document conversion is permitted.
- **Evidence boundary:** Every card remains machine-only and `human_confirmed: false`. The ledger retains its Phase 2A statuses and empty source locators; candidate metadata, machine locator extraction, human card confirmation, ledger promotion, and manuscript citation activation remain separate decisions.
- **Compute boundary:** Legacy CP2K execution exists, but no legacy calculation is manuscript-eligible. The CP2K paper supports software capability only, not validity of repository results.

## 2026-10-28 — import the external MatterGen Ca-P-O-H-C guidance-factor scan

- **Decision:** Import siqi's five-run Ca-P-O-H-C scan (gf 0.001, 0.5, 1, 1.5, 2; 3,200 structures each; MatterSim-5M relaxation) as `legacy_unverified` (external import). Raw outputs live untracked in `0.Generative_Discovery/results/`; the tracked summary, SHA-256 manifest, and shortlist live in `0.Generative_Discovery/runs/Ca-P-O-H-C_gfscan_RTX5090_20260329/`.
- **Use:** The scan completes the small-cell Ca-P-O-H-C phase space and convex hull for judging carbonation products. It is not a route to generating HAP itself: every structure has ≤20 atoms, and the HAP cell has 44.
- **Finding:** Higher guidance enforces the five-element system (4 → 1,066 per run) at the cost of stability (49.8 % → 8.7 %). Of 2,324 unique S.U.N. structures, 52 are charge-neutral Ca phases with isolated CO3 and PO4.
- **Boundary:** Stability labels are MatterSim energies against the MP2020 reference with a 0.1 eV/atom window. Every structure stays manuscript-ineligible until a second MLIP and a DFT hull at the downstream level of theory confirm it.
- **Revisit when:** Provenance gaps (MatterGen commit, seeds, gf1 metrics, empty gf2 trajectory archive) are answered, and DFT hull results exist for the shortlist.

## 2026-10-29 — pipeline v2: CPU CP2K labels, GPU MLIP stages, fixed labelling protocol

- **Decision:** Rebuild the CP2K → DeePMD → LAMMPS pipeline as described in `workflow/PIPELINE_V2.md`. The user asked for GPU acceleration and literature-level accuracy targets rather than maximal accuracy.
- **CP2K on GPU (measured):** Upstream publishes CUDA images only up to 2025.1 and only for P100. GitHub releases ship source and a CPU-only static `.ssmp`. A CP2K 2026.2 sm_89 image was built (`workflow/environment/Dockerfile.cp2k-cuda-sm89`; DBCSR needed an A40→sm_89 patch). On HAP bulk single points it took 251 s (88 atoms) and 902 s (352 atoms). The official CPU image with 24 MPI ranks took 97 s and 476 s, with identical energies. CP2K labelling therefore runs CPU-only, and the RTX 4090 runs the MLIP stages concurrently.
- **Protocol (measured):** PBE-D3(BJ)/DZVP-MOLOPT-SR-GTH/GTH-PBE, Γ point, 800/60 Ry. 800 Ry is the lowest tested cutoff within 1 meV/atom and 0.01 eV/Å of 1000 Ry; 600 Ry gave 0.0108 eV/Å.
- **Recorded deviation:** The 88-atom (8 Å) Γ cell is 2.67 meV/atom above the 528-atom result and fails the 1 meV/atom width gate. The 352-atom cell passes (0.04 meV/atom). Labels nevertheless use ≥ 8 Å cells, because 12 Å interface cells would multiply single-workstation labelling time by roughly 3–4×. Forces differ by only 0.017 eV/Å RMS. The energy bias is measured on paired 8/12 Å labels and reported with every energy gate.
- **Accuracy targets:** force RMSE ≤ 0.10 eV/Å, energy RMSE ≤ 2 meV/atom on a sealed test set from an independent exploration seed, plus lattice and stability gates. These follow the SCAN-ML CaCO3 Deep Potential (0.4 meV/atom, 101 meV/Å); that reference is not yet entered in `literature/references.csv`.
- **Boundary:** No v1 labels are reused. Nothing from v2 is manuscript-eligible until a gate report and run manifest are reviewed.

## 2026-11-04 — merge pipeline v2 into shawn_dev and publish a run snapshot

- **Decision:** Shawn asked for shawn_dev to carry the current workflow, README and results. `integration/shawn-4090-reproducible` (local v1 validation, MatterGen S1–S4 screening, pipeline-v2 code and records) is merged with the 17 October release already on shawn_dev, together with the uncommitted allocator fix, documentation and the 2026-11-03 novelty-check registry rows. The running orchestrator keeps its own checkout; the merge was done in a separate worktree.
- **Snapshot:** `results/pipeline-v2-20261104` is exported by `workflow/export_pipeline_v2_results.py` from the run directory. It holds stage records, training inputs, learning curves, sealed-test and exploration summaries, and exploration start structures, and `scripts/summarize.py --check` recomputes every table from them. Labels, data sets, models and trajectories are excluded and identified by SHA-256.
- **Inventory vocabulary:** Seven rows published on shawn_dev on 2026-10-17 used status and provenance labels that `make check` rejects (14 errors on the published branch). They are mapped to the checker vocabulary: all-candidates-rejected, failed-qualification, failed-learned-dynamics and unqualified become `failed`; development-only and conditional become `partial`; training-thresholds-met becomes `completed`; the hash-based provenance labels become `path_only`. Each row keeps its original label in its notes, and no finding changes. The two `P0-CAMPAIGN` rows are split: the bounded generation records keep the ID, and the imported scan becomes `P0-GFSCAN-20260329`.
- **Corrections:** The v1 force RMSE of 0.108 eV/Å is the retrospective test RMSE on 8 frames (32 frames in all; `workflow/RECOVERY.md`), not an RMSE on 32 frames. The SCAN-ML CaCO3 reference is now the unverified candidate `piaggi2025caco3`.
- **Boundary:** The snapshot is of an unfinished run. Only the sealed-test gate has been evaluated, on intermediate ensembles; nothing is manuscript-eligible. The biased-exploration plan drafted on 2026-11-04 is not included until it is reviewed.
