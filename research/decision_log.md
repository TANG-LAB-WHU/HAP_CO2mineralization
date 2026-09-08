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
