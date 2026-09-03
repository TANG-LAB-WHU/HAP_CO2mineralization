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
