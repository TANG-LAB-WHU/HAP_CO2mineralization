# Figure plan and provenance

For every final figure, record the generating script, input artifacts and run IDs, source checksums, Git commit, units, and every manual post-processing step. A figure is eligible for the manuscript only when all referenced artifacts are marked `manuscript_eligible=yes` in `research/artifact_inventory.csv` and all result claims are supported in `literature/evidence.csv`.

| Figure | Planned content | Evidence status | Artifact/run IDs | Source checksum(s) | Generator | Manual edits | Eligible |
|---|---|---|---|---|---|---|---|
| 1 | End-to-end multiscale workflow | schematic pending | n/a | n/a | pending | none | no |
| 2 | Generated candidate and filtering landscape | TODO-EVIDENCE R1 | pending | pending | pending | none | no |
| 3 | Composition-by-facet adsorption ranking | TODO-EVIDENCE R2 | pending | pending | pending | none | no |
| 4 | MLIP validation: energy/force parity and learning curves | TODO-EVIDENCE R3 | pending | pending | pending | none | no |
| 5 | Hydrated-interface structure and dynamics | TODO-EVIDENCE R4 | pending | pending | pending | none | no |
| 6 | Free-energy profile and proposed mechanism | TODO-EVIDENCE R4 | pending | pending | pending | none | no |

Prefer scripts under `5.Figure_Publication/`; copy only final publication assets here when their provenance is documented.

## Legacy figure exclusion

Figures under `4.PostProcess_Analysis/HAP112_DPA_NVTcsvr_Succeeded/` are preserved as legacy artifacts but are not publication figures. They must not be copied here or cited as evidence because the source data are absent, the recorded quality assessment is poor, and the historical RDF parser could substitute synthetic RDF-like data for invalid input. See artifact `POST-HAP112` and evidence record `X001`.
