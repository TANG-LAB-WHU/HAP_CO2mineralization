# Manuscript workspace

- `manuscript.qmd`: living paper draft.
- `references.bib`: verified bibliography used by Quarto and automatically exported by Better BibTeX from `HAP_CO2mineralization/10_References_Verified`.
- `figures/`: publication figures and provenance notes.
- `_output/`: generated HTML/PDF output; ignored by Git.

Build from the repository root:

```bash
make check
make paper
make preview
make paper-pdf
```

Do not remove `TODO-EVIDENCE` markers until the corresponding row in `literature/evidence.csv` identifies a verified publication or computational artifact.

Zotero 10 and Better BibTeX 9 use the item's native `Citation Key` field for registry-aligned keys; the older pin workflow is not required. Do not hand-edit the bibliography during normal operation. The auto-export omits `file,attachment,note,annote,abstract,keywords`; the regenerated 2026-09-05 export passed the repository's local-path and attachment-field checks.
