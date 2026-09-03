# Manuscript workspace

- `manuscript.qmd`: living paper draft.
- `references.bib`: verified bibliography used by Quarto.
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

