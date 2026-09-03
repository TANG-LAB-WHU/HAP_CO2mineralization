# Literature and evidence workflow

The bibliography answers “what may be cited”; the evidence ledger answers “which source supports this exact claim.”

## Source lifecycle

1. Search and record the query, source, and date in `search_log.md`.
2. Screen title and abstract for relevance.
3. Verify metadata from a publisher, DOI registry, or official proceedings.
4. Add the citation to `paper/references.bib`.
5. Add or update a row in `evidence.csv`.
6. Cite it in `paper/manuscript.qmd` only for claims it actually supports.

PDFs belong in `literature/pdfs/`, which is intentionally ignored. Do not commit copyrighted PDFs unless redistribution is clearly permitted.

## Evidence statuses

- `supported`: source/artifact checked and appropriate for the claim.
- `partial`: source supports only part of the claim.
- `planned`: method or experiment is proposed but not completed.
- `pending`: artifact exists or is expected but has not been verified.
- `rejected`: evidence is unsuitable; retain the reason to avoid repeating work.

