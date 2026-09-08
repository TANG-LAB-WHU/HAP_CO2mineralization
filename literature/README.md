# Literature and evidence workflow

The reference registry answers “has this publication's bibliographic identity been reviewed?” The evidence ledger answers “does this source support this exact claim?” The manuscript bibliography contains only references that have passed human bibliographic review.

## Files and row semantics

- `references.csv`: exactly one row per publication. `reference_status` is one of `candidate`, `metadata_partial`, `verified`, or `rejected`.
- `relevance.csv`: exactly one row per registry publication. It ranks project relevance and records actual repository-use maturity without copying bibliographic or evidence status.
- `evidence.csv`: exactly one row per claim-source relationship. `evidence_id` is unique; `claim_id` may repeat when several sources are assessed for one claim.
- `notes/<citation_key>.md`: a machine evidence card with legal-full-text provenance, checksum, candidate locator, scoped paraphrase, and an explicit human-review gate.
- `search_log.md`: exact databases/services, queries, dates, screening decisions, and DOI/title de-duplication results.
- `paper/references.bib`: the only bibliography available to Quarto. It must contain exactly the `verified` reference keys and no attachment paths.

## Reference lifecycle

1. Discover a source for an explicit outline evidence gap and add it as `candidate`.
2. Check title, complete ordered authors, issued year, venue, DOI status, and canonical publisher URL against Crossref plus a publisher or official-proceedings page.
3. An automated check may move the source only to `metadata_partial`. Record checked URLs, discrepancies, and the human fields still requiring confirmation; leave `verified_by` and `verified_at` blank.
4. A human reviewer compares the record with the authoritative pages. Only explicit approval may set `reference_status=verified`, `verified_by`, and an ISO `verified_at` date.
5. Read the legal full text and record a claim-specific page/section/paragraph/figure/table/equation locator in `evidence.csv`. Title or abstract screening alone cannot create supporting evidence.
6. Set a relationship to `supported` only when the exact claim is supported. Use `partial` with narrowed manuscript wording, keep inaccessible or unresolved material `pending`, and retain rejected sources with a reason.

The trace used by the checker is:

```text
<!-- CLAIM: Cxxx --> + [@citation_key]
  -> evidence.csv (claim_id, citation_key, supported|partial, source_locator)
  -> references.csv (citation_key, verified)
  -> paper/references.bib
```

Pending manuscript material must use `TODO-EVIDENCE <claim-id>` or `PLANNED`; candidate and metadata-partial sources are not active citations.

The checker requires a non-empty BibTeX `author` field but does not yet compare creator lists with the registry. Reliable comparison is deferred because consortium authors, LaTeX accent encoding, name particles, and ordered full creator lists require structured parsing; do not replace this with a lossy string comparison.

## Full-text evidence cards

Phase 2B keeps bibliographic identity, project relevance, machine extraction, human confirmation, ledger promotion, and citation activation as separate gates:

The four-source first wave validates this workflow only. It does not claim complete coverage of HAP–CO2 hydrated-interface literature.

```text
references.csv metadata status
  -> relevance.csv selection
  -> notes/<citation_key>.md machine extraction
  -> Shawn confirms the exact PDF SHA and locator
  -> separate evidence.csv promotion decision
  -> separate manuscript citation activation
```

Downloads fail closed. The final response URL must remain HTTPS; the response must be successful and must not be a login, error, or challenge page; bytes must start with `%PDF-`; `pdfinfo` must report at least one page; and `pdftotext` must yield non-empty text whose title and author identity match the registry. A DOI printed in the PDF is checked directly. When an accepted manuscript or preprint omits the journal DOI, the card must instead name the official landing metadata that associates that exact repository record with the DOI. Record SHA-256 only after all checks pass, and delete invalid bytes rather than converting HTML or XML into a PDF.

All Codex-generated Phase 2B-1 cards use `human_confirmed: false`. Their candidate locators are not copied into `evidence.csv`, and relationships remain `pending` until Shawn reviews the same checksum and explicitly authorizes a later reconciliation. The combined verbatim excerpt allowance is at most 25 words per card and excludes title, abstract, and search-snippet text.

`make check` validates source-controlled records and works in a clean clone without local PDFs. `make evidence-local-check` additionally requires Poppler and validates ignored PDFs through no-follow directory/file descriptors. It hashes the bytes copied from the opened descriptor, runs every Poppler command against that private stable snapshot, renders locator pages in the same atomically created system-temporary directory, and removes the directory on success or failure.

## Search boundary

Search only gaps present in the Introduction, Background, or Methods outline. Prefer 2–4 directly relevant primary sources or authoritative reviews per gap. Record every query and screening decision, de-duplicate by DOI first and normalized title plus first author/year second, and obey the active round's candidate cap. Do not add adjacent literature merely to fill a quota.

## Zotero and Better BibTeX configuration

Zotero's database and managed attachments must remain outside this repository. Use the default macOS data directory or another external location; never point Zotero at `literature/` or `paper/`.

Create these collections in `My Library`:

```text
HAP_CO2mineralization/
├── 00_Review_Queue
├── 10_References_Verified
└── 90_Rejected
```

With Zotero 10 and Better BibTeX 9:

1. Add/import registry items into `00_Review_Queue` and enter the exact registry key in Zotero's native `Citation Key` field. Zotero 10 plus Better BibTeX 9 does not require the legacy Better BibTeX pin action; the native field is the project authority.
2. Configure new keys as `auth.lower + year + shorttitle(1,1).lower`, plain text, unique within the library, with automatic key regeneration disabled. Explicit native `Citation Key` values must remain unchanged.
3. Move only human-verified items into `10_References_Verified`.
4. Export that collection as `Better BibTeX` to the absolute local path of `paper/references.bib`; select `Keep updated`.
5. Use on-change export with a 5-second delay, citation-key sorting, both DOI and URL, no journal abbreviation, no files, no notes, and no JabRef metadata. Omit `file,attachment,note,annote,abstract,keywords`.
6. Review every generated diff before commit and run `make check` plus `make paper`.

The Zotero GUI gate was completed on 2026-09-04 with Zotero 10.0.1 and Better BibTeX 9.0.63: the three collections exist, the eight seed items use native `Citation Key` values, and a real automatic export was observed. The export setting `Fields to omit from export` is `file,attachment,note,annote,abstract,keywords`. A regenerated export was validated on 2026-09-05 with no forbidden fields, local paths, or Zotero attachment identifiers. Better BibTeX Git integration remains disabled.

Long-term PDFs belong in Zotero-managed storage outside the repository. During an authorized evidence-extraction pass, local working copies may be kept only in ignored `literature/pdfs/`. Validator snapshots and renders use private system-temporary directories and must leave no residue; ignored `tmp/pdfs/` is only a manual scratch location. `literature/pdfs/`, `tmp/pdfs/`, `literature/attachments/`, `literature/zotero/`, `.zotero/`, and `zotero.sqlite*` must never be tracked.
