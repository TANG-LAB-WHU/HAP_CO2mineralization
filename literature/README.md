# Literature and evidence workflow

The reference registry answers “has this publication's bibliographic identity been reviewed?” The evidence ledger answers “does this source support this exact claim?” The manuscript bibliography contains only references that have passed human bibliographic review.

## Files and row semantics

- `references.csv`: exactly one row per publication. `reference_status` is one of `candidate`, `metadata_partial`, `verified`, or `rejected`.
- `evidence.csv`: exactly one row per claim-source relationship. `evidence_id` is unique; `claim_id` may repeat when several sources are assessed for one claim.
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

PDFs belong in Zotero-managed storage outside the repository. `literature/pdfs/`, `literature/attachments/`, `literature/zotero/`, `.zotero/`, and `zotero.sqlite*` are defense-in-depth ignore targets and must never be tracked.
