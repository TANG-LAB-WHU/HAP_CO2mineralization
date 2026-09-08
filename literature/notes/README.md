# Evidence cards

Each `<citation_key>.md` file contains one machine-readable JSON evidence card followed by optional reviewer notes. Cards record full-text provenance and candidate claim locators; they do not change `literature/evidence.csv` status and do not authorize a manuscript citation.

## Machine extraction gate

- Use only a legal publisher, official repository, author-manuscript repository, or user-supplied licensed copy.
- Validate an HTTPS response as a real PDF, run `pdfinfo` and `pdftotext`, match the publication identity, calculate SHA-256, and visually inspect each cited PDF page.
- Keep PDFs under ignored `literature/pdfs/`. Local validation opens that directory chain and each PDF without following symlinks, copies the already-open bytes into an atomically created private system-temporary directory, and runs Poppler only against that stable snapshot. Snapshot text and page renders are removed with the private directory on every exit path; `tmp/pdfs/` remains ignored only for manual extraction scratch files.
- Use `blocked_fulltext` when any acquisition or validation check fails. Never convert an HTML, XML, login, error, or challenge response into a PDF.
- Evidence-card excerpts must come from the body, not a title, abstract, or search snippet. All excerpts in one card have a combined limit of 25 words.
- Pin every extracted arXiv PDF URL to an explicit `v<number>` revision. Publisher and other official-repository URLs are unaffected.
- `selection_wave` uses `later` or `defer` for non-selected records and `phase_2b_<positive integer>` for extraction waves. Every phase-selected reference has exactly one card; cards cannot be orphaned from a phase-selected relevance row.
- Codex-generated cards remain `extracted_pending_human` with `human_confirmed: false` and null reviewer/date. A later `human_confirmed` card requires JSON boolean `true`, a reviewer, and an ISO-8601 UTC review timestamp; confirmation still does not promote `evidence.csv` or activate a citation.

## Human review gate

Shawn opens the PDF whose SHA-256 is recorded, then checks the source/version, section, PDF and printed page, locator, excerpt, paraphrase, limitations, and assessment. Human confirmation of a card and promotion of an evidence-ledger relationship are separate explicit decisions. Bibliographic verification alone is never claim evidence.

`make check` validates committed records without local PDFs. `make evidence-local-check` additionally requires Poppler, rechecks ignored PDFs and hashes, extracts non-empty text, validates the exact wrapped title and first-author surname from pages 1–2, and renders every locator page. A version of record must expose the registered DOI in those front-matter pages or trusted `pdfinfo` metadata; an accepted manuscript or preprint may omit it, but any DOI present there must agree with the registry. Full-document body or bibliography text cannot rescue a failed identity check.

Locator validation keeps a line-preserving source view for genuine section headings and Abstract/References/Bibliography boundaries. A separate exact excerpt view repairs conservative lowercase word continuations across line-ending hyphens while mapping every generated character back to its source line. The nearest genuine heading must precede the excerpt, with no forbidden heading line between them. Exact complete single-line headings are accepted. A schema-v1 two-line heading must begin with an Arabic section number, a canonical dotted Roman numeral, or a dotted alphabetic subsection label; contain no more than 160 normalized characters; and have no sentence-ending punctuation on either line except the dot in a label-only first line. Ambiguous unnumbered wrapping fails closed. Each Poppler subprocess has a 30-second fail-closed timeout and receives only the private stable snapshot.

For document identity, pages 1–2 are inputs from which the checker derives a bounded front-matter slice. Using the same complete-line heading recognizer, the slice ends immediately before the first Abstract, Introduction, References, or Bibliography heading, including numbered variants. Title, first-author, and page-text DOI checks use only that slice; later body, bibliography, and page-2 content after an earlier boundary cannot rescue identity. Trusted `pdfinfo` metadata may still supply DOI evidence under the version-specific rules.

## JSON schema v1

Required top-level fields are `schema_version`, `citation_key`, `doi`, `extraction_status`, `pdf_filename`, `pdf_sha256`, `landing_page_url`, `pdf_source_url`, `full_text_version`, `access_basis`, `relevance`, `coverage_area`, `blocker_reason`, `machine_extracted_at`, `human_confirmed`, `human_reviewer`, `human_reviewed_at`, and `claim_assessments`. Extracted cards also record `identity_checks` so a journal DOI absent from a repository/preprint PDF is not silently presented as a PDF-text match.

The only states are `blocked_fulltext`, `extracted_pending_human`, and `human_confirmed`. Blocked cards keep PDF fields and assessments empty and require no local PDF. Both extracted states require complete provenance, SHA-256, PDF, locator, excerpt, paraphrase, limitations, and assessment fields; `human_confirmed` additionally requires a true JSON boolean plus reviewer and UTC review timestamp. A confirmed card may remain linked to a pending ledger row because ledger promotion is a separate Shawn decision.
