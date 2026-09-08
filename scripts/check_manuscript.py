#!/usr/bin/env python3
"""Dependency-free checks for manuscript, evidence, and artifact integrity."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import unicodedata
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from urllib.parse import urlsplit


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_REL = Path("paper/manuscript.qmd")
BIBLIOGRAPHY_REL = Path("paper/references.bib")
EVIDENCE_REL = Path("literature/evidence.csv")
REFERENCES_REL = Path("literature/references.csv")
INVENTORY_REL = Path("research/artifact_inventory.csv")
RELEVANCE_REL = Path("literature/relevance.csv")
NOTES_REL = Path("literature/notes")
PDFS_REL = Path("literature/pdfs")
FIGURES_REL = Path("paper/figures")

BIB_ENTRY_START = re.compile(
    r"^\s*@([A-Za-z]+)\s*\{\s*([^,\s]+)\s*,", re.MULTILINE
)
CITE_KEY = re.compile(r"(?<![\w])@([A-Za-z0-9_:.-]+)")
TODO_ID = re.compile(r"\bTODO-EVIDENCE\s+([A-Z][A-Z0-9-]*)\b")
CLAIM_MARKER = re.compile(r"<!--\s*CLAIM:\s*([A-Z][A-Z0-9-]*)\s*-->")
CONFLICT_MARKER = re.compile(r"^(?:<<<<<<< .+|=======|>>>>>>> .+)$", re.MULTILINE)

EVIDENCE_COLUMNS = [
    "evidence_id",
    "claim_id",
    "status",
    "claim",
    "source_type",
    "evidence_path",
    "citation_key",
    "source_locator",
    "run_id",
    "owner",
    "verified_by",
    "verified_at",
    "notes",
]
REFERENCE_COLUMNS = [
    "citation_key",
    "reference_status",
    "title",
    "authors",
    "year",
    "venue",
    "doi",
    "doi_status",
    "publisher_url",
    "metadata_source_url",
    "verified_by",
    "verified_at",
    "notes",
]
INVENTORY_COLUMNS = [
    "artifact_id",
    "stage",
    "system",
    "artifact_type",
    "path",
    "status",
    "provenance_status",
    "manuscript_eligible",
    "notes",
]
RELEVANCE_COLUMNS = [
    "priority_rank",
    "citation_key",
    "relevance_category",
    "coverage_area",
    "repo_use_status",
    "selection_wave",
    "target_claim_ids",
    "rationale",
    "assessed_at",
]

EVIDENCE_STATUSES = {"supported", "partial", "planned", "pending", "rejected"}
REFERENCE_STATUSES = {"verified", "metadata_partial", "candidate", "rejected"}
DOI_STATUSES = {"verified", "not_assigned", "pending"}
RELEVANCE_CATEGORIES = {
    "core_direct",
    "method_support",
    "contextual",
    "conditional_on_actual_use",
    "reject_or_defer",
}
REPO_USE_STATUSES = {
    "not_applicable",
    "legacy_run_unqualified",
    "implemented_no_qualified_run",
    "placeholder_only",
}
NON_SELECTED_WAVES = {"later", "defer"}
EXTRACTION_WAVE = re.compile(r"phase_2b_[1-9]\d*")
PHASE_2B_1_FROZEN = {
    "qomi2022mineralization": "mineralization_interface_background",
    "astala2008hapwater": "hap_surface_hydration",
    "nowicki2024capture": "high_temperature_CO2_capture_boundary",
    "kuhne2020cp2k": "actual_compute_method",
}
PHASE_2B_1_KEYS = set(PHASE_2B_1_FROZEN)
PHASE_2B_1_COVERAGE = set(PHASE_2B_1_FROZEN.values())
CARD_EXTRACTION_STATUSES = {
    "blocked_fulltext",
    "extracted_pending_human",
    "human_confirmed",
}
CARD_ASSESSMENTS = {
    "supports",
    "partially_supports",
    "does_not_support",
    "unclear",
}
FULL_TEXT_VERSIONS = {
    "version_of_record",
    "author_accepted_manuscript",
    "preprint",
    "unknown_pending",
}
ACCESS_BASES = {
    "publisher_open_access",
    "official_repository",
    "author_manuscript_repository",
    "user_supplied_licensed_copy",
}
DOI_MATCH_BASES = {"pdf_text", "official_landing_metadata"}
FORBIDDEN_BIB_FIELDS = {
    "file",
    "attachment",
    "note",
    "annote",
    "abstract",
    "keywords",
}
LOCAL_PATH = re.compile(
    r"(?:file://|(?<![A-Za-z0-9])/(?:Users|Volumes|home|private|tmp)/|"
    r"(?<![A-Za-z0-9])~/|(?<![A-Za-z0-9])[A-Za-z]:[\\/])",
    re.IGNORECASE,
)
SOURCE_TYPES = {"review", "primary_publication", "repository", "computational_artifact"}
LITERATURE_SOURCE_TYPES = {"review", "primary_publication"}
INVENTORY_STATUSES = {
    "planned",
    "present_unreviewed",
    "completed",
    "partial",
    "failed",
    "invalid",
    "missing",
}
PROVENANCE_STATUSES = {"none", "path_only", "legacy_unverified", "manifest_complete"}
ELIGIBILITY_VALUES = {"yes", "no"}
PATH_REQUIRED_STATUSES = {"present_unreviewed", "completed", "partial", "failed", "invalid"}
MAX_TRACKED_FILE_BYTES = 10 * 1024 * 1024
POPPLER_TIMEOUT_SECONDS = 30
FORBIDDEN_CLAIM_SECTION_HEADINGS = {"abstract", "references", "bibliography"}
FRONT_MATTER_BOUNDARY_HEADINGS = {
    "abstract",
    "introduction",
    "references",
    "bibliography",
}
HEADING_LABEL = re.compile(
    r"^(?:(?:\d+(?:\.\d+)*|[ivxlcdm]+|[a-z])\.?)\s+(.+)$",
    flags=re.IGNORECASE,
)
SPLIT_HEADING_PREFIX_TOKEN = re.compile(
    r"^(?P<label>\d+(?:\.\d+)*\.?|[a-z]+\.)(?:\s+|$)",
    flags=re.IGNORECASE,
)
ROMAN_SECTION_NUMBER = re.compile(
    r"m{0,4}(?:cm|cd|d?c{0,3})(?:xc|xl|l?x{0,3})(?:ix|iv|v?i{0,3})",
    flags=re.IGNORECASE,
)
DOI_IN_DOCUMENT = re.compile(
    r"(?<![A-Za-z0-9])10\.\d{4,9}/[^\s<>{}\[\]\"']+",
    flags=re.IGNORECASE,
)
MAX_FRONT_MATTER_TITLE_START_LINE = 30
MAX_FRONT_MATTER_TITLE_LINES = 6
MAX_FRONT_MATTER_AUTHOR_LINES = 20
MAX_SPLIT_HEADING_CHARACTERS = 160
CARD_REQUIRED_FIELDS = {
    "schema_version",
    "citation_key",
    "doi",
    "extraction_status",
    "pdf_filename",
    "pdf_sha256",
    "landing_page_url",
    "pdf_source_url",
    "full_text_version",
    "access_basis",
    "relevance",
    "coverage_area",
    "blocker_reason",
    "machine_extracted_at",
    "human_confirmed",
    "human_reviewer",
    "human_reviewed_at",
    "claim_assessments",
}


@dataclass
class ValidationReport:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    counts: dict[str, int] = field(default_factory=dict)

    def error(self, message: str) -> None:
        self.errors.append(message)

    def warning(self, message: str) -> None:
        self.warnings.append(message)


@dataclass(frozen=True)
class PDFSourceLine:
    """One normalized pdftotext line with its original page-line position."""

    index: int
    text: str


@dataclass(frozen=True)
class PDFExcerptView:
    """Dehyphenated page text plus source-line provenance for every character."""

    text: str
    source_line_by_character: tuple[int, ...]


def _safe_repo_path(root: Path, value: str, label: str, report: ValidationReport) -> Path | None:
    """Resolve a non-empty repository-relative path without allowing traversal."""
    if not value:
        return None

    candidate = Path(value)
    if candidate.is_absolute() or ".." in candidate.parts:
        report.error(f"{label} must be a repository-relative path: {value}")
        return None

    resolved = (root / candidate).resolve()
    try:
        resolved.relative_to(root)
    except ValueError:
        report.error(f"{label} escapes the repository: {value}")
        return None
    return resolved


def _read_csv(
    path: Path,
    expected_columns: list[str],
    label: str,
    report: ValidationReport,
) -> list[dict[str, str]]:
    if not path.is_file():
        report.error(f"missing required file: {path}")
        return []

    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            actual_columns = reader.fieldnames or []
            if actual_columns != expected_columns:
                report.error(
                    f"{label} columns differ from the required schema: "
                    f"expected {expected_columns}, got {actual_columns}"
                )
            return list(reader)
    except (OSError, csv.Error, UnicodeError) as exc:
        report.error(f"could not read {label}: {exc}")
        return []


def _duplicates(values: list[str]) -> list[str]:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return sorted(duplicates)


def _matching_brace(text: str, opening: int) -> int | None:
    """Return the closing brace paired with *opening*, respecting escapes."""
    depth = 0
    escaped = False
    for index in range(opening, len(text)):
        character = text[index]
        if escaped:
            escaped = False
            continue
        if character == "\\":
            escaped = True
        elif character == "{":
            depth += 1
        elif character == "}":
            depth -= 1
            if depth == 0:
                return index
    return None


def _parse_bib_fields(
    body: str,
    citation_key: str,
    report: ValidationReport,
) -> dict[str, str]:
    """Parse the simple field forms emitted by Better BibTeX."""
    fields: dict[str, str] = {}
    position = 0
    while position < len(body):
        while position < len(body) and (body[position].isspace() or body[position] == ","):
            position += 1
        if position >= len(body):
            break

        name_match = re.match(r"[A-Za-z][A-Za-z0-9_-]*", body[position:])
        if name_match is None:
            report.error(f"{citation_key}: could not parse BibTeX field")
            break
        field_name = name_match.group(0).lower()
        position += name_match.end()
        while position < len(body) and body[position].isspace():
            position += 1
        if position >= len(body) or body[position] != "=":
            report.error(f"{citation_key}: BibTeX field '{field_name}' is missing '='")
            break
        position += 1
        while position < len(body) and body[position].isspace():
            position += 1
        if position >= len(body):
            report.error(f"{citation_key}: BibTeX field '{field_name}' has no value")
            break

        if body[position] == "{":
            closing = _matching_brace(body, position)
            if closing is None:
                report.error(
                    f"{citation_key}: BibTeX field '{field_name}' has unbalanced braces"
                )
                break
            value = body[position + 1 : closing]
            position = closing + 1
        elif body[position] == '"':
            start = position + 1
            position = start
            escaped = False
            while position < len(body):
                character = body[position]
                if escaped:
                    escaped = False
                elif character == "\\":
                    escaped = True
                elif character == '"':
                    break
                position += 1
            if position >= len(body):
                report.error(
                    f"{citation_key}: BibTeX field '{field_name}' has an unterminated quote"
                )
                break
            value = body[start:position]
            position += 1
        else:
            start = position
            while position < len(body) and body[position] != ",":
                position += 1
            value = body[start:position].strip()

        if field_name in fields:
            report.error(f"{citation_key}: duplicate BibTeX field '{field_name}'")
        else:
            fields[field_name] = value.strip()

    return fields


def _parse_bib_entries(
    bibliography_text: str,
    report: ValidationReport,
) -> list[tuple[str, dict[str, str]]]:
    """Build a citation-key-to-fields representation of Better BibTeX output."""
    entries: list[tuple[str, dict[str, str]]] = []
    position = 0
    while True:
        match = BIB_ENTRY_START.search(bibliography_text, position)
        if match is None:
            break
        entry_type, citation_key = match.groups()
        opening = bibliography_text.find("{", match.start(), match.end())
        closing = _matching_brace(bibliography_text, opening)
        if closing is None:
            report.error(f"{citation_key}: BibTeX entry has unbalanced braces")
            break
        position = closing + 1
        if entry_type.lower() in {"comment", "preamble", "string"}:
            continue
        entries.append(
            (
                citation_key,
                _parse_bib_fields(
                    bibliography_text[match.end() : closing],
                    citation_key,
                    report,
                ),
            )
        )
    return entries


def _normalize_doi(value: str) -> str:
    """Normalize DOI transport forms without changing identifier content."""
    normalized = value.strip()
    normalized = re.sub(
        r"^(?:https?://(?:dx\.)?doi\.org/|doi\s*:\s*)",
        "",
        normalized,
        count=1,
        flags=re.IGNORECASE,
    )
    return normalized.strip().lower().rstrip(".,;:").strip()


def _normalize_bibliographic_title(value: str) -> str:
    """Conservatively normalize display-only BibTeX title protection."""
    normalized = unicodedata.normalize("NFKC", value)
    normalized = normalized.replace("{", "").replace("}", "")
    return " ".join(normalized.split()).casefold()


def _normalize_pdf_identity(value: str) -> str:
    """Normalize PDF-extracted text without permitting fuzzy identity matches."""
    normalized = unicodedata.normalize("NFKC", value).casefold()
    normalized = re.sub(r"[^\w]+", " ", normalized, flags=re.UNICODE)
    return " ".join(normalized.split())


def _normalize_pdf_source_lines(value: str) -> list[PDFSourceLine]:
    """Normalize horizontal text while preserving every source line boundary."""
    normalized = unicodedata.normalize("NFKC", value).replace("\u00ad", "")
    normalized = normalized.replace("\r\n", "\n").replace("\r", "\n")
    return [
        PDFSourceLine(index=index, text=" ".join(line.split()))
        for index, line in enumerate(normalized.split("\n"))
    ]


def _is_genuine_cross_line_hyphenation(previous: str, following: str) -> bool:
    """Conservatively recognize a wrapped word without consuming heading lines."""
    if len(previous) < 2 or not previous.endswith("-") or not previous[-2].isalnum():
        return False
    if not following or not following[0].isalnum() or not following[0].islower():
        return False
    previous_word = re.search(r"[\w]+-$", previous, flags=re.UNICODE)
    following_word = re.match(r"[\w]+", following, flags=re.UNICODE)
    return previous_word is not None and following_word is not None


def _build_pdf_excerpt_view(lines: list[PDFSourceLine]) -> PDFExcerptView:
    """Flatten source lines for exact excerpts while retaining line provenance."""
    characters: list[str] = []
    source_lines: list[int] = []
    previous_text = ""
    for line in lines:
        if not line.text:
            continue
        if characters:
            if _is_genuine_cross_line_hyphenation(previous_text, line.text):
                characters.pop()
                source_lines.pop()
            else:
                characters.append(" ")
                source_lines.append(line.index)
        characters.extend(line.text)
        source_lines.extend([line.index] * len(line.text))
        previous_text = line.text
    return PDFExcerptView("".join(characters), tuple(source_lines))


def _normalize_pdf_locator_text(value: str) -> str:
    """Normalize and conservatively dehyphenate an exact locator value."""
    return _build_pdf_excerpt_view(_normalize_pdf_source_lines(value)).text


def _normalize_heading(value: str) -> str:
    return " ".join(unicodedata.normalize("NFKC", value).split()).casefold()


def _heading_without_label(value: str) -> str:
    normalized = _normalize_heading(value)
    match = HEADING_LABEL.fullmatch(normalized)
    return match.group(1) if match is not None else normalized


def _heading_values_match(candidate: str, declared: str) -> bool:
    candidate_normalized = _normalize_heading(candidate)
    declared_normalized = _normalize_heading(declared)
    return candidate_normalized == declared_normalized or (
        _heading_without_label(candidate_normalized)
        == _heading_without_label(declared_normalized)
    )


def _structured_split_heading_prefix_end(value: str) -> int | None:
    """Return a safe Arabic, canonical dotted Roman, or dotted-letter prefix."""
    normalized = " ".join(value.split())
    match = SPLIT_HEADING_PREFIX_TOKEN.match(normalized)
    if match is None:
        return None
    label = match.group("label")
    if label[0].isdigit():
        return match.end()
    stem = label[:-1]
    if len(stem) == 1 or ROMAN_SECTION_NUMBER.fullmatch(stem) is not None:
        return match.end()
    return None


def _is_sentence_like_split_heading_line(value: str) -> bool:
    """Reject sentence punctuation while allowing a label-only line such as A."""
    normalized = " ".join(value.split())
    if _structured_split_heading_prefix_end(normalized) == len(normalized):
        return False
    return re.search(r"[.!?]\s*$", normalized) is not None


def _genuine_heading_spans(
    lines: list[PDFSourceLine], declared: str
) -> list[tuple[int, int]]:
    """Return exact one-line or conservatively structured two-line headings.

    A two-line schema-v1 heading must start with an Arabic/Roman section
    number or a dotted alphabetic subsection label, contain at most 160
    normalized characters, and have no sentence-ending punctuation on either
    line. Ambiguous unnumbered wrapped headings fail closed.
    """
    spans: list[tuple[int, int]] = []
    for position, line in enumerate(lines):
        if not line.text:
            continue
        if _heading_values_match(line.text, declared):
            spans.append((line.index, line.index))
        if position + 1 >= len(lines):
            continue
        following = lines[position + 1]
        if not following.text:
            continue
        combined = f"{line.text} {following.text}"
        if (
            len(combined) <= MAX_SPLIT_HEADING_CHARACTERS
            and _structured_split_heading_prefix_end(line.text) is not None
            and not _is_sentence_like_split_heading_line(line.text)
            and not _is_sentence_like_split_heading_line(following.text)
            and _heading_values_match(combined, declared)
        ):
            spans.append((line.index, following.index))
    return spans


def _forbidden_heading_lines(lines: list[PDFSourceLine]) -> list[tuple[int, str]]:
    """Identify forbidden headings only from complete line-preserving records."""
    boundaries: list[tuple[int, str]] = []
    for line in lines:
        if line.text and _heading_without_label(line.text) in FORBIDDEN_CLAIM_SECTION_HEADINGS:
            boundaries.append((line.index, line.text))
    return boundaries


def _front_matter_source_lines(value: str) -> list[PDFSourceLine]:
    """Derive front matter from pages 1--2 using shared heading recognition."""
    lines = _normalize_pdf_source_lines(value)
    boundary_starts = [
        span[0]
        for heading in FRONT_MATTER_BOUNDARY_HEADINGS
        for span in _genuine_heading_spans(lines, heading)
    ]
    if not boundary_starts:
        return lines
    first_boundary = min(boundary_starts)
    return [line for line in lines if line.index < first_boundary]


def _front_matter_title_span(
    source_lines: list[PDFSourceLine], registered_title: str
) -> tuple[list[PDFSourceLine], tuple[int, int] | None]:
    """Find an exact, possibly wrapped title block near the document front."""
    lines = [line for line in source_lines if line.text]
    target = _normalize_pdf_identity(registered_title)
    search_limit = min(len(lines), MAX_FRONT_MATTER_TITLE_START_LINE)
    for start in range(search_limit):
        for width in range(1, MAX_FRONT_MATTER_TITLE_LINES + 1):
            end = start + width
            if end > len(lines):
                break
            candidate = _build_pdf_excerpt_view(lines[start:end]).text
            if _normalize_pdf_identity(candidate) == target:
                return lines, (start, end - 1)
    return lines, None


def _first_author_follows_title(
    lines: list[PDFSourceLine], title_span: tuple[int, int], first_author: str
) -> bool:
    """Require the registered first author in the author block after the title."""
    start = title_span[1] + 1
    end = min(len(lines), start + MAX_FRONT_MATTER_AUTHOR_LINES)
    author_text = _normalize_pdf_identity(
        _build_pdf_excerpt_view(lines[start:end]).text
    )
    normalized_author = _normalize_pdf_identity(first_author)
    return bool(
        normalized_author
        and f" {normalized_author} " in f" {author_text} "
    )


def _trim_document_doi(value: str) -> str:
    """Remove sentence punctuation and unmatched closing delimiters from a DOI."""
    trimmed = value.rstrip(".,;:")
    pairs = (("(", ")"), ("[", "]"), ("{", "}"))
    changed = True
    while changed and trimmed:
        changed = False
        for opening, closing in pairs:
            if trimmed.endswith(closing) and trimmed.count(opening) < trimmed.count(closing):
                trimmed = trimmed[:-1].rstrip(".,;:")
                changed = True
    return trimmed


def _document_dois(value: str) -> set[str]:
    return {
        _normalize_doi(_trim_document_doi(match.group(0)))
        for match in DOI_IN_DOCUMENT.finditer(value)
        if _normalize_doi(_trim_document_doi(match.group(0)))
    }


def _validate_front_matter_identity(
    front_matter: str,
    pdfinfo_text: str,
    card: dict[str, object],
    reference: dict[str, str],
    report: ValidationReport,
) -> None:
    """Validate document identity without allowing body or bibliography rescue."""
    key = str(card.get("citation_key", "")).strip()
    front_matter_lines = _front_matter_source_lines(front_matter)
    lines, title_span = _front_matter_title_span(
        front_matter_lines, reference.get("title", "")
    )
    if title_span is None:
        report.error(f"{key}: registry title was not found as a front-matter title block")
    else:
        first_author = reference.get("authors", "").split("|", 1)[0].split(",", 1)[0]
        if not _first_author_follows_title(lines, title_span, first_author):
            report.error(
                f"{key}: registry first author was not found in the front-matter author block"
            )

    registered_doi = _normalize_doi(reference.get("doi", ""))
    front_matter_text = "\n".join(line.text for line in front_matter_lines)
    observed_dois = _document_dois(front_matter_text) | _document_dois(pdfinfo_text)
    conflicting_dois = sorted(doi for doi in observed_dois if doi != registered_doi)
    if conflicting_dois:
        report.error(
            f"{key}: front matter or PDF metadata contains a conflicting DOI: "
            + ", ".join(conflicting_dois)
        )

    version = str(card.get("full_text_version", "")).strip()
    if version == "version_of_record" and registered_doi not in observed_dois:
        report.error(
            f"{key}: version_of_record PDF requires the registered DOI in front matter "
            "or PDF metadata"
        )


def _all_occurrences(text: str, value: str) -> list[int]:
    """Return every exact, potentially overlapping occurrence of *value*."""
    if not value:
        return []
    positions: list[int] = []
    start = 0
    while True:
        position = text.find(value, start)
        if position < 0:
            return positions
        positions.append(position)
        start = position + 1


def _validate_locator_on_page(
    page_text: str,
    section_heading: str,
    excerpt: str,
) -> tuple[str | None, str | None]:
    """Match an excerpt exactly using source-line-aware heading boundaries."""
    source_lines = _normalize_pdf_source_lines(page_text)
    excerpt_view = _build_pdf_excerpt_view(source_lines)
    heading_spans = _genuine_heading_spans(source_lines, section_heading)
    forbidden_boundaries = _forbidden_heading_lines(source_lines)
    normalized_excerpt = _normalize_pdf_locator_text(excerpt)
    excerpt_positions = _all_occurrences(excerpt_view.text, normalized_excerpt)

    if not heading_spans:
        return "heading_missing", None
    if not excerpt_positions:
        return "excerpt_missing", None

    saw_preceding_heading = False
    blocked_heading: str | None = None
    for excerpt_position in excerpt_positions:
        excerpt_line = excerpt_view.source_line_by_character[excerpt_position]
        preceding = [span for span in heading_spans if span[1] < excerpt_line]
        if not preceding:
            continue
        saw_preceding_heading = True
        nearest_heading = max(preceding, key=lambda span: (span[1], span[0]))
        heading_end_line = nearest_heading[1]
        intervening = [
            (line_index, heading)
            for line_index, heading in forbidden_boundaries
            if heading_end_line < line_index <= excerpt_line
        ]
        if not intervening:
            return None, None
        blocked_heading = min(intervening)[1]

    if not saw_preceding_heading:
        return "excerpt_precedes_heading", None
    return "forbidden_boundary", blocked_heading


def _is_forbidden_claim_section_heading(value: str) -> bool:
    normalized = unicodedata.normalize("NFKC", value).casefold()
    words = re.findall(r"[a-z]+", normalized)
    return bool(set(words) & FORBIDDEN_CLAIM_SECTION_HEADINGS)


def _is_versioned_arxiv_pdf_url(value: str) -> bool | None:
    """Return whether an arXiv PDF URL has an explicit version, or None otherwise."""
    parsed = urlsplit(value)
    hostname = (parsed.hostname or "").casefold()
    if hostname != "arxiv.org" and not hostname.endswith(".arxiv.org"):
        return None
    return bool(
        parsed.path.startswith("/pdf/")
        and re.search(r"v[1-9]\d*(?:\.pdf)?/?$", parsed.path, flags=re.IGNORECASE)
    )


def _card_word_count(value: str) -> int:
    return len(re.findall(r"\b[\w’'-]+\b", value, flags=re.UNICODE))


def _is_iso_utc_timestamp(value: object) -> bool:
    if not isinstance(value, str) or not re.fullmatch(
        r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z", value.strip()
    ):
        return False
    try:
        datetime.strptime(value.strip(), "%Y-%m-%dT%H:%M:%SZ")
    except ValueError:
        return False
    return True


def _read_evidence_card(
    path: Path,
    citation_key: str,
    report: ValidationReport,
) -> dict[str, object] | None:
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as exc:
        report.error(f"{citation_key}: could not read evidence card: {exc}")
        return None
    if LOCAL_PATH.search(text):
        report.error(f"{citation_key}: evidence card contains a machine-local path")

    blocks = re.findall(r"```json\s*(\{.*?\})\s*```", text, flags=re.DOTALL)
    if len(blocks) != 1:
        report.error(
            f"{citation_key}: evidence card must contain exactly one JSON code block"
        )
        return None
    try:
        card = json.loads(blocks[0])
    except json.JSONDecodeError as exc:
        report.error(f"{citation_key}: evidence card JSON is invalid: {exc}")
        return None
    if not isinstance(card, dict):
        report.error(f"{citation_key}: evidence card JSON must be an object")
        return None
    missing = sorted(CARD_REQUIRED_FIELDS - card.keys())
    if missing:
        report.error(
            f"{citation_key}: evidence card missing required field(s): "
            + ", ".join(missing)
        )
    return card


def _validate_evidence_card(
    card: dict[str, object],
    filename_key: str,
    reference_by_key: dict[str, dict[str, str]],
    relationships_by_claim_key: dict[tuple[str, str], dict[str, str]],
    report: ValidationReport,
) -> None:
    key = str(card.get("citation_key", "")).strip()
    label = key or filename_key
    if key != filename_key:
        report.error(
            f"{label}: evidence-card filename and citation_key must match"
        )
    if card.get("schema_version") != 1:
        report.error(f"{label}: evidence card schema_version must be 1")

    reference = reference_by_key.get(key)
    if reference is None:
        report.error(f"{label}: evidence card citation_key is missing from references.csv")
    elif _normalize_doi(str(card.get("doi", ""))) != _normalize_doi(
        reference.get("doi", "")
    ):
        report.error(f"{label}: evidence-card DOI does not match references.csv")

    status = str(card.get("extraction_status", "")).strip()
    if status not in CARD_EXTRACTION_STATUSES:
        report.error(f"{label}: unsupported evidence-card extraction_status '{status}'")

    for url_field in ("landing_page_url", "pdf_source_url"):
        value = str(card.get(url_field, "") or "").strip()
        if not value.startswith("https://"):
            report.error(f"{label}: {url_field} must use HTTPS")
    pdf_source_url = str(card.get("pdf_source_url", "") or "").strip()
    arxiv_versioned = _is_versioned_arxiv_pdf_url(pdf_source_url)
    if status in {"extracted_pending_human", "human_confirmed"} and arxiv_versioned is False:
        report.error(
            f"{label}: extracted arXiv pdf_source_url must include an explicit v<number> suffix"
        )

    version = str(card.get("full_text_version", "")).strip()
    if version not in FULL_TEXT_VERSIONS:
        report.error(f"{label}: unsupported full_text_version '{version}'")
    access_basis = str(card.get("access_basis", "")).strip()
    if access_basis not in ACCESS_BASES:
        report.error(f"{label}: unsupported access_basis '{access_basis}'")
    if not str(card.get("relevance", "")).strip():
        report.error(f"{label}: evidence card requires a relevance statement")
    if not str(card.get("coverage_area", "")).strip():
        report.error(f"{label}: evidence card requires coverage_area")
    timestamp = card.get("machine_extracted_at")
    if not _is_iso_utc_timestamp(timestamp):
        report.error(f"{label}: machine_extracted_at must be ISO-8601 UTC")

    human_confirmed = card.get("human_confirmed")
    human_reviewer = card.get("human_reviewer")
    human_reviewed_at = card.get("human_reviewed_at")
    if status in {"blocked_fulltext", "extracted_pending_human"}:
        if human_confirmed is not False:
            report.error(
                f"{label}: {status} requires human_confirmed to be JSON false"
            )
        if human_reviewer is not None or human_reviewed_at is not None:
            report.error(
                f"{label}: unconfirmed evidence card must leave human reviewer/date null"
            )
    elif status == "human_confirmed":
        if human_confirmed is not True:
            report.error(
                f"{label}: human_confirmed status requires human_confirmed to be JSON true"
            )
        if not isinstance(human_reviewer, str) or not human_reviewer.strip():
            report.error(f"{label}: human_confirmed status requires human_reviewer")
        if not _is_iso_utc_timestamp(human_reviewed_at):
            report.error(
                f"{label}: human_reviewed_at must be an ISO-8601 UTC timestamp"
            )

    assessments = card.get("claim_assessments")
    if not isinstance(assessments, list):
        report.error(f"{label}: claim_assessments must be a JSON array")
        assessments = []

    if status == "blocked_fulltext":
        if not str(card.get("blocker_reason", "") or "").strip():
            report.error(f"{label}: blocked_fulltext requires blocker_reason")
        for field_name in ("pdf_filename", "pdf_sha256"):
            if card.get(field_name) not in (None, ""):
                report.error(
                    f"{label}: blocked_fulltext must leave {field_name} empty"
                )
        if assessments:
            report.error(
                f"{label}: blocked_fulltext cannot contain claim assessments or locators"
            )
        return

    if card.get("blocker_reason") is not None:
        report.error(f"{label}: extracted card must leave blocker_reason null")
    expected_filename = f"{key}.pdf"
    pdf_filename = str(card.get("pdf_filename", "")).strip()
    if (
        pdf_filename != expected_filename
        or Path(pdf_filename).name != pdf_filename
        or "/" in pdf_filename
        or "\\" in pdf_filename
    ):
        report.error(
            f"{label}: pdf_filename must be the basename {expected_filename}"
        )
    checksum = str(card.get("pdf_sha256", "")).strip()
    if not re.fullmatch(r"[0-9a-f]{64}", checksum):
        report.error(f"{label}: extracted card requires a lowercase SHA-256")

    identity = card.get("identity_checks")
    if not isinstance(identity, dict):
        report.error(f"{label}: extracted card requires identity_checks")
    else:
        if identity.get("title_in_pdf_text") is not True:
            report.error(f"{label}: PDF title identity check must pass")
        if identity.get("author_in_pdf_text") is not True:
            report.error(f"{label}: PDF author identity check must pass")
        doi_basis = str(identity.get("doi_match_basis", "")).strip()
        if doi_basis not in DOI_MATCH_BASES:
            report.error(f"{label}: unsupported DOI match basis '{doi_basis}'")
        association_url = str(identity.get("doi_association_url", "") or "").strip()
        if doi_basis == "official_landing_metadata" and not association_url.startswith(
            "https://"
        ):
            report.error(
                f"{label}: official landing DOI match requires an HTTPS association URL"
            )

    if not assessments:
        report.error(f"{label}: extracted card requires at least one claim assessment")
    quote_words = 0
    seen_claims: set[str] = set()
    for index, item in enumerate(assessments, start=1):
        item_label = f"{label}: claim_assessments[{index}]"
        if not isinstance(item, dict):
            report.error(f"{item_label} must be a JSON object")
            continue
        claim_id = str(item.get("claim_id", "")).strip()
        if not claim_id:
            report.error(f"{item_label} requires claim_id")
        elif claim_id in seen_claims:
            report.error(f"{label}: duplicate claim assessment for {claim_id}")
        else:
            seen_claims.add(claim_id)
        if claim_id and (claim_id, key) not in relationships_by_claim_key:
            report.error(
                f"{label}: claim assessment has no evidence.csv relationship: {claim_id}"
            )
        assessment = str(item.get("assessment", "")).strip()
        if assessment not in CARD_ASSESSMENTS:
            report.error(f"{item_label} has unsupported assessment '{assessment}'")
        section_heading = str(item.get("section_heading", "")).strip()
        if not section_heading:
            report.error(f"{item_label} requires section_heading")
        elif _is_forbidden_claim_section_heading(section_heading):
            report.error(
                f"{item_label} cannot use Abstract, References, or Bibliography "
                "as a claim-bearing section heading"
            )
        pdf_page = item.get("pdf_page")
        if not isinstance(pdf_page, int) or isinstance(pdf_page, bool) or pdf_page < 1:
            report.error(f"{item_label} requires a positive integer pdf_page")
        for field_name in (
            "printed_page",
            "figure_table_equation_locator",
            "verbatim_excerpt",
            "faithful_paraphrase",
            "scope_limitations",
        ):
            if not str(item.get(field_name, "")).strip():
                report.error(f"{item_label} requires {field_name}")
        quote_words += _card_word_count(str(item.get("verbatim_excerpt", "")))
    if quote_words > 25:
        report.error(
            f"{label}: verbatim excerpts total {quote_words} words; maximum is 25"
        )


def _open_local_pdf_descriptor(
    root: Path,
    citation_key: str,
    filename: str,
    report: ValidationReport,
) -> int | None:
    """Open a PDF through a no-follow directory-descriptor chain."""
    expected_filename = f"{citation_key}.pdf"
    if (
        filename != expected_filename
        or Path(filename).name != filename
        or "/" in filename
        or "\\" in filename
    ):
        report.error(
            f"{citation_key}: local pdf_filename must be the basename {expected_filename}"
        )
        return None

    directory_flags = (
        os.O_RDONLY
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
        | getattr(os, "O_CLOEXEC", 0)
    )
    file_flags = (
        os.O_RDONLY
        | getattr(os, "O_NOFOLLOW", 0)
        | getattr(os, "O_NONBLOCK", 0)
        | getattr(os, "O_CLOEXEC", 0)
    )
    directory_descriptors: list[int] = []
    file_descriptor: int | None = None
    try:
        try:
            root_descriptor = os.open(root, directory_flags)
        except OSError as exc:
            report.error(
                f"{citation_key}: repository root could not be opened safely: {exc}"
            )
            return None
        directory_descriptors.append(root_descriptor)

        parent_descriptor = root_descriptor
        for component in ("literature", "pdfs"):
            try:
                descriptor = os.open(
                    component,
                    directory_flags,
                    dir_fd=parent_descriptor,
                )
            except OSError as exc:
                report.error(
                    f"{citation_key}: local evidence PDF is missing or its directory "
                    f"chain is unsafe; directories must be real and must not follow "
                    f"symlinks: {component}: {exc}"
                )
                return None
            directory_descriptors.append(descriptor)
            parent_descriptor = descriptor

        try:
            file_descriptor = os.open(
                filename,
                file_flags,
                dir_fd=parent_descriptor,
            )
        except OSError as exc:
            report.error(
                f"{citation_key}: local evidence PDF is missing, inaccessible, or a "
                f"symlink: {PDFS_REL / filename}: {exc}"
            )
            return None

        try:
            file_status = os.fstat(file_descriptor)
        except OSError as exc:
            report.error(
                f"{citation_key}: could not inspect local evidence PDF descriptor: {exc}"
            )
            return None
        if not stat.S_ISREG(file_status.st_mode):
            report.error(f"{citation_key}: local evidence PDF must be a regular file")
            return None

        result = file_descriptor
        file_descriptor = None
        return result
    finally:
        if file_descriptor is not None:
            os.close(file_descriptor)
        for descriptor in reversed(directory_descriptors):
            os.close(descriptor)


def _copy_pdf_snapshot(
    source_descriptor: int,
    snapshot_path: Path,
    citation_key: str,
    expected_hash: str,
    report: ValidationReport,
) -> Path | None:
    """Copy and hash exactly one already-open source into a private snapshot."""
    digest = hashlib.sha256()
    header = bytearray()
    destination_descriptor: int | None = None
    try:
        source = os.fdopen(source_descriptor, "rb", closefd=True)
    except OSError as exc:
        os.close(source_descriptor)
        report.error(f"{citation_key}: could not read opened PDF descriptor: {exc}")
        return None
    try:
        with source:
            destination_descriptor = os.open(
                snapshot_path,
                os.O_WRONLY
                | os.O_CREAT
                | os.O_EXCL
                | getattr(os, "O_NOFOLLOW", 0)
                | getattr(os, "O_CLOEXEC", 0),
                0o600,
            )
            with os.fdopen(destination_descriptor, "wb", closefd=True) as destination:
                destination_descriptor = None
                while True:
                    chunk = source.read(1024 * 1024)
                    if not chunk:
                        break
                    if len(header) < 5:
                        header.extend(chunk[: 5 - len(header)])
                    digest.update(chunk)
                    destination.write(chunk)
    except OSError as exc:
        report.error(f"{citation_key}: could not create stable PDF snapshot: {exc}")
        return None
    finally:
        if destination_descriptor is not None:
            os.close(destination_descriptor)

    if bytes(header) != b"%PDF-":
        report.error(f"{citation_key}: local evidence file does not begin with %PDF-")
        return None
    actual_hash = digest.hexdigest()
    if actual_hash != expected_hash:
        report.error(
            f"{citation_key}: local evidence PDF SHA-256 mismatch: expected "
            f"{expected_hash}, actual {actual_hash}"
        )
        return None
    return snapshot_path


def _run_poppler(
    arguments: list[str],
    operation: str,
    citation_key: str,
    report: ValidationReport,
    *,
    text: bool = True,
) -> subprocess.CompletedProcess[str] | subprocess.CompletedProcess[bytes] | None:
    """Run one Poppler command with bounded, deterministic failure handling."""
    try:
        result = subprocess.run(
            arguments,
            check=False,
            capture_output=True,
            text=text,
            timeout=POPPLER_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired:
        report.error(
            f"{citation_key}: {operation} timed out after "
            f"{POPPLER_TIMEOUT_SECONDS} seconds"
        )
        return None
    except FileNotFoundError:
        report.error(f"{citation_key}: {operation} executable was not found")
        return None
    except OSError:
        report.error(f"{citation_key}: {operation} could not be started")
        return None
    if result.returncode != 0:
        report.error(
            f"{citation_key}: {operation} failed with exit status {result.returncode}"
        )
        return None
    return result


def _poppler_executable(name: str) -> str:
    return shutil.which(name) or name


def _validate_local_pdf(
    root: Path,
    card: dict[str, object],
    reference: dict[str, str],
    report: ValidationReport,
) -> None:
    key = str(card.get("citation_key", "")).strip()
    filename = str(card.get("pdf_filename", "")).strip()
    try:
        with tempfile.TemporaryDirectory(prefix=f"{key}-evidence-") as directory:
            source_descriptor = _open_local_pdf_descriptor(
                root, key, filename, report
            )
            if source_descriptor is None:
                return
            snapshot_path = _copy_pdf_snapshot(
                source_descriptor,
                Path(directory) / "source.pdf",
                key,
                str(card.get("pdf_sha256", "")).strip(),
                report,
            )
            if snapshot_path is None:
                return
            _validate_pdf_snapshot(
                snapshot_path,
                Path(directory),
                card,
                reference,
                report,
            )
    except OSError as exc:
        report.error(f"{key}: private PDF validation directory failed: {exc}")


def _validate_pdf_snapshot(
    pdf_path: Path,
    private_directory: Path,
    card: dict[str, object],
    reference: dict[str, str],
    report: ValidationReport,
) -> None:
    """Run identity, locator, and rendering checks on one immutable snapshot."""
    key = str(card.get("citation_key", "")).strip()
    info = _run_poppler(
        [_poppler_executable("pdfinfo"), str(pdf_path)],
        "pdfinfo",
        key,
        report,
    )
    if info is None:
        return
    assert isinstance(info.stdout, str)
    match = re.search(r"^Pages:\s+(\d+)\s*$", info.stdout, flags=re.MULTILINE)
    page_count = int(match.group(1)) if match else 0
    if page_count < 1:
        report.error(f"{key}: local evidence PDF has no readable pages")
        return

    front_page_count = min(page_count, 2)
    front_matter = _run_poppler(
        [
            _poppler_executable("pdftotext"),
            "-f",
            "1",
            "-l",
            str(front_page_count),
            "-layout",
            str(pdf_path),
            "-",
        ],
        "pdftotext front matter",
        key,
        report,
    )
    if front_matter is None:
        return
    assert isinstance(front_matter.stdout, str)
    if not front_matter.stdout.strip():
        report.error(f"{key}: pdftotext returned no usable front-matter text")
        return
    assert isinstance(info.stdout, str)
    _validate_front_matter_identity(
        front_matter.stdout,
        info.stdout,
        card,
        reference,
        report,
    )

    extracted = _run_poppler(
        [_poppler_executable("pdftotext"), "-layout", str(pdf_path), "-"],
        "pdftotext full document",
        key,
        report,
    )
    if extracted is None:
        return
    assert isinstance(extracted.stdout, str)
    if not extracted.stdout.strip():
        report.error(f"{key}: pdftotext returned no usable text")
        return

    assessments = [
        item
        for item in card.get("claim_assessments", [])
        if isinstance(item, dict)
    ]
    pages = sorted(
        {
            item.get("pdf_page")
            for item in assessments
            if isinstance(item.get("pdf_page"), int)
        }
    )
    valid_pages: list[int] = []
    for page in pages:
        if page < 1 or page > page_count:
            report.error(
                f"{key}: locator page {page} is outside PDF page count {page_count}"
            )
            continue
        valid_pages.append(page)

    page_text_by_page: dict[int, str] = {}
    for page in valid_pages:
        page_text = _run_poppler(
            [
                _poppler_executable("pdftotext"),
                "-f",
                str(page),
                "-l",
                str(page),
                str(pdf_path),
                "-",
            ],
            f"pdftotext locator page {page}",
            key,
            report,
        )
        if page_text is not None:
            assert isinstance(page_text.stdout, str)
            if not page_text.stdout.strip():
                report.error(f"{key}: locator PDF page {page} has no usable text")
            else:
                page_text_by_page[page] = page_text.stdout

        output_prefix = private_directory / f"locator-{page}"
        rendered = _run_poppler(
            [
                _poppler_executable("pdftoppm"),
                "-f",
                str(page),
                "-l",
                str(page),
                "-singlefile",
                "-png",
                str(pdf_path),
                str(output_prefix),
            ],
            f"pdftoppm locator page {page}",
            key,
            report,
            text=False,
        )
        image_path = output_prefix.with_suffix(".png")
        if rendered is not None:
            try:
                valid_image = image_path.is_file() and image_path.stat().st_size > 0
            except OSError:
                valid_image = False
            if not valid_image:
                report.error(f"{key}: could not render locator PDF page {page}")

    for index, item in enumerate(assessments, start=1):
        page = item.get("pdf_page")
        if not isinstance(page, int) or page not in page_text_by_page:
            continue
        item_label = f"{key}: claim_assessments[{index}]"
        failure, boundary = _validate_locator_on_page(
            page_text_by_page[page],
            str(item.get("section_heading", "")).strip(),
            str(item.get("verbatim_excerpt", "")).strip(),
        )
        if failure == "heading_missing":
            report.error(
                f"{item_label}: section heading was not found on locator PDF page {page}"
            )
        elif failure == "excerpt_missing":
            report.error(
                f"{item_label}: verbatim excerpt was not found on locator PDF page {page}"
            )
        elif failure == "excerpt_precedes_heading":
            report.error(
                f"{item_label}: verbatim excerpt precedes the declared section heading "
                f"on locator PDF page {page}"
            )
        elif failure == "forbidden_boundary":
            report.error(
                f"{item_label}: locator crosses forbidden section heading '{boundary}' "
                f"on PDF page {page}"
            )


def _git_tracked_files(root: Path, report: ValidationReport) -> list[Path]:
    if not (root / ".git").exists():
        return []
    result = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=root,
        check=False,
        capture_output=True,
    )
    if result.returncode != 0:
        report.error("git ls-files failed while checking tracked-file policy")
        return []
    return [root / item.decode("utf-8") for item in result.stdout.split(b"\0") if item]


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _normalized_repository_identifier(value: str) -> str:
    normalized = value.strip().rstrip("/")
    if normalized.endswith(".git"):
        normalized = normalized[:-4]
    if "://" in normalized:
        normalized = normalized.split("://", 1)[1]
        normalized = normalized.split("/", 1)[1] if "/" in normalized else normalized
    elif "@" in normalized and ":" in normalized:
        normalized = normalized.split(":", 1)[1]
    return normalized


def _validate_manifest_file_reference(
    root: Path,
    item: object,
    label: str,
    claim_id: str,
    report: ValidationReport,
) -> None:
    """Verify one manifest file reference against repository-local bytes."""
    if not isinstance(item, dict):
        report.error(f"{claim_id}: {label} must be a JSON object")
        return

    checksum = str(item.get("sha256", "")).strip()
    if not re.fullmatch(r"[0-9a-fA-F]{64}", checksum):
        report.error(f"{claim_id}: {label} requires a SHA-256 checksum")

    path_value = str(item.get("path", "")).strip()
    if not path_value:
        report.error(f"{claim_id}: {label} requires a repository-relative path")
        return
    resolved = _safe_repo_path(root, path_value, f"{claim_id} {label} path", report)
    if resolved is None:
        return
    if not resolved.is_file():
        report.error(f"{claim_id}: {label} file does not exist: {path_value}")
        return
    if not re.fullmatch(r"[0-9a-fA-F]{64}", checksum):
        return

    try:
        actual_checksum = _sha256_file(resolved)
    except OSError as exc:
        report.error(f"{claim_id}: could not checksum {label} file {path_value}: {exc}")
        return
    if actual_checksum.lower() != checksum.lower():
        report.error(
            f"{claim_id}: {label} checksum mismatch for {path_value}: "
            f"manifest has {checksum.lower()}, actual is {actual_checksum}"
        )


def _validate_computational_manifest(
    root: Path,
    path: Path,
    expected_run_id: str,
    claim_id: str,
    report: ValidationReport,
) -> None:
    """Apply the minimum Phase-1 eligibility checks to a run manifest."""
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        report.error(f"{claim_id}: could not read run manifest: {exc}")
        return
    if not isinstance(manifest, dict):
        report.error(f"{claim_id}: run manifest must contain a JSON object")
        return

    if manifest.get("run_id") != expected_run_id:
        report.error(f"{claim_id}: run_id does not match its manifest")

    git_info = manifest.get("git")
    if not isinstance(git_info, dict):
        report.error(f"{claim_id}: run manifest requires git provenance")
    else:
        repository = str(git_info.get("repository", "")).strip()
        if not repository:
            report.error(f"{claim_id}: run manifest requires a Git repository identifier")
        branch = str(git_info.get("branch", "")).strip()
        if branch != "shawn_dev":
            report.error(f"{claim_id}: run manifest Git branch must be shawn_dev")
        commit = str(git_info.get("commit", ""))
        if not re.fullmatch(r"[0-9a-fA-F]{40}", commit):
            report.error(f"{claim_id}: run manifest requires a full 40-character Git SHA")
        elif not (root / ".git").exists():
            report.error(f"{claim_id}: cannot verify the manifest Git SHA outside a Git repository")
        else:
            commit_result = subprocess.run(
                ["git", "cat-file", "-e", f"{commit}^{{commit}}"],
                cwd=root,
                check=False,
                capture_output=True,
            )
            if commit_result.returncode != 0:
                report.error(
                    f"{claim_id}: manifest Git SHA does not identify a commit in this repository"
                )
            elif branch == "shawn_dev":
                branch_ref = None
                for candidate_ref in (
                    "refs/heads/shawn_dev",
                    "refs/remotes/origin/shawn_dev",
                ):
                    ref_result = subprocess.run(
                        ["git", "show-ref", "--verify", "--quiet", candidate_ref],
                        cwd=root,
                        check=False,
                        capture_output=True,
                    )
                    if ref_result.returncode == 0:
                        branch_ref = candidate_ref
                        break
                if branch_ref is None:
                    report.error(
                        f"{claim_id}: cannot verify manifest commit against shawn_dev"
                    )
                else:
                    ancestry_result = subprocess.run(
                        ["git", "merge-base", "--is-ancestor", commit, branch_ref],
                        cwd=root,
                        check=False,
                        capture_output=True,
                    )
                    if ancestry_result.returncode != 0:
                        report.error(
                            f"{claim_id}: manifest Git SHA is not reachable from shawn_dev"
                        )

            remote_result = subprocess.run(
                ["git", "config", "--get", "remote.origin.url"],
                cwd=root,
                check=False,
                capture_output=True,
                text=True,
            )
            if remote_result.returncode != 0 or not remote_result.stdout.strip():
                report.error(f"{claim_id}: cannot verify manifest repository without origin")
            elif repository and _normalized_repository_identifier(
                repository
            ) != _normalized_repository_identifier(remote_result.stdout):
                report.error(
                    f"{claim_id}: manifest repository does not match remote.origin.url"
                )
        if git_info.get("dirty") is not False:
            report.error(f"{claim_id}: manuscript evidence must come from a clean Git commit")

    configuration = manifest.get("configuration")
    _validate_manifest_file_reference(
        root, configuration, "configuration", claim_id, report
    )

    for collection_name in ("inputs", "outputs"):
        collection = manifest.get(collection_name)
        if not isinstance(collection, list) or not collection:
            report.error(f"{claim_id}: run manifest requires at least one {collection_name} entry")
            continue
        for index, item in enumerate(collection):
            _validate_manifest_file_reference(
                root,
                item,
                f"{collection_name}[{index}]",
                claim_id,
                report,
            )

    validation = manifest.get("validation")
    if not isinstance(validation, dict) or validation.get("passed") is not True:
        report.error(f"{claim_id}: run manifest requires a passed validation decision")
    elif not validation.get("reviewer") or not validation.get("reviewed_at"):
        report.error(f"{claim_id}: validation reviewer and date are required")

    metrics = manifest.get("metrics")
    if not isinstance(metrics, dict) or not metrics:
        report.error(f"{claim_id}: supported computational evidence requires metrics")
    else:
        for metric_name, metric in metrics.items():
            if not isinstance(metric, dict) or "value" not in metric or not metric.get("unit"):
                report.error(
                    f"{claim_id}: metric '{metric_name}' requires value and unit"
                )


def _validate_tracked_file_policy(root: Path, report: ValidationReport) -> list[Path]:
    tracked = _git_tracked_files(root, report)
    forbidden_prefixes = (
        "paper/_output/",
        "paper/.quarto/",
        "literature/pdfs/",
        "tmp/pdfs/",
        "research/private/",
        ".research-runs/",
    )
    zotero_prefixes = (
        ".zotero/",
        "literature/zotero/",
        "literature/attachments/",
    )

    for path in tracked:
        relative = path.relative_to(root).as_posix()
        if relative.startswith(forbidden_prefixes) or (
            relative.startswith("workflow/profiles/") and relative.endswith(".local.yaml")
        ):
            report.error(f"generated, private, or machine-local file is tracked: {relative}")
        if relative.startswith(zotero_prefixes) or path.name.startswith("zotero.sqlite"):
            report.error(
                f"Zotero database or attachment path must not be tracked: {relative}"
            )
        try:
            if path.is_file() and path.stat().st_size > MAX_TRACKED_FILE_BYTES:
                report.error(
                    f"tracked file exceeds 10 MiB and must use external storage: {relative}"
                )
        except OSError as exc:
            report.error(f"could not inspect tracked file {relative}: {exc}")
    return tracked


def validate_repository(
    root: Path = ROOT,
    check_git: bool = True,
    check_local_pdfs: bool = False,
) -> ValidationReport:
    """Validate the research-writing contract rooted at *root*."""
    report = ValidationReport()
    try:
        root = root.resolve(strict=True)
    except (FileNotFoundError, OSError) as exc:
        report.error(f"repository root is unavailable: {exc}")
        return report
    required = [
        MANUSCRIPT_REL,
        BIBLIOGRAPHY_REL,
        EVIDENCE_REL,
        REFERENCES_REL,
        INVENTORY_REL,
        RELEVANCE_REL,
    ]
    for relative in required:
        if not (root / relative).is_file():
            report.error(f"missing required file: {relative.as_posix()}")
    if report.errors:
        return report

    manuscript_text = (root / MANUSCRIPT_REL).read_text(encoding="utf-8")
    bibliography_text = (root / BIBLIOGRAPHY_REL).read_text(encoding="utf-8")
    bib_entries = _parse_bib_entries(bibliography_text, report)
    bib_key_list = [citation_key for citation_key, _ in bib_entries]
    bib_keys = set(bib_key_list)
    bib_by_key = {citation_key: fields for citation_key, fields in bib_entries}
    cite_keys = set(CITE_KEY.findall(manuscript_text))

    duplicate_bib_keys = _duplicates(bib_key_list)
    if duplicate_bib_keys:
        report.error("duplicate bibliography keys: " + ", ".join(duplicate_bib_keys))
    for citation_key, fields in bib_entries:
        forbidden_fields = sorted(FORBIDDEN_BIB_FIELDS & fields.keys())
        if forbidden_fields:
            report.error(
                f"{citation_key}: BibTeX contains forbidden field(s): "
                + ", ".join(forbidden_fields)
            )
        for field_name, value in fields.items():
            if LOCAL_PATH.search(value):
                report.error(
                    f"{citation_key}: BibTeX {field_name} contains a local path"
                )
    missing_citations = sorted(cite_keys - bib_keys)
    if missing_citations:
        report.error("undefined citation keys: " + ", ".join(missing_citations))

    evidence_rows = _read_csv(
        root / EVIDENCE_REL, EVIDENCE_COLUMNS, "evidence.csv", report
    )
    reference_rows = _read_csv(
        root / REFERENCES_REL, REFERENCE_COLUMNS, "references.csv", report
    )
    inventory_rows = _read_csv(
        root / INVENTORY_REL, INVENTORY_COLUMNS, "artifact_inventory.csv", report)
    relevance_rows = _read_csv(
        root / RELEVANCE_REL, RELEVANCE_COLUMNS, "relevance.csv", report
    )

    reference_keys = [
        row.get("citation_key", "").strip() for row in reference_rows
    ]
    if any(not item for item in reference_keys):
        report.error("references.csv contains an empty citation_key")
    duplicate_reference_keys = _duplicates([item for item in reference_keys if item])
    if duplicate_reference_keys:
        report.error(
            "duplicate reference citation keys: "
            + ", ".join(duplicate_reference_keys)
        )

    reference_by_key: dict[str, dict[str, str]] = {}
    verified_reference_keys: set[str] = set()
    seen_dois: dict[str, str] = {}
    for row_number, row in enumerate(reference_rows, start=2):
        citation_key = row.get("citation_key", "").strip() or f"row {row_number}"
        reference_status = row.get("reference_status", "").strip()
        doi_status = row.get("doi_status", "").strip()
        doi = row.get("doi", "").strip()
        authors = row.get("authors", "").strip()

        if citation_key != f"row {row_number}":
            reference_by_key[citation_key] = row
        if reference_status not in REFERENCE_STATUSES:
            report.error(
                f"{citation_key}: unsupported reference status '{reference_status}'"
            )
        if doi_status not in DOI_STATUSES:
            report.error(f"{citation_key}: unsupported doi_status '{doi_status}'")

        normalized_doi = doi.lower()
        if normalized_doi:
            previous_key = seen_dois.get(normalized_doi)
            if previous_key is not None:
                report.error(
                    f"duplicate reference DOI '{doi}' for {previous_key} and {citation_key}"
                )
            else:
                seen_dois[normalized_doi] = citation_key

        for url_field in ("publisher_url", "metadata_source_url"):
            value = row.get(url_field, "").strip()
            if value and not value.startswith("https://"):
                report.error(f"{citation_key}: {url_field} must use HTTPS")

        year = row.get("year", "").strip()
        if year and not re.fullmatch(r"\d{4}", year):
            report.error(f"{citation_key}: year must use YYYY format")
        verified_at = row.get("verified_at", "").strip()
        if verified_at and not re.fullmatch(r"\d{4}-\d{2}-\d{2}", verified_at):
            report.error(f"{citation_key}: verified_at must use YYYY-MM-DD format")

        if reference_status == "verified":
            verified_reference_keys.add(citation_key)
            for field_name in (
                "citation_key",
                "title",
                "authors",
                "year",
                "venue",
                "publisher_url",
                "metadata_source_url",
                "verified_by",
                "verified_at",
            ):
                if not row.get(field_name, "").strip():
                    report.error(
                        f"{citation_key}: {field_name} is required for a verified reference"
                    )
            if doi_status not in {"verified", "not_assigned"}:
                report.error(
                    f"{citation_key}: verified reference requires a resolved doi_status"
                )
            if doi_status == "verified" and not re.fullmatch(
                r"10\.\d{4,9}/\S+", doi, re.IGNORECASE
            ):
                report.error(
                    f"{citation_key}: doi_status verified requires a canonical DOI"
                )
            if doi_status == "not_assigned" and doi:
                report.error(
                    f"{citation_key}: doi_status not_assigned requires an empty DOI"
                )
            if re.search(r"\b(?:and\s+others|et\s+al\.?)\b", authors, re.IGNORECASE):
                report.error(
                    f"{citation_key}: verified reference requires a complete ordered author list"
                )
        elif reference_status == "metadata_partial":
            for field_name in ("metadata_source_url", "notes"):
                if not row.get(field_name, "").strip():
                    report.error(
                        f"{citation_key}: {field_name} is required for metadata_partial"
                    )
            if row.get("verified_by", "").strip() or row.get(
                "verified_at", ""
            ).strip():
                report.error(
                    f"{citation_key}: metadata_partial must leave verified_by and "
                    "verified_at blank until explicit human approval"
                )
        elif reference_status == "rejected" and not row.get("notes", "").strip():
            report.error(f"{citation_key}: rejected reference requires a reason in notes")

    if bib_keys != verified_reference_keys:
        missing_from_bib = sorted(verified_reference_keys - bib_keys)
        unverified_in_bib = sorted(bib_keys - verified_reference_keys)
        if missing_from_bib:
            report.error(
                "verified references missing from paper/references.bib: "
                + ", ".join(missing_from_bib)
            )
        if unverified_in_bib:
            report.error(
                "paper/references.bib contains non-verified references: "
                + ", ".join(unverified_in_bib)
            )

    for citation_key in sorted(verified_reference_keys & bib_keys):
        reference_row = reference_by_key[citation_key]
        bib_fields = bib_by_key[citation_key]
        for field_name in ("title", "author", "year"):
            if not bib_fields.get(field_name, "").strip():
                report.error(f"{citation_key}: BibTeX {field_name} is required")

        doi_status = reference_row.get("doi_status", "").strip()
        bib_doi = _normalize_doi(bib_fields.get("doi", ""))
        registry_doi = _normalize_doi(reference_row.get("doi", ""))
        if doi_status == "verified":
            if not bib_doi:
                report.error(f"{citation_key}: BibTeX DOI is required")
            elif bib_doi != registry_doi:
                report.error(
                    f"{citation_key}: BibTeX DOI does not match references.csv"
                )
        elif doi_status == "not_assigned" and bib_doi:
            report.error(
                f"{citation_key}: BibTeX contains a DOI not assigned in references.csv"
            )

        bib_year = bib_fields.get("year", "").strip()
        registry_year = reference_row.get("year", "").strip()
        if bib_year and bib_year != registry_year:
            report.error(
                f"{citation_key}: BibTeX year does not match references.csv"
            )

        bib_title = bib_fields.get("title", "").strip()
        registry_title = reference_row.get("title", "").strip()
        if bib_title and _normalize_bibliographic_title(
            bib_title
        ) != _normalize_bibliographic_title(registry_title):
            report.error(
                f"{citation_key}: BibTeX title does not match references.csv"
            )

    relevance_keys = [row.get("citation_key", "").strip() for row in relevance_rows]
    if any(not item for item in relevance_keys):
        report.error("relevance.csv contains an empty citation_key")
    duplicate_relevance_keys = _duplicates([item for item in relevance_keys if item])
    if duplicate_relevance_keys:
        report.error(
            "duplicate relevance citation keys: " + ", ".join(duplicate_relevance_keys)
        )
    if set(relevance_keys) != set(reference_keys):
        missing_relevance = sorted(set(reference_keys) - set(relevance_keys))
        unknown_relevance = sorted(set(relevance_keys) - set(reference_keys))
        if missing_relevance:
            report.error(
                "references missing from relevance.csv: " + ", ".join(missing_relevance)
            )
        if unknown_relevance:
            report.error(
                "relevance.csv contains unknown references: "
                + ", ".join(unknown_relevance)
            )

    known_claim_ids = {
        row.get("claim_id", "").strip()
        for row in evidence_rows
        if row.get("claim_id", "").strip()
    }
    relevance_by_key: dict[str, dict[str, str]] = {}
    priority_ranks: list[int] = []
    phase_selected_keys: set[str] = set()
    phase_2b_1_keys: set[str] = set()
    phase_2b_1_coverage: list[str] = []
    frozen_contract_active = PHASE_2B_1_KEYS <= set(reference_keys)
    for row_number, row in enumerate(relevance_rows, start=2):
        citation_key = row.get("citation_key", "").strip() or f"row {row_number}"
        if citation_key != f"row {row_number}":
            relevance_by_key[citation_key] = row
        try:
            rank = int(row.get("priority_rank", ""))
            if rank < 1:
                raise ValueError
            priority_ranks.append(rank)
        except ValueError:
            report.error(f"{citation_key}: priority_rank must be a positive integer")
        category = row.get("relevance_category", "").strip()
        if category not in RELEVANCE_CATEGORIES:
            report.error(f"{citation_key}: unsupported relevance_category '{category}'")
        repo_use_status = row.get("repo_use_status", "").strip()
        if repo_use_status not in REPO_USE_STATUSES:
            report.error(f"{citation_key}: unsupported repo_use_status '{repo_use_status}'")
        selection_wave = row.get("selection_wave", "").strip()
        is_extraction_wave = EXTRACTION_WAVE.fullmatch(selection_wave) is not None
        if selection_wave not in NON_SELECTED_WAVES and not is_extraction_wave:
            report.error(f"{citation_key}: unsupported selection_wave '{selection_wave}'")
        coverage_area = row.get("coverage_area", "").strip()
        if not coverage_area:
            report.error(f"{citation_key}: coverage_area is required")
        if not row.get("rationale", "").strip():
            report.error(f"{citation_key}: rationale is required")
        if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", row.get("assessed_at", "").strip()):
            report.error(f"{citation_key}: assessed_at must use YYYY-MM-DD")
        raw_target_claims = row.get("target_claim_ids", "")
        target_tokens = [item.strip() for item in raw_target_claims.split("|")]
        if not raw_target_claims.strip():
            report.error(f"{citation_key}: target_claim_ids is required")
            target_claims: list[str] = []
        else:
            if any(not item for item in target_tokens):
                report.error(
                    f"{citation_key}: target_claim_ids cannot contain empty tokens"
                )
            target_claims = [item for item in target_tokens if item]
        duplicate_target_claims = _duplicates(target_claims)
        if duplicate_target_claims:
            report.error(
                f"{citation_key}: duplicate target_claim_ids: "
                + ", ".join(duplicate_target_claims)
            )
        unknown_target_claims = sorted(set(target_claims) - known_claim_ids)
        if unknown_target_claims:
            report.error(
                f"{citation_key}: target_claim_ids not found in evidence.csv: "
                + ", ".join(unknown_target_claims)
            )
        if is_extraction_wave:
            phase_selected_keys.add(citation_key)
        if selection_wave == "phase_2b_1":
            phase_2b_1_keys.add(citation_key)
            phase_2b_1_coverage.append(coverage_area)
    if sorted(priority_ranks) != list(range(1, len(relevance_rows) + 1)):
        report.error("relevance priority_rank values must be unique and contiguous from 1")
    if frozen_contract_active:
        if phase_2b_1_keys != PHASE_2B_1_KEYS:
            report.error(
                "Phase 2B-1 selected keys must remain: "
                + ", ".join(sorted(PHASE_2B_1_KEYS))
            )
        for frozen_key, frozen_coverage in PHASE_2B_1_FROZEN.items():
            frozen_row = relevance_by_key.get(frozen_key)
            if frozen_row is None:
                report.error(
                    f"Phase 2B-1 frozen relevance row is missing: {frozen_key}"
                )
                continue
            if frozen_row.get("selection_wave", "").strip() != "phase_2b_1":
                report.error(
                    f"{frozen_key}: frozen selection_wave must remain phase_2b_1"
                )
            if frozen_row.get("coverage_area", "").strip() != frozen_coverage:
                report.error(
                    f"{frozen_key}: frozen coverage_area must remain {frozen_coverage}"
                )
    elif phase_2b_1_keys:
        if len(phase_2b_1_keys) != 4:
            report.error("Phase 2B-1 must select exactly four references")
        if set(phase_2b_1_coverage) != PHASE_2B_1_COVERAGE:
            report.error(
                "Phase 2B-1 selection must cover exactly: "
                + ", ".join(sorted(PHASE_2B_1_COVERAGE))
            )
        if phase_2b_1_keys != PHASE_2B_1_KEYS:
            report.error(
                "Phase 2B-1 selected keys must remain: "
                + ", ".join(sorted(PHASE_2B_1_KEYS))
            )

    evidence_ids = [row.get("evidence_id", "").strip() for row in evidence_rows]
    if any(not item for item in evidence_ids):
        report.error("evidence.csv contains an empty evidence_id")
    duplicate_evidence_ids = _duplicates([item for item in evidence_ids if item])
    if duplicate_evidence_ids:
        report.error("duplicate evidence IDs: " + ", ".join(duplicate_evidence_ids))

    claim_ids = [row.get("claim_id", "").strip() for row in evidence_rows]
    if any(not item for item in claim_ids):
        report.error("evidence.csv contains an empty claim_id")
    evidence_by_claim: dict[str, list[dict[str, str]]] = {}
    for row in evidence_rows:
        claim_id = row.get("claim_id", "").strip()
        if claim_id:
            evidence_by_claim.setdefault(claim_id, []).append(row)

    relationship_keys: list[str] = []
    relationships_by_claim_key: dict[tuple[str, str], dict[str, str]] = {}
    for row in evidence_rows:
        claim_id = row.get("claim_id", "").strip()
        citation_key = row.get("citation_key", "").strip()
        source_identity = (
            citation_key
            or row.get("evidence_path", "").strip()
            or row.get("run_id", "").strip()
            or row.get("source_type", "").strip()
        )
        if claim_id and source_identity:
            relationship_keys.append(f"{claim_id}\x1f{source_identity}")
        if claim_id and citation_key:
            relationships_by_claim_key[(claim_id, citation_key)] = row
    duplicate_relationships = _duplicates(relationship_keys)
    if duplicate_relationships:
        formatted = [item.replace("\x1f", " -> ") for item in duplicate_relationships]
        report.error("duplicate claim-source relationships: " + ", ".join(formatted))
    supported_computational_paths: list[tuple[str, str]] = []

    for row_number, row in enumerate(evidence_rows, start=2):
        claim_id = row.get("claim_id", "").strip() or f"row {row_number}"
        status = row.get("status", "").strip()
        source_type = row.get("source_type", "").strip()
        citation_key = row.get("citation_key", "").strip()
        evidence_path = row.get("evidence_path", "").strip()

        if status not in EVIDENCE_STATUSES:
            report.error(f"{claim_id}: unsupported evidence status '{status}'")
        if source_type not in SOURCE_TYPES:
            report.error(f"{claim_id}: unsupported source_type '{source_type}'")
        if not row.get("claim", "").strip():
            report.error(f"{claim_id}: claim text is required")
        if not row.get("owner", "").strip():
            report.error(f"{claim_id}: owner is required")
        reference_row = reference_by_key.get(citation_key) if citation_key else None
        if citation_key and reference_row is None:
            report.error(f"{claim_id}: citation key is missing from references.csv: '{citation_key}'")

        resolved_evidence = _safe_repo_path(
            root, evidence_path, f"{claim_id} evidence_path", report
        )
        if resolved_evidence is not None and not resolved_evidence.exists():
            report.error(f"{claim_id}: evidence_path does not exist: {evidence_path}")

        if status in {"supported", "partial"} and source_type in LITERATURE_SOURCE_TYPES:
            for field_name in ("citation_key", "source_locator", "verified_by", "verified_at"):
                if not row.get(field_name, "").strip():
                    report.error(
                        f"{claim_id}: {field_name} is required for {status} literature evidence"
                    )
            if reference_row is not None and reference_row.get(
                "reference_status", ""
            ).strip() != "verified":
                report.error(
                    f"{claim_id}: {status} literature evidence requires a verified reference"
                )
        if status in {"supported", "partial"} and source_type == "computational_artifact":
            for field_name in ("evidence_path", "run_id", "verified_by", "verified_at"):
                if not row.get(field_name, "").strip():
                    report.error(
                        f"{claim_id}: {field_name} is required for {status} computational evidence"
                    )
            if evidence_path and not evidence_path.endswith(".json"):
                report.error(
                    f"{claim_id}: computational evidence_path must identify a JSON run manifest"
                )
            if evidence_path:
                supported_computational_paths.append((claim_id, evidence_path))
            if resolved_evidence is not None and resolved_evidence.is_file():
                _validate_computational_manifest(
                    root,
                    resolved_evidence,
                    row.get("run_id", "").strip(),
                    claim_id,
                    report,
                )
        if status in {"supported", "partial"} and source_type == "repository":
            if not evidence_path:
                report.error(f"{claim_id}: repository evidence requires evidence_path")

    notes_dir = root / NOTES_REL
    card_paths = sorted(
        path for path in notes_dir.glob("*.md") if path.name != "README.md"
    ) if notes_dir.is_dir() else []
    cards_by_key: dict[str, dict[str, object]] = {}
    card_key_counts: dict[str, int] = {}
    for card_path in card_paths:
        filename_key = card_path.stem
        card = _read_evidence_card(card_path, filename_key, report)
        if card is None:
            continue
        key = str(card.get("citation_key", "")).strip()
        card_key_counts[key] = card_key_counts.get(key, 0) + 1
        if key in cards_by_key:
            report.error(f"duplicate evidence card for citation key: {key}")
        else:
            cards_by_key[key] = card
        _validate_evidence_card(
            card,
            filename_key,
            reference_by_key,
            relationships_by_claim_key,
            report,
        )

    if frozen_contract_active:
        for frozen_key in PHASE_2B_1_KEYS:
            if card_key_counts.get(frozen_key, 0) != 1:
                report.error(
                    f"{frozen_key}: Phase 2B-1 requires exactly one evidence card"
                )

    if set(cards_by_key) != phase_selected_keys:
        missing_cards = sorted(phase_selected_keys - set(cards_by_key))
        unselected_cards = sorted(set(cards_by_key) - phase_selected_keys)
        if missing_cards:
            report.error(
                "phase-selected references missing evidence cards: "
                + ", ".join(missing_cards)
            )
        if unselected_cards:
            report.error(
                "orphan evidence cards exist without a phase-selected relevance row: "
                + ", ".join(unselected_cards)
            )

    for key, card in cards_by_key.items():
        card_status = str(card.get("extraction_status", "")).strip()
        expected_path = (NOTES_REL / f"{key}.md").as_posix()
        matching_rows = [
            row
            for row in evidence_rows
            if row.get("citation_key", "").strip() == key
            and row.get("evidence_path", "").strip() == expected_path
        ]
        if not matching_rows:
            report.error(
                f"{key}: evidence.csv must link the evidence card using {expected_path}"
            )
        for row in matching_rows:
            claim_id = row.get("claim_id", "").strip()
            ledger_status = row.get("status", "").strip()
            if card_status != "human_confirmed" and ledger_status != "pending":
                report.error(
                    f"{claim_id}: machine evidence card must not upgrade evidence status"
                )
            if card_status != "human_confirmed" or ledger_status == "pending":
                for field_name in ("source_locator", "verified_by", "verified_at"):
                    if row.get(field_name, "").strip():
                        report.error(
                            f"{claim_id}: unpromoted evidence card must leave "
                            f"{field_name} blank"
                        )
        relevance = relevance_by_key.get(key, {})
        card_coverage = str(card.get("coverage_area", "")).strip()
        if card_coverage != relevance.get("coverage_area", "").strip():
            report.error(f"{key}: card coverage_area does not match relevance.csv")
        target_claims = {
            item.strip()
            for item in relevance.get("target_claim_ids", "").split("|")
            if item.strip()
        }
        relationship_claims = {
            row.get("claim_id", "").strip()
            for row in matching_rows
            if row.get("claim_id", "").strip()
        }
        if not relationship_claims <= target_claims:
            report.error(f"{key}: evidence-card relationship is outside target_claim_ids")
        if check_local_pdfs and card_status in {
            "extracted_pending_human",
            "human_confirmed",
        }:
            reference = reference_by_key.get(key)
            if reference is not None:
                _validate_local_pdf(root, card, reference, report)

    evidence_citations = {
        row.get("citation_key", "").strip() for row in evidence_rows if row.get("citation_key", "").strip()
    }
    unmapped_citations = sorted(cite_keys - evidence_citations)
    if unmapped_citations:
        report.error(
            "manuscript citations missing from evidence.csv: " + ", ".join(unmapped_citations)
        )

    for block in re.split(r"\n\s*\n", manuscript_text):
        block_citations = set(CITE_KEY.findall(block))
        markers = CLAIM_MARKER.findall(block)
        if block_citations and not markers:
            report.error(
                "citation block requires a CLAIM marker: "
                + ", ".join(sorted(block_citations))
            )
            continue
        if len(markers) > 1:
            report.error("a manuscript block may contain only one CLAIM marker")
            continue
        if not markers:
            continue

        claim_id = markers[0]
        claim_rows = evidence_by_claim.get(claim_id, [])
        if not claim_rows:
            report.error(f"CLAIM marker is missing from evidence.csv: {claim_id}")
            continue
        for citation_key in block_citations:
            relationship = relationships_by_claim_key.get((claim_id, citation_key))
            if relationship is None:
                report.error(
                    f"{claim_id} -> {citation_key}: claim-source relationship is missing"
                )
            elif relationship.get("status", "").strip() not in {"supported", "partial"}:
                report.error(
                    f"{claim_id} -> {citation_key}: citation does not have supported or partial evidence"
                )

        if not block_citations and not any(
            row.get("status", "").strip() in {"supported", "partial", "planned"}
            for row in claim_rows
        ):
            has_todo = re.search(
                rf"TODO-EVIDENCE\s+{re.escape(claim_id)}\b", block
            )
            has_planned = re.search(r"\bPLANNED\b", block)
            if not has_todo and not has_planned:
                report.error(
                    f"{claim_id}: pending claim block requires TODO-EVIDENCE or PLANNED"
                )

    todo_ids = TODO_ID.findall(manuscript_text)
    for todo_id in sorted(set(todo_ids)):
        rows = evidence_by_claim.get(todo_id, [])
        if not rows:
            report.error(f"TODO-EVIDENCE ID is missing from evidence.csv: {todo_id}")
        elif any(row.get("status", "").strip() == "supported" for row in rows):
            report.error(f"TODO-EVIDENCE {todo_id} is still marked supported")

    artifact_ids = [row.get("artifact_id", "").strip() for row in inventory_rows]
    if any(not item for item in artifact_ids):
        report.error("artifact_inventory.csv contains an empty artifact_id")
    duplicate_artifact_ids = _duplicates([item for item in artifact_ids if item])
    if duplicate_artifact_ids:
        report.error("duplicate artifact IDs: " + ", ".join(duplicate_artifact_ids))

    conflict_inventory_paths: set[str] = set()
    eligible_inventory_paths: set[str] = set()
    for row_number, row in enumerate(inventory_rows, start=2):
        artifact_id = row.get("artifact_id", "").strip() or f"row {row_number}"
        status = row.get("status", "").strip()
        provenance = row.get("provenance_status", "").strip()
        eligible = row.get("manuscript_eligible", "").strip()
        artifact_path = row.get("path", "").strip()

        if status not in INVENTORY_STATUSES:
            report.error(f"{artifact_id}: unsupported artifact status '{status}'")
        if provenance not in PROVENANCE_STATUSES:
            report.error(f"{artifact_id}: unsupported provenance_status '{provenance}'")
        if eligible not in ELIGIBILITY_VALUES:
            report.error(f"{artifact_id}: manuscript_eligible must be yes or no")
        if status in PATH_REQUIRED_STATUSES and not artifact_path:
            report.error(f"{artifact_id}: path is required for artifact status {status}")

        resolved_artifact = _safe_repo_path(
            root, artifact_path, f"{artifact_id} path", report
        )
        if resolved_artifact is not None and status in PATH_REQUIRED_STATUSES:
            if not resolved_artifact.exists():
                report.error(f"{artifact_id}: artifact path does not exist: {artifact_path}")
        if eligible == "yes" and not (
            status == "completed" and provenance == "manifest_complete"
        ):
            report.error(
                f"{artifact_id}: manuscript eligibility requires completed status and manifest_complete provenance"
            )
        if eligible == "yes" and artifact_path:
            eligible_inventory_paths.add(artifact_path)
        if row.get("artifact_type", "").strip() == "conflicted_input":
            if not (status == "invalid" and eligible == "no"):
                report.error(
                    f"{artifact_id}: conflicted inputs must be invalid and manuscript-ineligible"
                )
            if artifact_path:
                conflict_inventory_paths.add(artifact_path)

    tracked = _validate_tracked_file_policy(root, report) if check_git else []
    tracked_paths = {
        path.relative_to(root).as_posix()
        for path in tracked
        if path.exists()
    }
    for claim_id, manifest_path in supported_computational_paths:
        if manifest_path not in eligible_inventory_paths:
            report.error(
                f"{claim_id}: run manifest is not marked manuscript-eligible in artifact inventory"
            )
        if check_git and manifest_path not in tracked_paths:
            report.error(f"{claim_id}: run manifest must be tracked by Git: {manifest_path}")

    marker_paths: set[str] = set()
    files_to_scan = tracked or [path for path in root.rglob("*") if path.is_file()]
    for path in files_to_scan:
        try:
            if path.stat().st_size > MAX_TRACKED_FILE_BYTES:
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        if CONFLICT_MARKER.search(text):
            marker_paths.add(path.relative_to(root).as_posix())
    unquarantined_markers = sorted(marker_paths - conflict_inventory_paths)
    if unquarantined_markers:
        report.error(
            "files with conflict markers are not quarantined in artifact inventory: "
            + ", ".join(unquarantined_markers)
        )

    figure_dir = root / FIGURES_REL
    if figure_dir.is_dir():
        for figure_path in figure_dir.rglob("*"):
            if not figure_path.is_file() or figure_path.name == "README.md":
                continue
            relative = figure_path.relative_to(root).as_posix()
            if relative not in eligible_inventory_paths:
                report.error(
                    f"publication figure lacks an eligible artifact-inventory row: {relative}"
                )

    report.counts = {
        "citations_used": len(cite_keys),
        "bibliography_entries": len(bib_keys),
        "reference_rows": len(reference_rows),
        "relevance_rows": len(relevance_rows),
        "evidence_cards": len(cards_by_key),
        "evidence_rows": len(evidence_rows),
        "todo_evidence_ids": len(set(todo_ids)),
        "artifact_rows": len(inventory_rows),
        "eligible_artifacts": sum(
            row.get("manuscript_eligible", "").strip() == "yes" for row in inventory_rows
        ),
        "quarantined_conflict_files": len(marker_paths),
    }
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check-local-pdfs",
        action="store_true",
        help="validate ignored local PDFs, extracted text, checksums, and locator rendering",
    )
    arguments = parser.parse_args()
    report = validate_repository(ROOT, check_local_pdfs=arguments.check_local_pdfs)
    for warning in report.warnings:
        print(f"WARNING: {warning}", file=sys.stderr)
    for error in report.errors:
        print(f"ERROR: {error}", file=sys.stderr)

    if report.counts:
        print(
            f"Citations: {report.counts['citations_used']} used, "
            f"{report.counts['bibliography_entries']} available"
        )
        print(f"Evidence rows: {report.counts['evidence_rows']}")
        print(f"Reference rows: {report.counts['reference_rows']}")
        print(f"Relevance rows: {report.counts['relevance_rows']}")
        print(f"Evidence cards: {report.counts['evidence_cards']}")
        print(f"Tracked TODO-EVIDENCE IDs: {report.counts['todo_evidence_ids']}")
        print(
            f"Artifacts: {report.counts['artifact_rows']} inventoried, "
            f"{report.counts['eligible_artifacts']} manuscript-eligible"
        )
        print(
            "Quarantined conflict-marker files: "
            f"{report.counts['quarantined_conflict_files']}"
        )

    if report.errors:
        return 1
    print("Research manuscript checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
