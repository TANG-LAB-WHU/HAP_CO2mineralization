#!/usr/bin/env python3
"""Dependency-free checks for manuscript, evidence, and artifact integrity."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import sys
import unicodedata
from dataclasses import dataclass, field
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_REL = Path("paper/manuscript.qmd")
BIBLIOGRAPHY_REL = Path("paper/references.bib")
EVIDENCE_REL = Path("literature/evidence.csv")
REFERENCES_REL = Path("literature/references.csv")
INVENTORY_REL = Path("research/artifact_inventory.csv")
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

EVIDENCE_STATUSES = {"supported", "partial", "planned", "pending", "rejected"}
REFERENCE_STATUSES = {"verified", "metadata_partial", "candidate", "rejected"}
DOI_STATUSES = {"verified", "not_assigned", "pending"}
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


@dataclass
class ValidationReport:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    counts: dict[str, int] = field(default_factory=dict)

    def error(self, message: str) -> None:
        self.errors.append(message)

    def warning(self, message: str) -> None:
        self.warnings.append(message)


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
        resolved.relative_to(root.resolve())
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


def validate_repository(root: Path = ROOT, check_git: bool = True) -> ValidationReport:
    """Validate the research-writing contract rooted at *root*."""
    report = ValidationReport()
    required = [
        MANUSCRIPT_REL,
        BIBLIOGRAPHY_REL,
        EVIDENCE_REL,
        REFERENCES_REL,
        INVENTORY_REL,
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
    report = validate_repository(ROOT)
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
