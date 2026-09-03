#!/usr/bin/env python3
"""Dependency-free checks for manuscript, evidence, and artifact integrity."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_REL = Path("paper/manuscript.qmd")
BIBLIOGRAPHY_REL = Path("paper/references.bib")
EVIDENCE_REL = Path("literature/evidence.csv")
INVENTORY_REL = Path("research/artifact_inventory.csv")
FIGURES_REL = Path("paper/figures")

BIB_KEY = re.compile(r"@[A-Za-z]+\s*\{\s*([^,\s]+)\s*,")
CITE_KEY = re.compile(r"(?<![\w])@([A-Za-z0-9_:.-]+)")
TODO_ID = re.compile(r"^#{1,6}\s+TODO-EVIDENCE\s+([A-Z][A-Z0-9-]*)\b", re.MULTILINE)
CONFLICT_MARKER = re.compile(r"^(?:<<<<<<< .+|=======|>>>>>>> .+)$", re.MULTILINE)

EVIDENCE_COLUMNS = [
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

    for path in tracked:
        relative = path.relative_to(root).as_posix()
        if relative.startswith(forbidden_prefixes) or (
            relative.startswith("workflow/profiles/") and relative.endswith(".local.yaml")
        ):
            report.error(f"generated, private, or machine-local file is tracked: {relative}")
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
    required = [MANUSCRIPT_REL, BIBLIOGRAPHY_REL, EVIDENCE_REL, INVENTORY_REL]
    for relative in required:
        if not (root / relative).is_file():
            report.error(f"missing required file: {relative.as_posix()}")
    if report.errors:
        return report

    manuscript_text = (root / MANUSCRIPT_REL).read_text(encoding="utf-8")
    bibliography_text = (root / BIBLIOGRAPHY_REL).read_text(encoding="utf-8")
    bib_key_list = BIB_KEY.findall(bibliography_text)
    bib_keys = set(bib_key_list)
    cite_keys = set(CITE_KEY.findall(manuscript_text))

    duplicate_bib_keys = _duplicates(bib_key_list)
    if duplicate_bib_keys:
        report.error("duplicate bibliography keys: " + ", ".join(duplicate_bib_keys))
    missing_citations = sorted(cite_keys - bib_keys)
    if missing_citations:
        report.error("undefined citation keys: " + ", ".join(missing_citations))

    evidence_rows = _read_csv(
        root / EVIDENCE_REL, EVIDENCE_COLUMNS, "evidence.csv", report
    )
    inventory_rows = _read_csv(
        root / INVENTORY_REL, INVENTORY_COLUMNS, "artifact_inventory.csv", report)

    evidence_ids = [row.get("claim_id", "").strip() for row in evidence_rows]
    if any(not item for item in evidence_ids):
        report.error("evidence.csv contains an empty claim_id")
    duplicate_evidence_ids = _duplicates([item for item in evidence_ids if item])
    if duplicate_evidence_ids:
        report.error("duplicate claim IDs: " + ", ".join(duplicate_evidence_ids))
    evidence_by_id = {
        row.get("claim_id", "").strip(): row
        for row in evidence_rows
        if row.get("claim_id", "").strip()
    }
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
        if citation_key and citation_key not in bib_keys:
            report.error(f"{claim_id}: undefined citation key '{citation_key}'")

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

    todo_ids = TODO_ID.findall(manuscript_text)
    duplicate_todo_ids = _duplicates(todo_ids)
    if duplicate_todo_ids:
        report.error("duplicate TODO-EVIDENCE IDs: " + ", ".join(duplicate_todo_ids))
    for todo_id in todo_ids:
        row = evidence_by_id.get(todo_id)
        if row is None:
            report.error(f"TODO-EVIDENCE ID is missing from evidence.csv: {todo_id}")
        elif row.get("status", "").strip() == "supported":
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
        "evidence_rows": len(evidence_rows),
        "todo_evidence_ids": len(todo_ids),
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
