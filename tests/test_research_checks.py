from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

from scripts.check_manuscript import (
    EVIDENCE_COLUMNS,
    INVENTORY_COLUMNS,
    PHASE_2B_1_FROZEN,
    POPPLER_TIMEOUT_SECONDS,
    RELEVANCE_COLUMNS,
    ROOT as PROJECT_ROOT,
    validate_repository,
)
from scripts import check_manuscript as manuscript_checks


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


class ResearchValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self._write(
            "paper/manuscript.qmd",
            "# Draft\n\n<!-- CLAIM: C1 -->\nClaim [@seed].\n\n## TODO-EVIDENCE R1\n",
        )
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
            "  year = {2020},\n  doi = {10.1000/seed}\n}\n",
        )
        self.evidence_rows = [
            self._evidence_row(
                evidence_id="E001",
                claim_id="C1",
                status="supported",
                claim="Seed literature claim.",
                source_type="primary_publication",
                citation_key="seed",
                source_locator="p. 1, Results, paragraph 1",
                verified_by="Researcher",
                verified_at="2026-09-03",
            ),
            self._evidence_row(
                evidence_id="E002",
                claim_id="R1",
                status="pending",
                claim="Planned computational result.",
                source_type="computational_artifact",
            ),
        ]
        self.reference_rows = [
            self._reference_row(
                citation_key="seed",
                reference_status="verified",
                title="Seed",
                authors="Example, A.",
                year="2020",
                venue="Example Journal",
                doi="10.1000/seed",
                doi_status="verified",
                publisher_url="https://publisher.example/seed",
                metadata_source_url="https://doi.org/10.1000/seed",
                verified_by="Researcher",
                verified_at="2026-09-03",
            )
        ]
        self.relevance_rows = [
            self._relevance_row(
                priority_rank="1",
                citation_key="seed",
                coverage_area="test_claim",
                target_claim_ids="C1",
            )
        ]
        self.inventory_rows = [
            self._inventory_row(
                artifact_id="PLANNED",
                stage="publication",
                system="test",
                artifact_type="figure set",
                status="missing",
                provenance_status="none",
            )
        ]
        self._write_csv("literature/evidence.csv", EVIDENCE_COLUMNS, self.evidence_rows)
        self._write_csv(
            "literature/references.csv", REFERENCE_COLUMNS, self.reference_rows
        )
        self._write_csv(
            "literature/relevance.csv", RELEVANCE_COLUMNS, self.relevance_rows
        )
        self._write_csv(
            "research/artifact_inventory.csv", INVENTORY_COLUMNS, self.inventory_rows
        )
        (self.root / "paper/figures").mkdir(parents=True)
        self._write("paper/figures/README.md", "# Figures\n")

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def _write(self, relative: str, content: str) -> None:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    def _write_csv(
        self, relative: str, columns: list[str], rows: list[dict[str, str]]
    ) -> None:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def _evidence_row(**overrides: str) -> dict[str, str]:
        row = {column: "" for column in EVIDENCE_COLUMNS}
        row.update({"owner": "Researcher", "notes": "Test fixture."})
        row.update(overrides)
        return row

    @staticmethod
    def _inventory_row(**overrides: str) -> dict[str, str]:
        row = {column: "" for column in INVENTORY_COLUMNS}
        row.update({"manuscript_eligible": "no", "notes": "Test fixture."})
        row.update(overrides)
        return row

    @staticmethod
    def _reference_row(**overrides: str) -> dict[str, str]:
        row = {column: "" for column in REFERENCE_COLUMNS}
        row.update({"reference_status": "candidate", "doi_status": "pending"})
        row.update(overrides)
        return row

    @staticmethod
    def _relevance_row(**overrides: str) -> dict[str, str]:
        row = {column: "" for column in RELEVANCE_COLUMNS}
        row.update(
            {
                "relevance_category": "core_direct",
                "repo_use_status": "not_applicable",
                "selection_wave": "later",
                "rationale": "Test fixture relevance decision.",
                "assessed_at": "2026-09-05",
            }
        )
        row.update(overrides)
        return row

    def _validate(self):
        return validate_repository(self.root, check_git=False)

    def _git(self, *arguments: str) -> str:
        result = subprocess.run(
            ["git", *arguments],
            cwd=self.root,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(0, result.returncode, result.stderr)
        return result.stdout.strip()

    def _initialize_git(self) -> str:
        self._git("init", "-q", "-b", "shawn_dev")
        self._git("config", "user.name", "Research Check Tests")
        self._git("config", "user.email", "research-checks@example.invalid")
        self._git(
            "remote",
            "add",
            "origin",
            "https://github.com/TANG-LAB-WHU/HAP_CO2mineralization.git",
        )
        self._git("add", ".")
        self._git("commit", "-q", "-m", "test baseline")
        return self._git("rev-parse", "HEAD")

    def _sha256(self, relative: str) -> str:
        return hashlib.sha256((self.root / relative).read_bytes()).hexdigest()

    def _add_supported_manifest(self, commit: str) -> str:
        manifest_path = "workflow/runs/run-1.json"
        manifest = {
            "schema_version": 1,
            "run_id": "run-1",
            "status": "completed",
            "backend": "local",
            "stage": "analysis",
            "git": {
                "repository": "TANG-LAB-WHU/HAP_CO2mineralization",
                "branch": "shawn_dev",
                "commit": commit,
                "dirty": False,
            },
            "configuration": {
                "path": "run-data/config.txt",
                "sha256": self._sha256("run-data/config.txt"),
            },
            "inputs": [
                {
                    "path": "run-data/input.dat",
                    "sha256": self._sha256("run-data/input.dat"),
                }
            ],
            "outputs": [
                {
                    "path": "run-data/output.dat",
                    "sha256": self._sha256("run-data/output.dat"),
                }
            ],
            "metrics": {"energy": {"value": -1.0, "unit": "eV"}},
            "validation": {
                "passed": True,
                "reviewer": "Researcher",
                "reviewed_at": "2026-09-03",
            },
        }
        self._write(manifest_path, json.dumps(manifest, indent=2) + "\n")
        self.evidence_rows.append(
            self._evidence_row(
                evidence_id="E003",
                claim_id="RUN1",
                status="supported",
                claim="Test computational result.",
                source_type="computational_artifact",
                evidence_path=manifest_path,
                run_id="run-1",
                verified_by="Researcher",
                verified_at="2026-09-03",
            )
        )
        self.inventory_rows.append(
            self._inventory_row(
                artifact_id="RUN-1",
                stage="analysis",
                system="test",
                artifact_type="run manifest",
                path=manifest_path,
                status="completed",
                provenance_status="manifest_complete",
                manuscript_eligible="yes",
            )
        )
        self._rewrite_evidence()
        self._rewrite_inventory()
        return manifest_path

    def _rewrite_evidence(self) -> None:
        self._write_csv("literature/evidence.csv", EVIDENCE_COLUMNS, self.evidence_rows)

    def _rewrite_inventory(self) -> None:
        self._write_csv(
            "research/artifact_inventory.csv", INVENTORY_COLUMNS, self.inventory_rows
        )

    def _rewrite_references(self) -> None:
        known_keys = {row["citation_key"] for row in self.relevance_rows}
        for reference in self.reference_rows:
            key = reference.get("citation_key", "")
            if key and key not in known_keys:
                self.relevance_rows.append(
                    self._relevance_row(
                        priority_rank=str(len(self.relevance_rows) + 1),
                        citation_key=key,
                        coverage_area="test_claim",
                        target_claim_ids="C1",
                    )
                )
                known_keys.add(key)
        self._write_csv(
            "literature/references.csv", REFERENCE_COLUMNS, self.reference_rows
        )
        self._rewrite_relevance()

    def _rewrite_relevance(self) -> None:
        self._write_csv(
            "literature/relevance.csv", RELEVANCE_COLUMNS, self.relevance_rows
        )

    def _write_card(self, card: dict[str, object], filename_key: str = "seed") -> None:
        self._write(
            f"literature/notes/{filename_key}.md",
            f"# Evidence card: {filename_key}\n\n```json\n"
            + json.dumps(card, indent=2)
            + "\n```\n",
        )

    def _pending_card(
        self,
        *,
        selection_wave: str = "phase_2b_2",
        coverage_area: str = "test_claim",
        **overrides: object,
    ) -> dict[str, object]:
        self._write(
            "paper/manuscript.qmd",
            "# Draft\n\n<!-- CLAIM: C1 -->\n## TODO-EVIDENCE C1\n",
        )
        self.evidence_rows[0].update(
            {
                "status": "pending",
                "evidence_path": "literature/notes/seed.md",
                "source_locator": "",
                "verified_by": "",
                "verified_at": "",
            }
        )
        self._rewrite_evidence()
        self.relevance_rows[0].update(
            {
                "selection_wave": selection_wave,
                "coverage_area": coverage_area,
            }
        )
        self._rewrite_relevance()
        card: dict[str, object] = {
            "schema_version": 1,
            "citation_key": "seed",
            "doi": "10.1000/seed",
            "extraction_status": "extracted_pending_human",
            "pdf_filename": "seed.pdf",
            "pdf_sha256": "a" * 64,
            "landing_page_url": "https://publisher.example/seed",
            "pdf_source_url": "https://publisher.example/seed.pdf",
            "full_text_version": "version_of_record",
            "access_basis": "publisher_open_access",
            "relevance": "Test relevance.",
            "coverage_area": coverage_area,
            "blocker_reason": None,
            "machine_extracted_at": "2026-09-06T00:00:00Z",
            "human_confirmed": False,
            "human_reviewer": None,
            "human_reviewed_at": None,
            "identity_checks": {
                "title_in_pdf_text": True,
                "author_in_pdf_text": True,
                "doi_match_basis": "pdf_text",
                "doi_association_url": "https://publisher.example/seed",
            },
            "claim_assessments": [
                {
                    "claim_id": "C1",
                    "assessment": "supports",
                    "section_heading": "Results",
                    "pdf_page": 1,
                    "printed_page": "1",
                    "figure_table_equation_locator": "not_applicable",
                    "verbatim_excerpt": "Short source text.",
                    "faithful_paraphrase": "Narrow test paraphrase.",
                    "scope_limitations": "Test-only scope.",
                }
            ],
        }
        card.update(overrides)
        return card

    def _local_pdf_card(
        self, **overrides: object
    ) -> tuple[dict[str, object], Path]:
        card = self._pending_card(**overrides)
        pdf_path = self.root / "literature/pdfs/seed.pdf"
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        pdf_path.write_bytes(b"%PDF-1.4\nlocal evidence test fixture\n")
        card["pdf_sha256"] = hashlib.sha256(pdf_path.read_bytes()).hexdigest()
        self._write_card(card)
        return card, pdf_path

    def _install_frozen_phase_fixture(self) -> None:
        """Install the frozen four-card contract plus one future-wave card."""
        template = self._pending_card(selection_wave="phase_2b_2")
        self._write_card(template)
        claims = {
            "qomi2022mineralization": "C001",
            "astala2008hapwater": "C010",
            "nowicki2024capture": "C011",
            "kuhne2020cp2k": "C006",
        }
        for index, (key, coverage) in enumerate(
            PHASE_2B_1_FROZEN.items(), start=2
        ):
            doi = f"10.1000/{key}"
            claim_id = claims[key]
            self.reference_rows.append(
                self._reference_row(citation_key=key, doi=doi)
            )
            self.relevance_rows.append(
                self._relevance_row(
                    priority_rank=str(index),
                    citation_key=key,
                    coverage_area=coverage,
                    selection_wave="phase_2b_1",
                    target_claim_ids=claim_id,
                )
            )
            self.evidence_rows.append(
                self._evidence_row(
                    evidence_id=f"EF{index}",
                    claim_id=claim_id,
                    status="pending",
                    claim=f"Synthetic frozen claim for {key}.",
                    source_type="primary_publication",
                    evidence_path=f"literature/notes/{key}.md",
                    citation_key=key,
                )
            )
            card = json.loads(json.dumps(template))
            card.update(
                citation_key=key,
                doi=doi,
                pdf_filename=f"{key}.pdf",
                pdf_source_url=f"https://publisher.example/{key}.pdf",
                landing_page_url=f"https://publisher.example/{key}",
                coverage_area=coverage,
            )
            card["identity_checks"]["doi_association_url"] = (
                f"https://publisher.example/{key}"
            )
            card["claim_assessments"][0]["claim_id"] = claim_id
            self._write_card(card, filename_key=key)
        self._rewrite_references()
        self._rewrite_evidence()
        self._rewrite_relevance()

    @staticmethod
    def _assert_poppler_invocation(
        arguments: list[str], kwargs: dict[str, object]
    ) -> tuple[str, Path]:
        if not isinstance(arguments, list):
            raise AssertionError("Poppler command must be represented as a list")
        operation = Path(arguments[0]).name
        if operation not in {"pdfinfo", "pdftotext", "pdftoppm"}:
            raise AssertionError(f"unexpected command: {arguments}")
        if arguments[0] != f"/mock/{operation}":
            raise AssertionError(f"unexpected executable: {arguments[0]}")
        if kwargs.get("timeout") != POPPLER_TIMEOUT_SECONDS:
            raise AssertionError("Poppler timeout contract was not enforced")
        if kwargs.get("shell", False) is not False:
            raise AssertionError("Poppler must not invoke a shell")
        source_index = 1 if operation == "pdfinfo" else -2
        source_path = Path(arguments[source_index])
        if source_path.name != "source.pdf" or not source_path.is_file():
            raise AssertionError(f"Poppler did not receive a stable snapshot: {source_path}")
        if "literature/pdfs" in source_path.as_posix():
            raise AssertionError("Poppler received the original untrusted PDF pathname")
        if operation == "pdftoppm" and Path(arguments[-1]).parent != source_path.parent:
            raise AssertionError("render destination is outside the private snapshot directory")
        return operation, source_path

    @staticmethod
    def _successful_poppler(
        page_text_by_page: dict[int, str] | None = None,
        *,
        page_count: int = 2,
        front_matter_text: str | None = None,
        full_document_text: str | None = None,
        pdfinfo_metadata: str = "",
        expected_pdf_bytes: bytes | None = None,
        observed_snapshot_paths: list[Path] | None = None,
    ):
        if page_text_by_page is None:
            page_text_by_page = {1: "Results\nShort source text.\n"}
        if front_matter_text is None:
            front_matter_text = "Seed\nExample, A.\nDOI: 10.1000/seed\n"
        if full_document_text is None:
            full_document_text = (
                front_matter_text + "\n" + "\n".join(page_text_by_page.values())
            )

        def run(arguments: list[str], **kwargs: object):
            operation, source_path = ResearchValidationTests._assert_poppler_invocation(
                arguments, kwargs
            )
            if expected_pdf_bytes is not None and source_path.read_bytes() != expected_pdf_bytes:
                raise AssertionError("Poppler snapshot bytes differ from the opened source")
            if observed_snapshot_paths is not None:
                observed_snapshot_paths.append(source_path)
            if operation == "pdfinfo":
                return subprocess.CompletedProcess(
                    arguments,
                    0,
                    stdout=f"Pages: {page_count}\n{pdfinfo_metadata}",
                    stderr="",
                )
            if operation == "pdftotext":
                if "-f" in arguments:
                    if "-layout" in arguments:
                        output = front_matter_text
                    else:
                        page = int(arguments[arguments.index("-f") + 1])
                        output = page_text_by_page.get(page, "")
                else:
                    output = full_document_text
                return subprocess.CompletedProcess(
                    arguments, 0, stdout=output, stderr=""
                )
            if operation == "pdftoppm":
                Path(f"{arguments[-1]}.png").write_bytes(b"PNG test fixture")
                return subprocess.CompletedProcess(
                    arguments, 0, stdout=b"", stderr=b""
                )
            raise AssertionError(f"unexpected command: {arguments}")

        return run

    def _validate_local_with_poppler(self, side_effect) -> list[str]:
        with mock.patch(
            "scripts.check_manuscript.shutil.which",
            side_effect=lambda name: f"/mock/{name}",
        ), mock.patch(
            "scripts.check_manuscript.subprocess.run", side_effect=side_effect
        ):
            return validate_repository(
                self.root, check_git=False, check_local_pdfs=True
            ).errors

    @classmethod
    def _failing_poppler(
        cls,
        target: str,
        failure: str,
        observed_snapshot_paths: list[Path] | None = None,
    ):
        successful = cls._successful_poppler(
            observed_snapshot_paths=observed_snapshot_paths
        )

        def run(arguments: list[str], **kwargs: object):
            operation, source_path = cls._assert_poppler_invocation(arguments, kwargs)
            if observed_snapshot_paths is not None:
                observed_snapshot_paths.append(source_path)
            if operation != target:
                return successful(arguments, **kwargs)
            if failure == "timeout":
                raise subprocess.TimeoutExpired(
                    arguments, timeout=kwargs["timeout"]
                )
            if failure == "missing":
                raise FileNotFoundError(arguments[0])
            if failure == "oserror":
                raise OSError("test execution failure")
            if failure == "nonzero":
                output = "" if kwargs.get("text", True) else b""
                return subprocess.CompletedProcess(
                    arguments, 17, stdout=output, stderr=output
                )
            raise AssertionError(f"unknown failure: {failure}")

        return run

    def test_valid_pending_baseline_passes(self) -> None:
        self.assertEqual([], self._validate().errors)

    def test_missing_reference_registry_fails(self) -> None:
        (self.root / "literature/references.csv").unlink()
        errors = self._validate().errors
        self.assertTrue(any("literature/references.csv" in error for error in errors))

    def test_duplicate_reference_key_fails(self) -> None:
        self.reference_rows.append(dict(self.reference_rows[0]))
        self._rewrite_references()
        errors = self._validate().errors
        self.assertTrue(any("duplicate reference citation keys" in error for error in errors))

    def test_duplicate_bibtex_key_fails(self) -> None:
        entry = (
            "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
            "  year = {2020},\n  doi = {10.1000/seed}\n}\n"
        )
        self._write("paper/references.bib", entry + entry)
        errors = self._validate().errors
        self.assertTrue(any("duplicate bibliography keys" in error for error in errors))

    def test_unsupported_reference_status_fails(self) -> None:
        self.reference_rows[0]["reference_status"] = "machine_verified"
        self._rewrite_references()
        errors = self._validate().errors
        self.assertTrue(any("unsupported reference status" in error for error in errors))

    def test_verified_reference_requires_complete_metadata_and_review(self) -> None:
        for field_name in (
            "authors",
            "venue",
            "publisher_url",
            "metadata_source_url",
            "verified_by",
            "verified_at",
        ):
            with self.subTest(field_name=field_name):
                original = self.reference_rows[0][field_name]
                self.reference_rows[0][field_name] = ""
                self._rewrite_references()
                errors = self._validate().errors
                self.assertTrue(
                    any(field_name in error and "verified reference" in error for error in errors)
                )
                self.reference_rows[0][field_name] = original

    def test_verified_reference_rejects_truncated_authors(self) -> None:
        self.reference_rows[0]["authors"] = "Example, A. and others"
        self._rewrite_references()
        errors = self._validate().errors
        self.assertTrue(any("complete ordered author list" in error for error in errors))

    def test_verified_reference_doi_state_must_be_consistent(self) -> None:
        self.reference_rows[0]["doi_status"] = "not_assigned"
        self._rewrite_references()
        errors = self._validate().errors
        self.assertTrue(any("doi_status" in error for error in errors))

    def test_metadata_partial_cannot_claim_human_verification(self) -> None:
        self.reference_rows[0].update(
            {
                "reference_status": "metadata_partial",
                "verified_by": "Shawn",
                "verified_at": "2026-09-03",
                "notes": "Machine metadata check; human review pending.",
            }
        )
        self._rewrite_references()
        errors = self._validate().errors
        self.assertTrue(
            any("must leave verified_by and verified_at blank" in error for error in errors)
        )

    def test_nonverified_reference_cannot_enter_bibliography(self) -> None:
        self.reference_rows[0].update(
            {
                "reference_status": "metadata_partial",
                "verified_by": "",
                "verified_at": "",
                "notes": "Machine metadata check; human review pending.",
            }
        )
        self._rewrite_references()
        errors = self._validate().errors
        self.assertTrue(any("non-verified references" in error for error in errors))

    def test_verified_bibliography_requires_doi(self) -> None:
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
            "  year = {2020}\n}\n",
        )
        errors = self._validate().errors
        self.assertTrue(any("BibTeX DOI is required" in error for error in errors))

    def test_verified_bibliography_rejects_wrong_doi(self) -> None:
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
            "  year = {2020},\n  doi = {10.1000/wrong}\n}\n",
        )
        errors = self._validate().errors
        self.assertTrue(any("BibTeX DOI does not match references.csv" in error for error in errors))

    def test_verified_bibliography_rejects_swapped_dois(self) -> None:
        self.reference_rows.append(
            self._reference_row(
                citation_key="seed2",
                reference_status="verified",
                title="Second Seed",
                authors="Example, B.",
                year="2021",
                venue="Example Journal",
                doi="10.1000/seed2",
                doi_status="verified",
                publisher_url="https://publisher.example/seed2",
                metadata_source_url="https://doi.org/10.1000/seed2",
                verified_by="Researcher",
                verified_at="2026-09-03",
            )
        )
        self._rewrite_references()
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
            "  year = {2020},\n  doi = {10.1000/seed2}\n}\n"
            "@article{seed2,\n  author = {Example, B.},\n  title = {Second Seed},\n"
            "  year = {2021},\n  doi = {10.1000/seed}\n}\n",
        )
        errors = self._validate().errors
        mismatches = [
            error for error in errors if "BibTeX DOI does not match references.csv" in error
        ]
        self.assertEqual(2, len(mismatches))

    def test_bibliography_normalizes_legal_doi_variants(self) -> None:
        for doi in (
            "HTTPS://DOI.ORG/10.1000/SEED.",
            "http://doi.org/10.1000/SEED,",
            "DOI: 10.1000/SEED;",
        ):
            with self.subTest(doi=doi):
                self._write(
                    "paper/references.bib",
                    "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
                    f"  year = {{2020}},\n  doi = {{{doi}}}\n}}\n",
                )
                self.assertEqual([], self._validate().errors)

    def test_bibliography_rejects_doi_when_registry_says_not_assigned(self) -> None:
        self.reference_rows[0]["doi"] = ""
        self.reference_rows[0]["doi_status"] = "not_assigned"
        self._rewrite_references()
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
            "  year = {2020},\n  doi = {10.1000/unregistered}\n}\n",
        )
        errors = self._validate().errors
        self.assertTrue(any("DOI not assigned in references.csv" in error for error in errors))

    def test_bibliography_allows_missing_doi_when_not_assigned(self) -> None:
        self.reference_rows[0]["doi"] = ""
        self.reference_rows[0]["doi_status"] = "not_assigned"
        self._rewrite_references()
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
            "  year = {2020}\n}\n",
        )
        self.assertEqual([], self._validate().errors)

    def test_verified_bibliography_rejects_year_mismatch(self) -> None:
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
            "  year = {2021},\n  doi = {10.1000/seed}\n}\n",
        )
        errors = self._validate().errors
        self.assertTrue(any("BibTeX year does not match references.csv" in error for error in errors))

    def test_verified_bibliography_rejects_title_mismatch(self) -> None:
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A.},\n  title = {Different words},\n"
            "  year = {2020},\n  doi = {10.1000/seed}\n}\n",
        )
        errors = self._validate().errors
        self.assertTrue(any("BibTeX title does not match references.csv" in error for error in errors))

    def test_bibliography_normalizes_unicode_case_braces_and_whitespace_in_title(self) -> None:
        self.reference_rows[0]["title"] = "Café surface"
        self._rewrite_references()
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A.},\n"
            "  title = {{{CAFÉ}}   Surface},\n  year = {2020},\n"
            "  doi = {10.1000/seed}\n}\n",
        )
        self.assertEqual([], self._validate().errors)

    def test_verified_bibliography_requires_core_fields(self) -> None:
        fields = {
            "title": "  title = {Seed},\n",
            "author": "  author = {Example, A.},\n",
            "year": "  year = {2020},\n",
        }
        for missing in fields:
            with self.subTest(missing=missing):
                body = "".join(value for name, value in fields.items() if name != missing)
                self._write(
                    "paper/references.bib",
                    "@article{seed,\n" + body + "  doi = {10.1000/seed}\n}\n",
                )
                errors = self._validate().errors
                self.assertTrue(
                    any(f"BibTeX {missing} is required" in error for error in errors)
                )

    def test_real_verified_bibliography_passes_metadata_checks(self) -> None:
        self.assertEqual([], validate_repository(PROJECT_ROOT, check_git=True).errors)

    def test_tracked_zotero_database_fails(self) -> None:
        self._write("literature/zotero/zotero.sqlite", "not-a-real-database\n")
        self._initialize_git()
        errors = validate_repository(self.root, check_git=True).errors
        self.assertTrue(any("Zotero database or attachment" in error for error in errors))

    def test_undefined_citation_fails(self) -> None:
        self._write(
            "paper/manuscript.qmd",
            "# Draft\n\n<!-- CLAIM: C1 -->\nClaim [@missing].\n",
        )
        errors = self._validate().errors
        self.assertTrue(any("undefined citation keys" in error for error in errors))

    def test_duplicate_claim_source_relationship_fails(self) -> None:
        duplicate = dict(self.evidence_rows[0])
        duplicate["evidence_id"] = "E099"
        self.evidence_rows.append(duplicate)
        self._rewrite_evidence()
        errors = self._validate().errors
        self.assertTrue(any("duplicate claim-source relationships" in error for error in errors))

    def test_multiple_sources_can_share_a_claim_id(self) -> None:
        relationship_columns = [
            "evidence_id",
            *[column for column in EVIDENCE_COLUMNS if column != "evidence_id"],
        ]
        self.evidence_rows[0]["evidence_id"] = "E001"
        self.evidence_rows[1]["evidence_id"] = "E002"
        second_source = self._evidence_row(
            evidence_id="E003",
            claim_id="C1",
            status="pending",
            claim="Seed literature claim awaiting verification.",
            source_type="primary_publication",
            citation_key="seed2",
        )
        self.evidence_rows.append(second_source)
        self.reference_rows.append(
            self._reference_row(
                citation_key="seed2",
                reference_status="verified",
                title="Second Seed",
                authors="Example, B.",
                year="2021",
                venue="Example Journal",
                doi="10.1000/seed2",
                doi_status="verified",
                publisher_url="https://publisher.example/seed2",
                metadata_source_url="https://doi.org/10.1000/seed2",
                verified_by="Researcher",
                verified_at="2026-09-03",
            )
        )
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
            "  year = {2020},\n  doi = {10.1000/seed}\n}\n"
            "@article{seed2,\n  author = {Example, B.},\n  title = {Second Seed},\n"
            "  year = {2021},\n  doi = {10.1000/seed2}\n}\n",
        )
        self._write_csv(
            "literature/evidence.csv", relationship_columns, self.evidence_rows
        )
        self._rewrite_references()
        self.assertEqual([], self._validate().errors)

    def test_duplicate_evidence_id_fails(self) -> None:
        relationship_columns = [
            "evidence_id",
            *[column for column in EVIDENCE_COLUMNS if column != "evidence_id"],
        ]
        self.evidence_rows[0]["evidence_id"] = "E001"
        self.evidence_rows[1]["evidence_id"] = "E001"
        self._write_csv(
            "literature/evidence.csv", relationship_columns, self.evidence_rows
        )
        errors = self._validate().errors
        self.assertTrue(any("duplicate evidence IDs" in error for error in errors))

    def test_unsupported_evidence_status_fails(self) -> None:
        self.evidence_rows[0]["status"] = "unreviewed"
        self._rewrite_evidence()
        errors = self._validate().errors
        self.assertTrue(any("unsupported evidence status" in error for error in errors))

    def test_supported_literature_requires_locator_and_verification(self) -> None:
        self.evidence_rows[0].update(
            {"source_locator": "", "verified_by": "", "verified_at": ""}
        )
        self._rewrite_evidence()
        errors = self._validate().errors
        self.assertTrue(any("source_locator is required" in error for error in errors))
        self.assertTrue(any("verified_by is required" in error for error in errors))
        self.assertTrue(any("verified_at is required" in error for error in errors))

    def test_citation_requires_claim_marker(self) -> None:
        self._write("paper/manuscript.qmd", "# Draft\n\nClaim [@seed].\n")
        errors = self._validate().errors
        self.assertTrue(any("citation block requires a CLAIM marker" in error for error in errors))

    def test_pending_relationship_cannot_be_cited_as_fact(self) -> None:
        self.evidence_rows[0]["status"] = "pending"
        self._rewrite_evidence()
        errors = self._validate().errors
        self.assertTrue(any("does not have supported or partial evidence" in error for error in errors))

    def test_claim_marker_and_citation_must_match_relationship(self) -> None:
        self.reference_rows.append(
            self._reference_row(
                citation_key="seed2",
                reference_status="verified",
                title="Second Seed",
                authors="Example, B.",
                year="2021",
                venue="Example Journal",
                doi="10.1000/seed2",
                doi_status="verified",
                publisher_url="https://publisher.example/seed2",
                metadata_source_url="https://doi.org/10.1000/seed2",
                verified_by="Researcher",
                verified_at="2026-09-03",
            )
        )
        self._rewrite_references()
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
            "  year = {2020},\n  doi = {10.1000/seed}\n}\n"
            "@article{seed2,\n  author = {Example, B.},\n  title = {Second Seed},\n"
            "  year = {2021},\n  doi = {10.1000/seed2}\n}\n",
        )
        self._write(
            "paper/manuscript.qmd",
            "# Draft\n\n<!-- CLAIM: C1 -->\nClaim [@seed2].\n",
        )
        errors = self._validate().errors
        self.assertTrue(any("claim-source relationship is missing" in error for error in errors))

    def test_pending_claim_block_requires_todo_or_planned_marker(self) -> None:
        self.evidence_rows[0]["status"] = "pending"
        self._rewrite_evidence()
        self._write("paper/manuscript.qmd", "# Draft\n\n<!-- CLAIM: C1 -->\nUnqualified claim.\n")
        errors = self._validate().errors
        self.assertTrue(any("pending claim block requires TODO-EVIDENCE or PLANNED" in error for error in errors))

    def test_bibliography_rejects_each_forbidden_field(self) -> None:
        for field_name in ("file", "attachment", "note", "annote", "abstract", "keywords"):
            with self.subTest(field_name=field_name):
                self._write(
                    "paper/references.bib",
                    "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
                    "  year = {2020},\n  doi = {10.1000/seed},\n"
                    f"  {field_name} = {{private export content}}\n}}\n",
                )
                errors = self._validate().errors
                self.assertTrue(
                    any("forbidden field" in error for error in errors)
                )

    def test_bibliography_rejects_local_paths_on_supported_platforms(self) -> None:
        paths = {
            "file URL": "file:///Users/researcher/Zotero/item.pdf",
            "macOS user": "/Users/researcher/Zotero/item.pdf",
            "macOS volume": "/Volumes/Research/item.pdf",
            "Linux home": "/home/researcher/item.pdf",
            "private": "/private/var/item.pdf",
            "temporary": "/tmp/item.pdf",
            "home shorthand": "~/Zotero/item.pdf",
            "Windows backslash": "C:\\Users\\researcher\\Zotero\\item.pdf",
            "Windows slash": "C:/Users/researcher/Zotero/item.pdf",
        }
        for label, path in paths.items():
            with self.subTest(label=label):
                self._write(
                    "paper/references.bib",
                    "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n"
                    "  year = {2020},\n  doi = {10.1000/seed},\n"
                    f"  howpublished = {{{path}}}\n}}\n",
                )
                errors = self._validate().errors
                self.assertTrue(any("local path" in error for error in errors))

    def test_bibliography_allows_normal_web_urls_and_text_slashes(self) -> None:
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A./Example, B.},\n"
            "  title = {Seed with dry/wet comparison},\n  year = {2020},\n"
            "  doi = {https://doi.org/10.1000/SEED},\n"
            "  url = {https://publisher.example/articles/seed}\n}\n",
        )
        self.reference_rows[0]["title"] = "Seed with dry/wet comparison"
        self.reference_rows[0]["authors"] = "Example, A./Example, B."
        self._rewrite_references()
        self.assertEqual([], self._validate().errors)

    def test_supported_computational_claim_requires_manifest(self) -> None:
        self.evidence_rows[1].update(
            {"status": "supported", "run_id": "run-1", "verified_by": "R", "verified_at": "2026-09-03"}
        )
        self._rewrite_evidence()
        errors = self._validate().errors
        self.assertTrue(any("evidence_path is required" in error for error in errors))

    def test_invalid_manifest_checksums_fail(self) -> None:
        manifest_path = "workflow/runs/run-1.json"
        self._write(
            manifest_path,
            '{"run_id":"run-1","git":{"commit":"short","dirty":true},'
            '"configuration":{"sha256":"bad"},"inputs":[{}],"outputs":[{}],'
            '"metrics":{"energy":{"value":1}},"validation":{"passed":false}}',
        )
        self.evidence_rows[1].update(
            {
                "status": "supported",
                "evidence_path": manifest_path,
                "run_id": "run-1",
                "verified_by": "R",
                "verified_at": "2026-09-03",
            }
        )
        self.inventory_rows.append(
            self._inventory_row(
                artifact_id="RUN-1",
                stage="analysis",
                system="test",
                artifact_type="run manifest",
                path=manifest_path,
                status="completed",
                provenance_status="manifest_complete",
                manuscript_eligible="yes",
            )
        )
        self._rewrite_evidence()
        self._rewrite_inventory()
        errors = self._validate().errors
        self.assertTrue(any("SHA-256" in error for error in errors))
        self.assertTrue(any("passed validation decision" in error for error in errors))
        self.assertTrue(any("requires value and unit" in error for error in errors))

    def test_supported_manifest_verifies_real_files_commit_and_tracking(self) -> None:
        self._write("run-data/config.txt", "configuration\n")
        self._write("run-data/input.dat", "input\n")
        commit = self._initialize_git()
        self._write("run-data/output.dat", "output\n")
        self._add_supported_manifest(commit)
        self._git("add", ".")

        self.assertEqual([], validate_repository(self.root, check_git=True).errors)

    def test_manifest_checksum_mismatch_fails(self) -> None:
        self._write("run-data/config.txt", "configuration\n")
        self._write("run-data/input.dat", "input\n")
        commit = self._initialize_git()
        self._write("run-data/output.dat", "original output\n")
        self._add_supported_manifest(commit)
        self._write("run-data/output.dat", "tampered output\n")
        self._git("add", ".")

        errors = validate_repository(self.root, check_git=True).errors
        self.assertTrue(any("checksum mismatch" in error for error in errors))

    def test_manifest_commit_must_exist_in_repository(self) -> None:
        self._write("run-data/config.txt", "configuration\n")
        self._write("run-data/input.dat", "input\n")
        self._initialize_git()
        self._write("run-data/output.dat", "output\n")
        self._add_supported_manifest("f" * 40)
        self._git("add", ".")

        errors = validate_repository(self.root, check_git=True).errors
        self.assertTrue(any("does not identify a commit" in error for error in errors))

    def test_manifest_commit_must_be_reachable_from_shawn_dev(self) -> None:
        self._write("run-data/config.txt", "configuration\n")
        self._write("run-data/input.dat", "input\n")
        self._initialize_git()
        self._git("switch", "-q", "-c", "unrelated-run")
        self._write("unrelated.txt", "side branch only\n")
        self._git("add", "unrelated.txt")
        self._git("commit", "-q", "-m", "unrelated run commit")
        unrelated_commit = self._git("rev-parse", "HEAD")
        self._git("switch", "-q", "shawn_dev")
        self._write("run-data/output.dat", "output\n")
        self._add_supported_manifest(unrelated_commit)
        self._git("add", ".")

        errors = validate_repository(self.root, check_git=True).errors
        self.assertTrue(any("not reachable from shawn_dev" in error for error in errors))

    def test_supported_manifest_must_be_tracked(self) -> None:
        self._write("run-data/config.txt", "configuration\n")
        self._write("run-data/input.dat", "input\n")
        commit = self._initialize_git()
        self._write("run-data/output.dat", "output\n")
        manifest_path = self._add_supported_manifest(commit)

        errors = validate_repository(self.root, check_git=True).errors
        self.assertTrue(
            any(manifest_path in error and "must be tracked" in error for error in errors)
        )

    def test_run_manifest_example_contains_checker_required_sections(self) -> None:
        repository_root = Path(__file__).resolve().parents[1]
        example = json.loads(
            (repository_root / "workflow/schemas/run_manifest.example.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertTrue(example["metrics"])
        self.assertIn("validation", example)
        self.assertIn("passed", example["validation"])
        self.assertIn("reviewer", example["validation"])
        self.assertIn("reviewed_at", example["validation"])

    def test_unknown_todo_evidence_id_fails(self) -> None:
        self._write("paper/manuscript.qmd", "# Draft\n\n## TODO-EVIDENCE R9\n")
        errors = self._validate().errors
        self.assertTrue(any("R9" in error and "missing" in error for error in errors))

    def test_unregistered_publication_figure_fails(self) -> None:
        (self.root / "paper/figures/result.png").write_bytes(b"not-a-real-image")
        errors = self._validate().errors
        self.assertTrue(any("publication figure lacks" in error for error in errors))

    def test_conflict_marker_must_be_quarantined(self) -> None:
        marker = "<" * 7 + " branch\nvalue\n" + "=" * 7 + "\nother\n" + ">" * 7 + " branch\n"
        self._write("legacy/input.inp", marker)
        errors = self._validate().errors
        self.assertTrue(any("not quarantined" in error for error in errors))

        self.inventory_rows.append(
            self._inventory_row(
                artifact_id="CONFLICT",
                stage="legacy",
                system="test",
                artifact_type="conflicted_input",
                path="legacy/input.inp",
                status="invalid",
                provenance_status="path_only",
            )
        )
        self._rewrite_inventory()
        self.assertEqual([], self._validate().errors)

    def test_missing_relevance_registry_fails(self) -> None:
        (self.root / "literature/relevance.csv").unlink()
        errors = self._validate().errors
        self.assertTrue(any("literature/relevance.csv" in error for error in errors))

    def test_relevance_registry_must_match_reference_keys(self) -> None:
        self.relevance_rows[0]["citation_key"] = "unknown"
        self._rewrite_relevance()
        errors = self._validate().errors
        self.assertTrue(any("references missing from relevance.csv" in error for error in errors))
        self.assertTrue(any("unknown references" in error for error in errors))

    def test_relevance_ranks_must_be_unique_and_contiguous(self) -> None:
        self.relevance_rows[0]["priority_rank"] = "2"
        self._rewrite_relevance()
        errors = self._validate().errors
        self.assertTrue(any("unique and contiguous" in error for error in errors))

    def test_relevance_unknown_claim_fails_on_later_row(self) -> None:
        self.relevance_rows[0]["target_claim_ids"] = "C999"
        self._rewrite_relevance()
        errors = self._validate().errors
        self.assertTrue(any("not found in evidence.csv: C999" in error for error in errors))

    def test_relevance_multiple_valid_target_claims_pass(self) -> None:
        self.relevance_rows[0]["target_claim_ids"] = "C1|R1"
        self._rewrite_relevance()
        self.assertEqual([], self._validate().errors)

    def test_relevance_duplicate_target_claim_fails(self) -> None:
        self.relevance_rows[0]["target_claim_ids"] = "C1|C1"
        self._rewrite_relevance()
        errors = self._validate().errors
        self.assertTrue(any("duplicate target_claim_ids: C1" in error for error in errors))

    def test_relevance_empty_target_token_fails(self) -> None:
        for value in ("|C1", "C1|", "C1||R1", "C1|  |R1"):
            with self.subTest(value=value):
                self.relevance_rows[0]["target_claim_ids"] = value
                self._rewrite_relevance()
                errors = self._validate().errors
                self.assertTrue(
                    any("cannot contain empty tokens" in error for error in errors)
                )

    def test_future_extraction_wave_with_matching_card_passes(self) -> None:
        card = self._pending_card(selection_wave="phase_2b_2")
        self._write_card(card)
        self.assertEqual([], self._validate().errors)

    def test_phase_selected_reference_without_card_fails(self) -> None:
        self.relevance_rows[0]["selection_wave"] = "phase_2b_2"
        self._rewrite_relevance()
        errors = self._validate().errors
        self.assertTrue(any("phase-selected references missing" in error for error in errors))

    def test_orphan_card_without_phase_selected_reference_fails(self) -> None:
        card = self._pending_card(selection_wave="phase_2b_2")
        self.relevance_rows[0]["selection_wave"] = "later"
        self._rewrite_relevance()
        self._write_card(card)
        errors = self._validate().errors
        self.assertTrue(any("orphan evidence cards" in error for error in errors))

    def test_invalid_extraction_wave_identifier_fails(self) -> None:
        self.relevance_rows[0]["selection_wave"] = "phase_2b_next"
        self._rewrite_relevance()
        errors = self._validate().errors
        self.assertTrue(any("unsupported selection_wave" in error for error in errors))

    def test_phase_selection_requires_all_four_coverage_areas(self) -> None:
        card = self._pending_card(
            selection_wave="phase_2b_1",
            coverage_area="mineralization_interface_background",
        )
        self._write_card(card)
        errors = self._validate().errors
        self.assertTrue(any("selection must cover exactly" in error for error in errors))

    def test_phase_selection_cannot_exceed_four_references(self) -> None:
        coverage = [
            "hap_surface_hydration",
            "high_temperature_CO2_capture_boundary",
            "actual_compute_method",
            "mineralization_interface_background",
        ]
        self.relevance_rows[0].update(
            {
                "selection_wave": "phase_2b_1",
                "coverage_area": "mineralization_interface_background",
            }
        )
        for index, area in enumerate(coverage, start=2):
            key = f"candidate{index}"
            self.reference_rows.append(self._reference_row(citation_key=key))
            self.relevance_rows.append(
                self._relevance_row(
                    priority_rank=str(index),
                    citation_key=key,
                    coverage_area=area,
                    selection_wave="phase_2b_1",
                    target_claim_ids="C1",
                )
            )
        self._write_csv(
            "literature/references.csv", REFERENCE_COLUMNS, self.reference_rows
        )
        self._rewrite_relevance()
        errors = self._validate().errors
        self.assertTrue(any("exactly four" in error for error in errors))

    def test_phase_2b_1_frozen_key_set_fails_closed(self) -> None:
        card = self._pending_card(
            selection_wave="phase_2b_1",
            coverage_area="mineralization_interface_background",
        )
        self._write_card(card)
        errors = self._validate().errors
        self.assertTrue(any("selected keys must remain" in error for error in errors))

    def test_frozen_phase_rejects_relabeling_all_rows(self) -> None:
        self._install_frozen_phase_fixture()
        for row in self.relevance_rows:
            if row["citation_key"] in PHASE_2B_1_FROZEN:
                row["selection_wave"] = "phase_2b_2"
        self._rewrite_relevance()
        errors = self._validate().errors
        self.assertTrue(any("selected keys must remain" in error for error in errors))
        for key in PHASE_2B_1_FROZEN:
            self.assertTrue(any(key in error and "selection_wave" in error for error in errors))

    def test_frozen_phase_rejects_relabeling_one_row(self) -> None:
        self._install_frozen_phase_fixture()
        key = next(iter(PHASE_2B_1_FROZEN))
        row = next(row for row in self.relevance_rows if row["citation_key"] == key)
        row["selection_wave"] = "phase_2b_2"
        self._rewrite_relevance()
        errors = self._validate().errors
        self.assertTrue(any(key in error and "selection_wave" in error for error in errors))

    def test_frozen_phase_rejects_removing_all_four_cards(self) -> None:
        self._install_frozen_phase_fixture()
        for key in PHASE_2B_1_FROZEN:
            (self.root / f"literature/notes/{key}.md").unlink()
        errors = self._validate().errors
        for key in PHASE_2B_1_FROZEN:
            self.assertTrue(any(key in error and "exactly one evidence card" in error for error in errors))

    def test_frozen_phase_rejects_missing_relevance_rows(self) -> None:
        self._install_frozen_phase_fixture()
        self.relevance_rows = [
            row
            for row in self.relevance_rows
            if row["citation_key"] not in PHASE_2B_1_FROZEN
        ]
        self._rewrite_relevance()
        errors = self._validate().errors
        for key in PHASE_2B_1_FROZEN:
            self.assertTrue(any(key in error and "frozen relevance row is missing" in error for error in errors))

    def test_future_wave_coexists_with_frozen_phase_contract(self) -> None:
        self._install_frozen_phase_fixture()
        self.assertEqual([], self._validate().errors)

    def test_evidence_card_missing_required_field_fails(self) -> None:
        card = self._pending_card()
        card.pop("landing_page_url")
        self._write_card(card)
        errors = self._validate().errors
        self.assertTrue(any("missing required field" in error for error in errors))

    def test_evidence_card_doi_must_match_registry(self) -> None:
        card = self._pending_card(doi="10.1000/wrong")
        self._write_card(card)
        errors = self._validate().errors
        self.assertTrue(any("evidence-card DOI does not match" in error for error in errors))

    def test_evidence_card_requires_valid_sha256(self) -> None:
        card = self._pending_card(pdf_sha256="not-a-hash")
        self._write_card(card)
        errors = self._validate().errors
        self.assertTrue(any("requires a lowercase SHA-256" in error for error in errors))

    def test_versioned_arxiv_pdf_url_passes(self) -> None:
        card = self._pending_card(
            pdf_source_url="https://arxiv.org/pdf/2003.03868v2",
            full_text_version="preprint",
        )
        self._write_card(card)
        self.assertEqual([], self._validate().errors)

    def test_unversioned_arxiv_pdf_url_fails(self) -> None:
        card = self._pending_card(
            pdf_source_url="https://arxiv.org/pdf/2003.03868",
            full_text_version="preprint",
        )
        self._write_card(card)
        errors = self._validate().errors
        self.assertTrue(any("explicit v<number> suffix" in error for error in errors))

    def test_non_arxiv_pdf_url_is_unaffected(self) -> None:
        card = self._pending_card(
            pdf_source_url="https://publisher.example/fulltext.pdf"
        )
        self._write_card(card)
        self.assertEqual([], self._validate().errors)

    def test_pdf_filename_rejects_traversal_absolute_and_alternate_separators(self) -> None:
        for filename in ("../seed.pdf", "/outside/seed.pdf", "folder/seed.pdf", "folder\\seed.pdf"):
            with self.subTest(filename=filename):
                card = self._pending_card(pdf_filename=filename)
                self._write_card(card)
                errors = self._validate().errors
                self.assertTrue(any("must be the basename seed.pdf" in error for error in errors))

    def test_machine_evidence_card_cannot_claim_human_confirmation(self) -> None:
        card = self._pending_card(
            human_confirmed=True,
            human_reviewer="Reviewer",
            human_reviewed_at="2026-09-05",
        )
        self._write_card(card)
        errors = self._validate().errors
        self.assertTrue(any("requires human_confirmed to be JSON false" in error for error in errors))
        self.assertTrue(any("reviewer/date null" in error for error in errors))

    def test_valid_human_confirmed_card_with_pending_ledger_passes(self) -> None:
        card = self._pending_card(
            extraction_status="human_confirmed",
            human_confirmed=True,
            human_reviewer="Reviewer",
            human_reviewed_at="2026-09-06T12:34:56Z",
        )
        self._write_card(card)
        self.assertEqual([], self._validate().errors)
        self.assertEqual("pending", self.evidence_rows[0]["status"])

    def test_human_confirmed_status_requires_true_json_boolean(self) -> None:
        for value in (False, "true", "false", 1):
            with self.subTest(value=value):
                card = self._pending_card(
                    extraction_status="human_confirmed",
                    human_confirmed=value,
                    human_reviewer="Reviewer",
                    human_reviewed_at="2026-09-06T12:34:56Z",
                )
                self._write_card(card)
                errors = self._validate().errors
                self.assertTrue(any("requires human_confirmed to be JSON true" in error for error in errors))

    def test_unconfirmed_statuses_reject_true_human_boolean(self) -> None:
        for status in ("extracted_pending_human", "blocked_fulltext"):
            with self.subTest(status=status):
                overrides: dict[str, object] = {
                    "extraction_status": status,
                    "human_confirmed": True,
                }
                if status == "blocked_fulltext":
                    overrides.update(
                        pdf_filename=None,
                        pdf_sha256=None,
                        blocker_reason="Access blocked.",
                        claim_assessments=[],
                    )
                card = self._pending_card(**overrides)
                self._write_card(card)
                errors = self._validate().errors
                self.assertTrue(any("requires human_confirmed to be JSON false" in error for error in errors))

    def test_human_confirmed_card_requires_reviewer_and_valid_timestamp(self) -> None:
        cases = (
            ({"human_reviewer": None}, "requires human_reviewer"),
            ({"human_reviewer": ""}, "requires human_reviewer"),
            ({"human_reviewed_at": None}, "must be an ISO-8601 UTC timestamp"),
            ({"human_reviewed_at": "2026-09-06"}, "must be an ISO-8601 UTC timestamp"),
            ({"human_reviewed_at": "not-a-timestamp"}, "must be an ISO-8601 UTC timestamp"),
            ({"human_reviewed_at": "2026-99-99T12:34:56Z"}, "must be an ISO-8601 UTC timestamp"),
        )
        for override, expected in cases:
            with self.subTest(override=override):
                settings: dict[str, object] = {
                    "extraction_status": "human_confirmed",
                    "human_confirmed": True,
                    "human_reviewer": "Reviewer",
                    "human_reviewed_at": "2026-09-06T12:34:56Z",
                }
                settings.update(override)
                card = self._pending_card(**settings)
                self._write_card(card)
                self.assertTrue(any(expected in error for error in self._validate().errors))

    def test_human_confirmed_card_requires_pdf_sha_and_assessment(self) -> None:
        cases = (
            ({"pdf_filename": None}, "pdf_filename must be the basename"),
            ({"pdf_sha256": None}, "requires a lowercase SHA-256"),
            ({"claim_assessments": []}, "requires at least one claim assessment"),
        )
        for override, expected in cases:
            with self.subTest(override=override):
                card = self._pending_card(
                    extraction_status="human_confirmed",
                    human_confirmed=True,
                    human_reviewer="Reviewer",
                    human_reviewed_at="2026-09-06T12:34:56Z",
                    **override,
                )
                self._write_card(card)
                self.assertTrue(any(expected in error for error in self._validate().errors))

    def test_human_confirmed_card_hash_failure_is_rejected_locally(self) -> None:
        card, _ = self._local_pdf_card(
            extraction_status="human_confirmed",
            human_confirmed=True,
            human_reviewer="Reviewer",
            human_reviewed_at="2026-09-06T12:34:56Z",
        )
        card["pdf_sha256"] = "0" * 64
        self._write_card(card)
        errors = validate_repository(
            self.root, check_git=False, check_local_pdfs=True
        ).errors
        self.assertTrue(any("SHA-256 mismatch" in error for error in errors))

    def test_current_project_cards_remain_machine_unconfirmed(self) -> None:
        cards = sorted((PROJECT_ROOT / "literature/notes").glob("*.md"))
        cards = [path for path in cards if path.name != "README.md"]
        self.assertEqual(4, len(cards))
        for path in cards:
            block = re.search(
                r"```json\s*(\{.*?\})\s*```",
                path.read_text(encoding="utf-8"),
                flags=re.DOTALL,
            )
            self.assertIsNotNone(block)
            card = json.loads(block.group(1))
            self.assertIs(card["human_confirmed"], False)
            self.assertIsNone(card["human_reviewer"])
            self.assertIsNone(card["human_reviewed_at"])

    def test_evidence_card_rejects_more_than_25_quoted_words(self) -> None:
        card = self._pending_card()
        card["claim_assessments"][0]["verbatim_excerpt"] = " ".join(
            f"word{index}" for index in range(26)
        )
        self._write_card(card)
        errors = self._validate().errors
        self.assertTrue(any("maximum is 25" in error for error in errors))

    def test_evidence_card_rejects_unknown_claim(self) -> None:
        card = self._pending_card()
        card["claim_assessments"][0]["claim_id"] = "C999"
        self._write_card(card)
        errors = self._validate().errors
        self.assertTrue(any("no evidence.csv relationship" in error for error in errors))

    def test_blocked_card_cannot_contain_locator_assessments(self) -> None:
        card = self._pending_card(
            extraction_status="blocked_fulltext",
            pdf_filename=None,
            pdf_sha256=None,
            blocker_reason="Publisher returned HTTP 403.",
        )
        self._write_card(card)
        errors = self._validate().errors
        self.assertTrue(any("cannot contain claim assessments" in error for error in errors))

    def test_duplicate_evidence_card_is_rejected(self) -> None:
        card = self._pending_card()
        self._write_card(card)
        self._write_card(card, filename_key="duplicate")
        errors = self._validate().errors
        self.assertTrue(any("duplicate evidence card" in error for error in errors))

    def test_machine_card_does_not_upgrade_evidence_ledger(self) -> None:
        card = self._pending_card()
        self._write_card(card)
        self.evidence_rows[0]["status"] = "partial"
        self.evidence_rows[0]["source_locator"] = "p. 1"
        self._rewrite_evidence()
        errors = self._validate().errors
        self.assertTrue(any("must not upgrade evidence status" in error for error in errors))

    def test_local_pdf_check_rejects_missing_pdf(self) -> None:
        card = self._pending_card()
        self._write_card(card)
        errors = validate_repository(
            self.root, check_git=False, check_local_pdfs=True
        ).errors
        self.assertTrue(any("local evidence PDF is missing" in error for error in errors))

    def test_default_check_does_not_invoke_poppler_for_extracted_card(self) -> None:
        card = self._pending_card()
        self._write_card(card)
        with mock.patch(
            "scripts.check_manuscript.subprocess.run",
            side_effect=AssertionError("default validation invoked an external command"),
        ):
            self.assertEqual([], self._validate().errors)

    def test_blocked_fulltext_card_does_not_require_local_pdf(self) -> None:
        card = self._pending_card(
            extraction_status="blocked_fulltext",
            pdf_filename=None,
            pdf_sha256=None,
            blocker_reason="Publisher access failed closed.",
            claim_assessments=[],
        )
        self._write_card(card)
        errors = validate_repository(
            self.root, check_git=False, check_local_pdfs=True
        ).errors
        self.assertEqual([], errors)

    def test_local_pdf_check_rejects_hash_mismatch(self) -> None:
        card = self._pending_card()
        self._write_card(card)
        pdf = self.root / "literature/pdfs/seed.pdf"
        pdf.parent.mkdir(parents=True)
        pdf.write_bytes(b"%PDF-test fixture")
        errors = validate_repository(
            self.root, check_git=False, check_local_pdfs=True
        ).errors
        self.assertTrue(any("SHA-256 mismatch" in error for error in errors))

    def test_local_pdf_uses_clean_private_snapshot_directory(self) -> None:
        self._local_pdf_card()
        snapshots: list[Path] = []
        errors = self._validate_local_with_poppler(
            self._successful_poppler(observed_snapshot_paths=snapshots)
        )
        self.assertEqual([], errors)
        self.assertTrue(snapshots)
        self.assertTrue(all(not path.exists() for path in snapshots))
        self.assertFalse((self.root / "tmp/pdfs").exists())

    def test_local_pdf_file_symlink_escape_is_rejected(self) -> None:
        card = self._pending_card()
        with tempfile.TemporaryDirectory() as external_directory:
            sentinel = Path(external_directory) / "sentinel.pdf"
            sentinel_bytes = b"%PDF-1.4\nexternal sentinel\n"
            sentinel.write_bytes(sentinel_bytes)
            pdf_path = self.root / "literature/pdfs/seed.pdf"
            pdf_path.parent.mkdir(parents=True)
            pdf_path.symlink_to(sentinel)
            card["pdf_sha256"] = hashlib.sha256(sentinel_bytes).hexdigest()
            self._write_card(card)

            errors = validate_repository(
                self.root, check_git=False, check_local_pdfs=True
            ).errors

            self.assertTrue(any("or a symlink" in error for error in errors))
            self.assertEqual(sentinel_bytes, sentinel.read_bytes())

    def test_local_pdf_directory_symlink_escape_is_rejected(self) -> None:
        card = self._pending_card()
        with tempfile.TemporaryDirectory() as external_directory:
            external_root = Path(external_directory)
            sentinel = external_root / "sentinel.txt"
            sentinel.write_text("unchanged", encoding="utf-8")
            (external_root / "seed.pdf").write_bytes(b"%PDF-1.4\nexternal PDF\n")
            pdf_directory = self.root / "literature/pdfs"
            pdf_directory.parent.mkdir(parents=True, exist_ok=True)
            pdf_directory.symlink_to(external_root, target_is_directory=True)
            self._write_card(card)

            errors = validate_repository(
                self.root, check_git=False, check_local_pdfs=True
            ).errors

            self.assertTrue(any("must not follow symlinks" in error for error in errors))
            self.assertEqual("unchanged", sentinel.read_text(encoding="utf-8"))

    def test_local_pdf_source_is_opened_with_no_follow(self) -> None:
        self._local_pdf_card()
        real_open = os.open
        with mock.patch(
            "scripts.check_manuscript.os.open", wraps=real_open
        ) as open_spy:
            errors = self._validate_local_with_poppler(self._successful_poppler())
        self.assertEqual([], errors)
        source_calls = [
            call
            for call in open_spy.call_args_list
            if call.args and call.args[0] == "seed.pdf" and "dir_fd" in call.kwargs
        ]
        self.assertEqual(1, len(source_calls))
        if hasattr(os, "O_NOFOLLOW"):
            self.assertTrue(source_calls[0].args[1] & os.O_NOFOLLOW)

    def test_nonregular_local_pdf_source_is_rejected(self) -> None:
        card = self._pending_card()
        pdf_path = self.root / "literature/pdfs/seed.pdf"
        pdf_path.mkdir(parents=True)
        self._write_card(card)
        errors = validate_repository(
            self.root, check_git=False, check_local_pdfs=True
        ).errors
        self.assertTrue(any("must be a regular file" in error for error in errors))

    def test_pathname_replaced_by_symlink_after_open_uses_stable_snapshot(self) -> None:
        _, pdf_path = self._local_pdf_card()
        original_bytes = pdf_path.read_bytes()
        snapshot_paths: list[Path] = []
        real_open_source = manuscript_checks._open_local_pdf_descriptor
        with tempfile.TemporaryDirectory() as external_directory:
            sentinel = Path(external_directory) / "replacement.pdf"
            sentinel.write_bytes(b"%PDF-1.4\nattacker replacement\n")

            def open_then_replace(*args: object, **kwargs: object):
                descriptor = real_open_source(*args, **kwargs)
                if descriptor is not None:
                    pdf_path.unlink()
                    pdf_path.symlink_to(sentinel)
                return descriptor

            with mock.patch(
                "scripts.check_manuscript._open_local_pdf_descriptor",
                side_effect=open_then_replace,
            ):
                errors = self._validate_local_with_poppler(
                    self._successful_poppler(
                        expected_pdf_bytes=original_bytes,
                        observed_snapshot_paths=snapshot_paths,
                    )
                )
            self.assertEqual([], errors)
            self.assertTrue(pdf_path.is_symlink())
            self.assertEqual(b"%PDF-1.4\nattacker replacement\n", sentinel.read_bytes())
        self.assertTrue(snapshot_paths)
        self.assertTrue(all(not path.parent.exists() for path in snapshot_paths))

    def test_pathname_replaced_with_different_sha_keeps_snapshot_consistent(self) -> None:
        _, pdf_path = self._local_pdf_card()
        original_bytes = pdf_path.read_bytes()
        replacement_bytes = b"%PDF-1.4\ndifferent replacement hash\n"
        real_open_source = manuscript_checks._open_local_pdf_descriptor

        def open_then_replace(*args: object, **kwargs: object):
            descriptor = real_open_source(*args, **kwargs)
            if descriptor is not None:
                pdf_path.unlink()
                pdf_path.write_bytes(replacement_bytes)
            return descriptor

        with mock.patch(
            "scripts.check_manuscript._open_local_pdf_descriptor",
            side_effect=open_then_replace,
        ):
            errors = self._validate_local_with_poppler(
                self._successful_poppler(expected_pdf_bytes=original_bytes)
            )
        self.assertEqual([], errors)
        self.assertEqual(replacement_bytes, pdf_path.read_bytes())
        self.assertNotEqual(
            hashlib.sha256(original_bytes).hexdigest(),
            hashlib.sha256(replacement_bytes).hexdigest(),
        )

    def test_repository_render_symlink_is_never_followed(self) -> None:
        self._local_pdf_card()
        with tempfile.TemporaryDirectory() as external_directory:
            external_root = Path(external_directory)
            sentinel = external_root / "sentinel.txt"
            sentinel.write_text("unchanged", encoding="utf-8")
            tmp_root = self.root / "tmp"
            tmp_root.mkdir()
            (tmp_root / "pdfs").symlink_to(external_root, target_is_directory=True)

            errors = self._validate_local_with_poppler(self._successful_poppler())

            self.assertEqual([], errors)
            self.assertEqual("unchanged", sentinel.read_text(encoding="utf-8"))
            self.assertEqual([sentinel], list(external_root.iterdir()))

    def test_poppler_timeouts_fail_closed_for_every_operation(self) -> None:
        self._local_pdf_card()
        for operation in ("pdfinfo", "pdftotext", "pdftoppm"):
            with self.subTest(operation=operation):
                snapshots: list[Path] = []
                errors = self._validate_local_with_poppler(
                    self._failing_poppler(operation, "timeout", snapshots)
                )
                self.assertTrue(
                    any(operation in error and "timed out" in error for error in errors)
                )
                self.assertTrue(snapshots)
                self.assertTrue(all(not path.parent.exists() for path in snapshots))

    def test_missing_poppler_executables_fail_closed_for_every_operation(self) -> None:
        self._local_pdf_card()
        for operation in ("pdfinfo", "pdftotext", "pdftoppm"):
            with self.subTest(operation=operation):
                errors = self._validate_local_with_poppler(
                    self._failing_poppler(operation, "missing")
                )
                self.assertTrue(
                    any(
                        operation in error and "executable was not found" in error
                        for error in errors
                    )
                )

    def test_poppler_oserrors_fail_closed_for_every_operation(self) -> None:
        self._local_pdf_card()
        for operation in ("pdfinfo", "pdftotext", "pdftoppm"):
            with self.subTest(operation=operation):
                snapshots: list[Path] = []
                errors = self._validate_local_with_poppler(
                    self._failing_poppler(operation, "oserror", snapshots)
                )
                self.assertTrue(
                    any(
                        operation in error and "could not be started" in error
                        for error in errors
                    )
                )
                self.assertTrue(snapshots)
                self.assertTrue(all(not path.parent.exists() for path in snapshots))

    def test_poppler_nonzero_status_fails_closed_for_every_operation(self) -> None:
        self._local_pdf_card()
        for operation in ("pdfinfo", "pdftotext", "pdftoppm"):
            with self.subTest(operation=operation):
                snapshots: list[Path] = []
                errors = self._validate_local_with_poppler(
                    self._failing_poppler(operation, "nonzero", snapshots)
                )
                self.assertTrue(
                    any(
                        operation in error and "exit status 17" in error
                        for error in errors
                    )
                )
                self.assertTrue(snapshots)
                self.assertTrue(all(not path.parent.exists() for path in snapshots))

    def test_failed_render_leaves_no_partial_output(self) -> None:
        self._local_pdf_card()
        successful = self._successful_poppler()
        snapshot_paths: list[Path] = []

        def fail_render(arguments: list[str], **kwargs: object):
            operation, source_path = self._assert_poppler_invocation(arguments, kwargs)
            snapshot_paths.append(source_path)
            if operation == "pdftoppm":
                Path(f"{arguments[-1]}.png").write_bytes(b"partial")
                return subprocess.CompletedProcess(
                    arguments, 9, stdout=b"", stderr=b""
                )
            return successful(arguments, **kwargs)

        errors = self._validate_local_with_poppler(fail_render)
        self.assertTrue(any("pdftoppm" in error for error in errors))
        self.assertTrue(snapshot_paths)
        self.assertTrue(all(not path.parent.exists() for path in snapshot_paths))

    def test_locator_excerpt_must_exist_on_declared_page(self) -> None:
        self._local_pdf_card()
        page_text = {
            1: "Results\nDifferent source text.\n",
            2: "Results\nShort source text.\n",
        }
        errors = self._validate_local_with_poppler(
            self._successful_poppler(page_text)
        )
        self.assertTrue(any("excerpt was not found" in error for error in errors))

    def test_abstract_cannot_be_a_claim_bearing_section(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = "Abstract"
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler({1: "Abstract\nShort source text.\n"})
        )
        self.assertTrue(any("cannot use Abstract" in error for error in errors))

    def test_reference_sections_cannot_be_claim_bearing(self) -> None:
        for heading in ("References", "Bibliography"):
            with self.subTest(heading=heading):
                card, _ = self._local_pdf_card()
                card["claim_assessments"][0]["section_heading"] = heading
                self._write_card(card)
                errors = self._validate_local_with_poppler(
                    self._successful_poppler(
                        {1: f"{heading}\nShort source text.\n"}
                    )
                )
                self.assertTrue(any("cannot use Abstract" in error for error in errors))

    def test_forbidden_heading_boundary_after_declared_heading_is_rejected(self) -> None:
        cases = (
            ("Introduction", "I. ABSTRACT"),
            ("Introduction", "7 REFERENCES"),
            ("Introduction", "VII. REFERENCES"),
            ("Acknowledgements", "Bibliography"),
        )
        for declared_heading, forbidden_heading in cases:
            with self.subTest(forbidden_heading=forbidden_heading):
                card, _ = self._local_pdf_card()
                card["claim_assessments"][0]["section_heading"] = declared_heading
                self._write_card(card)
                page_text = {
                    1: f"{declared_heading}\nBody text.\n{forbidden_heading}\n"
                    "Short source text.\n"
                }
                errors = self._validate_local_with_poppler(
                    self._successful_poppler(page_text)
                )
                self.assertTrue(
                    any("crosses forbidden section heading" in error for error in errors)
                )

    def test_hyphenated_body_line_cannot_hide_forbidden_heading_boundary(self) -> None:
        for forbidden_heading in ("ABSTRACT", "REFERENCES", "BIBLIOGRAPHY"):
            with self.subTest(forbidden_heading=forbidden_heading):
                card, _ = self._local_pdf_card()
                card["claim_assessments"][0]["section_heading"] = "Introduction"
                self._write_card(card)
                page_text = {
                    1: "Introduction\nA body line ends-\n"
                    f"{forbidden_heading}\nShort source text.\n"
                }
                errors = self._validate_local_with_poppler(
                    self._successful_poppler(page_text)
                )
                self.assertTrue(
                    any("crosses forbidden section heading" in error for error in errors)
                )

    def test_genuine_line_ending_hyphenation_still_matches_exact_excerpt(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = "Results"
        card["claim_assessments"][0]["verbatim_excerpt"] = (
            "The microstructure remains stable."
        )
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                {1: "Results\nThe micro-\nstructure remains stable.\n"}
            )
        )
        self.assertEqual([], errors)

    def test_dehyphenated_excerpt_view_retains_source_line_provenance(self) -> None:
        lines = manuscript_checks._normalize_pdf_source_lines(
            "Results\nThe micro-\nstructure remains stable.\n"
        )
        view = manuscript_checks._build_pdf_excerpt_view(lines)
        start = view.text.index("microstructure")
        continuation = start + len("micro")
        self.assertEqual(1, view.source_line_by_character[start])
        self.assertEqual(2, view.source_line_by_character[continuation])

    def test_heading_phrase_inside_abstract_prose_is_not_a_heading(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = "Introduction"
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                {
                    1: "ABSTRACT\nThis Introduction summarizes the work.\n"
                    "Short source text.\n"
                }
            )
        )
        self.assertTrue(any("section heading was not found" in error for error in errors))

    def test_exact_standalone_heading_line_passes(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = "Introduction"
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler({1: "INTRODUCTION\nShort source text.\n"})
        )
        self.assertEqual([], errors)

    def test_numbered_standalone_heading_line_passes(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = "Introduction"
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler({1: "I. INTRODUCTION\nShort source text.\n"})
        )
        self.assertEqual([], errors)

    def test_genuine_two_line_split_heading_passes(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = (
            "IX. AB-INITIO MOLECULAR DYNAMICS"
        )
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                {1: "IX. AB-INITIO\nMOLECULAR DYNAMICS\nShort source text.\n"}
            )
        )
        self.assertEqual([], errors)

    def test_complete_two_line_prose_cannot_be_declared_as_heading(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = (
            "This is ordinary prose split across two lines"
        )
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                {
                    1: "This is ordinary prose split\nacross two lines\n"
                    "Short source text.\n"
                }
            )
        )
        self.assertTrue(any("section heading was not found" in error for error in errors))

    def test_split_prose_containing_heading_phrase_is_rejected(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = "Introduction"
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                {
                    1: "ABSTRACT\nThis Introduction\nsummarizes the work.\n"
                    "Short source text.\n"
                }
            )
        )
        self.assertTrue(any("section heading was not found" in error for error in errors))

    def test_two_line_heading_requires_exact_combined_value(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = (
            "IX. AB-INITIO MOLECULAR DYNAMICS"
        )
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                {
                    1: "IX. AB-INITIO\nMOLECULAR DYNAMICS AND METHODS\n"
                    "Short source text.\n"
                }
            )
        )
        self.assertTrue(any("section heading was not found" in error for error in errors))

    def test_heading_cannot_span_more_than_two_lines(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = (
            "IX. AB-INITIO MOLECULAR DYNAMICS METHODS"
        )
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                {
                    1: "IX. AB-INITIO\nMOLECULAR DYNAMICS\nMETHODS\n"
                    "Short source text.\n"
                }
            )
        )
        self.assertTrue(any("section heading was not found" in error for error in errors))

    def test_duplicate_heading_uses_nearest_genuine_preceding_heading(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = "Introduction"
        self._write_card(card)
        page_text = {
            1: "Introduction\nEarlier body.\nREFERENCES\n"
            "I. INTRODUCTION\nShort source text.\n"
        }
        errors = self._validate_local_with_poppler(
            self._successful_poppler(page_text)
        )
        self.assertEqual([], errors)

    def test_valid_introduction_after_abstract_boundary_passes(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = "Introduction"
        self._write_card(card)
        page_text = {
            1: "ABSTRACT\nAbstract text.\nI. Introduction\nShort source text.\n"
        }
        errors = self._validate_local_with_poppler(
            self._successful_poppler(page_text)
        )
        self.assertEqual([], errors)

    def test_prose_references_mention_is_not_a_section_boundary(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["section_heading"] = "Introduction"
        self._write_card(card)
        page_text = {
            1: "Introduction\nThis prose references earlier work in context.\n"
            "Short source text.\n"
        }
        errors = self._validate_local_with_poppler(
            self._successful_poppler(page_text)
        )
        self.assertEqual([], errors)

    def test_current_extracted_card_locator_shapes_pass_boundary_validation(self) -> None:
        cases = (
            (
                "[H1] Aqueous Carbonation Reaction Mechanisms\n"
                "It begins with CO2 solubilization in water, i.e., solubility trapping, "
                "and CO2 speciation to carbonic acid, bicarbonate, and carbonate",
                "[H1] Aqueous Carbonation Reaction Mechanisms",
                "It begins with CO2 solubilization in water, i.e., solubility trapping, "
                "and CO2 speciation to carbonic acid, bicarbonate, and carbonate",
            ),
            (
                "IV. DISCUSSION\nDetails of the atomic structure of the surfaces lead "
                "to different strengths of adsorption",
                "IV. DISCUSSION",
                "Details of the atomic structure of the surfaces lead to different "
                "strengths of adsorption",
            ),
            (
                "A.\nBorn-Oppenheimer molecular dynamics\nIn Born-Oppenheimer MD "
                "the potential energy is minimized at every AIMD step",
                "A. Born-Oppenheimer molecular dynamics",
                "is minimized at every AIMD step",
            ),
        )
        for page_text, heading, excerpt in cases:
            with self.subTest(heading=heading):
                self.assertEqual(
                    (None, None),
                    manuscript_checks._validate_locator_on_page(
                        page_text, heading, excerpt
                    ),
                )

    def test_locator_section_heading_must_exist_on_declared_page(self) -> None:
        self._local_pdf_card()
        errors = self._validate_local_with_poppler(
            self._successful_poppler({1: "Discussion\nShort source text.\n"})
        )
        self.assertTrue(any("section heading was not found" in error for error in errors))

    def test_locator_excerpt_must_follow_section_heading(self) -> None:
        self._local_pdf_card()
        errors = self._validate_local_with_poppler(
            self._successful_poppler({1: "Short source text.\nResults\n"})
        )
        self.assertTrue(any("excerpt precedes" in error for error in errors))

    def test_locator_normalization_handles_unicode_soft_hyphen_and_line_breaks(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["verbatim_excerpt"] = (
            "The microstructure uses coffee analysis."
        )
        self._write_card(card)
        page_text = {1: "Results\nThe micro-\nstructure uses coﬀee\u00ad analysis.\n"}
        errors = self._validate_local_with_poppler(
            self._successful_poppler(page_text)
        )
        self.assertEqual([], errors)

    def test_locator_page_outside_pdf_is_rejected(self) -> None:
        card, _ = self._local_pdf_card()
        card["claim_assessments"][0]["pdf_page"] = 3
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler(page_count=2)
        )
        self.assertTrue(any("outside PDF page count" in error for error in errors))

    def test_invalid_pdf_header_is_rejected_before_poppler(self) -> None:
        card, pdf_path = self._local_pdf_card()
        pdf_path.write_bytes(b"not a PDF")
        card["pdf_sha256"] = hashlib.sha256(pdf_path.read_bytes()).hexdigest()
        self._write_card(card)
        errors = validate_repository(
            self.root, check_git=False, check_local_pdfs=True
        ).errors
        self.assertTrue(any("does not begin with %PDF-" in error for error in errors))

    def test_empty_full_document_text_is_rejected(self) -> None:
        self._local_pdf_card()
        successful = self._successful_poppler()

        def empty_text(arguments: list[str], **kwargs: object):
            if Path(arguments[0]).name == "pdftotext" and "-f" not in arguments:
                return subprocess.CompletedProcess(arguments, 0, stdout="", stderr="")
            return successful(arguments, **kwargs)

        errors = self._validate_local_with_poppler(empty_text)
        self.assertTrue(any("returned no usable text" in error for error in errors))

    def test_page_two_bibliography_identity_cannot_rescue_wrong_front_matter(self) -> None:
        self._local_pdf_card()
        pages_one_and_two = (
            "A Different Document\nDifferent, B.\n\fREFERENCES\n"
            "Seed\nExample, A.\nDOI: 10.1000/seed\n"
        )
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                front_matter_text=pages_one_and_two,
            )
        )
        self.assertTrue(any("front-matter title block" in error for error in errors))

    def test_version_of_record_doi_after_introduction_does_not_count(self) -> None:
        self._local_pdf_card()
        pages_one_and_two = (
            "Seed\nExample, A.\nI. INTRODUCTION\n"
            "The body cites DOI: 10.1000/seed.\n"
        )
        errors = self._validate_local_with_poppler(
            self._successful_poppler(front_matter_text=pages_one_and_two)
        )
        self.assertTrue(any("version_of_record" in error and "DOI" in error for error in errors))

    def test_exact_target_citation_after_abstract_cannot_rescue_identity(self) -> None:
        self._local_pdf_card()
        pages_one_and_two = (
            "A Different Document\nDifferent, B.\nABSTRACT\n"
            "Seed\nExample, A.\nDOI: 10.1000/seed\n"
        )
        errors = self._validate_local_with_poppler(
            self._successful_poppler(front_matter_text=pages_one_and_two)
        )
        self.assertTrue(any("front-matter title block" in error for error in errors))

    def test_page_one_boundary_prevents_page_two_identity_rescue(self) -> None:
        self._local_pdf_card()
        pages_one_and_two = (
            "A Different Document\nDifferent, B.\nABSTRACT\nUnrelated text.\n\f"
            "Seed\nExample, A.\nDOI: 10.1000/seed\n"
        )
        errors = self._validate_local_with_poppler(
            self._successful_poppler(front_matter_text=pages_one_and_two)
        )
        self.assertTrue(any("front-matter title block" in error for error in errors))

    def test_page_one_prose_citation_cannot_replace_document_title(self) -> None:
        self._local_pdf_card()
        front = (
            "A Different Document\nDifferent, B.\nABSTRACT\n"
            "This paper discusses Seed by Example and cites DOI 10.1000/seed.\n"
        )
        errors = self._validate_local_with_poppler(
            self._successful_poppler(front_matter_text=front)
        )
        self.assertTrue(any("front-matter title block" in error for error in errors))

    def test_version_of_record_with_matching_front_matter_doi_passes(self) -> None:
        self._local_pdf_card()
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                front_matter_text=(
                    "Seed\nExample, A.\nDOI: 10.1000/seed\n1 Introduction\n"
                    "Body text.\n"
                )
            )
        )
        self.assertEqual([], errors)

    def test_version_of_record_without_front_matter_or_metadata_doi_fails(self) -> None:
        self._local_pdf_card()
        errors = self._validate_local_with_poppler(
            self._successful_poppler(front_matter_text="Seed\nExample, A.\n")
        )
        self.assertTrue(any("version_of_record" in error and "DOI" in error for error in errors))

    def test_version_of_record_accepts_matching_pdfinfo_doi(self) -> None:
        self._local_pdf_card()
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                front_matter_text="Seed\nExample, A.\n",
                pdfinfo_metadata="Subject: DOI 10.1000/seed\n",
            )
        )
        self.assertEqual([], errors)

    def test_version_of_record_with_conflicting_front_matter_doi_fails(self) -> None:
        self._local_pdf_card()
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                front_matter_text="Seed\nExample, A.\nDOI: 10.9999/different\n"
            )
        )
        self.assertTrue(any("conflicting DOI" in error for error in errors))

    def test_preprint_without_publisher_doi_passes_front_matter_identity(self) -> None:
        card, _ = self._local_pdf_card(
            full_text_version="preprint",
            pdf_source_url="https://arxiv.org/pdf/2003.03868v2",
        )
        card["identity_checks"]["doi_match_basis"] = "official_landing_metadata"
        card["identity_checks"]["doi_association_url"] = "https://arxiv.org/abs/2003.03868"
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                front_matter_text="Seed\nExample, A.\nI. INTRODUCTION\nBody.\n"
            )
        )
        self.assertEqual([], errors)

    def test_accepted_manuscript_without_embedded_doi_passes(self) -> None:
        card, _ = self._local_pdf_card(
            full_text_version="author_accepted_manuscript"
        )
        card["identity_checks"]["doi_match_basis"] = "official_landing_metadata"
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                front_matter_text="Seed\nExample, A.\nABSTRACT\nBody.\n"
            )
        )
        self.assertEqual([], errors)

    def test_preprint_with_conflicting_front_matter_doi_fails(self) -> None:
        card, _ = self._local_pdf_card(
            full_text_version="preprint",
            pdf_source_url="https://arxiv.org/pdf/2003.03868v2",
        )
        card["identity_checks"]["doi_match_basis"] = "official_landing_metadata"
        card["identity_checks"]["doi_association_url"] = "https://arxiv.org/abs/2003.03868"
        self._write_card(card)
        errors = self._validate_local_with_poppler(
            self._successful_poppler(
                front_matter_text="Seed\nExample, A.\nDOI: 10.9999/different\n"
            )
        )
        self.assertTrue(any("conflicting DOI" in error for error in errors))

    def test_legacy_postprocessor_has_no_synthetic_fallback(self) -> None:
        repository_root = Path(__file__).resolve().parents[1]
        source_path = (
            repository_root
            / "4.PostProcess_Analysis/HAP112_DPA_NVTcsvr_Succeeded/lammps_postprocessing.py"
        )
        source = source_path.read_text(encoding="utf-8")
        compile(source, str(source_path), "exec")
        self.assertNotIn("generate_sample_rdf_data", source)
        self.assertNotIn("np.random", source)
        self.assertIn("raise DataValidationError", source)


class _FakeBooleanArray:
    def __init__(self, value: bool):
        self.value = value

    def all(self) -> bool:
        return self.value


class _FakeArray:
    def __init__(self, rows):
        self.values = [value for row in rows for value in row]
        self.size = len(self.values)

    def __eq__(self, other):
        return _FakeBooleanArray(all(value == other for value in self.values))


class _FakeSeries:
    def __init__(self, values):
        self.values = values

    def min(self):
        return min(self.values)

    def max(self):
        return max(self.values)


class _FakeDataFrame:
    def __init__(self, rows, columns=None):
        self.rows = rows
        self.columns = columns or [f"column_{index}" for index in range(len(rows[0]))]
        self.shape = (len(rows), len(self.columns))

    def __len__(self):
        return len(self.rows)

    def __getitem__(self, column):
        index = self.columns.index(column)
        return _FakeSeries([row[index] for row in self.rows])

    def select_dtypes(self, include=None):
        return self

    def to_numpy(self, dtype=float, copy=False):
        return _FakeArray(self.rows)


def _load_legacy_postprocessor():
    repository_root = Path(__file__).resolve().parents[1]
    source_path = (
        repository_root
        / "4.PostProcess_Analysis/HAP112_DPA_NVTcsvr_Succeeded/lammps_postprocessing.py"
    )

    fake_numpy = types.ModuleType("numpy")
    fake_numpy.number = float
    fake_numpy.isfinite = lambda values: _FakeBooleanArray(
        all(math.isfinite(value) for value in values.values)
    )
    fake_numpy.all = lambda values: values.all()

    fake_pandas = types.ModuleType("pandas")
    fake_pandas.DataFrame = _FakeDataFrame
    fake_pandas.read_csv = mock.Mock(side_effect=ValueError("unparseable test input"))

    fake_pyplot = types.ModuleType("matplotlib.pyplot")
    fake_pyplot.rcParams = {}
    fake_matplotlib = types.ModuleType("matplotlib")
    fake_matplotlib.pyplot = fake_pyplot
    fake_seaborn = types.ModuleType("seaborn")
    fake_stats = types.ModuleType("scipy.stats")
    fake_signal = types.ModuleType("scipy.signal")
    fake_signal.find_peaks = lambda *args, **kwargs: ([], {})
    fake_scipy = types.ModuleType("scipy")
    fake_scipy.stats = fake_stats

    fake_modules = {
        "numpy": fake_numpy,
        "pandas": fake_pandas,
        "matplotlib": fake_matplotlib,
        "matplotlib.pyplot": fake_pyplot,
        "seaborn": fake_seaborn,
        "scipy": fake_scipy,
        "scipy.stats": fake_stats,
        "scipy.signal": fake_signal,
    }
    spec = importlib.util.spec_from_file_location("legacy_lammps_postprocessor_test", source_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load legacy postprocessor for testing")
    module = importlib.util.module_from_spec(spec)
    with mock.patch.dict(sys.modules, fake_modules):
        spec.loader.exec_module(module)
    return module


class LegacyPostprocessorBehaviorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.module = _load_legacy_postprocessor()

    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def _processor(self):
        processor = self.module.LAMMPSPostProcessor.__new__(
            self.module.LAMMPSPostProcessor
        )
        processor.input_dir = self.root
        processor.output_dir = self.root / "output"
        processor.data = {}
        processor.results = {}
        processor.anomalies = {}
        processor.data_debug = {}
        return processor

    def _write(self, name: str, content: str = "placeholder\n") -> Path:
        path = self.root / name
        path.write_text(content, encoding="utf-8")
        return path

    def test_rdf_parser_accepts_complete_lammps_blocks(self) -> None:
        rdf_path = self._write(
            "rdf.dat",
            "# TimeStep Number-of-rows\n"
            "0 2\n"
            "1 0.5 0.0 0.0\n"
            "2 1.0 1.5 0.2\n"
            "100 2\n"
            "1 0.5 0.1 0.0\n"
            "2 1.0 1.4 0.3\n",
        )
        data = self._processor().load_rdf_data_enhanced(rdf_path)
        self.assertEqual(4, len(data))

    def test_rdf_parser_rejects_empty_zero_malformed_and_truncated_data(self) -> None:
        cases = {
            "empty": "",
            "zero": "0 2\n1 0.5 0.0 0.0\n2 1.0 0.0 0.0\n",
            "malformed": "0 1\n1 0.5 not-a-number 0.0\n",
            "truncated": "0 2\n1 0.5 1.0 0.0\n",
        }
        for name, content in cases.items():
            with self.subTest(name=name):
                path = self._write(f"{name}.dat", content)
                with self.assertRaises(self.module.DataValidationError):
                    self._processor().load_rdf_data_enhanced(path)

    def test_required_main_data_parse_failure_is_fatal(self) -> None:
        self._write("comprehensive_thermodynamics.dat")
        processor = self._processor()
        with mock.patch.object(
            processor, "load_data_with_enhanced_validation", return_value=None
        ):
            with self.assertRaises(self.module.DataValidationError):
                processor.load_data()

    def test_all_zero_main_data_is_fatal(self) -> None:
        path = self._write("all-zero.dat", "0 0\n0 0\n")
        self.module.pd.read_csv = mock.Mock(
            return_value=_FakeDataFrame([[0.0, 0.0], [0.0, 0.0]], ["a", "b"])
        )
        with self.assertRaisesRegex(self.module.DataValidationError, "All numeric values are zero"):
            self._processor().load_data_with_enhanced_validation(path)

    def test_missing_optional_rdfs_do_not_reject_current_producer_output(self) -> None:
        main_files = [
            "comprehensive_thermodynamics.dat",
            "co2_chemical_transformation.dat",
            "proton_transfer_water_formation.dat",
            "energy_fluctuations.dat",
            "msd_analysis.dat",
            "co2_orientation.dat",
            "temperature_evolution.log",
            "production_temperature.log",
        ]
        required_rdfs = [
            "rdf_CO_detailed.dat",
            "rdf_C_Ca.dat",
            "rdf_C_P.dat",
            "rdf_C_H_detailed.dat",
            "rdf_H_O_detailed.dat",
        ]
        for filename in main_files + required_rdfs:
            self._write(filename)

        processor = self._processor()
        with (
            mock.patch.object(
                processor, "load_data_with_enhanced_validation", return_value=object()
            ),
            mock.patch.object(processor, "analyze_log_file"),
            mock.patch.object(processor, "save_debug_information"),
        ):
            processor.load_data()
        self.assertEqual(5, len(processor.data["rdf"]))

    def test_present_invalid_optional_rdf_is_fatal(self) -> None:
        main_files = [
            "comprehensive_thermodynamics.dat",
            "co2_chemical_transformation.dat",
            "proton_transfer_water_formation.dat",
            "energy_fluctuations.dat",
            "msd_analysis.dat",
            "co2_orientation.dat",
            "temperature_evolution.log",
            "production_temperature.log",
        ]
        required_rdfs = [
            "rdf_CO_detailed.dat",
            "rdf_C_Ca.dat",
            "rdf_C_P.dat",
            "rdf_C_H_detailed.dat",
            "rdf_H_O_detailed.dat",
        ]
        for filename in main_files + required_rdfs + ["rdf_Ca_O.dat"]:
            self._write(filename)

        processor = self._processor()

        def load(path):
            if path.name == "rdf_Ca_O.dat":
                raise self.module.DataValidationError("invalid optional RDF")
            return object()

        with (
            mock.patch.object(
                processor, "load_data_with_enhanced_validation", side_effect=load
            ),
            mock.patch.object(processor, "analyze_log_file"),
            mock.patch.object(processor, "save_debug_information"),
        ):
            with self.assertRaises(self.module.DataValidationError):
                processor.load_data()

    def test_empty_or_incomplete_log_is_fatal_and_complete_log_passes(self) -> None:
        processor = self._processor()
        log_path = self._write("log.lammps", "")
        with self.assertRaises(self.module.DataValidationError):
            processor.analyze_log_file()

        log_path.write_text("LAMMPS test build\nrun 100\n", encoding="utf-8")
        with self.assertRaises(self.module.DataValidationError):
            processor.analyze_log_file()

        log_path.write_text(
            "LAMMPS test build\nLoop time of 1.0 on 1 procs for 100 steps\n",
            encoding="utf-8",
        )
        processor.analyze_log_file()
        self.assertTrue(processor.data["log_analysis"]["simulation_progress"]["completed"])

    def test_analysis_failure_maps_to_nonzero_cli_status(self) -> None:
        processor = mock.Mock()
        processor.run_complete_enhanced_analysis.return_value = False
        existing_input = mock.Mock()
        existing_input.exists.return_value = True
        with (
            mock.patch.object(self.module, "Path", return_value=existing_input),
            mock.patch.object(
                self.module, "LAMMPSPostProcessor", return_value=processor
            ),
        ):
            self.assertEqual(1, self.module.main())


if __name__ == "__main__":
    unittest.main()
