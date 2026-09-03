#!/usr/bin/env python3
"""Dependency-free checks for manuscript citations and evidence discipline."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT = ROOT / "paper" / "manuscript.qmd"
BIBLIOGRAPHY = ROOT / "paper" / "references.bib"
EVIDENCE = ROOT / "literature" / "evidence.csv"

BIB_KEY = re.compile(r"@[A-Za-z]+\s*\{\s*([^,\s]+)\s*,")
CITE_KEY = re.compile(r"(?<![\w])@([A-Za-z0-9_:.-]+)")


def fail(message: str) -> None:
    print(f"ERROR: {message}", file=sys.stderr)


def main() -> int:
    errors = 0
    for path in (MANUSCRIPT, BIBLIOGRAPHY, EVIDENCE):
        if not path.is_file():
            fail(f"missing required file: {path.relative_to(ROOT)}")
            errors += 1
    if errors:
        return 1

    manuscript_text = MANUSCRIPT.read_text(encoding="utf-8")
    bibliography_text = BIBLIOGRAPHY.read_text(encoding="utf-8")

    bib_keys = set(BIB_KEY.findall(bibliography_text))
    cite_keys = set(CITE_KEY.findall(manuscript_text))
    missing = sorted(cite_keys - bib_keys)
    unused = sorted(bib_keys - cite_keys)

    if missing:
        fail("undefined citation keys: " + ", ".join(missing))
        errors += 1

    with EVIDENCE.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        actual_columns = set(reader.fieldnames or [])

    required_columns = {
        "claim_id",
        "status",
        "claim",
        "source_type",
        "evidence_path",
        "citation_key",
        "owner",
        "notes",
    }
    if required_columns != actual_columns:
        fail("evidence.csv columns differ from the required schema")
        errors += 1

    claim_ids = [row.get("claim_id", "") for row in rows]
    duplicate_ids = sorted({item for item in claim_ids if item and claim_ids.count(item) > 1})
    if duplicate_ids:
        fail("duplicate claim IDs: " + ", ".join(duplicate_ids))
        errors += 1

    evidence_citations = {row.get("citation_key", "") for row in rows}
    missing_evidence_bib = sorted((evidence_citations - {""}) - bib_keys)
    if missing_evidence_bib:
        fail("evidence ledger uses undefined citations: " + ", ".join(missing_evidence_bib))
        errors += 1

    todo_count = manuscript_text.count("TODO-EVIDENCE")
    print(f"Citations: {len(cite_keys)} used, {len(bib_keys)} available")
    if unused:
        print("Unused bibliography entries: " + ", ".join(unused))
    print(f"Evidence rows: {len(rows)}")
    print(f"Open TODO-EVIDENCE markers: {todo_count}")

    if errors:
        return 1
    print("Manuscript checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

