"""Read-only package hashes and XYZ element counts; no scientific calculations."""
from pathlib import Path
from collections import Counter
import hashlib
import json

ROOT = Path(__file__).resolve().parent

def require(ok, message):
    if not ok:
        raise ValueError(message)

def local(name):
    p = (ROOT / name).resolve()
    require(p.is_relative_to(ROOT), "Path outside package")
    return p

def digest(name):
    return hashlib.sha256(local(name).read_bytes()).hexdigest()

def read(name):
    return json.loads(local(name).read_text(encoding="utf-8-sig"))

def elements(name):
    lines = local(name).read_text(encoding="utf-8").splitlines()
    n = int(lines[0])
    require(len(lines) >= n + 2 and not any(x.strip() for x in lines[n + 2:]), "Expected exactly one XYZ frame")
    counts = Counter(line.split()[0] for line in lines[2:n + 2])
    require(sum(counts.values()) == n, "XYZ atom count mismatch")
    return n, dict(counts)

def main():
    manifest = read("manifest.json")
    for entry in manifest["files"]:
        require(local(entry["path"]).stat().st_size == entry["bytes"], "Size mismatch: " + entry["path"])
        require(digest(entry["path"]) == entry["sha256"], "Hash mismatch: " + entry["path"])
    for entry in read("source-provenance.json")["sources"]:
        require(digest(entry["publication_path"]) == entry["publication_sha256"], "Publication provenance mismatch")
        if entry["transformation"] == "byte-for-byte copy":
            require(entry["publication_sha256"] == entry["source_sha256"], "Copied original hash mismatch")
    filt = read("hap44/filter_summary.json")
    rule = read("hap44/hap44-screening.json")
    require(len(filt["candidates"]) == 4, "Candidate count mismatch")
    for c in filt["candidates"]:
        name = "hap44/" + c["candidate_id"] + "/generated_crystals.extxyz"
        n, counts = elements(name)
        require(n == rule["hard_gates"]["atom_count_equals"], "HAP44 frame length mismatch")
        require(counts == c["element_counts"], "Recorded element counts mismatch")
        require(digest(name) == c["raw_structure_sha256"], "Candidate raw structure hash mismatch")
        require(counts != rule["hard_gates"]["element_counts_equal"], "Expected recorded composition rejection")
        require(c["reason_codes"] == ["exact_composition_mismatch"], "Reason record changed")
    sf = read("smoke/filter_summary.json")
    ng, cg = elements("smoke/generated_crystals.extxyz")
    nr, cr = elements("smoke/relaxed_structures.extxyz")
    require(ng == nr == 8 and cg == cr, "Smoke element counts changed on recorded relaxation")
    require("P" not in cg and "P" in sf["required_elements"], "Smoke missing-P rejection mismatch")
    require(sf["input_count"] == 1 and sf["saved_count"] == 0 and sf["chemistry_rejected"] == 1, "Smoke filter counts mismatch")
    require(digest("smoke/relaxed_structures.extxyz") == sf["structures_sha256"], "Smoke filter structure binding mismatch")
    print("PASS: package hashes and recorded XYZ element counts only; no scientific evaluation")

if __name__ == "__main__":
    main()
