"""Verify the published file inventory; never start scientific calculations."""
from pathlib import Path
import hashlib
import json


def main():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / "results/release-20261017.manifest.json").read_text(encoding="utf-8"))
    checked = 0
    for item in manifest["files"]:
        path = (root / item["path"]).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            raise ValueError("Missing or out-of-scope file: " + item["path"])
        data = path.read_bytes()
        if len(data) != item["bytes"] or hashlib.sha256(data).hexdigest() != item["sha256"]:
            raise ValueError("File identity differs: " + item["path"])
        checked += 1
    assert manifest["base_commit"] == "96d3967eb932495607095bf38a272fab4ebc4c08"
    assert manifest["new_joint40_executed"] is False
    assert manifest["scientific_calculations_run_for_publication"] is False
    assert manifest["manuscript_eligible"] is False
    print(f"PASS: {checked} release files verified; no scientific calculation performed.")


if __name__ == "__main__":
    main()
