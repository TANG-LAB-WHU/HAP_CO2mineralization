"""Evaluate cp2k_protocol.py runs: cutoff and Gamma-width convergence against the tightest case.

Gates (PIPELINE_V2.md): |dE| <= 1 meV/atom and force-RMS difference <= 0.01 eV/A versus
cut1000_rel60; |dE| <= 1 meV/atom versus the widest Gamma cell. The smallest passing
setting is reported; nothing is chosen if the reference run itself failed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cp2k_protocol import parse_output  # noqa: E402

HA = 27.211386245988
HA_BOHR = 51.42208619083232


def main():
    parser = argparse.ArgumentParser(description="CP2K protocol convergence report")
    parser.add_argument("--dir", type=Path, required=True)
    parser.add_argument("--run-subdir", default="run", help="subdirectory of each case holding output.out")
    args = parser.parse_args()
    plan = json.loads((args.dir / "plan.json").read_text())
    res = {}
    for case in plan:
        out = args.dir / case["name"] / args.run_subdir / "output.out"
        if out.exists():
            r = parse_output(out)
            res[case["name"]] = {**case, **r}
    report = {"cutoff": [], "width": []}
    ref = res.get("cut1000_rel60")
    if ref and ref["scf_converged"]:
        f_ref = np.array(ref["forces_Ha_bohr"]) * HA_BOHR
        for name, r in res.items():
            if not name.startswith("cut") or not r["scf_converged"]:
                continue
            de = (r["energy_Ha"] - ref["energy_Ha"]) * HA / r["natoms"] * 1000
            df = float(np.sqrt(np.mean((np.array(r["forces_Ha_bohr"]) * HA_BOHR - f_ref) ** 2)))
            report["cutoff"].append(dict(name=name, dE_meV_atom=round(de, 3), dF_rms_eV_A=round(df, 4),
                                         wall_s=r["wall_s"], passes=abs(de) <= 1.0 and df <= 0.01))
    widths = sorted((r for n, r in res.items() if n.startswith("width") and r["scf_converged"]), key=lambda r: r["natoms"])
    if widths:
        wref = widths[-1]
        for r in widths:
            de = (r["energy_Ha"] / r["natoms"] - wref["energy_Ha"] / wref["natoms"]) * HA * 1000
            report["width"].append(dict(name=r["name"], natoms=r["natoms"], dE_meV_atom=round(de, 3),
                                        wall_s=r["wall_s"], passes=abs(de) <= 1.0))
    passing = [c for c in report["cutoff"] if c["passes"]]
    report["chosen_cutoff"] = min(passing, key=lambda c: int(c["name"][3:].split("_")[0]))["name"] if passing else None
    wpass = [w for w in report["width"] if w["passes"]]
    report["chosen_width"] = min(wpass, key=lambda w: w["natoms"])["name"] if wpass else None
    (args.dir / "protocol_report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
