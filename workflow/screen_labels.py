"""Energy screen of finished CP2K labels against a trained 4-model DeePMD ensemble (run inside hap-mlip).

A label whose SCF converged to a wrong electronic state passes every check in cp2k_label.py collect: the SCF
converged, the forces parse and stay below the outlier cap. Six such labels were found on 2026-11-05 (5 of
HAP(001) + water + CO2, 1 portlandite), 50-86 meV/atom above the energy the ensemble predicts, while the other
labels of the same system scatter by about 1 meV/atom (workflow/PIPELINE_V2.md, "Wrong SCF states").

For each system, the residual (E_model - E_DFT) / N is compared with the system median. A frame is flagged when
it deviates by more than --threshold while the 4 models agree to within --max-spread, so frames the ensemble has
not learned (new systems, new chemistry) are not flagged. Systems with fewer than --min-frames frames are skipped.
A frame relabelled from the atomic guess (cp2k_label.py relabel keeps the first result as output.suspect.out) is
not flagged again when both runs agree to within --reproduce-tol: the energy is then reproducible and the
ensemble is what is wrong.

  python screen_labels.py --label-dir iter-004/labels --models iter-003/explore/models/model-*.pb
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cp2k_protocol import HA_EV, parse_output  # noqa: E402

TYPE_MAP = ["Ca", "P", "O", "H", "C"]


def flag(rows, threshold, max_spread, min_frames, reproduce_tol=1.0):
    """Mark rows (dicts with system, residual, spread and, after a relabel, relabel_change) deviating from their
    system median; returns the medians."""
    medians = {}
    for system in sorted({r["system"] for r in rows}):
        group = [r for r in rows if r["system"] == system]
        if len(group) < min_frames:
            continue
        med = float(np.median([r["residual_meV_atom"] for r in group]))
        medians[system] = med
        for r in group:
            r["deviation_meV_atom"] = r["residual_meV_atom"] - med
            reproduced = abs(r.get("relabel_change_meV_atom", np.inf)) <= reproduce_tol
            r["flagged"] = bool(abs(r["deviation_meV_atom"]) > threshold and r["spread_meV_atom"] < max_spread
                                and not reproduced)
    return medians


def main():
    parser = argparse.ArgumentParser(description="Flag CP2K labels whose energy disagrees with a trained ensemble")
    parser.add_argument("--label-dir", type=Path, required=True)
    parser.add_argument("--models", nargs=4, type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=10.0, help="meV/atom from the system median")
    parser.add_argument("--max-spread", type=float, default=5.0, help="meV/atom, largest - smallest model energy")
    parser.add_argument("--min-frames", type=int, default=5)
    parser.add_argument("--reproduce-tol", type=float, default=1.0,
                        help="meV/atom; a relabel this close to the first result confirms it")
    parser.add_argument("--out", default="energy_screen.json", help="written inside --label-dir")
    args = parser.parse_args()
    from ase.io import read
    from deepmd.infer import DeepPot
    pots = [DeepPot(str(m)) for m in args.models]
    rows = []
    for d in sorted(p for p in args.label_dir.iterdir() if (p / "done").exists()):
        res = parse_output(d / "output.out")
        if not res["scf_converged"] or res["energy_Ha"] is None:
            continue
        a = read(d / "frame.extxyz")
        types = [TYPE_MAP.index(s) for s in a.get_chemical_symbols()]
        es = [float(p.eval(a.positions.reshape(1, -1), a.cell.array.reshape(1, -1), types)[0].item()) for p in pots]
        n = len(a)
        row = dict(name=d.name, system=a.info.get("system", "?"), natoms=n, energy_dft_eV=res["energy_Ha"] * HA_EV,
                   residual_meV_atom=1000 * (np.mean(es) - res["energy_Ha"] * HA_EV) / n,
                   spread_meV_atom=1000 * (max(es) - min(es)) / n, flagged=False)
        if (d / "output.suspect.out").exists():
            first = parse_output(d / "output.suspect.out")["energy_Ha"]
            row["suspect_energy_dft_eV"] = first * HA_EV
            row["relabel_change_meV_atom"] = 1000 * (res["energy_Ha"] - first) * HA_EV / n
        rows.append(row)
    medians = flag(rows, args.threshold, args.max_spread, args.min_frames, args.reproduce_tol)
    flagged = [r["name"] for r in rows if r["flagged"]]
    (args.label_dir / args.out).write_text(json.dumps(dict(
        models=[str(m) for m in args.models], threshold_meV_atom=args.threshold, max_spread_meV_atom=args.max_spread,
        min_frames=args.min_frames, system_median_meV_atom=medians, flagged=flagged, rows=rows), indent=2) + "\n")
    for r in rows:
        if r["flagged"]:
            print(f"FLAGGED {r['name']} {r['system']} deviation {r['deviation_meV_atom']:+.1f} meV/atom")
    print(f"screened {len(rows)} labels in {args.label_dir}: {len(flagged)} flagged")


if __name__ == "__main__":
    main()
