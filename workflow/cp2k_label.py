"""Stage 4: CP2K labelling with the fixed v2 protocol.

  prepare: frames (extxyz) -> one directory per frame with input.inp (+ frame.extxyz)
  run:     shell loop for the hap-cp2k image; skips frames that already have a 'done' marker
  collect: CP2K outputs -> DeePMD npy systems grouped by (system, formula), with a label
           audit (SCF convergence, max force) so failed labels are excluded, not hidden.
           CP2K 2026.x prints "FORCES|" blocks that dpdata 0.2.24 cannot parse, so parsing uses
           cp2k_protocol.parse_output; virial = stress[GPa] * V / 160.2177 (dpdata's sign).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cp2k_protocol import HA_EV, HA_BOHR_EV_A, cp2k_input, parse_output  # noqa: E402

TYPE_MAP = ["Ca", "P", "O", "H", "C"]

RUN_SCRIPT = """#!/bin/bash
# Usage inside hap-cp2k: run_labels.sh <label_dir> <mpi_ranks> <omp_threads>
cd "$1" || exit 1
for d in */; do
  d=${d%/}
  [ -f "$d/done" ] && continue
  ( cd "$d" && OMP_NUM_THREADS=$3 mpirun -np $2 cp2k.psmp -i input.inp -o output.out \\
      && grep -q "PROGRAM ENDED AT" output.out && touch done ) || echo "FAILED $d"
done
echo "ALL_LABELS_FINISHED"
"""


def prepare(args):
    from ase.io import read, write
    args.label_dir.mkdir(parents=True, exist_ok=True)
    index = []
    for path in args.frames:
        for k, atoms in enumerate(read(path, ":")):
            name = f"{path.stem}_{k:05d}"
            d = args.label_dir / name
            d.mkdir(exist_ok=True)
            if (d / "done").exists():
                # Finished label: keep its frame, whose orientation matches output.out
                atoms = read(d / "frame.extxyz")
                index.append(dict(name=name, source=str(path), frame=k, natoms=len(atoms),
                                  system=atoms.info.get("system", path.stem),
                                  stress="STRESS_TENSOR" in (d / "input.inp").read_text()))
                continue
            # CP2K 2026.2 rejects some valid cells (e.g. ASE fcc primitive vectors) with "No meaningful
            # cell information"; the lower-triangular standard form is accepted and is only a rotation
            atoms.set_cell(atoms.cell.standard_form()[0], scale_atoms=True)
            # Slabs (vacuum along c) get no stress; every other cell is a periodic bulk or liquid
            periodic_bulk = args.stress or not atoms.info.get("system", path.stem).startswith("hap001")
            text = cp2k_input(atoms, "label", args.cutoff, args.rel_cutoff)
            if periodic_bulk:
                text = text.replace("  METHOD Quickstep\n", "  METHOD Quickstep\n  STRESS_TENSOR ANALYTICAL\n")
                text = text.replace("    &FORCES ON\n    &END FORCES\n",
                                    "    &FORCES ON\n    &END FORCES\n    &STRESS_TENSOR ON\n    &END STRESS_TENSOR\n")
            (d / "input.inp").write_text(text, newline="\n")
            write(d / "frame.extxyz", atoms)
            index.append(dict(name=name, source=str(path), frame=k, natoms=len(atoms),
                              system=atoms.info.get("system", path.stem), stress=periodic_bulk))
    (args.label_dir / "run_labels.sh").write_text(RUN_SCRIPT, newline="\n")  # bash rejects CRLF
    with open(args.label_dir / "index.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(index[0].keys()))
        writer.writeheader()
        writer.writerows(index)
    print(f"prepared {len(index)} frames in {args.label_dir}")


def collect(args):
    from ase.io import read
    rows = list(csv.DictReader(open(args.label_dir / "index.csv")))
    audit, groups = [], {}
    for r in rows:
        d = args.label_dir / r["name"]
        status = "missing"
        if (d / "done").exists():
            res = parse_output(d / "output.out")
            atoms = read(d / "frame.extxyz")
            forces = np.array(res["forces_Ha_bohr"]) * HA_BOHR_EV_A
            if not res["scf_converged"]:
                status = "scf_not_converged"
            elif forces.shape != (len(atoms), 3):
                status = "parse_error:forces"
            elif r["stress"] == "True" and res["stress_GPa"] is None:
                status = "parse_error:stress"
            elif np.linalg.norm(forces, axis=1).max() > args.max_force:
                status = "force_outlier"
            else:
                status = "ok"
                key = (r["system"], atoms.get_chemical_formula())
                g = groups.setdefault(key, dict(types=[TYPE_MAP.index(s) for s in atoms.get_chemical_symbols()],
                                                coord=[], box=[], energy=[], force=[], virial=[]))
                g["coord"].append(atoms.positions.reshape(-1))
                g["box"].append(atoms.cell.array.reshape(-1))
                g["energy"].append(res["energy_Ha"] * HA_EV)
                g["force"].append(forces.reshape(-1))
                g["virial"].append(None if res["stress_GPa"] is None
                                   else (np.array(res["stress_GPa"]) * atoms.get_volume() / 160.21766208).reshape(-1))
        audit.append({**r, "status": status})
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    for (system, formula), g in groups.items():
        n = len(g["energy"])
        valid = set(rng.choice(n, size=int(round(n * args.valid_fraction)), replace=False)) if n >= 10 else set()
        for split, idx in (("train", [i for i in range(n) if i not in valid]), ("valid", sorted(valid))):
            if not idx:
                continue
            sd = args.out / split / f"{system}__{formula}"
            (sd / "set.000").mkdir(parents=True, exist_ok=True)
            (sd / "type.raw").write_text("\n".join(map(str, g["types"])) + "\n")
            (sd / "type_map.raw").write_text("\n".join(TYPE_MAP) + "\n")
            for k in ("coord", "box", "energy", "force"):
                np.save(sd / "set.000" / f"{k}.npy", np.array([g[k][i] for i in idx]))
            if all(g["virial"][i] is not None for i in idx):
                np.save(sd / "set.000" / "virial.npy", np.array([g["virial"][i] for i in idx]))
    with open(args.out / "label_audit.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(audit[0].keys()))
        writer.writeheader()
        writer.writerows(audit)
    counts = {}
    for a in audit:
        counts[a["status"]] = counts.get(a["status"], 0) + 1
    (args.out / "label_summary.json").write_text(json.dumps(
        dict(counts=counts, systems={f"{s}__{f}": len(g["energy"]) for (s, f), g in groups.items()}), indent=2) + "\n")
    print(counts)


def main():
    parser = argparse.ArgumentParser(description="CP2K labelling")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--frames", nargs="+", type=Path, required=True)
    p.add_argument("--label-dir", type=Path, required=True)
    p.add_argument("--cutoff", type=int, required=True)
    p.add_argument("--rel-cutoff", type=int, required=True)
    p.add_argument("--stress", action="store_true", help="compute stress for every frame")
    c = sub.add_parser("collect")
    c.add_argument("--label-dir", type=Path, required=True)
    c.add_argument("--out", type=Path, required=True)
    c.add_argument("--max-force", type=float, default=15.0, help="eV/A; larger frames are excluded as outliers")
    c.add_argument("--valid-fraction", type=float, default=0.1, help="0 for a sealed test set")
    c.add_argument("--seed", type=int, default=20261029)
    args = parser.parse_args()
    prepare(args) if args.cmd == "prepare" else collect(args)


if __name__ == "__main__":
    main()
