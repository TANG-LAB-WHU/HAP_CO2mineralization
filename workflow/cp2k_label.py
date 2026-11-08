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
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cp2k_protocol import HA_EV, HA_BOHR_EV_A, cp2k_input, parse_output  # noqa: E402

TYPE_MAP = ["Ca", "P", "O", "H", "C"]

RUN_SCRIPT = """#!/bin/bash
# Usage inside hap-cp2k: run_labels.sh <label_dir> <mpi_ranks> <omp_threads> [<shard> <shards>]
# <mpi_ranks> may name a file holding the rank count, read again before each frame (seed containers).
# <shard> <shards> (from 2026-11-08): only every <shards>-th frame directory in sorted order, starting at <shard>
# (0-based), so separate jobs (e.g. a Slurm array) label disjoint frames of one directory.
# SCF fallback (from 2026-11-05): a frame whose OT DIIS SCF fails is run once more from its last wavefunction
# with OT CG + FULL_ALL. Functional, basis, cutoffs and EPS_SCF are unchanged; the first output is kept as
# output.diis.out. A frame that fails both gets a 'failed' marker and is not run again.
cd "$1" || exit 1
ranks() { if [ -f "$1" ]; then cat "$1"; else echo "$1"; fi; }
cp2k() { OMP_NUM_THREADS=$OMP mpirun -np "$(ranks "$RANKS")" cp2k.psmp -i "$1" -o output.out \\
           && grep -q "PROGRAM ENDED AT" output.out && ! grep -q "SCF run NOT converged" output.out; }
RANKS=$2 OMP=$3 SHARD=${4:-0} SHARDS=${5:-1}
i=-1
for d in */; do
  d=${d%/}
  i=$((i + 1))
  [ $((i % SHARDS)) -eq "$SHARD" ] || continue
  { [ -f "$d/done" ] || [ -f "$d/failed" ]; } && continue
  cd "$d" || continue
  if cp2k input.inp; then
    touch done
  elif [ -f label-RESTART.wfn ] && grep -q "MINIMIZER DIIS" input.inp; then
    mv output.out output.diis.out
    cp label-RESTART.wfn start.wfn
    sed -e 's/SCF_GUESS ATOMIC/SCF_GUESS RESTART/' \\
        -e 's/MINIMIZER DIIS/MINIMIZER CG\\n        LINESEARCH 2PNT/' \\
        -e 's/PRECONDITIONER FULL_SINGLE_INVERSE/PRECONDITIONER FULL_ALL/' \\
        -e 's/^  &DFT$/  \\&DFT\\n    WFN_RESTART_FILE_NAME start.wfn/' input.inp > input.cg.inp
    if cp2k input.cg.inp; then touch done; echo "RECOVERED $d"; else touch failed; echo "FAILED $d"; fi
  else
    touch failed; echo "FAILED $d"
  fi
  cd ..
done
echo "ALL_LABELS_FINISHED"
"""

# Labels flagged by screen_labels.py are run again from the atomic guess with OT CG + FULL_ALL
RELABEL_SCRIPT = """#!/bin/bash
# Usage inside hap-cp2k: relabel.sh <label_dir> <mpi_ranks> <omp_threads> [list]; frames listed in [list]
# (default relabel.txt)
# Not a 'while read' loop: mpirun reads stdin and swallowed the rest of the list (2026-11-05, iter-000)
cd "$1" || exit 1
for d in $(cat "${4:-relabel.txt}"); do
  [ -f "$d/done" ] && continue
  ( cd "$d" && OMP_NUM_THREADS=$3 mpirun -np $2 cp2k.psmp -i input.atomcg.inp -o output.out < /dev/null \\
      && grep -q "PROGRAM ENDED AT" output.out && ! grep -q "SCF run NOT converged" output.out \\
      && touch done && echo "RELABELLED $d" ) || { touch "$d/failed"; echo "FAILED $d"; }
done
echo "RELABEL_FINISHED"
"""


# Labels only: converged labels of iterations 2-3 needed 2 outer SCF loops (455 frames) or 3 (2 frames). A frame
# whose OT DIIS is still unconverged after OUTER_MAX_SCF loops is stuck (iter-004 candidates_00065 sat at the same
# energy for 8 loops on 2026-11-05) and goes to the CG fallback, instead of using all 20 loops of cp2k_protocol.
OUTER_MAX_SCF = 6


def cap_outer_scf(text):
    old = "        MAX_SCF 20\n      &END OUTER_SCF\n"
    if text.count(old) != 1:
        raise SystemExit("cannot cap the outer SCF: 'MAX_SCF 20' not found once in &OUTER_SCF")
    return text.replace(old, f"        MAX_SCF {OUTER_MAX_SCF}\n      &END OUTER_SCF\n")


def robust_scf(text):
    """The same input with OT CG + FULL_ALL in place of OT DIIS + FULL_SINGLE_INVERSE."""
    for old, new in (("MINIMIZER DIIS\n", "MINIMIZER CG\n        LINESEARCH 2PNT\n"),
                     ("PRECONDITIONER FULL_SINGLE_INVERSE\n", "PRECONDITIONER FULL_ALL\n")):
        if text.count(old) != 1:
            raise SystemExit(f"cannot switch the SCF: '{old.strip()}' not found once")
        text = text.replace(old, new)
    return text


def relabel(args):
    """Set up relabelling of the frames screen_labels.py flagged in --screen; their first output is kept as
    output.suspect.out and the screen as <screen>.first.json. A frame is relabelled once (a relabel that was
    interrupted, neither done nor failed, is listed again); if it is still flagged afterwards, collect
    --exclude-screen leaves it out."""
    screen = json.loads((args.label_dir / args.screen).read_text())
    (args.label_dir / f"{Path(args.screen).stem}.first.json").write_text(json.dumps(screen, indent=2) + "\n")
    todo = []
    for name in screen["flagged"]:
        d = args.label_dir / name
        if (d / "output.suspect.out").exists():
            if not ((d / "done").exists() or (d / "failed").exists()):
                todo.append(name)
            continue
        (d / "output.out").rename(d / "output.suspect.out")
        (d / "done").unlink()
        (d / "input.atomcg.inp").write_text(robust_scf((d / "input.inp").read_text()), newline="\n")
        todo.append(name)
    (args.label_dir / args.list).write_text("".join(f"{n}\n" for n in todo), newline="\n")
    (args.label_dir / "relabel.sh").write_text(RELABEL_SCRIPT, newline="\n")
    print(f"{len(todo)} frames to relabel in {args.label_dir}")


def scripts(args):
    """Write the current run_labels.sh. Replaced by a rename: bash reads a running script as it goes, and a seed
    container or another array task may still be running the old one."""
    tmp = args.label_dir / f"run_labels.sh.{os.getpid()}.new"
    tmp.write_text(RUN_SCRIPT, newline="\n")  # bash rejects CRLF
    tmp.replace(args.label_dir / "run_labels.sh")


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
            text = cap_outer_scf(cp2k_input(atoms, "label", args.cutoff, args.rel_cutoff))
            if periodic_bulk:
                text = text.replace("  METHOD Quickstep\n", "  METHOD Quickstep\n  STRESS_TENSOR ANALYTICAL\n")
                text = text.replace("    &FORCES ON\n    &END FORCES\n",
                                    "    &FORCES ON\n    &END FORCES\n    &STRESS_TENSOR ON\n    &END STRESS_TENSOR\n")
            (d / "input.inp").write_text(text, newline="\n")
            write(d / "frame.extxyz", atoms)
            index.append(dict(name=name, source=str(path), frame=k, natoms=len(atoms),
                              system=atoms.info.get("system", path.stem), stress=periodic_bulk))
    scripts(args)
    with open(args.label_dir / "index.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(index[0].keys()))
        writer.writeheader()
        writer.writerows(index)
    print(f"prepared {len(index)} frames in {args.label_dir}")


def collect(args):
    from ase.io import read
    rows = list(csv.DictReader(open(args.label_dir / "index.csv")))
    screens = [args.label_dir / s for s in args.screen] if args.exclude_screen else []
    outliers = {n for s in screens if s.exists() for n in json.loads(s.read_text())["flagged"]}
    audit, groups = [], {}
    for r in rows:
        d = args.label_dir / r["name"]
        status = "missing"
        if r["name"] in outliers:
            status = "energy_outlier"
        elif (d / "done").exists():
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
    c.add_argument("--exclude-screen", action="store_true",
                   help="leave out the frames still flagged in any of the --screen files")
    c.add_argument("--screen", nargs="+", default=["energy_screen.json"], help="screen results inside --label-dir")
    r = sub.add_parser("relabel")
    r.add_argument("--label-dir", type=Path, required=True)
    r.add_argument("--screen", default="energy_screen.json", help="screen result inside --label-dir")
    r.add_argument("--list", default="relabel.txt", help="frame list written for relabel.sh, inside --label-dir")
    s = sub.add_parser("scripts", help="rewrite run_labels.sh of a prepared directory with the current version")
    s.add_argument("--label-dir", type=Path, required=True)
    args = parser.parse_args()
    {"prepare": prepare, "collect": collect, "relabel": relabel, "scripts": scripts}[args.cmd](args)


if __name__ == "__main__":
    main()
