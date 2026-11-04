#!/usr/bin/env python3
"""
S4 gate: CP2K convex hull for the S3 survivors with the v2 labelling protocol (workflow/PIPELINE_V2.md).

Structures start from the MACE-MPA-0 geometries of S2. Each gets three jobs:
  opt    CELL_OPT in a Gamma supercell with every perpendicular width >= --opt-width (8 A, as the labels)
  final  single point at the relaxed geometry (removes the fixed-grid error of the cell optimisation)
  fs12   single point at the MACE geometry with widths >= --fs-width (12 A); fs12 - first opt step
         is the Gamma finite-size correction, applied when available and reported separately
References are every Alex-MP phase in the CaO-P2O5-CO2-H2O space (formally neutral with Ca2+ P5+ C4+ H+ O2-)
with DFT E_hull <= --ref-ehull: tier 1 is the S2 decomposition products of the survivors,
tier 2 the rest (at most --polymorphs per formula, smallest first).

  prepare: python cp2k_hull.py prepare --screening <run>/screening --out <run>/dft-hull
  run:     python cp2k_hull.py run --out <run>/dft-hull         (host-side scheduler, see below)
  collect: python cp2k_hull.py collect --out <run>/dft-hull

The scheduler shares the CPU with the active-learning pipeline: it uses --idle-ranks while no pipeline
CP2K label or DeePMD training container runs and --busy-ranks otherwise. Excess jobs are stopped and
later resume from their CP2K restart file, so a long cell optimisation never blocks the pipeline.
"""

import argparse
import csv
import json
import re
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "workflow"))
from cp2k_protocol import HA_EV, cp2k_input, parse_output, repeats_for  # noqa: E402

CP2K_IMAGE = "cp2k/cp2k:2026.2_mpich_generic_psmp"
SURVIVORS = "s3"  # screening_summary.json key holding the passed candidate keys
OX = {"Ca": 2, "P": 5, "C": 4, "H": 1, "O": -2}

RUN_JOB = """#!/bin/bash
# Usage inside the CP2K image: run_job.sh <job_dir> <mpi_ranks>
cd "$1" || exit 1
inp=input.inp
# A stopped cell optimisation resumes from its last step
# (a restart file cut off by the stop falls back to the backup of the step before)
for f in opt-1.restart opt-1.restart.bak-1; do
  [ -f $f ] && grep -q "&END FORCE_EVAL" $f && cp $f resume.inp && inp=resume.inp && break
done
n=$(ls output_*.out 2>/dev/null | wc -l)
OMP_NUM_THREADS=1 mpirun -np $2 cp2k.psmp -i $inp -o output_$n.out && grep -q "PROGRAM ENDED AT" output_$n.out && touch done
"""

CELL_OPT = """&MOTION
  &CELL_OPT
    # CP2K 2026 removed CELL_OPT/TYPE; cell optimisations are always DIRECT_CELL_OPT
    OPTIMIZER BFGS
    MAX_ITER 400
    EXTERNAL_PRESSURE 1.01325
    PRESSURE_TOLERANCE 100.0
  &END CELL_OPT
  &PRINT
    &RESTART
      BACKUP_COPIES 2
      &EACH
        CELL_OPT 1
      &END EACH
    &END RESTART
  &END PRINT
&END MOTION
"""


def with_stress(text):
    text = text.replace("  METHOD Quickstep\n", "  METHOD Quickstep\n  STRESS_TENSOR ANALYTICAL\n")
    return text.replace("    &FORCES ON\n    &END FORCES\n",
                        "    &FORCES ON\n    &END FORCES\n    &STRESS_TENSOR ON\n    &END STRESS_TENSOR\n")


def supercell(atoms, width):
    sc = atoms.repeat(repeats_for(atoms, width))
    # CP2K 2026.2 rejects some valid cells ("No meaningful cell information"); the standard form is a rotation
    sc.set_cell(sc.cell.standard_form()[0], scale_atoms=True)
    return sc


def write_job(d, atoms, project, args, kind):
    from ase.io import write
    d.mkdir(parents=True, exist_ok=True)
    text = with_stress(cp2k_input(atoms, project, args.cutoff, args.rel_cutoff))
    if kind == "opt":
        text = text.replace("  RUN_TYPE ENERGY_FORCE\n", "  RUN_TYPE CELL_OPT\n").replace("SCF_GUESS ATOMIC", "SCF_GUESS RESTART")
        text = text.replace("&FORCE_EVAL\n", CELL_OPT + "&FORCE_EVAL\n", 1)
    (d / "input.inp").write_text(text, newline="\n")
    write(d / "start.extxyz", atoms)


def select(screening, args):
    from ase.io import read
    refs = read(screening / "s2/inputs/reference_phases.extxyz", ":")
    relaxed = {a.info["key"]: a for a in read(screening / "s2/mace-mpa-0/relaxed.extxyz", ":")}
    survivors = json.loads((ROOT / args.summary).read_text())[SURVIVORS]["passed_keys"]

    # Tier 1: phases on the S2 decomposition of any survivor, for any of the three models
    from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
    from pymatgen.core import Composition
    tier1 = set()
    for model_dir in sorted((screening / "s2").glob("*/energies.csv")):
        rows = {r["key"]: r for r in csv.DictReader(open(model_dir))}
        pd = PhaseDiagram([PDEntry(Composition(r["formula"]), float(r["e_per_atom"]) * int(r["natoms"]), name=r["key"])
                           for r in rows.values() if r["source"] == "reference_phases"])
        for k in survivors:
            r = rows[k]
            dec, _ = pd.get_decomp_and_e_above_hull(PDEntry(Composition(r["formula"]), float(r["e_per_atom"]) * int(r["natoms"])),
                                                    allow_negative=True)
            tier1 |= {e.name for e in dec}

    by_formula = defaultdict(list)
    for i, a in enumerate(refs):
        s = a.get_chemical_symbols()
        if a.info["dft_ehull"] <= args.ref_ehull and "O" in s and sum(OX[x] for x in s) == 0:
            by_formula[a.info["formula"]].append(i)
    picked = []
    for formula, ids in by_formula.items():
        # Polymorphs within 10 meV/atom of the lowest count as equal; the smallest supercell is computed first
        low = min(refs[i].info["dft_ehull"] for i in ids)
        ids.sort(key=lambda i: (refs[i].info["dft_ehull"] > low + 0.01, len(supercell(relaxed[f"reference_phases:{i}"], args.opt_width))))
        keep = [i for i in ids if f"reference_phases:{i}" in tier1]
        keep += [i for i in ids if i not in keep][:max(0, args.polymorphs - len(keep))]
        picked += keep
    picked += [i for i in args.include if i not in picked]
    out = [dict(key=k, tier=0, entry_id="", atoms=relaxed[k]) for k in survivors]
    out += [dict(key=f"reference_phases:{i}", tier=1 if f"reference_phases:{i}" in tier1 else 2,
                 entry_id=refs[i].info["entry_id"], atoms=relaxed[f"reference_phases:{i}"]) for i in sorted(picked)]
    return out


def prepare(args):
    from pymatgen.core import Composition
    args.out.mkdir(parents=True, exist_ok=True)
    index = []
    for s in select(args.screening, args):
        atoms = s["atoms"]
        tag = "cand" if s["tier"] == 0 else "ref"
        comp = Composition(atoms.get_chemical_formula())
        name = f"{tag}-{int(s['key'].split(':')[1]):03d}_{comp.reduced_composition.alphabetical_formula.replace(' ', '')}"
        opt, fs = supercell(atoms, args.opt_width), supercell(atoms, args.fs_width)
        jobs = {"opt": len(opt)}
        write_job(args.out / name / "opt", opt, "opt", args, "opt")
        # The finite-size check is skipped for cells too large for a single workstation
        if len(fs) <= args.fs_max_atoms and s["tier"] < 2:
            write_job(args.out / name / "fs12", fs, "fs12", args, "sp")
            jobs["fs12"] = len(fs)
        index.append(dict(name=name, key=s["key"], tier=s["tier"], entry_id=s["entry_id"],
                          formula=comp.reduced_formula, natoms_cell=len(atoms), jobs=jobs))
    (args.out / "run_job.sh").write_text(RUN_JOB, newline="\n")  # bash rejects CRLF
    (args.out / "index.json").write_text(json.dumps(index, indent=1) + "\n")
    counts = defaultdict(int)
    for r in index:
        counts[r["tier"]] += 1
    print(f"prepared {len(index)} structures (tiers {dict(counts)}) in {args.out}")


def last_frame(d):
    """Relaxed supercell from the CELL_OPT restart file (cell and coordinates of the last step, in A).

    CP2K 2026.2 writes no .cell trajectory for CELL_OPT by default; the restart file holds both together.
    """
    from ase import Atoms
    text = (d / "opt-1.restart").read_text()
    subsys = text[text.index("&SUBSYS"):]
    cell_block = subsys[subsys.index("&CELL"):subsys.index("&END CELL")]
    cell = [list(map(float, re.search(rf"^\s*{v}\s+(\S+)\s+(\S+)\s+(\S+)", cell_block, re.M).groups())) for v in "ABC"]
    coords = [ln.split() for ln in subsys[subsys.index("&COORD"):subsys.index("&END COORD")].splitlines()[1:] if ln.strip()]
    return Atoms([c[0] for c in coords], positions=[list(map(float, c[1:4])) for c in coords], cell=cell, pbc=True)


def outputs(d):
    return sorted(d.glob("output_*.out"), key=lambda p: int(p.stem.split("_")[1]))


# ---------------------------------------------------------------- scheduler

def docker(*a):
    return subprocess.run(["docker", *a], capture_output=True, text=True)


def pipeline_busy():
    """True while the pipeline labels with CP2K or trains DeePMD (both need the CPU)."""
    for line in docker("ps", "--format", "{{.Names}} {{.Image}}").stdout.splitlines():
        name, image = line.split()[:2]
        if name.startswith("hull-"):
            continue
        if "cp2k" in image or name.startswith("dp-"):
            return True
    return False


def job_list(out):
    index = json.loads((out / "index.json").read_text())
    jobs = []
    for r in index:
        for kind, natoms in r["jobs"].items():
            # Priority: final single points, then survivor/tier-1 optimisations, tier-2, finite-size checks
            prio = {"opt": 1 if r["tier"] < 2 else 2, "fs12": 3}[kind]
            jobs.append(dict(name=r["name"], kind=kind, natoms=natoms, prio=(prio, r["tier"], natoms)))
        jobs.append(dict(name=r["name"], kind="final", natoms=r["jobs"]["opt"], prio=(0, r["tier"], r["jobs"]["opt"])))
    return sorted(jobs, key=lambda j: j["prio"])


def run(args):
    out = args.out.resolve()
    state_path = out / "scheduler_state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else {"failures": {}, "events": []}
    log = open(ROOT / "logs/dft_hull.log", "a")

    def event(msg, **kw):
        state["events"].append(dict(t=time.strftime("%Y-%m-%d %H:%M:%S"), msg=msg, **kw))
        state_path.write_text(json.dumps(state, indent=1) + "\n")
        print(time.strftime("%H:%M:%S"), msg, kw, file=log, flush=True)

    running = {}  # container -> (job, ranks, start)
    stopped = set()
    while True:
        jobs = job_list(out)
        # Reap finished containers
        for cname in list(running):
            st = docker("inspect", "-f", "{{.State.Running}} {{.State.ExitCode}}", cname).stdout.split()
            if st and st[0] == "true":
                continue
            job, ranks, _ = running.pop(cname)
            d = out / job["name"] / job["kind"]
            code = int(st[1]) if st else -1
            docker("rm", "-f", cname)
            if cname in stopped:
                stopped.discard(cname)
                event("preempted", job=cname)
            elif (d / "done").exists():
                event("done", job=cname, ranks=ranks)
            else:
                state["failures"][cname] = state["failures"].get(cname, 0) + 1
                event("failed", job=cname, exit_code=code, attempt=state["failures"][cname])
        # Final single points are written once their optimisation is done
        for j in jobs:
            d = out / j["name"]
            if j["kind"] == "final" and (d / "opt/done").exists() and not (d / "final/input.inp").exists():
                atoms = last_frame(d / "opt")
                atoms.set_cell(atoms.cell.standard_form()[0], scale_atoms=True)
                write_job(d / "final", atoms, "final", args, "sp")
        pending = [j for j in jobs
                   if (out / j["name"] / j["kind"] / "input.inp").exists()
                   and not (out / j["name"] / j["kind"] / "done").exists()
                   and state["failures"].get(f"hull-{j['name']}-{j['kind']}", 0) < args.max_failures
                   and f"hull-{j['name']}-{j['kind']}" not in running]
        if not pending and not running:
            event("all jobs finished or failed")
            return
        budget = args.busy_ranks if pipeline_busy() else args.idle_ranks
        used = sum(r for _, r, _ in running.values())
        # Over budget: stop the most recently started jobs; optimisations resume from their restart file
        for cname in sorted(running, key=lambda c: -running[c][2]):
            if used <= budget:
                break
            docker("stop", "-t", "10", cname)
            stopped.add(cname)
            used -= running[cname][1]
        while pending and used + args.job_ranks <= budget:
            j = pending.pop(0)
            cname = f"hull-{j['name']}-{j['kind']}"
            docker("rm", "-f", cname)
            res = docker("run", "-d", "--name", cname, "--shm-size=2g", "-v", f"{out.as_posix()}:/data", CP2K_IMAGE,
                         "bash", "/data/run_job.sh", f"/data/{j['name']}/{j['kind']}", str(args.job_ranks))
            if res.returncode != 0:
                event("start failed", job=cname, stderr=res.stderr.strip()[-300:])
                break
            running[cname] = (j, args.job_ranks, time.time())
            used += args.job_ranks
            event("started", job=cname, ranks=args.job_ranks, natoms=j["natoms"], budget=budget)
        time.sleep(args.poll)


# ---------------------------------------------------------------- collect

def hull_energy(formula, ref_formulas, ref_energies):
    """Lower convex envelope at one composition: min sum w_i e_i with sum w_i x_i = x, w >= 0.

    x are atomic fractions and e energies per atom. The references need not include the elements:
    every phase here lies in the CaO-P2O5-CO2-H2O space, so only compounds span the hull.
    Returns (energy per atom, [(ref index, atom fraction)]) or None when x is outside the references.
    """
    from pymatgen.core import Composition
    from scipy.optimize import linprog
    els = sorted({e.symbol for f in [formula, *ref_formulas] for e in Composition(f)})
    frac = lambda f: [Composition(f).get_atomic_fraction(e) for e in els]  # noqa: E731
    res = linprog(ref_energies, A_eq=np.array([frac(f) for f in ref_formulas]).T, b_eq=frac(formula),
                  bounds=(0, None), method="highs")
    if res.status != 0:
        return None
    return res.fun, [(i, w) for i, w in enumerate(res.x) if w > 1e-6]


def collect(args):
    index = json.loads((args.out / "index.json").read_text())
    rows = []
    for r in index:
        d = args.out / r["name"]
        row = dict(name=r["name"], key=r["key"], tier=r["tier"], entry_id=r["entry_id"], formula=r["formula"],
                   natoms_opt=r["jobs"]["opt"], status="pending")
        if (d / "opt/done").exists():
            texts = "".join(p.read_text(errors="replace") for p in outputs(d / "opt"))
            row["opt_converged"] = "GEOMETRY OPTIMIZATION COMPLETED" in texts
            e0 = re.search(r"ENERGY\| Total FORCE_EVAL \( QS \) energy \[(?:a\.u\.|hartree)\]:?\s+(-?\d+\.\d+)", texts)
            row["e_opt_first_eV_atom"] = float(e0.group(1)) * HA_EV / r["jobs"]["opt"] if e0 else None
            row["opt_wall_h"] = round(sum((parse_output(p)["wall_s"] or 0) for p in outputs(d / "opt")) / 3600, 2)
        if (d / "final/done").exists():
            res = parse_output(outputs(d / "final")[-1])
            if res["scf_converged"]:
                row["e_final_eV_atom"] = res["energy_Ha"] * HA_EV / r["jobs"]["opt"]
                p = np.trace(np.array(res["stress_GPa"])) / 3 if res["stress_GPa"] else None
                row["final_pressure_GPa"] = round(p, 3) if p is not None else None
                row["status"] = "ok" if row.get("opt_converged") else "opt_not_converged"
            else:
                row["status"] = "final_scf_not_converged"
        if (d / "fs12/done").exists() and row.get("e_opt_first_eV_atom") is not None:
            res = parse_output(outputs(d / "fs12")[-1])
            if res["scf_converged"]:
                row["fs_correction_meV_atom"] = (res["energy_Ha"] * HA_EV / r["jobs"]["fs12"] - row["e_opt_first_eV_atom"]) * 1000
        rows.append(row)

    def energy(row, corrected):
        e = row["e_final_eV_atom"]
        return e + row["fs_correction_meV_atom"] / 1000 if corrected and "fs_correction_meV_atom" in row else e

    ok = [r for r in rows if r["status"] in ("ok", "opt_not_converged")]
    refs = [r for r in ok if r["tier"] > 0]
    for corrected, col in ((False, "ehull_meV_atom"), (True, "ehull_fs_meV_atom")):
        if not refs:
            break
        for r in ok:
            point = hull_energy(r["formula"], [x["formula"] for x in refs if x is not r],
                                [energy(x, corrected) for x in refs if x is not r])
            if point is None:  # composition outside the computed references so far
                continue
            r[col] = round((energy(r, corrected) - point[0]) * 1000, 1)
            if not corrected:
                names = [x["name"] for x in refs if x is not r]
                r["decomposition"] = " + ".join(f"{w:.2f} {names[i]}" for i, w in point[1])
    fields = ["name", "key", "tier", "entry_id", "formula", "natoms_opt", "status", "opt_converged", "opt_wall_h",
              "e_opt_first_eV_atom", "e_final_eV_atom", "final_pressure_GPa", "fs_correction_meV_atom",
              "ehull_meV_atom", "ehull_fs_meV_atom", "decomposition"]
    with open(args.out / "dft_hull.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fields)
        writer.writeheader()
        writer.writerows(rows)
    counts = defaultdict(int)
    for r in rows:
        counts[(r["tier"], r["status"])] += 1
    print({f"tier{t}:{s}": n for (t, s), n in sorted(counts.items())})
    for r in rows:
        if r["tier"] == 0 and "ehull_meV_atom" in r:
            print(f"{r['name']:28s} E_hull {r['ehull_meV_atom']:7.1f} meV/atom  (fs {r.get('ehull_fs_meV_atom', 'n/a')})  -> {r['decomposition']}")


def cells(args):
    """CP2K-relaxed cells of the hull structures that match the given structures (lattice references, dp_validate).

    A match needs the same composition and atom count, lattice vectors in the same order (lengths within 5 %,
    angles within 3 degrees of the hull start cell) and a pymatgen StructureMatcher fit.
    """
    from ase.io import read, write
    from pymatgen.analysis.structure_matcher import StructureMatcher
    from pymatgen.io.ase import AseAtomsAdaptor
    index = json.loads((args.out / "index.json").read_text())
    starts = {r["name"]: read(args.out / r["name"] / "opt/start.extxyz") for r in index}
    matcher = StructureMatcher(ltol=0.2, stol=0.3, angle_tol=5)
    found = []
    for path in args.match:
        a = read(path)
        name = path.stem.replace("_initial", "")
        hit = None
        for r in index:
            b = starts[r["name"]]
            if len(b) != len(a) or b.get_chemical_formula() != a.get_chemical_formula():
                continue
            pa, pb = a.cell.cellpar(), b.cell.cellpar()
            if np.any(np.abs(pa[:3] - pb[:3]) > 0.05 * pb[:3]) or np.any(np.abs(pa[3:] - pb[3:]) > 3.0):
                continue
            if matcher.fit(AseAtomsAdaptor.get_structure(a), AseAtomsAdaptor.get_structure(b)):
                hit = r
                break
        if hit is None:
            print(f"{name}: no matching hull structure")
            continue
        d = args.out / hit["name"] / "opt"
        if not (d / "done").exists():
            print(f"{name}: {hit['name']} cell optimisation not finished")
            continue
        converged = "GEOMETRY OPTIMIZATION COMPLETED" in "".join(o.read_text(errors="replace") for o in outputs(d))
        relaxed = last_frame(d)
        relaxed.info = dict(system=name, source=hit["name"], entry_id=hit["entry_id"], opt_converged=converged)
        found.append(relaxed)
        print(f"{name}: {hit['name']} ({hit['entry_id']}), converged {converged}, cell {relaxed.cell.cellpar().round(3).tolist()}")
    write(args.write, found)
    print(f"wrote {len(found)} reference cells to {args.write}")


def main():
    parser = argparse.ArgumentParser(description="CP2K hull for S3 survivors")
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--screening", type=Path, required=True, help="pipeline-v2/screening (S2 inputs and relaxations)")
    p.add_argument("--summary", default="0.Generative_Discovery/runs/Ca-P-O-H-C_gfscan_RTX5090_20260329/screening/screening_summary.json")
    p.add_argument("--ref-ehull", type=float, default=0.03, help="Alex-MP DFT E_hull limit for references (eV/atom)")
    p.add_argument("--polymorphs", type=int, default=1, help="tier-2 polymorphs per formula")
    p.add_argument("--include", type=int, nargs="*", default=[], help="extra reference_phases indices (e.g. aragonite)")
    p.add_argument("--opt-width", type=float, default=8.0)
    p.add_argument("--fs-width", type=float, default=12.0)
    p.add_argument("--fs-max-atoms", type=int, default=450)
    r = sub.add_parser("run")
    r.add_argument("--idle-ranks", type=int, default=24)
    r.add_argument("--busy-ranks", type=int, default=8)
    r.add_argument("--job-ranks", type=int, default=8)
    r.add_argument("--max-failures", type=int, default=2)
    r.add_argument("--poll", type=int, default=60)
    for q in (p, r):
        q.add_argument("--cutoff", type=int, default=800)
        q.add_argument("--rel-cutoff", type=int, default=60)
    c = sub.add_parser("cells")
    c.add_argument("--match", nargs="+", type=Path, required=True, help="structures to find in the hull (extxyz)")
    c.add_argument("--write", type=Path, required=True)
    for q in (p, r, c, sub.add_parser("collect")):
        q.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    {"prepare": prepare, "run": run, "collect": collect, "cells": cells}[args.cmd](args)


if __name__ == "__main__":
    main()
