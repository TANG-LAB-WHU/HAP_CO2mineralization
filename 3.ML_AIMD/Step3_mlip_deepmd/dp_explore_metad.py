"""Biased exploration of CO2 conversion at the hydrated HAP(001) surface (well-tempered metadynamics).

Implements the "Biased exploration of CO2 conversion" plan in workflow/PIPELINE_V2.md (option A, approved
2026-11-05). LAMMPS + DeePMD with `fix plumed` on the TF models that dp_explore.py converted. One well-tempered
METAD per carbon biases its coordination by O atoms other than its own two (R_0 1.6 A, NN 8, MM 16), walls keep
both C between 2 and 10.5 A above the top slab Ca, and "surface" runs also pull C1 down to 3.5 A over 10 ps.

Frames are selected as in dp_explore.py (trust_lo <= max_devi_f < trust_hi, at most --max-per-run, evenly
spread in time); "test" runs instead give --test-per-run frames evenly spread in time, whatever their deviation,
for the sealed test. Every selected frame passes a geometry check and carries its CV values. No free energies
come from these runs: the bias is deposited on a model that is not converged. Run inside hap-mlip.

  python dp_explore_metad.py --models iter-003/explore/models/model-*.pb \
      --structure exploration/w8/hap001_water_co2_initial.extxyz \
      --runs free_T300_s1:free:300:1:train test_T300:free:300:777:test --steps 100000 --out metad/iter-003 \
      --selected exploration/metad/selected --test-selected exploration/metad-test/selected
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write

from dp_explore import TYPE_MAP, python_devi

SYSTEM = "hap001_water_co2_metad"  # starts with hap001, so cp2k_label.py treats it as a slab (no stress)

LAMMPS_IN = """units metal
boundary p p p
atom_style atomic
read_data conf.lmp
{masses}
pair_style deepmd {models} out_freq {stride} out_file model_devi.out
pair_coeff * *
# Skin 0, as in dp_explore.py (see "LAMMPS neighbour skin" in PIPELINE_V2.md)
neighbor 0.0 bin
neigh_modify every 1 delay 0 check yes
velocity all create {temp} {seed} mom yes rot yes dist gaussian
fix 1 all nvt temp {temp} {temp} 0.05
fix pl all plumed plumedfile plumed.dat outfile plumed.log
timestep 0.0005
thermo {stride}
thermo_style custom step temp pe etotal press vol f_pl
dump 1 all custom {stride} traj.dump id type x y z
dump_modify 1 sort id
run {steps}
"""


def groups(a):
    """0-based atom groups: the two C, their O, all other O, all H, and the top slab Ca (height reference)."""
    s = np.array(a.get_chemical_symbols())
    carbons = [int(k) for k in np.where(s == "C")[0]]
    oc = sorted(int(j) for c in carbons for j in np.where(s == "O")[0] if a.get_distance(c, j, mic=True) < 1.4)
    onon = [int(k) for k in np.where(s == "O")[0] if k not in oc]
    hydrogens = [int(k) for k in np.where(s == "H")[0]]
    ref = int(np.where(s == "Ca")[0][np.argmax(a.positions[s == "Ca", 2])])
    if len(carbons) != 2 or len(oc) != 4:
        raise SystemExit(f"expected 2 CO2, found C {carbons} with O {oc}")
    return carbons, oc, onon, hydrogens, ref


def plumed_input(a, kind, temp, steps, args):
    carbons, oc, onon, hydrogens, ref = groups(a)
    one = lambda ix: ",".join(str(k + 1) for k in ix)  # PLUMED indices are 1-based
    h1 = a.positions[carbons[0], 2] - a.positions[ref, 2]
    metad = "".join(
        f"m{k}: METAD ARG=cn{k} SIGMA={args.sigma} HEIGHT={args.height} PACE={args.pace} BIASFACTOR={args.biasfactor} "
        f"TEMP={temp:g} GRID_MIN=-0.5 GRID_MAX=3.0 GRID_BIN=350 FILE=HILLS{k}\n" for k in (1, 2))
    pull, printed = "", "m1.bias,m2.bias,uw.bias,lw.bias"
    if kind == "surface":
        # Bring C1 into the first water layer over 10 ps, then hold it there; its METAD keeps acting
        pull = (f"pull: MOVINGRESTRAINT ARG=d1.z STEP0=0 AT0={h1:.3f} KAPPA0={args.pull_kappa} "
                f"STEP1={args.pull_steps} AT1={args.pull_to} KAPPA1={args.pull_kappa} "
                f"STEP2={steps} AT2={args.pull_to} KAPPA2={args.pull_kappa}\n")
        printed += ",pull.bias"
    return f"""UNITS LENGTH=A ENERGY=eV TIME=ps
# New C-O bonds of each C to O other than its own two: about 0.01 for CO2 in water, 0.7-0.9 for a bonded O
cn1: COORDINATION GROUPA={carbons[0] + 1} GROUPB={one(onon)} R_0=1.6 NN=8 MM=16
cn2: COORDINATION GROUPA={carbons[1] + 1} GROUPB={one(onon)} R_0=1.6 NN=8 MM=16
# H on the four CO2 oxygens (HCO3- or H2CO3), printed only
ch: COORDINATION GROUPA={one(oc)} GROUPB={one(hydrogens)} R_0=1.2 NN=6 MM=12
# Height of each C above the top slab Ca (atom {ref + 1}); DISTANCE gives the second atom minus the first
d1: DISTANCE ATOMS={ref + 1},{carbons[0] + 1} COMPONENTS
d2: DISTANCE ATOMS={ref + 1},{carbons[1] + 1} COMPONENTS
uw: UPPER_WALLS ARG=d1.z,d2.z AT={args.wall_hi},{args.wall_hi} KAPPA=2.0,2.0 EXP=2,2
lw: LOWER_WALLS ARG=d1.z,d2.z AT={args.wall_lo},{args.wall_lo} KAPPA=2.0,2.0 EXP=2,2
{metad}{pull}PRINT ARG=cn1,cn2,ch,d1.z,d2.z,{printed} STRIDE={args.stride} FILE=COLVAR
"""


def selection_gate(a):
    """Failures that make a frame unphysical. Its H limit (1.45 A) keeps a proton in transfer; stage 7 used 1.2 A
    until 2026-11-06 and now uses the same limit. The stage-7 result is recorded separately."""
    from ase.neighborlist import neighbor_list
    sym = np.array(a.get_chemical_symbols())
    i, j, d = neighbor_list("ijd", a, 2.0)
    fails = []
    if len(d) and d.min() <= 0.75:
        fails.append(f"min distance {d.min():.2f}")
    for k in np.where(sym == "P")[0]:
        n = int(np.sum((i == k) & (sym[j] == "O") & (d < 1.9)))
        if n != 4:
            fails.append(f"P{k} {n} O")
    for k in np.where(sym == "C")[0]:
        m = d[(i == k) & (sym[j] == "O")]
        if np.sum(m < 1.6) not in (2, 3) or (len(m) and m.min() < 1.1):
            fails.append(f"C{k} O {np.round(np.sort(m)[:4], 2).tolist()}")
    for k in np.where(sym == "H")[0]:
        m = d[(i == k) & (sym[j] == "O")]
        if not len(m) or not 0.8 < m.min() < 1.45:
            fails.append(f"H{k} nearest O {m.min() if len(m) else float('inf'):.2f}")
    return fails


def read_colvar(path):
    names = open(path).readline().split()[2:]
    return names, np.loadtxt(path, comments="#", ndmin=2)


def main():
    parser = argparse.ArgumentParser(description="Biased (WT-metadynamics) exploration of the hydrated interface")
    parser.add_argument("--models", nargs=4, type=Path, required=True, help="TF models (.pb) from dp_explore.py")
    parser.add_argument("--structure", type=Path, required=True)
    parser.add_argument("--runs", nargs="+", required=True, help="name:kind:temperature:seed:role, kind free|surface, "
                                                                  "role train|test")
    parser.add_argument("--steps", type=int, default=100_000, help="0.5 fs steps")
    parser.add_argument("--stride", type=int, default=50)
    parser.add_argument("--sigma", type=float, default=0.04)
    parser.add_argument("--height", type=float, default=0.02, help="eV")
    parser.add_argument("--pace", type=int, default=50)
    parser.add_argument("--biasfactor", type=float, default=20)
    parser.add_argument("--wall-lo", type=float, default=2.0, help="A above the top slab Ca")
    parser.add_argument("--wall-hi", type=float, default=10.5)
    parser.add_argument("--pull-to", type=float, default=3.5)
    parser.add_argument("--pull-steps", type=int, default=20_000)
    parser.add_argument("--pull-kappa", type=float, default=0.5, help="eV/A^2")
    parser.add_argument("--trust-lo", type=float, default=0.10)
    parser.add_argument("--trust-hi", type=float, default=0.30)
    parser.add_argument("--max-per-run", type=int, default=30)
    parser.add_argument("--test-per-run", type=int, default=10)
    parser.add_argument("--n-check", type=int, default=5)
    parser.add_argument("--check-tol", type=float, default=0.005, help="eV/A")
    parser.add_argument("--parallel", type=int, default=3)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--selected", type=Path, required=True, help="directory for the training frames")
    parser.add_argument("--test-selected", type=Path, required=True, help="directory for the sealed-test frames")
    args = parser.parse_args()
    from deepmd.infer import DeepPot
    from dp_validate import geometry_gate
    model_files = [m.resolve() for m in args.models]
    pots = [DeepPot(str(m)) for m in model_files]
    atoms0 = read(args.structure)
    masses = "\n".join(f"mass {i + 1} {atomic_masses[atomic_numbers[s]]:.4f}" for i, s in enumerate(TYPE_MAP))
    runs = []
    for spec in args.runs:
        name, kind, temp, seed, role = spec.split(":")
        if kind not in ("free", "surface") or role not in ("train", "test"):
            raise SystemExit(f"bad run spec {spec}")
        d = (args.out / name).resolve()
        d.mkdir(parents=True, exist_ok=True)
        write(d / "conf.lmp", atoms0, format="lammps-data", specorder=TYPE_MAP, masses=False)
        (d / "plumed.dat").write_text(plumed_input(atoms0, kind, float(temp), args.steps, args))
        (d / "in.lmp").write_text(LAMMPS_IN.format(masses=masses, models=" ".join(map(str, model_files)),
                                                   stride=args.stride, temp=float(temp),
                                                   seed=20261105 + 1000 * int(seed) + int(temp), steps=args.steps))
        runs.append((name, kind, float(temp), int(seed), role, d))
    env = {"DP_INTRA_OP_PARALLELISM_THREADS": "2", "DP_INTER_OP_PARALLELISM_THREADS": "1", **os.environ}
    done = lambda d: (d / "log.lammps").exists() and "Loop time" in (d / "log.lammps").read_text(errors="replace")
    with ThreadPoolExecutor(args.parallel) as pool:
        list(pool.map(lambda d: subprocess.run(["lmp", "-in", "in.lmp", "-log", "log.lammps"], cwd=d, check=False,
                                               env=env, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT),
                      [r[-1] for r in runs if not done(r[-1])]))
    symbols = atoms0.get_chemical_symbols()
    summary, train, test = {}, [], []
    for name, kind, temp, seed, role, d in runs:
        if not (d / "model_devi.out").exists() or not (d / "COLVAR").exists():
            summary[name] = {"status": "lammps_failed"}
            continue
        dev = np.loadtxt(d / "model_devi.out", ndmin=2)
        cvn, cv = read_colvar(d / "COLVAR")
        frames = read(d / "traj.dump", ":", format="lammps-dump-text")
        n = min(len(frames), len(dev), len(cv))
        frames, step, fmax, cv = frames[:n], dev[:n, 0].astype(int), dev[:n, 4], cv[:n]
        for f in frames:
            f.set_chemical_symbols(symbols)
        check = np.unique(np.linspace(0, n - 1, args.n_check).round().astype(int))
        mismatch = float(np.abs(python_devi(pots, [frames[k] for k in check]) - fmax[check]).max())
        if mismatch > args.check_tol:
            raise SystemExit(f"{name}: LAMMPS max_devi_f differs from Python by {mismatch:.4f} eV/A")
        col = {c: cv[:, k] for k, c in enumerate(cvn)}
        gate = [selection_gate(f) for f in frames]
        ok = np.array([not g for g in gate])
        if role == "train":
            pool_idx = np.where((fmax >= args.trust_lo) & (fmax < args.trust_hi) & ok)[0]
            cap = args.max_per_run
        else:
            pool_idx = np.where(ok)[0]
            cap = args.test_per_run
        pick = pool_idx[np.unique(np.linspace(0, len(pool_idx) - 1, min(cap, len(pool_idx))).round().astype(int))] \
            if len(pool_idx) else []
        for k in pick:
            f = frames[k]
            f.calc = None
            f.info = {"system": SYSTEM, "temperature_K": temp, "step": int(step[k]), "max_devi_f": float(fmax[k]),
                      "run": name, "kind": kind, "stage7_gate_pass": not geometry_gate(f),
                      **{c.replace(".", "_"): float(col[c][k]) for c in cvn if c != "time"}}
            (train if role == "train" else test).append(f)
        cn = np.maximum(col["cn1"], col["cn2"])
        complete = done(d)
        summary[name] = dict(kind=kind, temperature_K=temp, seed=seed, role=role, frames=n, lammps_complete=complete,
                             accurate=float((fmax < args.trust_lo).mean()),
                             candidate=float(((fmax >= args.trust_lo) & (fmax < args.trust_hi)).mean()),
                             failed=float((fmax >= args.trust_hi).mean()), max_devi_f_max=float(fmax.max()),
                             geometry_fail=float((~ok).mean()), picked=len(pick),
                             cn_max=float(cn.max()), frames_cn_above_0_5=int((cn > 0.5).sum()),
                             ch_max=float(col["ch"].max()),
                             height_range=[float(min(col["d1.z"].min(), col["d2.z"].min())),
                                           float(max(col["d1.z"].max(), col["d2.z"].max()))],
                             final_bias_eV=float(col["m1.bias"][-1] + col["m2.bias"][-1]),
                             python_check_max_abs_diff=mismatch)
        print(name, summary[name], flush=True)
    for frames, out in ((train, args.selected), (test, args.test_selected)):
        out.mkdir(parents=True, exist_ok=True)
        if frames:
            write(out / f"{SYSTEM}_selected.extxyz", frames)
    (args.out / "metad_summary.json").write_text(json.dumps(dict(
        runs=summary, n_train=len(train), n_test=len(test), models=[str(m) for m in model_files],
        settings={k: v for k, v in vars(args).items() if k not in ("models", "runs")}), indent=2, default=str) + "\n")
    print(f"{len(train)} training frames, {len(test)} sealed-test frames")


if __name__ == "__main__":
    main()
