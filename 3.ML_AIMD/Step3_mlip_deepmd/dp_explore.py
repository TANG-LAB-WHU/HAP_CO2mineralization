"""Active-learning exploration with LAMMPS + DeePMD on the GPU (DP-GEN style).

For each starting structure and temperature, LAMMPS runs NVT with the 4-model ensemble and
writes model_devi.out. Frames with trust_lo <= max_devi_f < trust_hi become labelling
candidates (capped per run, evenly spread in time); frames above trust_hi are counted as
failures and never labelled. Run inside hap-mlip.

  python dp_explore.py --models iter-000/model-*/graph-compress.pth --structures init/*.extxyz \
      --temperatures 300 500 --steps 20000 --out iter-000/explore
"""
from __future__ import annotations

import argparse
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write

TYPE_MAP = ["Ca", "P", "O", "H", "C"]


def python_devi(models, frames):
    """DP-GEN max_devi_f evaluated with the Python API, to cross-check the LAMMPS plugin."""
    from deepmd.infer import DeepPot
    pots = [DeepPot(str(m)) for m in models]
    out = []
    for a in frames:
        types = [TYPE_MAP.index(s) for s in a.get_chemical_symbols()]
        f = np.array([p.eval(a.positions.reshape(1, -1), a.cell.array.reshape(1, -1), types)[1].reshape(-1, 3)
                      for p in pots])
        out.append(float(np.sqrt(((f - f.mean(0)) ** 2).sum(-1).mean(0)).max()))
    return np.array(out)

LAMMPS_IN = """units metal
boundary p p p
atom_style atomic
read_data conf.lmp
{masses}
pair_style deepmd {models} out_freq {stride} out_file model_devi.out
pair_coeff * *
# With any neighbour skin the DeePMD-kit 3.1.3 PyTorch plugin also counts neighbours between rcut and
# rcut + skin: forces are off by 1-3 eV/A at the metal default of 2 A (compressed models even at 0.3 A).
# Skin 0 matches the Python API to 1e-4 eV/A; the list is then rebuilt every step
neighbor 0.0 bin
neigh_modify every 1 delay 0 check yes
velocity all create {temp} {seed} mom yes rot yes dist gaussian
fix 1 all nvt temp {temp} {temp} 0.05
timestep 0.0005
thermo {stride}
thermo_style custom step temp pe etotal press vol
dump 1 all custom {stride} traj.dump id type x y z
dump_modify 1 sort id
run {steps}
"""


def main():
    parser = argparse.ArgumentParser(description="LAMMPS-DeePMD model-deviation exploration")
    parser.add_argument("--models", nargs=4, type=Path, required=True)
    parser.add_argument("--structures", nargs="+", type=Path, required=True)
    parser.add_argument("--temperatures", nargs="+", type=float, default=[300, 500])
    parser.add_argument("--steps", type=int, default=20000, help="0.5 fs steps")
    parser.add_argument("--stride", type=int, default=50)
    parser.add_argument("--trust-lo", type=float, default=0.10)
    parser.add_argument("--trust-hi", type=float, default=0.30)
    parser.add_argument("--max-per-run", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--n-check", type=int, default=5, help="frames re-evaluated in Python per run")
    parser.add_argument("--check-tol", type=float, default=0.005, help="eV/A")
    parser.add_argument("--parallel", type=int, default=1, help="LAMMPS runs at once on the GPU")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    models = " ".join(str(m.resolve()) for m in args.models)
    masses = "\n".join(f"mass {i + 1} {atomic_masses[atomic_numbers[s]]:.4f}" for i, s in enumerate(TYPE_MAP))
    summary, candidates = {}, []
    runs = []
    for path in args.structures:
        atoms0 = read(path)
        name = path.stem.replace("_initial", "")
        for temp in args.temperatures:
            d = (args.out / f"{name}_T{int(temp)}").resolve()
            d.mkdir(parents=True, exist_ok=True)
            write(d / "conf.lmp", atoms0, format="lammps-data", specorder=TYPE_MAP, masses=False)
            (d / "in.lmp").write_text(LAMMPS_IN.format(masses=masses, models=models, stride=args.stride,
                                                       temp=temp, seed=args.seed + int(temp), steps=args.steps))
            runs.append((name, temp, d, atoms0))
    # The runs are independent (each has its own seed); one leaves the GPU mostly idle, so --parallel share it.
    # Analysis follows in input order
    with ThreadPoolExecutor(args.parallel) as pool:
        list(pool.map(lambda d: subprocess.run(["lmp", "-in", "in.lmp", "-log", "log.lammps"], cwd=d, check=False,
                                               stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT),
                      [d for _, _, d, _ in runs if not (d / "model_devi.out").exists()]))
    for name, temp, d, atoms0 in runs:
        if not (d / "model_devi.out").exists():
            summary[d.name] = {"status": "lammps_failed"}
            continue
        dev = np.loadtxt(d / "model_devi.out", ndmin=2)
        step, fmax = dev[:, 0].astype(int), dev[:, 4]
        frames = read(d / "traj.dump", ":", format="lammps-dump-text")
        n = min(len(frames), len(fmax))  # a crashed run can leave one extra line in either file
        frames, step, fmax = frames[:n], step[:n], fmax[:n]
        symbols = atoms0.get_chemical_symbols()
        # Re-evaluate a few frames in Python; a mismatch means the LAMMPS forces are wrong
        check = np.unique(np.linspace(0, len(frames) - 1, args.n_check).round().astype(int))
        for f in frames:
            f.set_chemical_symbols(symbols)
        py = python_devi(args.models, [frames[k] for k in check])
        mismatch = float(np.abs(py - fmax[check]).max())
        if mismatch > args.check_tol:
            raise SystemExit(f"{d.name}: LAMMPS max_devi_f differs from Python by {mismatch:.4f} eV/A")
        accurate = fmax < args.trust_lo
        cand = np.where((fmax >= args.trust_lo) & (fmax < args.trust_hi))[0]
        failed = fmax >= args.trust_hi
        pick = cand[np.unique(np.linspace(0, len(cand) - 1, min(args.max_per_run, len(cand))).round().astype(int))] if len(cand) else []
        for k in pick:
            f = frames[k]
            f.calc = None
            f.info = {"system": name, "temperature_K": temp, "step": int(step[k]), "max_devi_f": float(fmax[k])}
            candidates.append(f)
        summary[d.name] = dict(frames=len(fmax), accurate=float(accurate.mean()), candidate=float(len(cand) / len(fmax)),
                               failed=float(failed.mean()), picked=len(pick), max_devi_f_max=float(fmax.max()),
                               python_check_max_abs_diff=mismatch)
        print(d.name, summary[d.name], flush=True)
    write(args.out / "candidates.extxyz", candidates)
    acc = np.mean([s["accurate"] for s in summary.values() if "accurate" in s]) if summary else 0.0
    (args.out / "explore_summary.json").write_text(json.dumps(dict(runs=summary, mean_accurate=acc,
                                                                    n_candidates=len(candidates)), indent=2) + "\n")
    print(f"accurate fraction {acc:.3f}; {len(candidates)} candidates")


if __name__ == "__main__":
    main()
