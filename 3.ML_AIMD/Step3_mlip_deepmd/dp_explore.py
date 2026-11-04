"""Active-learning exploration with LAMMPS + DeePMD on the GPU (DP-GEN style).

For each starting structure and temperature, LAMMPS runs NVT with the 4-model ensemble and
writes model_devi.out. Frames with trust_lo <= max_devi_f < trust_hi become labelling
candidates (capped per run, evenly spread in time); frames above trust_hi are counted as
failures and never labelled. Run inside hap-mlip.

With --backend tf (the default) each PyTorch model is converted to the TensorFlow backend
(`dp convert-backend` of the uncompressed graph.pth next to it) and checked against the original
before use: the 3.1.3 PyTorch backend spends most of its time in the TorchScript interpreter
(about 90 ms per step for 88 atoms on the RTX 4090; TensorFlow about 25 ms).

  python dp_explore.py --models iter-000/model-*/graph-compress.pth --structures init/*.extxyz \
      --temperatures 300 500 --steps 20000 --out iter-000/explore
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

TYPE_MAP = ["Ca", "P", "O", "H", "C"]


def eval_ef(pot, a):
    types = [TYPE_MAP.index(s) for s in a.get_chemical_symbols()]
    e, f = pot.eval(a.positions.reshape(1, -1), a.cell.array.reshape(1, -1), types)[:2]
    return float(e.item()), f.reshape(-1, 3)


def tf_models(models, out, structures, tol):
    """Convert the PyTorch models to TensorFlow and check that they give the same energies and forces."""
    from deepmd.infer import DeepPot
    out.mkdir(parents=True, exist_ok=True)
    pbs, checks = [], []
    for k, m in enumerate(models):
        # convert-backend needs the uncompressed model, which dp_ensemble writes next to the compressed one
        src = m.with_name("graph.pth") if m.name == "graph-compress.pth" else m
        pb = out / f"model-{k}.pb"
        if not pb.exists():
            subprocess.run(["dp", "convert-backend", str(src), str(pb)], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        pt, tf = DeepPot(str(src)), DeepPot(str(pb))
        de, df = 0.0, 0.0
        for a in structures:
            (e1, f1), (e2, f2) = eval_ef(pt, a), eval_ef(tf, a)
            de, df = max(de, abs(e1 - e2) / len(a)), max(df, float(np.abs(f1 - f2).max()))
        checks.append(dict(source=str(src), converted=str(pb), max_abs_dE_eV_per_atom=de, max_abs_dF_eV_A=df))
        if de > tol or df > tol:
            raise SystemExit(f"{pb}: TensorFlow model differs from {src} by {de:.2e} eV/atom, {df:.2e} eV/A")
        pbs.append(pb)
    return pbs, checks


def python_devi(pots, frames):
    """DP-GEN max_devi_f evaluated with the Python API, to cross-check the LAMMPS plugin."""
    out = []
    for a in frames:
        f = np.array([eval_ef(p, a)[1] for p in pots])
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
    parser.add_argument("--backend", choices=["tf", "pt"], default="tf", help="DeePMD backend for LAMMPS and the check")
    parser.add_argument("--convert-tol", type=float, default=1e-6, help="eV/atom and eV/A, TF vs PyTorch")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    from deepmd.infer import DeepPot
    model_files, convert = [m.resolve() for m in args.models], []
    if args.backend == "tf":
        model_files, convert = tf_models(model_files, args.out.resolve() / "models",
                                         [read(p) for p in args.structures], args.convert_tol)
    pots = [DeepPot(str(m)) for m in model_files]
    models = " ".join(map(str, model_files))
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
    # TensorFlow otherwise starts a thread per core in each job; fewer threads leave the CPU to CP2K (same speed)
    env = {"DP_INTRA_OP_PARALLELISM_THREADS": "2", "DP_INTER_OP_PARALLELISM_THREADS": "1", **os.environ}
    # A run counts as done once LAMMPS has finished it ("Loop time"); output left by an interrupted run (for
    # example a stopped orchestrator) is kept under *.interrupted-<k> and the run starts again
    for _, _, d, _ in runs:
        log = d / "log.lammps"
        if (d / "model_devi.out").exists() and not (log.exists() and "Loop time" in log.read_text(errors="replace")):
            k = 1
            while d.with_name(f"{d.name}.interrupted-{k}").exists():
                k += 1
            keep = d.with_name(f"{d.name}.interrupted-{k}")
            keep.mkdir()
            for f in ("model_devi.out", "traj.dump", "log.lammps"):
                if (d / f).exists():
                    (d / f).rename(keep / f)
            print(f"{d.name}: interrupted output kept in {keep.name}; rerunning", flush=True)
    with ThreadPoolExecutor(args.parallel) as pool:
        list(pool.map(lambda d: subprocess.run(["lmp", "-in", "in.lmp", "-log", "log.lammps"], cwd=d, check=False,
                                               env=env, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT),
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
        py = python_devi(pots, [frames[k] for k in check])
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
        complete = "Loop time" in (d / "log.lammps").read_text(errors="replace") if (d / "log.lammps").exists() else False
        summary[d.name] = dict(frames=len(fmax), lammps_complete=complete,
                               accurate=float(accurate.mean()), candidate=float(len(cand) / len(fmax)),
                               failed=float(failed.mean()), picked=len(pick), max_devi_f_max=float(fmax.max()),
                               python_check_max_abs_diff=mismatch)
        print(d.name, summary[d.name], flush=True)
    write(args.out / "candidates.extxyz", candidates)
    acc = np.mean([s["accurate"] for s in summary.values() if "accurate" in s]) if summary else 0.0
    (args.out / "explore_summary.json").write_text(json.dumps(dict(runs=summary, mean_accurate=acc,
                                                                    n_candidates=len(candidates), backend=args.backend,
                                                                    conversion_check=convert), indent=2) + "\n")
    print(f"accurate fraction {acc:.3f}; {len(candidates)} candidates")


if __name__ == "__main__":
    main()
