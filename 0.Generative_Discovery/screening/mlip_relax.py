#!/usr/bin/env python3
"""
S2 MLIP relaxation of candidates and reference phases with one foundation model.
Run once per model inside its venv of the hap-mlip image, e.g.
  /opt/venv-mace/bin/python mlip_relax.py --model mace-mpa-0 --inputs s1_candidates.extxyz refs.extxyz --out-dir s2/mace
Output: <out-dir>/relaxed.extxyz and <out-dir>/energies.csv (one row per input structure).
Energies are raw model energies; no MP2020/TRI corrections are applied, so a hull must
use one model's energies for candidates and references alike.
"""

import argparse
import csv
import time
from pathlib import Path

import numpy as np
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import FIRE


def load_calculator(model, device, dtype="float64"):
    """float64 for relaxations and phonons; exploration MD may pass float32."""
    if model.startswith("mace"):
        from mace.calculators import mace_mp
        name = {"mace-mpa-0": "medium-mpa-0", "mace-mp-0": "medium"}[model]
        try:
            return mace_mp(model=name, default_dtype=dtype, device=device, enable_cueq=True)
        except Exception:
            return mace_mp(model=name, default_dtype=dtype, device=device)
    if model == "mattersim-5m":
        from mattersim.forcefield import MatterSimCalculator
        return MatterSimCalculator(load_path="MatterSim-v1.0.0-5M.pth", device=device)
    if model == "sevennet-mf-ompa":
        from sevenn.calculator import SevenNetCalculator
        try:  # cuEquivariance kernels are ~2x faster on the RTX 4090
            return SevenNetCalculator(model="7net-mf-ompa", modal="mpa", device=device, enable_cueq=True)
        except Exception:
            return SevenNetCalculator(model="7net-mf-ompa", modal="mpa", device=device)
    raise ValueError(model)


def main():
    parser = argparse.ArgumentParser(description="Relax structures with a foundation MLIP")
    parser.add_argument("--model", required=True, choices=["mace-mpa-0", "mace-mp-0", "mattersim-5m", "sevennet-mf-ompa"])
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--fmax", type=float, default=0.02, help="eV/A, applied to atoms and cell")
    parser.add_argument("--steps", type=int, default=600)
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    calc = load_calculator(args.model, args.device)
    done = set()
    energies_path = args.out_dir / "energies.csv"
    if energies_path.exists():
        done = {r["key"] for r in csv.DictReader(open(energies_path))}
    fields = ["key", "source", "index", "formula", "natoms", "e_initial_per_atom", "e_per_atom",
              "fmax_final", "steps", "converged", "volume_ratio", "rmsd_A", "seconds"]
    new_file = not energies_path.exists()
    with open(energies_path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fields)
        if new_file:
            writer.writeheader()
        for path in args.inputs:
            for i, atoms in enumerate(read(path, ":")):
                key = f"{path.stem}:{i}"
                if key in done:
                    continue
                t0 = time.time()
                start = atoms.copy()
                atoms.calc = calc
                e0 = atoms.get_potential_energy() / len(atoms)
                opt = FIRE(FrechetCellFilter(atoms), logfile=None)
                converged = bool(opt.run(fmax=args.fmax, steps=args.steps))
                forces = atoms.get_forces()
                frac = atoms.get_scaled_positions(wrap=False) - start.get_scaled_positions(wrap=False)
                frac -= np.round(frac)
                row = dict(
                    key=key, source=path.stem, index=i, formula=atoms.get_chemical_formula(mode="metal"),
                    natoms=len(atoms), e_initial_per_atom=e0,
                    e_per_atom=atoms.get_potential_energy() / len(atoms),
                    fmax_final=float(np.sqrt((forces ** 2).sum(1)).max()), steps=opt.nsteps,
                    converged=converged, volume_ratio=atoms.get_volume() / start.get_volume(),
                    rmsd_A=float(np.sqrt(np.mean(np.sum((frac @ atoms.cell.array) ** 2, axis=1)))),
                    seconds=round(time.time() - t0, 2),
                )
                writer.writerow(row)
                fh.flush()
                out = atoms.copy()
                out.info = {**start.info, "key": key, "model": args.model, "e_per_atom": row["e_per_atom"]}
                write(args.out_dir / "relaxed.extxyz", out, format="extxyz", append=True)
    print(f"{args.model}: done -> {energies_path}")


if __name__ == "__main__":
    main()
