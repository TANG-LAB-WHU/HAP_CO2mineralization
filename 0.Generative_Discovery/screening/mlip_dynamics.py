#!/usr/bin/env python3
"""
S3 gate: harmonic phonons and a short 300 K NPT run with one MLIP.
Input:  s2_pass.extxyz
Output: s3_dynamics.csv and s3_pass.extxyz

Pass when the lowest phonon frequency on a uniform q mesh is above --min-freq-thz
(small negative values near Gamma are acoustic noise) and the NPT run keeps every
CO3/PO4 unit intact with a mean volume change within --max-volume-change.
"""

import argparse
import csv
from pathlib import Path

import numpy as np
from ase import units
from ase.io import read, write
from ase.md.nptberendsen import NPTBerendsen
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.neighborlist import neighbor_list
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

from mlip_relax import load_calculator


def supercell_matrix(atoms, min_length=9.0):
    lengths = atoms.cell.lengths()
    return np.diag([max(1, int(np.ceil(min_length / l))) for l in lengths])


def min_phonon_frequency(atoms, calc, mesh=(8, 8, 8)):
    unit = PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.cell.array,
                        scaled_positions=atoms.get_scaled_positions())
    ph = Phonopy(unit, supercell_matrix=supercell_matrix(atoms))
    ph.generate_displacements(distance=0.01)
    forces = []
    for sc in ph.supercells_with_displacements:
        from ase import Atoms
        a = Atoms(sc.symbols, cell=sc.cell, scaled_positions=sc.scaled_positions, pbc=True)
        a.calc = calc
        forces.append(a.get_forces())
    ph.forces = forces
    ph.produce_force_constants()
    ph.run_mesh(mesh)
    return float(ph.get_mesh_dict()["frequencies"].min())


def motif_counts(atoms):
    i, j, d = neighbor_list("ijd", atoms, 1.8)
    sym = np.array(atoms.get_chemical_symbols())
    co3 = sum(1 for k in np.where(sym == "C")[0] if np.sum((sym[j[i == k]] == "O") & (d[i == k] < 1.45)) == 3)
    po4 = sum(1 for k in np.where(sym == "P")[0] if np.sum((sym[j[i == k]] == "O") & (d[i == k] < 1.75)) == 4)
    return co3, po4


def npt_check(atoms, calc, temperature=300.0, steps=10000, dt_fs=1.0, seed=7):
    md_atoms = atoms.repeat([int(n) for n in np.diag(supercell_matrix(atoms))])
    md_atoms.calc = calc
    MaxwellBoltzmannDistribution(md_atoms, temperature_K=temperature, rng=np.random.default_rng(seed))
    ref = motif_counts(md_atoms)
    v0 = md_atoms.get_volume()
    dyn = NPTBerendsen(md_atoms, timestep=dt_fs * units.fs, temperature_K=temperature,
                       pressure_au=1.01325 * units.bar, taut=100 * units.fs, taup=1000 * units.fs,
                       compressibility_au=4.57e-5 / units.bar)
    vols, broken = [], 0
    for block in range(steps // 500):
        dyn.run(500)
        vols.append(md_atoms.get_volume() / v0)
        broken += int(motif_counts(md_atoms) != ref)
    return float(np.mean(vols[len(vols) // 2:]) - 1), broken, float(md_atoms.get_temperature())


def main():
    parser = argparse.ArgumentParser(description="Phonon and NPT stability gate")
    parser.add_argument("--model", required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--min-freq-thz", type=float, default=-0.3)
    parser.add_argument("--max-volume-change", type=float, default=0.05)
    parser.add_argument("--md-steps", type=int, default=10000, help="1 fs steps; default 10 ps")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--shard", default="0/1", help="k/n: process every n-th structure starting at k")
    args = parser.parse_args()
    k_shard, n_shard = map(int, args.shard.split("/"))
    args.out_dir.mkdir(parents=True, exist_ok=True)

    calc = load_calculator(args.model, args.device)
    rows, passed = [], []
    for atoms in read(args.input, ":")[k_shard::n_shard]:
        atoms.calc = calc
        fmin = min_phonon_frequency(atoms, calc)
        # NPT only for dynamically stable cells; a phonon failure already decides the gate
        dvol, broken, t_end = (npt_check(atoms, calc, steps=args.md_steps) if fmin >= args.min_freq_thz
                               else (float("nan"), -1, float("nan")))
        ok = fmin >= args.min_freq_thz and broken == 0 and abs(dvol) <= args.max_volume_change
        rows.append(dict(key=atoms.info["key"], formula=atoms.get_chemical_formula(mode="metal"),
                         min_freq_thz=round(fmin, 3), npt_volume_change=round(dvol, 4),
                         motif_breaks=broken, t_final_K=round(t_end, 1), s3_pass=ok))
        print(rows[-1], flush=True)
        if ok:
            keep = atoms.copy()  # drop the calculator: its last results belong to the MD supercell
            keep.calc = None
            passed.append(keep)
    tag = "" if n_shard == 1 else f"_{k_shard}of{n_shard}"
    with open(args.out_dir / f"s3_dynamics{tag}.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    write(args.out_dir / f"s3_pass{tag}.extxyz", passed, format="extxyz")
    print(f"S3: {len(passed)}/{len(rows)} pass")


if __name__ == "__main__":
    main()
