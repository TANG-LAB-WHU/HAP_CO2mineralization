"""Stage 3: build exploration systems and sample them with a foundation MLIP on the GPU.

Each system is run with Langevin NVT at several temperatures plus strained/rattled copies;
frames are written every --stride steps. Output is raw configurations only (no labels);
CP2K labels are added in stage 4.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase import Atoms, units
from ase.build import molecule, surface
from ase.io import read, write
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from cp2k_protocol import hap_bulk, repeats_for  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "0.Generative_Discovery" / "screening"))
from mlip_relax import load_calculator  # noqa: E402


def water_layer(cell_xy, z0, thickness, rng, n_co2=0, spacing=3.1):
    """Fill a layer above z0 with H2O on a ~3.1 A grid (about 1 g/cm3); n_co2 sites in the top
    layer hold flat CO2 instead. Clashes left by random orientations are removed by the
    pre-relaxation in main()."""
    na, nb = (max(1, int(round(np.linalg.norm(v) / spacing))) for v in cell_xy)
    nz = max(1, int(round(thickness / spacing)))
    top_sites = [(i, j) for i in range(na) for j in range(nb)]
    co2_sites = {top_sites[k] for k in rng.choice(len(top_sites), size=n_co2, replace=False)}
    out = Atoms()
    for k in range(nz):
        for i in range(na):
            for j in range(nb):
                is_co2 = k == nz - 1 and (i, j) in co2_sites
                m = molecule("CO2" if is_co2 else "H2O")
                m.rotate(rng.uniform(0, 360), "z" if is_co2 else rng.normal(size=3))
                frac = (np.array([i, j]) + 0.5 + rng.uniform(-0.1, 0.1, 2)) / [na, nb]
                pos = frac[0] * cell_xy[0] + frac[1] * cell_xy[1] + [0, 0, z0 + (k + 0.5) * spacing]
                m.translate(pos - m.get_center_of_mass())
                out += m
    return out


def build_systems(min_width, rng):
    bulk = hap_bulk()
    systems = {"hap_bulk": bulk.repeat(repeats_for(bulk, min_width))}
    slab = surface(bulk, (0, 0, 1), 2, vacuum=0.0, periodic=True)
    rep = [max(1, int(np.ceil(min_width / l))) for l in slab.cell.lengths()[:2]] + [1]
    slab = slab.repeat(rep)
    top = slab.positions[:, 2].max()
    dry = slab.copy()
    dry.cell[2, 2] = top + 15.0
    co2 = molecule("CO2")
    co2.rotate(90, "y")
    ca_top = max((i for i, a in enumerate(slab) if a.symbol == "Ca"), key=lambda i: slab.positions[i, 2])
    co2.translate(slab.positions[ca_top] + [0, 0, 2.6] - co2.positions[0])
    systems["hap001_co2"] = dry + co2
    wet = slab.copy()
    wet.cell[2, 2] = top + 8.0 + 12.0
    systems["hap001_water_co2"] = wet + water_layer(wet.cell.array[:2], top + 1.6, 9.3, rng, n_co2=2)
    for a in systems.values():
        a.pbc = True
        a.wrap()
    return systems


def sample(atoms, calc, temperatures, steps, stride, seed, dt_fs=0.5):
    frames = []
    for t in temperatures:
        a = atoms.copy()
        a.calc = calc
        MaxwellBoltzmannDistribution(a, temperature_K=t, rng=np.random.default_rng(seed + int(t)))
        dyn = Langevin(a, dt_fs * units.fs, temperature_K=t, friction=0.02, rng=np.random.default_rng(seed))
        for _ in range(steps // stride):
            dyn.run(stride)
            f = a.copy()
            f.calc = None
            f.info = {"temperature_K": t, "seed": seed}
            frames.append(f)
    return frames


def main():
    parser = argparse.ArgumentParser(description="Foundation-MLIP exploration")
    parser.add_argument("--model", default="mace-mpa-0")
    parser.add_argument("--min-width", type=float, required=True, help="converged Gamma-point width from cp2k_protocol")
    parser.add_argument("--extra", nargs="*", type=Path, default=[], help="extra bulk extxyz (calcite, S3 survivors)")
    parser.add_argument("--only-extra", action="store_true", help="skip the built-in HAP systems")
    parser.add_argument("--temperatures", nargs="+", type=float, default=[300, 500, 700])
    parser.add_argument("--steps", type=int, default=4000)
    parser.add_argument("--stride", type=int, default=40)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    calc = load_calculator(args.model, "cuda", dtype="float32")
    systems = {} if args.only_extra else build_systems(args.min_width, rng)
    for path in args.extra:
        for k, a in enumerate(read(path, ":")):
            systems[a.info.get("system", f"{path.stem}_{k}")] = a.repeat(repeats_for(a, args.min_width))
    summary = {}
    from ase.optimize import FIRE
    for name, atoms in systems.items():
        atoms.calc = calc
        FIRE(atoms, logfile=None).run(fmax=0.1, steps=300)  # remove packing clashes before MD
        atoms.calc = None
        write(args.out / f"{name}_initial.extxyz", atoms)
        frames = sample(atoms, calc, args.temperatures, args.steps, args.stride, args.seed)
        for f in frames:
            f.info["system"] = name
        write(args.out / f"{name}_frames.extxyz", frames)
        summary[name] = dict(natoms=len(atoms), frames=len(frames), formula=atoms.get_chemical_formula())
        print(name, summary[name], flush=True)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
