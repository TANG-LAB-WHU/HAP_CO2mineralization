"""Compare the 8 A (88-atom) and 12 A (352-atom) Gamma labels of the same HAP bulk frames.

Each 12 A frame is its 8 A frame repeated 2 x 2 x 1 (`finite-size/` in the pipeline-v2 run directory), so the two
labels differ only in Gamma sampling. Every 12 A atom is mapped to its 8 A image before forces are compared.

  python workflow/finite_size_compare.py --labels .research-runs/pipeline-v2/finite-size/labels --out finite_size_pairs.csv
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from ase.io import read

from cp2k_protocol import parse_output

HA = 27.211386245988
BOHR = 0.529177210903


def compare(lab, i):
    a8 = read(lab / f"hapbulk_w8_{i:05d}/frame.extxyz")
    a12 = read(lab / f"hapbulk_w12_{i:05d}/frame.extxyz")
    r8 = parse_output(lab / f"hapbulk_w8_{i:05d}/output.out")
    r12 = parse_output(lab / f"hapbulk_w12_{i:05d}/output.out")
    if not (r8["scf_converged"] and r12["scf_converged"] and r8["finished"] and r12["finished"]):
        raise SystemExit(f"frame {i}: a label did not converge or finish")
    # Fractional coordinate of each 12 A atom in the 8 A cell, then its nearest 8 A atom
    frac8 = a8.get_scaled_positions()
    fr = np.linalg.solve(a8.cell.array.T, a12.positions.T).T % 1.0
    d = fr[:, None, :] - frac8[None, :, :]
    d -= np.rint(d)
    dist = np.linalg.norm(d @ a8.cell.array, axis=2)
    idx = dist.argmin(1)
    if a12.get_chemical_symbols() != [a8.get_chemical_symbols()[j] for j in idx]:
        raise SystemExit(f"frame {i}: 12 A atoms do not map onto 8 A atoms of the same element")
    df = np.array(r12["forces_Ha_bohr"]) * HA / BOHR - (np.array(r8["forces_Ha_bohr"]) * HA / BOHR)[idx]
    e8 = r8["energy_Ha"] * HA / len(a8)
    e12 = r12["energy_Ha"] * HA / len(a12)
    s8, s12 = np.array(r8["stress_GPa"]), np.array(r12["stress_GPa"])
    return dict(frame=i, natoms_8A=len(a8), natoms_12A=len(a12),
                dE_8A_minus_12A_meV_atom=(e8 - e12) * 1000, dF_rms_eV_A=float(np.sqrt((df ** 2).mean())),
                dF_max_abs_eV_A=float(np.abs(df).max()), dstress_rms_GPa=float(np.sqrt(((s8 - s12) ** 2).mean())),
                dpressure_GPa=float(np.trace(s8 - s12) / 3), max_mapping_distance_A=float(dist.min(1).max()))


def main():
    parser = argparse.ArgumentParser(description="8 A vs 12 A Gamma-label comparison")
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--frames", type=int, default=10)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    rows = [compare(args.labels, i) for i in range(args.frames)]
    with args.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6g}" if isinstance(v, float) else v) for k, v in r.items()})
    de = np.array([r["dE_8A_minus_12A_meV_atom"] for r in rows])
    print(f"dE 8-12 A: mean {de.mean():.2f} SD {de.std(ddof=1):.2f} range {de.min():.2f}-{de.max():.2f} meV/atom; "
          f"dF rms mean {np.mean([r['dF_rms_eV_A'] for r in rows]):.3f} eV/A")


if __name__ == "__main__":
    main()
