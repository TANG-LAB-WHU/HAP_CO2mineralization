"""CP2K labelling-protocol tests: cutoff convergence, Gamma-supercell size and GPU/CPU timing.

Every input is a single-point ENERGY_FORCE with PBE-D3(BJ), DZVP-MOLOPT-SR-GTH and GTH-PBE
pseudopotentials, Gamma point only, OT. Bulk cells are replicated until every perpendicular
width reaches the requested length, which is how small MatterGen cells and AIMD frames share
one protocol.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
VALENCE = {"Ca": 10, "P": 5, "O": 6, "H": 1, "C": 4}


def hap_bulk():
    from pymatgen.core import Structure
    from pymatgen.io.ase import AseAtomsAdaptor
    s = Structure.from_file(ROOT / "3.ML_AIMD/Step4_lammps_scaleup/HAP_ConventionalCell.cif")
    allowed = {"H": [0.0608, 0.5608], "O": [0.1978, 0.6978]}
    keep = [x for x in s if x.is_ordered or any(abs(x.frac_coords[2] - z) < 1e-4 for z in allowed[next(iter(x.species)).symbol])]
    ordered = Structure(s.lattice, [next(iter(x.species)).symbol for x in keep], [x.frac_coords for x in keep])
    atoms = AseAtomsAdaptor.get_atoms(ordered)
    assert len(atoms) == 44 and atoms.get_chemical_formula() == "H2Ca10O26P6"
    return atoms


def repeats_for(atoms, min_width):
    cell = atoms.cell.array
    vol = abs(np.linalg.det(cell))
    widths = [vol / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3)]
    return [max(1, int(np.ceil(min_width / w - 1e-6))) for w in widths]


def cp2k_input(atoms, project, cutoff=600, rel_cutoff=60, eps_scf="1e-7"):
    kinds = "\n".join(
        f"    &KIND {el}\n      BASIS_SET DZVP-MOLOPT-SR-GTH\n      POTENTIAL GTH-PBE-q{VALENCE[el]}\n    &END KIND"
        for el in sorted(set(atoms.get_chemical_symbols())))
    cell = "\n".join(f"      {v} {x:.10f} {y:.10f} {z:.10f}" for v, (x, y, z) in zip("ABC", atoms.cell.array))
    coords = "\n".join(f"      {a.symbol} {a.position[0]:.10f} {a.position[1]:.10f} {a.position[2]:.10f}" for a in atoms)
    return f"""&GLOBAL
  PROJECT {project}
  RUN_TYPE ENERGY_FORCE
  PRINT_LEVEL LOW
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    &MGRID
      CUTOFF {cutoff}
      REL_CUTOFF {rel_cutoff}
      NGRIDS 5
    &END MGRID
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF {eps_scf}
      MAX_SCF 50
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
      &OUTER_SCF
        EPS_SCF {eps_scf}
        MAX_SCF 20
      &END OUTER_SCF
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
          R_CUTOFF 15.0
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
  &END DFT
  &SUBSYS
    &CELL
{cell}
      PERIODIC XYZ
    &END CELL
    &COORD
{coords}
    &END COORD
{kinds}
  &END SUBSYS
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
&END FORCE_EVAL
"""


NUM = r"(-?\d+\.\d+(?:[EeDd][-+]?\d+)?)"
HA_EV = 27.211386245988
HA_BOHR_EV_A = 51.42208619083232


def parse_output(path):
    """Energy (Ha), forces (Ha/bohr), stress (GPa) and wall time from a CP2K 2026.x output."""
    text = Path(path).read_text(errors="replace")
    energy = re.findall(r"ENERGY\| Total FORCE_EVAL \( QS \) energy \[(?:a\.u\.|hartree)\]:?\s+" + NUM, text)
    scf_fail = "SCF run NOT converged" in text
    forces = [list(map(float, m.groups())) for m in re.finditer(
        r"^ FORCES\|\s+\d+\s+" + NUM + r"\s+" + NUM + r"\s+" + NUM + r"\s+" + NUM + r"\s*$", text, re.M)]
    forces = [f[:3] for f in forces]
    stress = None
    # CP2K 2026.2 prints the tensor in bar by default; older inputs may request GPa
    block = re.search(r"STRESS\| Analytical stress tensor \[(GPa|bar)\]\s*\n.*?\n((?: STRESS\|\s+[xyz].*\n){3})", text)
    if block:
        scale = 1.0 if block.group(1) == "GPa" else 1e-4
        stress = [[v * scale for v in map(float, re.findall(NUM, line)[-3:])] for line in block.group(2).splitlines()]
    timing = re.search(r"^ CP2K\s+1\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+[\d.]+\s+([\d.]+)", text, re.M)
    return dict(energy_Ha=float(energy[-1]) if energy else None, scf_converged=not scf_fail and bool(energy),
                forces_Ha_bohr=forces, stress_GPa=stress, wall_s=float(timing.group(1)) if timing else None,
                finished="PROGRAM ENDED AT" in text)


def main():
    parser = argparse.ArgumentParser(description="Write CP2K protocol-test inputs")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    from ase.io import write
    bulk = hap_bulk()
    plan = []
    # Cutoff series on the 88-atom 1x1x2 cell (all widths >= ~8 A).
    cell88 = bulk.repeat([1, 1, 2])
    for cutoff in (400, 500, 600, 800, 1000):
        for rel in (60,) if cutoff != 600 else (40, 60, 80):
            plan.append((f"cut{cutoff}_rel{rel}", cell88, cutoff, rel))
    # Gamma-only supercell convergence at the default cutoff.
    for width in (8, 12, 16):
        rep = repeats_for(bulk, width)
        plan.append((f"width{width}_{'x'.join(map(str, rep))}", bulk.repeat(rep), 600, 60))
    manifest = []
    for name, atoms, cutoff, rel in plan:
        d = args.out / name
        d.mkdir(parents=True, exist_ok=True)
        (d / "input.inp").write_text(cp2k_input(atoms, name, cutoff, rel))
        write(d / "structure.extxyz", atoms)
        manifest.append(dict(name=name, natoms=len(atoms), cutoff=cutoff, rel_cutoff=rel))
    (args.out / "plan.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=1))


if __name__ == "__main__":
    main()
