"""2026-11-08: in explore-4 (iteration-4 ensemble, NPT 300 K, 1 atm) s3_055 (MatterGen CaCO3, 80 atoms) went from
6- to about 9-coordinated Ca within 25 ps and 7.5 % below its CP2K-relaxed volume; calcite in the same conditions
stayed 6-coordinated. Relax the last frame with iteration-4 model 0 (graph.pth, CPU) at 0 K, cell free, and write it
for comparison with the DFT-hull CaCO3 references (compare.py). Run in hap-mlip:
  CUDA_VISIBLE_DEVICES= OMP_NUM_THREADS=1 python relax.py"""
from pathlib import Path

import numpy as np
from ase.filters import FrechetCellFilter
from ase.io import read, write
from ase.optimize import FIRE
from deepmd.calculator import DP

run = Path("/work/.research-runs/pipeline-v2")
here = Path(__file__).resolve().parent
sym = {1: "Ca", 2: "P", 3: "O", 4: "H", 5: "C"}
calc = DP(model=str(run / "iter-004/train/model-0/graph.pth"), type_dict={"Ca": 0, "P": 1, "O": 2, "H": 3, "C": 4})
for tag, path in (("transformed", run / "iter-004/explore/s3_055_C1Ca1O3_T300_npt/traj.dump"),
                  ("start", run / "iter-004/explore/s3_055_C1Ca1O3_T300_npt/traj.dump")):
    a = read(path, index=-1 if tag == "transformed" else 0, format="lammps-dump-text")
    a.set_chemical_symbols([sym[t] for t in a.numbers])
    a.calc = calc
    FIRE(FrechetCellFilter(a), logfile=str(here / f"relax_{tag}.log")).run(fmax=0.005, steps=3000)
    e = a.get_potential_energy() / len(a)
    print(f"{tag}: E={e:.5f} eV/atom, V/fu={a.get_volume() / 16:.2f}, max|F|={np.abs(a.get_forces()).max():.4f}")
    write(here / f"s3_055_{tag}_dp0K.extxyz", a)
