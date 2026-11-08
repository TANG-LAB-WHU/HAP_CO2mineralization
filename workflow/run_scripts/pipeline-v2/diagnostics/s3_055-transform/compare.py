"""2026-11-08: which DFT-hull CaCO3 reference (ref-186/187/188, CP2K-relaxed, our protocol) does the DP-relaxed
transformed s3_055 match? Space groups with spglib at several tolerances, StructureMatcher with default tolerances,
Ca-O coordination (< 3.0 A). Run on the host: python -I compare.py (from the repository root)."""
import sys
from pathlib import Path

import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

root = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(root / "0.Generative_Discovery/screening"))
import cp2k_hull  # noqa: E402

here = Path(__file__).resolve().parent
hull = here.parents[1] / "dft-hull"
structs = {t: read(here / f"s3_055_{t}_dp0K.extxyz") for t in ("transformed", "start")}
for n in ("ref-186_C1Ca1O3", "ref-187_C1Ca1O3", "ref-188_C1Ca1O3", "cand-055_C1Ca1O3"):
    structs[n] = cp2k_hull.last_frame(hull / n / "opt")


def cn(a):
    i, j = neighbor_list("ij", a, 3.0)
    s = a.get_chemical_symbols()
    sel = np.array([s[x] == "Ca" and s[y] == "O" for x, y in zip(i, j)])
    c = np.bincount(i[sel], minlength=len(a))[[k for k, x in enumerate(s) if x == "Ca"]]
    return f"{c.mean():.2f} [{c.min()}-{c.max()}]"


pm = {k: AseAtomsAdaptor.get_structure(a) for k, a in structs.items()}
for k, a in structs.items():
    nfu = a.get_chemical_symbols().count("Ca")
    sg = []
    for tol in (0.01, 0.1, 0.3):
        try:
            sg.append(SpacegroupAnalyzer(pm[k], symprec=tol).get_space_group_symbol())
        except Exception as e:  # noqa: BLE001
            sg.append(f"err {e}")
    print(f"{k:18s} atoms={len(a):3d} V/fu={a.get_volume() / nfu:6.2f} CN(Ca-O)={cn(a)}  spacegroup@0.01/0.1/0.3={sg}")
sm = StructureMatcher()
for k in ("transformed", "start"):
    for n in ("ref-186_C1Ca1O3", "ref-187_C1Ca1O3", "ref-188_C1Ca1O3", "cand-055_C1Ca1O3"):
        rms = sm.get_rms_dist(pm[k], pm[n])
        print(f"match {k:12s} vs {n:18s}: {sm.fit(pm[k], pm[n])}  rms={rms}")
