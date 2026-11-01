"""Unit checks for pipeline v2 helpers (no CP2K, GPU or network needed)."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow"))
sys.path.insert(0, str(ROOT / "0.Generative_Discovery" / "screening"))

from cp2k_protocol import cp2k_input, parse_output, repeats_for  # noqa: E402

CP2K_2026_SNIPPET = """\
 ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]          -1655.482851690643884

 FORCES| Atomic forces [hartree/bohr]
 FORCES|   Atom     x               y               z               |f|
 FORCES|      1 -1.87665900E-09 -8.94238071E-08  2.31777305E-03   2.31777305E-03
 FORCES|      2  4.81418568E-05 -1.38211328E-03 -4.79254860E-06   1.38295977E-03
 FORCES| Sum   0.00000000E+00  0.00000000E+00  0.00000000E+00   0.00000000E+00

 STRESS| Analytical stress tensor [bar]
 STRESS|                        x                   y                   z
 STRESS|      x        2.30135599210E+04  -9.35568486540E+03  -1.76539366405E+03
 STRESS|      y       -9.35568486540E+03   3.37196029454E+04  -2.44056312001E+03
 STRESS|      z       -1.76539366405E+03  -2.44056312001E+03   1.70286626981E+04
 STRESS| 1/3 Trace                                             2.45872751881E+04

 CP2K                                 1  1.0    0.018    0.021   95.026   95.026
  **** **** ******  **  PROGRAM ENDED AT                 2026-10-29 07:00:00.000
"""


class ParseOutputTests(unittest.TestCase):
    def test_cp2k_2026_forces_energy_stress(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "output.out"
            p.write_text(CP2K_2026_SNIPPET)
            r = parse_output(p)
        self.assertAlmostEqual(r["energy_Ha"], -1655.482851690643884)
        self.assertTrue(r["scf_converged"] and r["finished"])
        self.assertEqual(len(r["forces_Ha_bohr"]), 2)  # the "Sum" row is not an atom
        self.assertAlmostEqual(r["forces_Ha_bohr"][1][1], -1.38211328e-03)
        # copied from a real CP2K 2026.2 label (bar), converted to GPa
        np.testing.assert_allclose(r["stress_GPa"][0], [2.30135599210, -0.935568486540, -0.176539366405])
        np.testing.assert_allclose(np.trace(r["stress_GPa"]) / 3, 2.45872751881, rtol=1e-9)
        self.assertAlmostEqual(r["wall_s"], 95.026)

    def test_unconverged_scf_is_flagged(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "output.out"
            p.write_text(CP2K_2026_SNIPPET + "\n *** SCF run NOT converged ***\n")
            self.assertFalse(parse_output(p)["scf_converged"])


class CellAndInputTests(unittest.TestCase):
    def test_repeats_reach_min_width_for_hexagonal_cell(self):
        from ase.build import bulk
        a = bulk("Mg", "hcp", a=9.42, c=6.88)
        rep = repeats_for(a, 12.0)
        b = a.repeat(rep)
        cell = b.cell.array
        vol = abs(np.linalg.det(cell))
        widths = [vol / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3)]
        self.assertTrue(all(w >= 12.0 - 1e-6 for w in widths))
        self.assertEqual(rep, [2, 2, 2])

    def test_input_has_protocol_settings(self):
        from ase import Atoms
        text = cp2k_input(Atoms("CO2", positions=[[0, 0, 0], [1.16, 0, 0], [-1.16, 0, 0]], cell=[10, 10, 10], pbc=True),
                          "t", cutoff=800, rel_cutoff=60)
        for token in ("CUTOFF 800", "REL_CUTOFF 60", "DFTD3(BJ)", "DZVP-MOLOPT-SR-GTH", "GTH-PBE-q4", "EPS_SCF 1e-7"):
            self.assertIn(token, text)


class PrepareTests(unittest.TestCase):
    def test_cell_standard_form_and_done_frames_kept(self):
        from argparse import Namespace
        from ase.build import bulk
        from ase.io import read, write
        from cp2k_label import prepare
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            frames = tmp / "lime_selected.extxyz"
            write(frames, [bulk("CaO", "rocksalt", a=4.81).repeat(3)] * 2)
            done = tmp / "labels/lime_selected_00001"
            done.mkdir(parents=True)
            write(done / "frame.extxyz", read(frames, "1"))
            (done / "input.inp").write_text("old")
            (done / "done").touch()
            prepare(Namespace(frames=[frames], label_dir=tmp / "labels", cutoff=800, rel_cutoff=60, stress=False))
            cell = read(tmp / "labels/lime_selected_00000/frame.extxyz").cell.array
            self.assertTrue(np.allclose(np.triu(cell, 1), 0) and np.all(np.diag(cell) > 0))
            self.assertEqual((done / "input.inp").read_text(), "old")
            self.assertTrue(np.allclose(read(done / "frame.extxyz").cell.array, read(frames, "1").cell.array))


class WaterLayerTests(unittest.TestCase):
    def test_density_and_composition(self):
        from explore_systems import water_layer
        cell_xy = np.array([[9.42, 0, 0], [-4.71, 8.158, 0]])
        layer = water_layer(cell_xy, 0.0, 9.3, np.random.default_rng(0), n_co2=2)
        symbols = layer.get_chemical_symbols()
        self.assertEqual(symbols.count("C"), 2)
        n_water = symbols.count("H") // 2
        density = n_water / (abs(np.cross(cell_xy[0], cell_xy[1])[2]) * 9.3)
        self.assertGreater(density, 0.025)
        self.assertLess(density, 0.040)


class HullTests(unittest.TestCase):
    def test_hull_energy_without_elemental_references(self):
        from cp2k_hull import hull_energy
        # Ca2CO4 = 2 CaO + CO2 by atoms: (4/7)(-1) + (3/7)(-2)
        e, dec = hull_energy("Ca2CO4", ["CaO", "CO2", "CaCO3"], [-1.0, -2.0, -1.0])
        self.assertAlmostEqual(e, -10 / 7)
        self.assertEqual([i for i, _ in dec], [0, 1])
        # A stable CaCO3 (below the CaO-CO2 line) becomes a hull vertex
        e, dec = hull_energy("Ca2CO4", ["CaO", "CO2", "CaCO3"], [-1.0, -2.0, -2.0])
        self.assertEqual(sorted(i for i, _ in dec), [0, 2])
        self.assertAlmostEqual(e, (2 * -1.0 + 5 * -2.0) / 7)
        # Composition outside the references
        self.assertIsNone(hull_energy("Ca2CO5", ["CaO", "CO2"], [-1.0, -2.0]))


if __name__ == "__main__":
    unittest.main()
