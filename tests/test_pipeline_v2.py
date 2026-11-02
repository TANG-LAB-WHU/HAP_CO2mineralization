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


class EnsembleConfigTests(unittest.TestCase):
    def test_every_loss_term_is_active_in_pytorch(self):
        sys.path.insert(0, str(ROOT / "3.ML_AIMD/Step3_mlip_deepmd"))
        from dp_ensemble import config
        loss = config(["a"], ["b"], 11, 1000)["loss"]
        # DeePMD-kit 3.1.3 PyTorch keeps a term only when both prefactors are nonzero
        for term in "efv":
            self.assertTrue(loss[f"start_pref_{term}"] != 0 and loss[f"limit_pref_{term}"] != 0, term)


class ValidateTests(unittest.TestCase):
    def setUp(self):
        sys.path.insert(0, str(ROOT / "3.ML_AIMD/Step3_mlip_deepmd"))

    def test_thermo_blocks_and_drift(self):
        from dp_validate import drift, thermo_blocks
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "log.lammps"
            rows = "\n".join(f"{s} 300.0 {-100.0 + 0.002 * s:.12g}" for s in range(0, 2001, 20))
            log.write_text(f"Step Temp TotEng\n{rows}\nLoop time of 1.0 on 1 procs\n"
                           "WARNING: x\nStep Temp TotEng\n0 1 2\nLoop time of 1.0\n")
            blocks = thermo_blocks(log)
            self.assertEqual(len(blocks), 2)
            self.assertEqual(len(blocks[0]["Step"]), 101)
            t = blocks[0]["Step"] * 0.5 / 1000  # ps
            # 0.002 eV per step = 4 eV/ps; over 100 atoms that is 40 meV/atom/ps
            self.assertAlmostEqual(drift(t, blocks[0]["TotEng"], 100), 40.0, places=6)

    def test_geometry_gate(self):
        from ase import Atoms
        from dp_validate import geometry_gate
        po4 = [(0, 0, 0)] + [tuple(1.55 * np.array(v) / np.sqrt(3)) for v in ([1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1])]
        good = Atoms("PO4", positions=np.array(po4) + 5, cell=[10, 10, 10], pbc=True)
        good += Atoms("OH", positions=[(2, 2, 2), (2.97, 2, 2)])
        good += Atoms("CO2", positions=[(8, 2, 2), (9.16, 2, 2), (6.84, 2, 2)])
        self.assertEqual(geometry_gate(good), [])
        bad = good.copy()
        bad.positions[1] += [0.6, 0.6, 0.6]  # one P-O bond broken
        bad.positions[6] = (3.4, 2, 2)  # H 1.4 A from its O
        fails = geometry_gate(bad)
        self.assertTrue(any(f.startswith("P0 has 3 O") for f in fails))
        self.assertTrue(any(f.startswith("H6 nearest O") for f in fails))

    def test_repeats_and_volume(self):
        from ase.build import bulk
        from dp_validate import lattice_compare, repeats
        a = bulk("CaO", "rocksalt", a=4.8)
        rep = repeats(a, 14.0)
        widths = [a.cell.volume / np.linalg.norm(np.cross(a.cell[(k + 1) % 3], a.cell[(k + 2) % 3])) for k in range(3)]
        self.assertTrue(all(r * w >= 14.0 for r, w in zip(rep, widths)))
        ref = a.cell.cellpar()
        md = np.r_[ref[:3] * 1.01 * np.array(rep), ref[3:]]
        cmp_ = lattice_compare(md, ref, rep)
        self.assertTrue(np.allclose(cmp_["rel_length_dev"], 0.01, atol=1e-6))
        self.assertAlmostEqual(cmp_["rel_volume_dev"], 1.01 ** 3 - 1, places=4)


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
