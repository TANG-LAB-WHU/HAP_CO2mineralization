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
            self.assertIn("        MAX_SCF 6\n      &END OUTER_SCF",
                          (tmp / "labels/lime_selected_00000/input.inp").read_text())
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
        shared = good.copy()
        shared.positions[6] = (3.25, 2, 2)  # H 1.25 A from its O, as in a short hydrogen bond
        self.assertEqual(geometry_gate(shared), [])
        bad = good.copy()
        bad.positions[1] += [0.6, 0.6, 0.6]  # one P-O bond broken
        bad.positions[6] = (3.5, 2, 2)  # H 1.5 A from its O
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


class LabelScreenTests(unittest.TestCase):
    def test_flag_needs_agreeing_models_and_no_reproduced_relabel(self):
        from screen_labels import flag
        rows = [dict(system="w", residual_meV_atom=1.0 + 0.1 * k, spread_meV_atom=1.0) for k in range(6)]
        rows.append(dict(system="w", residual_meV_atom=-55.0, spread_meV_atom=1.0))  # wrong SCF state
        rows.append(dict(system="w", residual_meV_atom=-40.0, spread_meV_atom=12.0))  # models disagree
        rows.append(dict(system="w", residual_meV_atom=-30.0, spread_meV_atom=1.0, relabel_change_meV_atom=0.4))
        rows.append(dict(system="s", residual_meV_atom=-50.0, spread_meV_atom=1.0))  # too few frames
        flag(rows, threshold=10.0, max_spread=5.0, min_frames=5)
        self.assertEqual([r.get("flagged", False) for r in rows], [False] * 6 + [True, False, False, False])

    def test_robust_scf_and_collect_exclusion(self):
        import json
        from argparse import Namespace
        from ase import Atoms
        from cp2k_label import collect, relabel, robust_scf
        text = robust_scf(cp2k_input(Atoms("CO2", positions=[[0, 0, 0], [1.16, 0, 0], [-1.16, 0, 0]],
                                           cell=[10, 10, 10], pbc=True), "label", cutoff=800, rel_cutoff=60))
        self.assertIn("MINIMIZER CG\n        LINESEARCH 2PNT\n", text)
        self.assertIn("PRECONDITIONER FULL_ALL\n", text)
        self.assertIn("SCF_GUESS ATOMIC", text)
        with tempfile.TemporaryDirectory() as tmp:
            lab = Path(tmp) / "labels"
            from ase.io import write
            h2 = Atoms("H2", positions=[[0, 0, 0], [0.74, 0, 0]], cell=[8, 8, 8], pbc=True)
            for name in ("f0", "f1"):
                (lab / name).mkdir(parents=True)
                (lab / name / "output.out").write_text(CP2K_2026_SNIPPET)
                (lab / name / "input.inp").write_text(cp2k_input(h2, "label"))
                write(lab / name / "frame.extxyz", h2)
                (lab / name / "done").touch()
            (lab / "index.csv").write_text("name,source,frame,natoms,system,stress\n"
                                           "f0,x,0,2,hap001_t,False\nf1,x,1,2,hap001_t,False\n")
            (lab / "energy_screen.json").write_text(json.dumps({"flagged": ["f1"]}))
            first = Namespace(label_dir=lab, screen="energy_screen.json", list="relabel.txt")
            relabel(first)
            self.assertTrue((lab / "f1/output.suspect.out").exists() and not (lab / "f1/done").exists())
            self.assertTrue((lab / "energy_screen.first.json").exists())
            self.assertEqual((lab / "relabel.txt").read_text(), "f1\n")
            relabel(first)  # interrupted (neither done nor failed): listed again
            self.assertEqual((lab / "relabel.txt").read_text(), "f1\n")
            (lab / "f1/output.out").write_text(CP2K_2026_SNIPPET)
            (lab / "f1/done").touch()
            relabel(first)  # a frame is relabelled once
            self.assertEqual((lab / "relabel.txt").read_text(), "")
            # A later screen (recheck) keeps its own first copy and list; f1 stays relabelled once
            (lab / "energy_screen.recheck-5.json").write_text(json.dumps({"flagged": ["f0", "f1"]}))
            relabel(Namespace(label_dir=lab, screen="energy_screen.recheck-5.json", list="relabel.recheck-5.txt"))
            self.assertTrue((lab / "energy_screen.recheck-5.first.json").exists())
            self.assertEqual((lab / "relabel.recheck-5.txt").read_text(), "f0\n")
            self.assertEqual((lab / "relabel.txt").read_text(), "")
            for f in ("f0/output.out", "f0/done"):
                (lab / f).write_text(CP2K_2026_SNIPPET)
            collect(Namespace(label_dir=lab, out=Path(tmp) / "data", max_force=15.0, valid_fraction=0.0,
                              seed=1, exclude_screen=True, screen=["energy_screen.json"]))
            audit = (Path(tmp) / "data/label_audit.csv").read_text()
            self.assertIn("f0,x,0,2,hap001_t,False,ok", audit)
            self.assertIn("f1,x,1,2,hap001_t,False,energy_outlier", audit)
            collect(Namespace(label_dir=lab, out=Path(tmp) / "data2", max_force=15.0, valid_fraction=0.0,
                              seed=1, exclude_screen=True, screen=["energy_screen.json", "energy_screen.recheck-5.json"]))
            self.assertIn("f0,x,0,2,hap001_t,False,energy_outlier", (Path(tmp) / "data2/label_audit.csv").read_text())

    def test_recheck_acts_only_on_newly_flagged_frames(self):
        import json
        from argparse import Namespace
        from unittest import mock
        import pipeline_v2
        with tempfile.TemporaryDirectory() as tmp:
            lab, data = Path(tmp) / "labels", Path(tmp) / "data"
            lab.mkdir()
            data.mkdir()
            (lab / "energy_screen.json").write_text(json.dumps({"flagged": ["f1"]}))
            calls = []
            flagged = [["f1"]]
            with mock.patch.object(pipeline_v2, "screen", side_effect=lambda *a: flagged[0]), \
                    mock.patch.object(pipeline_v2, "sh", side_effect=lambda cmd, **kw: calls.append(cmd)), \
                    mock.patch.object(pipeline_v2, "cp2k", side_effect=lambda *a: calls.append(a)), \
                    mock.patch.object(pipeline_v2, "seed_containers", return_value=[]):
                args = Namespace(ranks=20)
                # Already left out by the screen at labelling time: nothing to do
                self.assertFalse(pipeline_v2.recheck(lab, data, ["m"] * 4, args, "recheck-5"))
                self.assertEqual(calls, [])
                flagged[0] = ["f1", "f2"]
                self.assertTrue(pipeline_v2.recheck(lab, data, ["m"] * 4, args, "recheck-5"))
            self.assertTrue((Path(tmp) / "data.pre-recheck-5").is_dir())
            collect = [c for c in calls if "collect" in map(str, c)][0]
            self.assertEqual(list(map(str, collect[-3:])), ["--screen", "energy_screen.json",
                                                            "energy_screen.recheck-5.json"])
            self.assertIn("relabel.recheck-5.txt", [str(x) for c in calls for x in c])


class HpcTests(unittest.TestCase):
    def tearDown(self):
        import pipeline_v2
        pipeline_v2.RUNTIME.update(runtime="docker", cp2k=pipeline_v2.CP2K_IMAGE, mlip=pipeline_v2.MLIP_IMAGE)
        del pipeline_v2.MOUNTS[:-1]

    def test_apptainer_commands_and_run_mount(self):
        from argparse import Namespace
        from unittest import mock
        import pipeline_v2
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp).resolve()
            pipeline_v2.RUNTIME.update(runtime="apptainer", cp2k="/sif/cp2k.sif", mlip="/sif/mlip.sif")
            pipeline_v2.MOUNTS.insert(0, (run, "/hap-run"))
            self.assertEqual(pipeline_v2.in_container(run / "iter-005/labels"), "/hap-run/iter-005/labels")
            self.assertEqual(pipeline_v2.in_container(ROOT / "workflow/x.py"), "/work/workflow/x.py")
            calls = []
            with mock.patch.object(pipeline_v2, "sh", side_effect=lambda cmd, **kw: calls.append(list(map(str, cmd)))):
                pipeline_v2.mlip(["python", "x.py"])
                pipeline_v2.cp2k(run / "labels", Namespace(ranks=64), "run_labels.sh", 2, 8)
            m, c = calls
            self.assertEqual(m[:3], ["apptainer", "exec", "--nv"])
            self.assertIn(f"{run.as_posix()}:/hap-run", m)
            self.assertEqual(m[-3:], ["/sif/mlip.sif", "python", "x.py"])
            self.assertIn("HYDRA_LAUNCHER=fork,HYDRA_RMK=user", c)
            self.assertEqual(c[c.index("/sif/cp2k.sif"):], ["/sif/cp2k.sif", pipeline_v2.CP2K_LAUNCH, "bash",
                                                            "/data/run_labels.sh", "/data", "64", "1", "2", "8"])
            self.assertEqual(pipeline_v2.seed_containers("running"), [])

    def test_phases_stop_at_the_other_phase(self):
        import json
        from unittest import mock
        import pipeline_v2
        self.assertEqual([pipeline_v2.step_phase(s) for s in ("label-5", "label-extra-metad", "recheck-5", "sealed_test",
                                                              "train-5", "test-5", "explore-5", "validate-w8")],
                         ["cpu"] * 4 + ["gpu"] * 4)
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp).resolve()
            done = ["sealed_test"] + [f"{s}-{k}" for k in range(5) for s in ("label", "train", "test", "explore")]
            (run / "pipeline_state.json").write_text(json.dumps({"done": done}))
            (run / "iter-004/explore/models").mkdir(parents=True)
            for i in range(4):
                (run / f"iter-004/explore/models/model-{i}.pb").touch()
            labelled, mlips = [], []

            def fake_mlip(a, **kw):
                mlips.append(a)
                if "dp_explore.py" in a[1]:
                    out = Path(a[-1].replace("/hap-run", str(run)))
                    out.mkdir(parents=True, exist_ok=True)
                    (out / "explore_summary.json").write_text(json.dumps({"mean_accurate": 0.8, "n_candidates": 9}))

            def main(phase):
                argv = ["pipeline_v2.py", "--run", str(run), "--cutoff", "800", "--rel-cutoff", "60", "--iters", "8",
                        "--runtime", "apptainer", "--cp2k-image", "c.sif", "--mlip-image", "m.sif", "--phase", phase]
                with mock.patch.object(sys, "argv", argv), \
                        mock.patch.object(pipeline_v2, "label", side_effect=lambda r, d, *a, **kw: labelled.append(d)), \
                        mock.patch.object(pipeline_v2, "recheck", return_value=False), \
                        mock.patch.object(pipeline_v2, "sh"), mock.patch.object(pipeline_v2, "mlip", side_effect=fake_mlip), \
                        mock.patch.object(pipeline_v2, "mlip_start"), mock.patch.object(pipeline_v2, "mlip_wait", return_value=0):
                    with self.assertRaises(SystemExit) as stop:
                        pipeline_v2.main()
                del pipeline_v2.MOUNTS[:-1]
                self.assertEqual(stop.exception.code, 0)
                return json.loads((run / "pipeline_state.json").read_text())["done"], (run / "next_phase").read_text()

            state, nxt = main("cpu")
            self.assertEqual(state[-2:], ["label-5", "recheck-5"])
            self.assertEqual(nxt, "gpu\n")
            self.assertEqual(labelled, [run / "iter-005/labels"])
            state, nxt = main("gpu")
            self.assertEqual(state[-3:], ["train-5", "test-5", "explore-5"])
            self.assertEqual(nxt, "cpu\n")
            self.assertEqual(len(labelled), 1)


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
