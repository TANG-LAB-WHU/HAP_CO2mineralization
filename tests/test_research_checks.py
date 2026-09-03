from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest import mock

from scripts.check_manuscript import (
    EVIDENCE_COLUMNS,
    INVENTORY_COLUMNS,
    validate_repository,
)


class ResearchValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self._write("paper/manuscript.qmd", "# Draft\n\nClaim [@seed].\n\n## TODO-EVIDENCE R1\n")
        self._write(
            "paper/references.bib",
            "@article{seed,\n  author = {Example, A.},\n  title = {Seed},\n  year = {2020}\n}\n",
        )
        self.evidence_rows = [
            self._evidence_row(
                claim_id="C1",
                status="pending",
                claim="Seed literature claim awaiting verification.",
                source_type="primary_publication",
                citation_key="seed",
            ),
            self._evidence_row(
                claim_id="R1",
                status="pending",
                claim="Planned computational result.",
                source_type="computational_artifact",
            ),
        ]
        self.inventory_rows = [
            self._inventory_row(
                artifact_id="PLANNED",
                stage="publication",
                system="test",
                artifact_type="figure set",
                status="missing",
                provenance_status="none",
            )
        ]
        self._write_csv("literature/evidence.csv", EVIDENCE_COLUMNS, self.evidence_rows)
        self._write_csv(
            "research/artifact_inventory.csv", INVENTORY_COLUMNS, self.inventory_rows
        )
        (self.root / "paper/figures").mkdir(parents=True)
        self._write("paper/figures/README.md", "# Figures\n")

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def _write(self, relative: str, content: str) -> None:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")

    def _write_csv(
        self, relative: str, columns: list[str], rows: list[dict[str, str]]
    ) -> None:
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def _evidence_row(**overrides: str) -> dict[str, str]:
        row = {column: "" for column in EVIDENCE_COLUMNS}
        row.update({"owner": "Researcher", "notes": "Test fixture."})
        row.update(overrides)
        return row

    @staticmethod
    def _inventory_row(**overrides: str) -> dict[str, str]:
        row = {column: "" for column in INVENTORY_COLUMNS}
        row.update({"manuscript_eligible": "no", "notes": "Test fixture."})
        row.update(overrides)
        return row

    def _validate(self):
        return validate_repository(self.root, check_git=False)

    def _git(self, *arguments: str) -> str:
        result = subprocess.run(
            ["git", *arguments],
            cwd=self.root,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(0, result.returncode, result.stderr)
        return result.stdout.strip()

    def _initialize_git(self) -> str:
        self._git("init", "-q", "-b", "shawn_dev")
        self._git("config", "user.name", "Research Check Tests")
        self._git("config", "user.email", "research-checks@example.invalid")
        self._git(
            "remote",
            "add",
            "origin",
            "https://github.com/TANG-LAB-WHU/HAP_CO2mineralization.git",
        )
        self._git("add", ".")
        self._git("commit", "-q", "-m", "test baseline")
        return self._git("rev-parse", "HEAD")

    def _sha256(self, relative: str) -> str:
        return hashlib.sha256((self.root / relative).read_bytes()).hexdigest()

    def _add_supported_manifest(self, commit: str) -> str:
        manifest_path = "workflow/runs/run-1.json"
        manifest = {
            "schema_version": 1,
            "run_id": "run-1",
            "status": "completed",
            "backend": "local",
            "stage": "analysis",
            "git": {
                "repository": "TANG-LAB-WHU/HAP_CO2mineralization",
                "branch": "shawn_dev",
                "commit": commit,
                "dirty": False,
            },
            "configuration": {
                "path": "run-data/config.txt",
                "sha256": self._sha256("run-data/config.txt"),
            },
            "inputs": [
                {
                    "path": "run-data/input.dat",
                    "sha256": self._sha256("run-data/input.dat"),
                }
            ],
            "outputs": [
                {
                    "path": "run-data/output.dat",
                    "sha256": self._sha256("run-data/output.dat"),
                }
            ],
            "metrics": {"energy": {"value": -1.0, "unit": "eV"}},
            "validation": {
                "passed": True,
                "reviewer": "Researcher",
                "reviewed_at": "2026-09-03",
            },
        }
        self._write(manifest_path, json.dumps(manifest, indent=2) + "\n")
        self.evidence_rows.append(
            self._evidence_row(
                claim_id="RUN1",
                status="supported",
                claim="Test computational result.",
                source_type="computational_artifact",
                evidence_path=manifest_path,
                run_id="run-1",
                verified_by="Researcher",
                verified_at="2026-09-03",
            )
        )
        self.inventory_rows.append(
            self._inventory_row(
                artifact_id="RUN-1",
                stage="analysis",
                system="test",
                artifact_type="run manifest",
                path=manifest_path,
                status="completed",
                provenance_status="manifest_complete",
                manuscript_eligible="yes",
            )
        )
        self._rewrite_evidence()
        self._rewrite_inventory()
        return manifest_path

    def _rewrite_evidence(self) -> None:
        self._write_csv("literature/evidence.csv", EVIDENCE_COLUMNS, self.evidence_rows)

    def _rewrite_inventory(self) -> None:
        self._write_csv(
            "research/artifact_inventory.csv", INVENTORY_COLUMNS, self.inventory_rows
        )

    def test_valid_pending_baseline_passes(self) -> None:
        self.assertEqual([], self._validate().errors)

    def test_undefined_citation_fails(self) -> None:
        self._write("paper/manuscript.qmd", "# Draft\n\nClaim [@missing].\n")
        errors = self._validate().errors
        self.assertTrue(any("undefined citation keys" in error for error in errors))

    def test_duplicate_claim_id_fails(self) -> None:
        self.evidence_rows.append(dict(self.evidence_rows[0]))
        self._rewrite_evidence()
        errors = self._validate().errors
        self.assertTrue(any("duplicate claim IDs" in error for error in errors))

    def test_unsupported_evidence_status_fails(self) -> None:
        self.evidence_rows[0]["status"] = "unreviewed"
        self._rewrite_evidence()
        errors = self._validate().errors
        self.assertTrue(any("unsupported evidence status" in error for error in errors))

    def test_supported_literature_requires_locator_and_verification(self) -> None:
        self.evidence_rows[0]["status"] = "supported"
        self._rewrite_evidence()
        errors = self._validate().errors
        self.assertTrue(any("source_locator is required" in error for error in errors))
        self.assertTrue(any("verified_by is required" in error for error in errors))
        self.assertTrue(any("verified_at is required" in error for error in errors))

    def test_supported_computational_claim_requires_manifest(self) -> None:
        self.evidence_rows[1].update(
            {"status": "supported", "run_id": "run-1", "verified_by": "R", "verified_at": "2026-09-03"}
        )
        self._rewrite_evidence()
        errors = self._validate().errors
        self.assertTrue(any("evidence_path is required" in error for error in errors))

    def test_invalid_manifest_checksums_fail(self) -> None:
        manifest_path = "workflow/runs/run-1.json"
        self._write(
            manifest_path,
            '{"run_id":"run-1","git":{"commit":"short","dirty":true},'
            '"configuration":{"sha256":"bad"},"inputs":[{}],"outputs":[{}],'
            '"metrics":{"energy":{"value":1}},"validation":{"passed":false}}',
        )
        self.evidence_rows[1].update(
            {
                "status": "supported",
                "evidence_path": manifest_path,
                "run_id": "run-1",
                "verified_by": "R",
                "verified_at": "2026-09-03",
            }
        )
        self.inventory_rows.append(
            self._inventory_row(
                artifact_id="RUN-1",
                stage="analysis",
                system="test",
                artifact_type="run manifest",
                path=manifest_path,
                status="completed",
                provenance_status="manifest_complete",
                manuscript_eligible="yes",
            )
        )
        self._rewrite_evidence()
        self._rewrite_inventory()
        errors = self._validate().errors
        self.assertTrue(any("SHA-256" in error for error in errors))
        self.assertTrue(any("passed validation decision" in error for error in errors))
        self.assertTrue(any("requires value and unit" in error for error in errors))

    def test_supported_manifest_verifies_real_files_commit_and_tracking(self) -> None:
        self._write("run-data/config.txt", "configuration\n")
        self._write("run-data/input.dat", "input\n")
        commit = self._initialize_git()
        self._write("run-data/output.dat", "output\n")
        self._add_supported_manifest(commit)
        self._git("add", ".")

        self.assertEqual([], validate_repository(self.root, check_git=True).errors)

    def test_manifest_checksum_mismatch_fails(self) -> None:
        self._write("run-data/config.txt", "configuration\n")
        self._write("run-data/input.dat", "input\n")
        commit = self._initialize_git()
        self._write("run-data/output.dat", "original output\n")
        self._add_supported_manifest(commit)
        self._write("run-data/output.dat", "tampered output\n")
        self._git("add", ".")

        errors = validate_repository(self.root, check_git=True).errors
        self.assertTrue(any("checksum mismatch" in error for error in errors))

    def test_manifest_commit_must_exist_in_repository(self) -> None:
        self._write("run-data/config.txt", "configuration\n")
        self._write("run-data/input.dat", "input\n")
        self._initialize_git()
        self._write("run-data/output.dat", "output\n")
        self._add_supported_manifest("f" * 40)
        self._git("add", ".")

        errors = validate_repository(self.root, check_git=True).errors
        self.assertTrue(any("does not identify a commit" in error for error in errors))

    def test_manifest_commit_must_be_reachable_from_shawn_dev(self) -> None:
        self._write("run-data/config.txt", "configuration\n")
        self._write("run-data/input.dat", "input\n")
        self._initialize_git()
        self._git("switch", "-q", "-c", "unrelated-run")
        self._write("unrelated.txt", "side branch only\n")
        self._git("add", "unrelated.txt")
        self._git("commit", "-q", "-m", "unrelated run commit")
        unrelated_commit = self._git("rev-parse", "HEAD")
        self._git("switch", "-q", "shawn_dev")
        self._write("run-data/output.dat", "output\n")
        self._add_supported_manifest(unrelated_commit)
        self._git("add", ".")

        errors = validate_repository(self.root, check_git=True).errors
        self.assertTrue(any("not reachable from shawn_dev" in error for error in errors))

    def test_supported_manifest_must_be_tracked(self) -> None:
        self._write("run-data/config.txt", "configuration\n")
        self._write("run-data/input.dat", "input\n")
        commit = self._initialize_git()
        self._write("run-data/output.dat", "output\n")
        manifest_path = self._add_supported_manifest(commit)

        errors = validate_repository(self.root, check_git=True).errors
        self.assertTrue(
            any(manifest_path in error and "must be tracked" in error for error in errors)
        )

    def test_run_manifest_example_contains_checker_required_sections(self) -> None:
        repository_root = Path(__file__).resolve().parents[1]
        example = json.loads(
            (repository_root / "workflow/schemas/run_manifest.example.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertTrue(example["metrics"])
        self.assertIn("validation", example)
        self.assertIn("passed", example["validation"])
        self.assertIn("reviewer", example["validation"])
        self.assertIn("reviewed_at", example["validation"])

    def test_unknown_todo_evidence_id_fails(self) -> None:
        self._write("paper/manuscript.qmd", "# Draft\n\n## TODO-EVIDENCE R9\n")
        errors = self._validate().errors
        self.assertTrue(any("R9" in error and "missing" in error for error in errors))

    def test_unregistered_publication_figure_fails(self) -> None:
        (self.root / "paper/figures/result.png").write_bytes(b"not-a-real-image")
        errors = self._validate().errors
        self.assertTrue(any("publication figure lacks" in error for error in errors))

    def test_conflict_marker_must_be_quarantined(self) -> None:
        marker = "<" * 7 + " branch\nvalue\n" + "=" * 7 + "\nother\n" + ">" * 7 + " branch\n"
        self._write("legacy/input.inp", marker)
        errors = self._validate().errors
        self.assertTrue(any("not quarantined" in error for error in errors))

        self.inventory_rows.append(
            self._inventory_row(
                artifact_id="CONFLICT",
                stage="legacy",
                system="test",
                artifact_type="conflicted_input",
                path="legacy/input.inp",
                status="invalid",
                provenance_status="path_only",
            )
        )
        self._rewrite_inventory()
        self.assertEqual([], self._validate().errors)

    def test_legacy_postprocessor_has_no_synthetic_fallback(self) -> None:
        repository_root = Path(__file__).resolve().parents[1]
        source_path = (
            repository_root
            / "4.PostProcess_Analysis/HAP112_DPA_NVTcsvr_Succeeded/lammps_postprocessing.py"
        )
        source = source_path.read_text(encoding="utf-8")
        compile(source, str(source_path), "exec")
        self.assertNotIn("generate_sample_rdf_data", source)
        self.assertNotIn("np.random", source)
        self.assertIn("raise DataValidationError", source)


class _FakeBooleanArray:
    def __init__(self, value: bool):
        self.value = value

    def all(self) -> bool:
        return self.value


class _FakeArray:
    def __init__(self, rows):
        self.values = [value for row in rows for value in row]
        self.size = len(self.values)

    def __eq__(self, other):
        return _FakeBooleanArray(all(value == other for value in self.values))


class _FakeSeries:
    def __init__(self, values):
        self.values = values

    def min(self):
        return min(self.values)

    def max(self):
        return max(self.values)


class _FakeDataFrame:
    def __init__(self, rows, columns=None):
        self.rows = rows
        self.columns = columns or [f"column_{index}" for index in range(len(rows[0]))]
        self.shape = (len(rows), len(self.columns))

    def __len__(self):
        return len(self.rows)

    def __getitem__(self, column):
        index = self.columns.index(column)
        return _FakeSeries([row[index] for row in self.rows])

    def select_dtypes(self, include=None):
        return self

    def to_numpy(self, dtype=float, copy=False):
        return _FakeArray(self.rows)


def _load_legacy_postprocessor():
    repository_root = Path(__file__).resolve().parents[1]
    source_path = (
        repository_root
        / "4.PostProcess_Analysis/HAP112_DPA_NVTcsvr_Succeeded/lammps_postprocessing.py"
    )

    fake_numpy = types.ModuleType("numpy")
    fake_numpy.number = float
    fake_numpy.isfinite = lambda values: _FakeBooleanArray(
        all(math.isfinite(value) for value in values.values)
    )
    fake_numpy.all = lambda values: values.all()

    fake_pandas = types.ModuleType("pandas")
    fake_pandas.DataFrame = _FakeDataFrame
    fake_pandas.read_csv = mock.Mock(side_effect=ValueError("unparseable test input"))

    fake_pyplot = types.ModuleType("matplotlib.pyplot")
    fake_pyplot.rcParams = {}
    fake_matplotlib = types.ModuleType("matplotlib")
    fake_matplotlib.pyplot = fake_pyplot
    fake_seaborn = types.ModuleType("seaborn")
    fake_stats = types.ModuleType("scipy.stats")
    fake_signal = types.ModuleType("scipy.signal")
    fake_signal.find_peaks = lambda *args, **kwargs: ([], {})
    fake_scipy = types.ModuleType("scipy")
    fake_scipy.stats = fake_stats

    fake_modules = {
        "numpy": fake_numpy,
        "pandas": fake_pandas,
        "matplotlib": fake_matplotlib,
        "matplotlib.pyplot": fake_pyplot,
        "seaborn": fake_seaborn,
        "scipy": fake_scipy,
        "scipy.stats": fake_stats,
        "scipy.signal": fake_signal,
    }
    spec = importlib.util.spec_from_file_location("legacy_lammps_postprocessor_test", source_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load legacy postprocessor for testing")
    module = importlib.util.module_from_spec(spec)
    with mock.patch.dict(sys.modules, fake_modules):
        spec.loader.exec_module(module)
    return module


class LegacyPostprocessorBehaviorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.module = _load_legacy_postprocessor()

    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def _processor(self):
        processor = self.module.LAMMPSPostProcessor.__new__(
            self.module.LAMMPSPostProcessor
        )
        processor.input_dir = self.root
        processor.output_dir = self.root / "output"
        processor.data = {}
        processor.results = {}
        processor.anomalies = {}
        processor.data_debug = {}
        return processor

    def _write(self, name: str, content: str = "placeholder\n") -> Path:
        path = self.root / name
        path.write_text(content, encoding="utf-8")
        return path

    def test_rdf_parser_accepts_complete_lammps_blocks(self) -> None:
        rdf_path = self._write(
            "rdf.dat",
            "# TimeStep Number-of-rows\n"
            "0 2\n"
            "1 0.5 0.0 0.0\n"
            "2 1.0 1.5 0.2\n"
            "100 2\n"
            "1 0.5 0.1 0.0\n"
            "2 1.0 1.4 0.3\n",
        )
        data = self._processor().load_rdf_data_enhanced(rdf_path)
        self.assertEqual(4, len(data))

    def test_rdf_parser_rejects_empty_zero_malformed_and_truncated_data(self) -> None:
        cases = {
            "empty": "",
            "zero": "0 2\n1 0.5 0.0 0.0\n2 1.0 0.0 0.0\n",
            "malformed": "0 1\n1 0.5 not-a-number 0.0\n",
            "truncated": "0 2\n1 0.5 1.0 0.0\n",
        }
        for name, content in cases.items():
            with self.subTest(name=name):
                path = self._write(f"{name}.dat", content)
                with self.assertRaises(self.module.DataValidationError):
                    self._processor().load_rdf_data_enhanced(path)

    def test_required_main_data_parse_failure_is_fatal(self) -> None:
        self._write("comprehensive_thermodynamics.dat")
        processor = self._processor()
        with mock.patch.object(
            processor, "load_data_with_enhanced_validation", return_value=None
        ):
            with self.assertRaises(self.module.DataValidationError):
                processor.load_data()

    def test_all_zero_main_data_is_fatal(self) -> None:
        path = self._write("all-zero.dat", "0 0\n0 0\n")
        self.module.pd.read_csv = mock.Mock(
            return_value=_FakeDataFrame([[0.0, 0.0], [0.0, 0.0]], ["a", "b"])
        )
        with self.assertRaisesRegex(self.module.DataValidationError, "All numeric values are zero"):
            self._processor().load_data_with_enhanced_validation(path)

    def test_missing_optional_rdfs_do_not_reject_current_producer_output(self) -> None:
        main_files = [
            "comprehensive_thermodynamics.dat",
            "co2_chemical_transformation.dat",
            "proton_transfer_water_formation.dat",
            "energy_fluctuations.dat",
            "msd_analysis.dat",
            "co2_orientation.dat",
            "temperature_evolution.log",
            "production_temperature.log",
        ]
        required_rdfs = [
            "rdf_CO_detailed.dat",
            "rdf_C_Ca.dat",
            "rdf_C_P.dat",
            "rdf_C_H_detailed.dat",
            "rdf_H_O_detailed.dat",
        ]
        for filename in main_files + required_rdfs:
            self._write(filename)

        processor = self._processor()
        with (
            mock.patch.object(
                processor, "load_data_with_enhanced_validation", return_value=object()
            ),
            mock.patch.object(processor, "analyze_log_file"),
            mock.patch.object(processor, "save_debug_information"),
        ):
            processor.load_data()
        self.assertEqual(5, len(processor.data["rdf"]))

    def test_present_invalid_optional_rdf_is_fatal(self) -> None:
        main_files = [
            "comprehensive_thermodynamics.dat",
            "co2_chemical_transformation.dat",
            "proton_transfer_water_formation.dat",
            "energy_fluctuations.dat",
            "msd_analysis.dat",
            "co2_orientation.dat",
            "temperature_evolution.log",
            "production_temperature.log",
        ]
        required_rdfs = [
            "rdf_CO_detailed.dat",
            "rdf_C_Ca.dat",
            "rdf_C_P.dat",
            "rdf_C_H_detailed.dat",
            "rdf_H_O_detailed.dat",
        ]
        for filename in main_files + required_rdfs + ["rdf_Ca_O.dat"]:
            self._write(filename)

        processor = self._processor()

        def load(path):
            if path.name == "rdf_Ca_O.dat":
                raise self.module.DataValidationError("invalid optional RDF")
            return object()

        with (
            mock.patch.object(
                processor, "load_data_with_enhanced_validation", side_effect=load
            ),
            mock.patch.object(processor, "analyze_log_file"),
            mock.patch.object(processor, "save_debug_information"),
        ):
            with self.assertRaises(self.module.DataValidationError):
                processor.load_data()

    def test_empty_or_incomplete_log_is_fatal_and_complete_log_passes(self) -> None:
        processor = self._processor()
        log_path = self._write("log.lammps", "")
        with self.assertRaises(self.module.DataValidationError):
            processor.analyze_log_file()

        log_path.write_text("LAMMPS test build\nrun 100\n", encoding="utf-8")
        with self.assertRaises(self.module.DataValidationError):
            processor.analyze_log_file()

        log_path.write_text(
            "LAMMPS test build\nLoop time of 1.0 on 1 procs for 100 steps\n",
            encoding="utf-8",
        )
        processor.analyze_log_file()
        self.assertTrue(processor.data["log_analysis"]["simulation_progress"]["completed"])

    def test_analysis_failure_maps_to_nonzero_cli_status(self) -> None:
        processor = mock.Mock()
        processor.run_complete_enhanced_analysis.return_value = False
        existing_input = mock.Mock()
        existing_input.exists.return_value = True
        with (
            mock.patch.object(self.module, "Path", return_value=existing_input),
            mock.patch.object(
                self.module, "LAMMPSPostProcessor", return_value=processor
            ),
        ):
            self.assertEqual(1, self.module.main())


if __name__ == "__main__":
    unittest.main()
