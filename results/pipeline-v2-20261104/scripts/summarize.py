"""Recompute the pipeline-v2 snapshot tables from the bundled records; standard library only.

  python results/pipeline-v2-20261104/scripts/summarize.py           # write tables/*.csv and print the headline numbers
  python results/pipeline-v2-20261104/scripts/summarize.py --check   # verify file hashes and that tables/ is unchanged

No model is loaded and no calculation is run. Aggregation rules:
- sealed-test force RMSE "all" and "interfaces" is atom-weighted over frames: sqrt(sum n_frames n_atoms RMSE^2 /
  sum n_frames n_atoms), from the per-system `dp test` RMSEs; "interfaces" are the hap001* slabs;
- sealed-test energy RMSE is frame-weighted: sqrt(sum n_frames RMSE^2 / sum n_frames);
- the mean accurate fraction of an exploration is the unweighted mean over its runs, as the orchestrator uses it.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import statistics
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parents[1]
GATES = dict(force_all=0.10, force_interfaces=0.12, energy_meV_atom=2.0)
HOST_PATH = re.compile(r"[A-Za-z]:[\\/]|/work/|\\Users\\|/Users/")


def natoms(system):
    formula = system.split("__")[-1]
    return sum(int(n or 1) for _, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula))


def iterations():
    return sorted(p for p in PKG.glob("iter-*") if p.is_dir())


def table_csv(rows):
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=list(rows[0]), lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return buf.getvalue()


def sealed_tables():
    overall, by_system = [], []
    for d in iterations():
        if not (d / "sealed_test.json").exists():
            continue
        report = json.loads((d / "sealed_test.json").read_text())
        for model, systems in report.items():
            f_all = [0.0, 0]
            f_int = [0.0, 0]
            e = [0.0, 0]
            for system, r in sorted(systems.items()):
                n, m = r["frames"], natoms(system)
                f_all[0] += n * m * r["force_rmse_eV_A"] ** 2
                f_all[1] += n * m
                if system.startswith("hap001"):
                    f_int[0] += n * m * r["force_rmse_eV_A"] ** 2
                    f_int[1] += n * m
                e[0] += n * r["energy_rmse_eV_atom"] ** 2
                e[1] += n
                by_system.append(dict(iteration=d.name, model=model, system=system, frames=n, natoms=m,
                                      force_rmse_eV_A=f"{r['force_rmse_eV_A']:.4f}",
                                      energy_rmse_meV_atom=f"{1000 * r['energy_rmse_eV_atom']:.3f}"))
            fa, fi, en = (math.sqrt(s / w) for s, w in (f_all, f_int, e))
            overall.append(dict(iteration=d.name, model=model, frames=e[1],
                                force_rmse_all_eV_A=f"{fa:.4f}", force_rmse_interfaces_eV_A=f"{fi:.4f}",
                                energy_rmse_meV_atom=f"{1000 * en:.3f}",
                                gate_passed=fa <= GATES["force_all"] and fi <= GATES["force_interfaces"]
                                and 1000 * en <= GATES["energy_meV_atom"]))
    return overall, by_system


def exploration_table():
    rows = []
    for d in iterations():
        for name in sorted(d.glob("explore*_summary.json")):
            s = json.loads(name.read_text())
            valid = name.name == "explore_summary.json"
            for run, r in s["runs"].items():
                rows.append(dict(iteration=d.name, valid=valid, run=run, frames=r.get("frames", ""),
                                 lammps_complete=r.get("lammps_complete", ""), accurate=f"{r['accurate']:.4f}",
                                 candidate=f"{r['candidate']:.4f}", failed=f"{r['failed']:.4f}",
                                 picked=r.get("picked", ""), max_devi_f_max_eV_A=f"{r['max_devi_f_max']:.4f}"))
    return rows


def label_table():
    rows = []
    for d in iterations():
        if not (d / "label_audit.csv").exists():
            continue
        with (d / "label_audit.csv").open(newline="") as f:
            audit = list(csv.DictReader(f))
        for system in sorted({r["system"] for r in audit}):
            sel = [r for r in audit if r["system"] == system]
            rows.append(dict(iteration=d.name, system=system, frames=len(sel),
                             ok=sum(r["status"] == "ok" for r in sel),
                             not_ok=";".join(f"{r['name']}:{r['status']}" for r in sel if r["status"] != "ok")))
    return rows


def finite_size():
    p = PKG / "finite-size/finite_size_pairs.csv"
    with p.open(newline="") as f:
        rows = list(csv.DictReader(f))
    de = [float(r["dE_8A_minus_12A_meV_atom"]) for r in rows]
    return dict(frames=len(rows), dE_mean=statistics.mean(de), dE_sd=statistics.stdev(de), dE_min=min(de),
                dE_max=max(de), dF_rms_mean=statistics.mean(float(r["dF_rms_eV_A"]) for r in rows),
                dF_max=max(float(r["dF_max_abs_eV_A"]) for r in rows),
                dstress_rms_mean=statistics.mean(float(r["dstress_rms_GPa"]) for r in rows))


def check_files():
    manifest = json.loads((PKG / "manifest.json").read_text())
    assert manifest["scientific_calculations_run_for_export"] is False
    assert manifest["manuscript_eligible"] is False and manifest["production_qualified"] is False
    listed = set()
    for item in manifest["files"]:
        p = PKG / item["path"]
        data = p.read_bytes()
        if len(data) != item["bytes"] or hashlib.sha256(data).hexdigest() != item["sha256"]:
            raise SystemExit(f"FAIL: file identity differs: {item['path']}")
        if p.suffix in (".json", ".csv", ".out") and HOST_PATH.search(data.decode("utf-8")):
            raise SystemExit(f"FAIL: host path left in {item['path']}")
        listed.add(item["path"])
    extra = {p.relative_to(PKG).as_posix() for p in PKG.rglob("*") if p.is_file()} - listed - {"manifest.json"}
    unexpected = sorted(p for p in extra if not (p.startswith(("tables/", "scripts/")) or p in ("README.md", "REPRODUCE.md")))
    if unexpected:
        raise SystemExit("FAIL: files not in the manifest: " + ", ".join(unexpected))
    return len(listed)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    overall, by_system = sealed_tables()
    tables = {"sealed_test_overall.csv": overall, "sealed_test_by_system.csv": by_system,
              "exploration_runs.csv": exploration_table(), "labels_by_system.csv": label_table()}
    if args.check:
        n = check_files()
        for name, rows in tables.items():
            if (PKG / "tables" / name).read_text(encoding="utf-8") != table_csv(rows):
                raise SystemExit(f"FAIL: tables/{name} differs from the recomputed table")
        print(f"PASS: {n} files match manifest.json; {len(tables)} tables recomputed identically; no calculation run")
    else:
        (PKG / "tables").mkdir(exist_ok=True)
        for name, rows in tables.items():
            (PKG / "tables" / name).write_text(table_csv(rows), encoding="utf-8", newline="")
    for r in overall:
        print(f"{r['iteration']} {r['model']}: force all {r['force_rmse_all_eV_A']}, interfaces "
              f"{r['force_rmse_interfaces_eV_A']} eV/A, energy {r['energy_rmse_meV_atom']} meV/atom, gate {r['gate_passed']}")
    for d in iterations():
        if (d / "explore_summary.json").exists():
            s = json.loads((d / "explore_summary.json").read_text())
            print(f"{d.name} exploration: mean accurate {s['mean_accurate']:.3f}, {s['n_candidates']} candidates")
    fs = finite_size()
    print("finite size (8 A minus 12 A, {frames} frames): dE {dE_mean:.2f} (SD {dE_sd:.2f}, {dE_min:.2f}-{dE_max:.2f}) "
          "meV/atom; dF rms {dF_rms_mean:.3f}, max {dF_max:.3f} eV/A; stress rms {dstress_rms_mean:.3f} GPa".format(**fs))
    progress = json.loads((PKG / "state/progress.json").read_text())
    print(f"snapshot {progress['snapshot_utc']}: " + "; ".join(
        f"{k} {v['converged_done']}/{v['frames']} converged, SCF failures {v['scf_not_converged'] or 'none'}"
        for k, v in progress["uncollected_label_steps"].items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
