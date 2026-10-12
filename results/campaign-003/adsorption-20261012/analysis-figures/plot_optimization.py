#!/usr/bin/env python3
"""Figure 2: redraw archived optimization evidence, without scientific execution.

Default (portable; reads only the compact CSV/JSON beside this script)::
    python plot_optimization.py

Optional, read-only re-extraction from the explicitly allowlisted raw records::
    python plot_optimization.py --extract /path/to/campaign-003

Only data/optimization-* and figures/c003-optimization.png are written.
Requires Python, NumPy and Matplotlib. No CP2K, model inference or network calls.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter

HERE = Path(__file__).resolve().parent
HA_EV = 27.211386245988
BOHR_A = 0.529177210903
FORCE_FACTOR = HA_EV / BOHR_A
SPECS = [
    ("B91", "B-initial", "first-batch-run/jobs/B91/native", 91, 71, 0),
    ("B91", "B-recovery", "first-batch-recovery1/jobs/B91/native", 91, 7, 71),
    ("B91", "B-continuation", "B-continuation1-run/jobs/B91/native", 91, 20, 78),
    ("A91", "A-initial", "first-batch-run/jobs/A91/native", 91, 1, 0),
    ("clean88", "clean-reference", "references-run/jobs/clean88/native", 88, 33, 0),
    ("gas3", "gas-reference", "references-run/jobs/gas3/native", 3, 2, 0),
]
FILES = ("cp2k.log", "HAP91-relax-pos.xyz", "HAP91-relax-frc.xyz", "result.json")
FLOAT = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"
FIELDNAMES = [
    "system", "segment", "raw_step", "unique_update", "attempted_unique_update",
    "record_status", "complete_force_evaluation", "include_unique_curve",
    "restart_duplicate", "energy_Ha", "relative_energy_eV", "force_unit",
    "maximum_force_Ha_per_bohr", "RMS_force_Ha_per_bohr",
    "maximum_force_eV_per_A", "RMS_force_eV_per_A",
    "maximum_step_bohr", "RMS_step_bohr", "maximum_step_A", "RMS_step_A",
    "max_step_pass", "rms_step_pass", "max_force_pass", "rms_force_pass",
    "geometry_source", "source_path", "source_sha256", "position_source_sha256",
    "force_xyz_source_sha256", "initial_json_source_sha256",
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def xyz_frames(path, n):
    lines = path.read_text(encoding="utf-8").splitlines()
    frames = {}
    for begin in range(0, len(lines), n + 2):
        chunk = lines[begin:begin + n + 2]
        assert len(chunk) == n + 2 and int(chunk[0]) == n, path
        match = re.fullmatch(r"\s*i\s*=\s*(\d+),\s*E\s*=\s*(" + FLOAT + r")\s*", chunk[1])
        assert match, (path, chunk[1])
        step = int(match[1])
        assert step not in frames, (path, step)
        atoms = [line.split() for line in chunk[2:]]
        frames[step] = {
            "energy": float(match[2]), "species": [a[0] for a in atoms],
            "array": np.array([[float(x) for x in a[1:4]] for a in atoms]),
        }
    return frames


def scalar(block, label):
    match = re.search(r"(?m)^\s*OPT\| " + re.escape(label) + r"\s+(" + FLOAT + r")\s*$", block)
    return float(match[1]) if match else None


def gate(block, label):
    match = re.search(r"(?m)^\s*OPT\| " + re.escape(label) + r" is converged\s+(YES|NO)\s*$", block)
    return match[1] if match else ""


def extract(root):
    """Read only the 24 static whitelist entries; rebuild compact evidence."""
    root = root.resolve()
    out = HERE / "data"
    out.mkdir(exist_ok=True)
    sources, rows, segment_info, restarts, energy_only = [], [], [], [], []
    prior_b_positions = None
    prior_b_energy = None
    initial_energies = {}
    for system, segment, relative, n, expected_updates, offset in SPECS:
        paths = {name: root / relative / name for name in FILES}
        hashes = {name: sha(path) for name, path in paths.items()}
        for name in FILES:
            sources.append({"path": f"{relative}/{name}", "sha256": hashes[name], "bytes": paths[name].stat().st_size})
        doc = json.loads(paths["result.json"].read_text(encoding="utf-8"))
        initial = doc["evaluations"][0]
        assert initial["complete"] and initial["optimizer_step"] == 0
        cell = np.array(initial["cell_A"])
        initial_positions = np.array(initial["positions_A"])
        assert initial_positions.shape == (n, 3)
        pos = xyz_frames(paths["HAP91-relax-pos.xyz"], n)
        frc = xyz_frames(paths["HAP91-relax-frc.xyz"], n)
        log = paths["cp2k.log"].read_text(encoding="utf-8")
        pieces = re.split(r"(?m)^\s*OPTIMIZATION STEP:\s*(\d+)\s*$", log)
        blocks = [(0, pieces[0])] + [(int(pieces[i]), pieces[i + 1]) for i in range(1, len(pieces), 2)]
        completed, incomplete = [], []
        max_force_xyz_difference = 0.0
        for raw_step, block in blocks:
            er = re.findall(r"ENERGY\| Total FORCE_EVAL.*?energy \[hartree\]\s+(" + FLOAT + r")", block)
            fr = re.findall(r"(?m)^\s*FORCES\|\s+(\d+)\s+(" + FLOAT + r")\s+(" + FLOAT + r")\s+(" + FLOAT + r")\s+" + FLOAT, block)
            duplicate = system == "B91" and offset > 0 and raw_step == 0
            row = dict.fromkeys(FIELDNAMES, "")
            row.update(system=system, segment=segment, raw_step=raw_step,
                       force_unit="eV/angstrom", restart_duplicate=int(duplicate),
                       source_path=f"{relative}/cp2k.log", source_sha256=hashes["cp2k.log"],
                       position_source_sha256=hashes["HAP91-relax-pos.xyz"],
                       force_xyz_source_sha256=hashes["HAP91-relax-frc.xyz"],
                       initial_json_source_sha256=hashes["result.json"])
            if not er or len(fr) != n:
                assert segment == "B-initial" and raw_step == 72 and not er and not fr
                assert "SCF WAVEFUNCTION OPTIMIZATION" in block
                row.update(attempted_unique_update=offset + raw_step,
                           record_status="incomplete_SCF_no_energy_or_force", complete_force_evaluation=0,
                           include_unique_curve=0, geometry_source="none_committed")
                incomplete.append(raw_step)
                rows.append(row)
                continue
            assert [int(f[0]) for f in fr] == list(range(1, n + 1))
            forces = np.array([[float(x) for x in f[1:]] for f in fr])
            energy = float(er[0])
            if len(er) > 1:
                assert len(er) == 2 and system in ("A91", "clean88", "gas3")
                energy_only.append({"system": system, "segment": segment, "local_final_frame": raw_step + 1,
                                    "energy_Ha": float(er[1]), "excluded_from_force_evaluations": True,
                                    "reason": "Final energy-only calculation; trajectory force frame can retain old forces."})
            if raw_step == 0:
                coordinates = initial_positions
                assert abs(energy - initial["energy_Ha"]) < 1e-10
                assert np.max(np.abs(forces - np.array(initial["forces_Ha_per_bohr"]))) < 1e-14
                if duplicate:
                    delta = (coordinates - prior_b_positions) @ np.linalg.inv(cell)
                    delta -= np.rint(delta)
                    cart = delta @ cell
                    maximum_difference = float(np.max(np.abs(cart)))
                    assert maximum_difference < 1e-8
                    assert abs(energy - prior_b_energy) < 1e-8
                    restarts.append({"segment": segment, "unique_update": offset,
                                     "maximum_periodic_component_difference_A": maximum_difference,
                                     "energy_reevaluation_difference_Ha": energy - prior_b_energy,
                                     "preserved_in_csv": True, "excluded_from_unique_curve": True})
            else:
                assert raw_step in pos and raw_step in frc
                coordinates = pos[raw_step]["array"]
                assert pos[raw_step]["species"] == initial["species"]
                assert frc[raw_step]["species"] == initial["species"]
                assert abs(pos[raw_step]["energy"] - energy) < 6e-10
                assert abs(frc[raw_step]["energy"] - energy) < 6e-10
                difference = float(np.max(np.abs(forces - frc[raw_step]["array"])))
                assert difference < 5.1e-11
                max_force_xyz_difference = max(max_force_xyz_difference, difference)
            max_force = float(np.max(np.abs(forces)))
            rms_force = float(np.sqrt(np.mean(forces ** 2)))
            max_step = scalar(block, "Maximum step size")
            rms_step = scalar(block, "RMS step size")
            if raw_step:
                assert max_step is not None and rms_step is not None
                assert abs(scalar(block, "Maximum gradient") - max_force) < 5.1e-11
                assert abs(scalar(block, "RMS gradient") - rms_force) < 5.1e-11
                assert scalar(block, "Convergence limit for maximum gradient") == 0.0001944690
                assert scalar(block, "Convergence limit for RMS gradient") == 0.0000972345
                assert scalar(block, "Convergence limit for maximum step size") == 0.003
                assert scalar(block, "Convergence limit for RMS step size") == 0.0015
            initial_energies.setdefault(system, energy)
            row.update(unique_update=offset + raw_step, complete_force_evaluation=1,
                       include_unique_curve=int(not duplicate),
                       record_status="restart_geometry_reevaluation" if duplicate else "complete_force_evaluation",
                       energy_Ha=energy, relative_energy_eV=(energy - initial_energies[system]) * HA_EV,
                       maximum_force_Ha_per_bohr=max_force, RMS_force_Ha_per_bohr=rms_force,
                       maximum_force_eV_per_A=max_force * FORCE_FACTOR, RMS_force_eV_per_A=rms_force * FORCE_FACTOR,
                       maximum_step_bohr=max_step if max_step is not None else "",
                       RMS_step_bohr=rms_step if rms_step is not None else "",
                       maximum_step_A=max_step * BOHR_A if max_step is not None else "",
                       RMS_step_A=rms_step * BOHR_A if rms_step is not None else "",
                       max_step_pass=gate(block, "Maximum step size"), rms_step_pass=gate(block, "RMS step size"),
                       max_force_pass=gate(block, "Maximum gradient"), rms_force_pass=gate(block, "RMS gradient"),
                       geometry_source="initial_result_json" if raw_step == 0 else "motion_xyz")
            rows.append(row)
            completed.append(raw_step)
            if system == "B91":
                prior_b_positions, prior_b_energy = coordinates, energy
        assert completed == list(range(expected_updates + 1)), (segment, completed)
        expected_frames = expected_updates + (system != "B91")
        assert sorted(pos) == list(range(1, expected_frames + 1))
        assert sorted(frc) == list(range(1, expected_frames + 1))
        segment_info.append({"system": system, "segment": segment, "raw_directory": relative,
                             "complete_force_evaluations": len(completed), "completed_updates": expected_updates,
                             "unique_update_start": offset, "unique_update_end": offset + expected_updates,
                             "incomplete_local_steps": incomplete,
                             "max_force_xyz_vs_log_difference_Ha_per_bohr": max_force_xyz_difference,
                             "status": "engineering_truncation_then_recovered" if incomplete else doc["status"],
                             "native_optimization_converged": doc["native_optimization_converged"] if not incomplete else False,
                             "result_json_used_for_initial_evaluation_only": True})
    csv_path = out / "optimization-history.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "schema": "campaign003-optimization-figure/v1", "raw_root_at_extraction": str(root),
        "source_access_policy": "Only the 24 listed raw files are read. Raw files are never written.",
        "source_allowlist": sources,
        "constants": {"eV_per_Ha": HA_EV, "angstrom_per_bohr": BOHR_A, "eV_per_A_per_Ha_per_bohr": FORCE_FACTOR},
        "definitions": {"maximum_force": "max(abs(F_i,alpha)) over Cartesian force components, not atom-vector norm",
                        "RMS_force": "sqrt(sum(F_i,alpha**2)/(3*N)); no force correction",
                        "energy_reference": "Each system has its own initial energy; only B91 energy is plotted.",
                        "unique_update": "Initial geometry=0; completed B segments add 71+7+20 updates. Restart initial re-evaluations do not add updates.",
                        "lines": "Recorded points joined within each engineering segment only; no line bridges restart boundaries.",
                        "missing_values": "Empty CSV numeric cells are missing, never zero-filled or interpolated."},
        "thresholds": {"maximum_force_eV_per_A": 0.010, "RMS_force_eV_per_A": 0.005,
                       "maximum_step_bohr": 0.003, "RMS_step_bohr": 0.0015,
                       "maximum_force_native_printed_Ha_per_bohr": 0.0001944690,
                       "RMS_force_native_printed_Ha_per_bohr": 0.0000972345},
        "segments": segment_info, "restart_duplicate_checks": restarts,
        "excluded_final_energy_only": energy_only,
        "B_counts": {"complete_force_evaluations_including_restarts": 101, "restart_initial_reevaluations": 2,
                     "unique_evaluated_geometries": 99, "completed_updates": 98, "unfinished_scf_attempts": 1},
        "first_B_snapshot_warning": "The first-run result.json is an obsolete 53-complete-evaluation snapshot; its totals and termination label are not used. Its complete initial evaluation alone supplies initial coordinates, cross-checked against the full raw log.",
        "scope": "Historical optimization evidence. B is budget-limited and unconverged; no relaxed A/B ranking. Optimizer convergence is not vibrational or method validation.",
        "scientifically_validated": False, "manuscript_eligible": False,
        "compact_data": {"path": "optimization-history.csv", "sha256": sha(csv_path), "rows": len(rows)},
    }
    (out / "optimization-provenance.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Extracted {len(rows)} rows from {len(sources)} allowlisted files.")


def load_compact():
    manifest = json.loads((HERE / "data/optimization-provenance.json").read_text(encoding="utf-8"))
    csv_path = HERE / "data" / manifest["compact_data"]["path"]
    assert sha(csv_path) == manifest["compact_data"]["sha256"], "Compact CSV hash mismatch"
    with csv_path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == manifest["compact_data"]["rows"]
    b = [r for r in rows if r["system"] == "B91"]
    complete = [r for r in b if r["complete_force_evaluation"] == "1"]
    unique = [r for r in b if r["include_unique_curve"] == "1"]
    assert len(complete) == 101 and len(unique) == 99
    assert [int(r["unique_update"]) for r in unique] == list(range(99))
    assert sum(r["restart_duplicate"] == "1" for r in complete) == 2
    incomplete = [r for r in b if r["complete_force_evaluation"] == "0"]
    assert len(incomplete) == 1 and not incomplete[0]["energy_Ha"] and not incomplete[0]["maximum_force_eV_per_A"]
    assert all(unique[-1][key] == "NO" for key in ("max_step_pass", "rms_step_pass", "max_force_pass", "rms_force_pass"))
    assert abs(float(unique[-1]["maximum_force_eV_per_A"]) - 0.04825949182) < 1e-10
    assert abs(float(unique[-1]["RMS_force_eV_per_A"]) - 0.00620483577) < 1e-10
    return rows, manifest


def draw(rows, manifest):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.labelsize": 10, "axes.titlesize": 10.5,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.7, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(8.4, 7.15))
    grid = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.05, 0.85],
                          left=0.11, right=0.965, bottom=0.065, top=0.925, hspace=0.39)
    ax_e, ax_f, ax_t = [fig.add_subplot(grid[i]) for i in range(3)]
    ax_f.sharex(ax_e)
    colors = ("#275D8C", "#218777", "#825BA0")
    shades = ("#F0F5FA", "#EDF7F2", "#F6F0FA")
    ranges = ((0, 71), (71, 78), (78, 98))
    labels = ("I  Initial: 71 updates", "II  Recovery: 7 updates", "III  Continuation: 20 updates")
    fig.legend([Line2D([0], [0], color=c, lw=2, marker="o", ms=3) for c in colors],
               labels, loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.54, 0.995), columnspacing=1.8, fontsize=9)
    for ax in (ax_e, ax_f):
        for (start, end), shade in zip(ranges, shades):
            ax.axvspan(start, end, facecolor=shade, zorder=0)
        for boundary in (71, 78):
            ax.axvline(boundary, color="#79838D", ls=(0, (3, 3)), lw=0.9, zorder=1)
        ax.set_xlim(-1.5, 100)
        ax.grid(axis="y", color="#DCE1E6", lw=0.5)
        ax.set_axisbelow(True)
        ax.xaxis.set_major_locator(FixedLocator([0, 20, 40, 60, 71, 78, 98]))
    for (system, segment, *_), color in zip(SPECS[:3], colors):
        samples = [r for r in rows if r["segment"] == segment and r["include_unique_curve"] == "1"]
        x = np.array([int(r["unique_update"]) for r in samples])
        ax_e.plot(x, [float(r["relative_energy_eV"]) for r in samples],
                  color=color, lw=1.3, marker="o", ms=2.6, mew=0, zorder=3)
        for field, style, metric_color, marker in [
            ("maximum_force_eV_per_A", "-", "#263B50", "o"),
            ("RMS_force_eV_per_A", "-", "#C27029", "s"),
        ]:
            ax_f.plot(x, [float(r[field]) for r in samples], color=metric_color,
                      lw=1.1, ls=style, marker=marker, ms=2.2, mew=0, zorder=3)
    ax_e.set_title("(a)  B91 electronic energy", loc="left", fontweight="bold", pad=8)
    ax_e.set_ylabel(r"$E - E_0$ (eV)")
    ax_e.set_ylim(-0.70, 0.05)
    ax_e.tick_params(labelbottom=False)
    final = [r for r in rows if r["system"] == "B91" and r["include_unique_curve"] == "1"][-1]
    ax_e.annotate(f"98 updates; $\\Delta E$ = {float(final['relative_energy_eV']):.3f} eV",
                  xy=(98, float(final["relative_energy_eV"])), xytext=(45, -0.37),
                  arrowprops={"arrowstyle": "-", "color": "#68737E", "lw": 0.8},
                  fontsize=9, color="#283643")
    ax_e.text(0.02, 0.11, "Only B91 energy is shown; $E_0$ is its initial energy.",
              transform=ax_e.transAxes, fontsize=8.3, color="#5C6670")
    for x, tag in ((35, "I"), (74.5, "II"), (88, "III")):
        ax_e.text(x, 0.98, tag, transform=ax_e.get_xaxis_transform(), va="top", ha="center", color="#697580", fontsize=8)
    ax_f.set_title("(b)  B91 force criteria", loc="left", fontweight="bold", pad=8)
    ax_f.set_ylabel(r"Force (eV $\mathrm{\AA}^{-1}$)")
    ax_f.set_xlabel("Completed geometry updates (restart initial evaluations excluded)")
    ax_f.set_yscale("log")
    ax_f.set_ylim(0.0035, 1.0)
    ax_f.yaxis.set_major_locator(FixedLocator([0.005, 0.01, 0.05, 0.1, 0.5]))
    ax_f.yaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:g}"))
    for value, color in ((0.010, "#263B50"), (0.005, "#C27029")):
        ax_f.axhline(value, lw=0.9, color=color, ls=(0, (5, 3)), zorder=2)
    force_handles = [Line2D([0], [0], color="#263B50", marker="o", ms=3, lw=1.1),
                     Line2D([0], [0], color="#C27029", marker="s", ms=3, lw=1.1),
                     Line2D([0], [0], color="#65717C", lw=0.9, ls=(0, (5, 3)))]
    ax_f.legend(force_handles, ["Max. component", "RMS component", "Limits: 0.010 / 0.005"],
                loc="upper right", ncol=3, columnspacing=1.3, handlelength=2.0,
                fontsize=8, frameon=True, facecolor="white", edgecolor="none", framealpha=0.88)
    ax_t.axis("off")
    ax_t.set_title("(c)  Final force-evaluated endpoints", loc="left", fontweight="bold", pad=7)
    table_rows = []
    for system, label in [("A91", "A91"), ("clean88", "clean88"), ("gas3", "gas3 (periodic)"), ("B91", "B91 (3 segments)")]:
        selected = [r for r in rows if r["system"] == system and r["complete_force_evaluation"] == "1"]
        last = selected[-1]
        passes = sum(last[key] == "YES" for key in ("max_step_pass", "rms_step_pass", "max_force_pass", "rms_force_pass"))
        table_rows.append([label, last["unique_update"], str(len(selected)),
                           f"{float(last['maximum_force_eV_per_A']):.5f}",
                           f"{float(last['RMS_force_eV_per_A']):.5f}",
                           f"{passes}/4", "Converged" if passes == 4 else "Budget partial"])
    table = ax_t.table(cellText=table_rows,
                       colLabels=["System", "Updates", "Force evals", "Max. F", "RMS F", "Tests passed", "Endpoint status"],
                       colWidths=[0.215, 0.09, 0.12, 0.11, 0.11, 0.135, 0.22],
                       bbox=[0, 0.23, 1, 0.70], cellLoc="center", colLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8.3)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#D6DEE5")
        cell.set_linewidth(0.5)
        cell.set_facecolor("#EAF0F5" if r == 0 else ("#F7F1F8" if r == 4 else "white"))
        if r == 0:
            cell.set_text_props(fontweight="bold", color="#263B50")
        if c == 0:
            cell.set_text_props(ha="left")
    ax_t.text(0, 0.11, r"Forces in eV $\mathrm{\AA}^{-1}$. Four tests: max./RMS force and max./RMS optimizer step.",
              transform=ax_t.transAxes, fontsize=8, color="#53606A")
    ax_t.text(0, -0.035, "B91: 101 force evaluations = 99 unique geometries + 2 restart re-evaluations; all four final tests fail.",
              transform=ax_t.transAxes, fontsize=8, color="#53606A")
    fig.text(0.11, 0.027, "Boundaries at 71 and 78 retain the engineering segments. Initial-run local step 72 ended during SCF; no point is added.",
             fontsize=7.8, color="#53606A")
    dest = HERE / "figures"
    dest.mkdir(exist_ok=True)
    for ext in ("png",):
        path = dest / f"c003-optimization.{ext}"
        fig.savefig(path, dpi=300, metadata={"Creator": "plot_optimization.py"})
        print(path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--extract", type=Path, help="Optional allowlisted raw root; extraction reads files only.")
    args = parser.parse_args()
    if args.extract is not None:
        extract(args.extract)
    rows, manifest = load_compact()
    draw(rows, manifest)


if __name__ == "__main__":
    main()
