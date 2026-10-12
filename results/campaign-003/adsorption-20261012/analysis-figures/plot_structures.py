#!/usr/bin/env python3
"""Reproduce Figure 1 from compact, provenance-bearing structure data.

Default: python plot_structures.py
Refresh data (read-only raw inputs):
    python plot_structures.py --extract ROOT_OF_CAMPAIGN_003

Only coordinate analysis and plotting are performed; no electronic-structure
or geometry-optimization calculation is launched. NumPy/Matplotlib are required.
All public atom IDs are one based. Internal index0 values are zero based.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
FIGURES = HERE / "figures"
JSON_PATH = DATA / "structures-source.json"
CELL_INDEX = (0, 1, 2)
ADSORBATE = (88, 89, 90)
SUBSTRATE_COUNT = 88
X_LIMIT = (-4.9, 4.9)
Z_LIMIT = (21.3, 31.2)
Y_HALF_WIDTH = 3.0
COLORS = {"Ca": "#0072B2", "P": "#CC79A7", "O": "#D55E00",
          "H": "#009E73", "C": "#20262D"}
SIZES = {"Ca": 130, "P": 84, "O": 68, "H": 30, "C": 130}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def minimum_image(displacement, cell, pbc):
    """Exact closest lattice image, including skew cells.

    Vectors of the cell are rows; Cartesian = fractional @ cell. A rounded
    fractional trial gives an upper bound R. For every shorter candidate,
    |(f+n)_i| <= R * ||inv(cell)[:, i]||. Enumerating these finite integer
    intervals therefore includes every possible improving image. Fractional
    rounding alone is NOT used as the minimum-image answer.
    """
    displacement = np.asarray(displacement, dtype=float)
    cell = np.asarray(cell, dtype=float)
    inverse = np.linalg.inv(cell)
    fractional = displacement @ inverse
    trial = np.where(pbc, -np.rint(fractional).astype(int), 0)
    radius = np.linalg.norm(displacement + trial @ cell)
    bounds = radius * np.linalg.norm(inverse, axis=0) + 1e-10
    intervals = [range(int(np.ceil(-fractional[i] - bounds[i])),
                       int(np.floor(-fractional[i] + bounds[i])) + 1)
                 if pbc[i] else [0] for i in CELL_INDEX]
    shifts = np.asarray(list(itertools.product(*intervals)), dtype=int)
    vectors = displacement + shifts @ cell
    distances = np.linalg.norm(vectors, axis=1)
    winner = int(np.argmin(distances))
    return vectors[winner], shifts[winner], float(distances[winner])


def read_extxyz(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    count = int(lines[0])
    header = lines[1]
    lattice = re.search(r'Lattice="([^"]+)"', header).group(1)
    pbc = re.search(r'pbc="([^"]+)"', header).group(1)
    atoms = [line.split() for line in lines[2:2 + count]]
    return dict(species=[a[0] for a in atoms],
                positions_A=[[float(x) for x in a[1:4]] for a in atoms],
                cell_A=np.asarray([float(x) for x in lattice.split()]).reshape(3, 3).tolist(),
                pbc=[value == "T" for value in pbc.split()],
                extxyz_header=header)


def extract(root: Path):
    specifications = [
        ("A_initial", "A", "Initial", "first-batch-run/jobs/A91/native/result.json", None),
        ("A_converged", "A", "Converged", "first-batch-run/jobs/A91/native/endpoint.extxyz",
         "first-batch-run/jobs/A91/native/result.json"),
        ("B_initial", "B", "Initial", "first-batch-run/jobs/B91/native/result.json", None),
        ("B_partial", "B", "Partial endpoint", "B-continuation1-run/jobs/B91/native/endpoint.extxyz",
         "B-continuation1-run/jobs/B91/native/result.json"),
    ]
    frames = []
    for key, family, label, relative, result_relative in specifications:
        source = root / relative
        source_info = dict(path=str(source), relative_path=relative, sha256=sha256(source))
        if result_relative is None:
            result = json.loads(source.read_text(encoding="utf-8"))
            evaluation = result["evaluations"][0]
            assert evaluation["evaluation_index"] == 0
            geometry = {k: evaluation[k] for k in ("species", "positions_A", "cell_A", "pbc")}
            source_info["selector"] = "evaluations[0]"
            state = dict(status="initial", evaluation_index=0,
                         source_run="first-batch-run", converged=None,
                         note="Initial means first-batch evaluation 0; later records are not used.")
            sources = [source_info]
        else:
            result_path = root / result_relative
            result = json.loads(result_path.read_text(encoding="utf-8"))
            geometry = read_extxyz(source)
            retained = result["retained_endpoint"]
            assert retained["sha256"] == source_info["sha256"], "Endpoint checksum mismatch"
            source_info["selector"] = "single retained endpoint frame"
            state = dict(status=result["termination_classification"],
                         evaluation_index=retained["evaluation_index"],
                         source_run=relative.split("/")[0],
                         converged=bool(result["optimization_converged"]),
                         partial_budget=bool(result["partial_budget"]),
                         note="Retained endpoint; B continuation endpoint is unconverged and not ranked.")
            sources = [source_info, dict(path=str(result_path), relative_path=result_relative,
                                        sha256=sha256(result_path),
                                        selector="termination_classification, optimization_converged, retained_endpoint")]
        assert len(geometry["species"]) == 91
        assert geometry["species"][88:91] == ["C", "O", "O"]
        frames.append(dict(frame_id=key, family=family, label=label,
                           **geometry, state=state, sources=sources))
    payload = dict(
        schema="campaign003.structure-figure.v1",
        units=dict(positions="angstrom", cell="angstrom", distance="angstrom", angle="degree"),
        indexing=dict(index0="zero-based array index", atom_id="one-based atom ID = index0 + 1",
                      substrate_index0=[0, 87], adsorbate_index0=[88, 89, 90]),
        cell_convention="Cell vectors are rows; Cartesian coordinate = fractional coordinate @ cell_A.",
        provenance_scope="Read-only extraction of existing coordinates; no new simulation.",
        representation=dict(view="orthographic x-z projection; sight line parallel to Cartesian y",
                            horizontal_axis="x minus the instantaneous CO2 carbon x coordinate",
                            vertical_axis="original Cartesian z", x_limits_A=list(X_LIMIT),
                            z_limits_A=list(Z_LIMIT),
                            substrate_slice="absolute(y_image - y_C89) <= 3.0 angstrom",
                            periodic_images="All translated substrate atoms intersecting the local crop; hollow circles mark nonzero lattice translations.",
                            geometry_distances="Full three-dimensional exact minimum-image distances, not lengths measured on the projection.",
                            contact_selection="A: nearest substrate atom to adsorbate O90. B: nearest substrate atom to C89, and nearest substrate calcium to O90. This fixed selection is used for each family's initial and retained frame.",
                            interpretation="Distance guides do not assign chemical bonds, identify products, establish mechanisms, or rank A versus the unconverged B."),
        frames=frames)
    DATA.mkdir(exist_ok=True)
    JSON_PATH.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def pair_measure(frame, origin, target, kind):
    xyz = np.asarray(frame["positions_A"])
    vector, shift, distance = minimum_image(xyz[target] - xyz[origin], frame["cell_A"], frame["pbc"])
    return dict(frame_id=frame["frame_id"], kind=kind,
                origin_index0=origin, origin_atom_id=origin + 1,
                origin_species=frame["species"][origin], target_index0=target,
                target_atom_id=target + 1, target_species=frame["species"][target],
                target_image_na=int(shift[0]), target_image_nb=int(shift[1]), target_image_nc=int(shift[2]),
                dx_A=float(vector[0]), dy_A=float(vector[1]), dz_A=float(vector[2]),
                distance_A=distance)


def nearest(frame, origin, species=None, kind="nearest_substrate"):
    candidates = [pair_measure(frame, origin, j, kind) for j in range(SUBSTRATE_COUNT)
                  if species is None or frame["species"][j] == species]
    return min(candidates, key=lambda x: x["distance_A"])


def measure(frame):
    bonds = [pair_measure(frame, 88, j, "intramolecular_distance") for j in (89, 90)]
    vectors = [np.array([b["dx_A"], b["dy_A"], b["dz_A"]]) for b in bonds]
    cosine = vectors[0] @ vectors[1] / np.linalg.norm(vectors[0]) / np.linalg.norm(vectors[1])
    angle = float(np.degrees(np.arccos(np.clip(cosine, -1, 1))))
    nearest_all = [nearest(frame, i) for i in ADSORBATE]
    nearest_ca = [nearest(frame, i, "Ca", "nearest_substrate_Ca") for i in ADSORBATE]
    contacts = [nearest_all[1]] if frame["family"] == "A" else [nearest_all[0], nearest_ca[1]]
    return dict(bonds=bonds, angle_deg=angle, nearest_substrate=nearest_all,
                nearest_calcium=nearest_ca, displayed_contacts=contacts)


def cropped_images(frame):
    xyz, cell = np.asarray(frame["positions_A"]), np.asarray(frame["cell_A"])
    carbon = xyz[88]
    low = np.array([carbon[0] + X_LIMIT[0], carbon[1] - Y_HALF_WIDTH, Z_LIMIT[0]])
    high = np.array([carbon[0] + X_LIMIT[1], carbon[1] + Y_HALF_WIDTH, Z_LIMIT[1]])
    corners = np.array(list(itertools.product(*zip(low, high))))
    inverse = np.linalg.inv(cell)
    images = []
    for i in range(SUBSTRATE_COUNT):
        fractional_shifts = (corners - xyz[i]) @ inverse
        ranges = [range(int(np.floor(fractional_shifts[:, j].min())) - 1,
                        int(np.ceil(fractional_shifts[:, j].max())) + 2)
                  if frame["pbc"][j] else [0] for j in CELL_INDEX]
        for shift in itertools.product(*ranges):
            point = xyz[i] + np.asarray(shift) @ cell
            if np.all(point >= low - 1e-10) and np.all(point <= high + 1e-10):
                images.append(dict(frame_id=frame["frame_id"], index0=i, atom_id=i + 1,
                                   species=frame["species"][i], image_na=shift[0], image_nb=shift[1], image_nc=shift[2],
                                   x_A=float(point[0]), y_A=float(point[1]), z_A=float(point[2]),
                                   plot_x_minus_carbon_A=float(point[0] - carbon[0]),
                                   depth_y_minus_carbon_A=float(point[1] - carbon[1])))
    return images


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def label(species, atom_id):
    return rf"$\mathrm{{{species}}}_{{{atom_id}}}$"


def plot_panel(ax, frame, metrics, images, letter):
    xyz, cell = np.asarray(frame["positions_A"]), np.asarray(frame["cell_A"])
    carbon = xyz[88]
    project = lambda p: np.array([p[0] - carbon[0], p[2]])
    ax.set(xlim=X_LIMIT, ylim=Z_LIMIT, aspect="equal", xlabel=r"$x-x_{\mathrm{C89}}$ (Å)", ylabel=r"$z$ (Å)")
    ax.set_xticks([-4, -2, 0, 2, 4])
    ax.set_yticks([22, 24, 26, 28, 30])
    ax.set_facecolor("#FBFCFD")
    ax.grid(color="#E4E9EF", linewidth=.55, zorder=0)
    for spine in ax.spines.values():
        spine.set_color("#B7C2CD")
    for atom in sorted(images, key=lambda p: p["depth_y_minus_carbon_A"]):
        periodic = any(atom[k] for k in ("image_na", "image_nb", "image_nc"))
        ax.scatter(atom["plot_x_minus_carbon_A"], atom["z_A"], s=SIZES[atom["species"]],
                   facecolors="none" if periodic else COLORS[atom["species"]],
                   edgecolors=COLORS[atom["species"]], alpha=.26, linewidths=1.2,
                   zorder=2)
    # The adsorbate is unwrapped around its carbon with the exact lattice metric.
    display_xyz = {88: carbon}
    for bond in metrics["bonds"]:
        endpoint = carbon + np.array([bond["dx_A"], bond["dy_A"], bond["dz_A"]])
        display_xyz[bond["target_index0"]] = endpoint
        a, b = project(carbon), project(endpoint)
        ax.plot([a[0], b[0]], [a[1], b[1]], color="#38424D", linewidth=2.0, zorder=6)
    for contact in metrics["displayed_contacts"]:
        origin = display_xyz[contact["origin_index0"]]
        endpoint = origin + np.array([contact["dx_A"], contact["dy_A"], contact["dz_A"]])
        a, b = project(origin), project(endpoint)
        ax.plot([a[0], b[0]], [a[1], b[1]], color="#4B5969", linewidth=1.5,
                linestyle=(0, (3, 3)), zorder=4)
        image_shift = np.array([contact[k] for k in ("target_image_na", "target_image_nb", "target_image_nc")])
        periodic = bool(np.any(image_shift))
        symbol = contact["target_species"]
        ax.scatter(*b, s=SIZES[symbol] * 1.25,
                   facecolors="white" if periodic else COLORS[symbol],
                   edgecolors=COLORS[symbol], linewidths=2, zorder=7)
        atom_text = label(symbol, contact["target_atom_id"]) + (r" ($-\mathbf{a}$)" if periodic else "")
        if frame["family"] == "A":
            atom_offset, distance_offset = (12, -3), (12, 0)
        elif symbol == "Ca":
            atom_offset, distance_offset = (-9, -22), (-20, 8)
        else:
            atom_offset, distance_offset = (12, -6), (18, 2)
        ax.annotate(atom_text, b, xytext=atom_offset, textcoords="offset points",
                    fontsize=10.5, color="#26323D",
                    ha="right" if atom_offset[0] < 0 else "left", zorder=10)
        midpoint = (a + b) / 2
        ax.annotate(f"{contact['distance_A']:.3f} Å", midpoint, xytext=distance_offset,
                    textcoords="offset points", ha="right" if distance_offset[0] < 0 else "left",
                    va="center", fontsize=10.5, color="#26323D", zorder=10,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=.94))
    for i in ADSORBATE:
        point = project(display_xyz[i])
        symbol = frame["species"][i]
        ax.scatter(*point, s=185 if symbol == "C" else 165, c=COLORS[symbol],
                   edgecolors="white", linewidths=.9, zorder=8)
        if frame["family"] == "A":
            offset, ha = (-12, 0), "right"
        else:
            offset, ha = ((0, 16) if i == 88 else ((-7, 14) if i == 89 else (8, 14))), "center"
        ax.annotate(label(symbol, i + 1), point, xytext=offset, textcoords="offset points",
                    ha=ha, va="center", fontsize=11, color="#18232D", zorder=10)
    title = f"{letter}   {frame['family']} · {frame['label']}"
    ax.set_title(title, loc="left", fontsize=14, weight="bold", pad=35)
    if frame["state"]["status"] == "partial_budget":
        state_text = "Continuation 1 · retained evaluation 20 · unconverged"
    elif frame["state"]["status"] == "converged":
        state_text = "Refined endpoint · evaluation 1"
    else:
        state_text = "First-batch evaluation 0"
    ax.text(0, 1.033, state_text, transform=ax.transAxes, fontsize=9.5, color="#596778")
    ax.text(.025, .975, "CO₂ geometry (3D)", transform=ax.transAxes, va="top",
            fontsize=9.4, color="#344558", weight="bold")
    distances = [b["distance_A"] for b in metrics["bonds"]]
    ax.text(.025, .931,
            f"C89–O90  {distances[0]:.3f} Å\nC89–O91  {distances[1]:.3f} Å\nO90–C89–O91  {metrics['angle_deg']:.1f}°",
            transform=ax.transAxes, va="top", fontsize=9.5, color="#253341",
            linespacing=1.5, zorder=12,
            bbox=dict(boxstyle="round,pad=.35", facecolor="white", edgecolor="#DEE5EB", alpha=.98))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract", type=Path, help="Read raw run coordinates and refresh compact inputs before rendering")
    args = parser.parse_args()
    if args.extract:
        extract(args.extract)
    payload = json.loads(JSON_PATH.read_text(encoding="utf-8"))
    frames = payload["frames"]
    measurements = {f["frame_id"]: measure(f) for f in frames}
    display_images = {f["frame_id"]: cropped_images(f) for f in frames}
    DATA.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    atom_rows, distance_rows, angle_rows, image_rows, cell_rows = [], [], [], [], []
    for f in frames:
        for i, (symbol, xyz) in enumerate(zip(f["species"], f["positions_A"])):
            atom_rows.append(dict(frame_id=f["frame_id"], index0=i, atom_id=i + 1,
                                  species=symbol, role="adsorbate" if i in ADSORBATE else "substrate",
                                  x_A=xyz[0], y_A=xyz[1], z_A=xyz[2]))
        for i, vector in enumerate(f["cell_A"]):
            cell_rows.append(dict(frame_id=f["frame_id"], vector="abc"[i], x_A=vector[0],
                                  y_A=vector[1], z_A=vector[2], periodic=f["pbc"][i]))
        m = measurements[f["frame_id"]]
        distance_rows.extend(m["bonds"] + m["nearest_substrate"] + m["nearest_calcium"])
        angle_rows.append(dict(frame_id=f["frame_id"], outer_atom_id_1=90, center_atom_id=89,
                               outer_atom_id_2=91, angle_deg=m["angle_deg"]))
        image_rows.extend(display_images[f["frame_id"]])
    write_csv(DATA / "structures-atoms.csv", atom_rows)
    write_csv(DATA / "structures-cell.csv", cell_rows)
    write_csv(DATA / "structures-distances.csv", distance_rows)
    write_csv(DATA / "structures-angles.csv", angle_rows)
    write_csv(DATA / "structures-display-images.csv", image_rows)
    (DATA / "structures-measurements.json").write_text(json.dumps(dict(
        source_data="structures-source.json", source_data_sha256=sha256(JSON_PATH),
        units=dict(distance="angstrom", angle="degree"), measurements=measurements), indent=2) + "\n", encoding="utf-8")
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.labelsize": 11, "xtick.labelsize": 9, "ytick.labelsize": 9,
                         "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 14.0))
    fig.subplots_adjust(left=.075, right=.955, top=.89, bottom=.165, wspace=.18, hspace=.40)
    for ax, frame, letter in zip(axes.flat, frames, "abcd"):
        plot_panel(ax, frame, measurements[frame["frame_id"]], display_images[frame["frame_id"]], letter)
    fig.suptitle("CO₂ near the HAP surface", x=.075, y=.982, ha="left", fontsize=21, weight="bold", color="#162B3C")
    fig.text(.075, .953, "Initial and retained coordinates · matched local projections · measured geometry only",
             fontsize=11.5, color="#546577")
    handles = [Line2D([], [], marker="o", color="none", markerfacecolor=COLORS[s],
                      markeredgecolor=COLORS[s], markersize=7, label=s) for s in ["Ca", "P", "O", "H", "C"]]
    handles += [Line2D([], [], marker="o", color="none", markerfacecolor="none",
                       markeredgecolor="#6A7787", markersize=7, label="Periodic image"),
                Line2D([], [], color="#4B5969", linestyle=(0, (3, 3)), label="3D distance guide")]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(.075, .110), ncol=7,
               frameon=False, fontsize=10, handletextpad=.5, columnspacing=1.3, borderaxespad=0)
    caption = (
        "Cartesian x–z projection; all panels share the same scale. The crop is recentered laterally on C89; z remains absolute.\n"
        "Substrate slice: |y − y(C89)| ≤ 3.0 Å; pale atoms provide context. Hollow circles are translated periodic images.\n"
        "Numbers are full 3D minimum-image measurements in the triclinic cell. B uses the Ca46 image translated by −a.\n"
        "Atom labels are one based. Dashed guides indicate proximity only. B is a partial-budget endpoint, not a converged result."
    )
    fig.text(.075, .092, caption, va="top", fontsize=9.3, color="#475869", linespacing=1.6)
    fig.savefig(FIGURES / "c003-structures.png", dpi=220)
    plt.close(fig)
    print(json.dumps(dict(figure_png=str(FIGURES / "c003-structures.png"),
                          frames=len(frames), atoms_per_frame=91,
                          source_json_sha256=sha256(JSON_PATH)), indent=2))


if __name__ == "__main__":
    main()
