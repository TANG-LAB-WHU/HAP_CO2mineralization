#!/usr/bin/env python3
"""
Summarize a MatterGen guidance-factor scan without modifying its outputs.
Input:  <scan-dir>/<run>/ with generated_crystals.extxyz, relaxed_structures.extxyz,
        metrics.json and output/{stable,sun}/*_structures.extxyz
Output: scan_summary.json, sun_all_runs.csv, sun_unique_cross_run.csv,
        priority_Ca_CO3_PO4_shortlist.csv, MANIFEST.json

"stable" and "S.U.N." are inherited from mattergen-evaluate (MatterSim energies,
E_hull < 0.1 eV/atom vs reference_MP2020correction); nothing here is DFT-validated.
"""

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from ase.io import read
from ase.neighborlist import neighbor_list
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Composition
from pymatgen.io.ase import AseAtomsAdaptor

# Fixed formal charges; a nonzero sum flags peroxide/hydride/reduced-carbon compositions
FORMAL_OX = {"Ca": 2, "P": 5, "C": 4, "H": 1, "O": -2}


def guidance_from_name(name):
    return float(name.split("gf")[-1].replace("dot", "."))


def count_motifs(atoms):
    """Count isolated CO3 (C-O < 1.45 A) and PO4 (P-O < 1.75 A) units."""
    i, j, d = neighbor_list("ijd", atoms, 1.8)
    sym = np.array(atoms.get_chemical_symbols())
    co3 = po4 = 0
    for k, s in enumerate(sym):
        if s not in ("C", "P"):
            continue
        mask = i == k
        nb, dd = sym[j[mask]], d[mask]
        if s == "C":
            n_o = np.sum((nb == "O") & (dd < 1.45))
            co3 += int(n_o == 3 and not np.any((nb != "O") & (dd < 1.6)))
        else:
            n_o = np.sum((nb == "O") & (dd < 1.75))
            po4 += int(n_o == 4 and not np.any((nb != "O") & (dd < 1.8)))
    return co3, po4


def classify(atoms):
    elements = set(atoms.get_chemical_symbols())
    co3, po4 = count_motifs(atoms)
    if "Ca" in elements and co3 and po4:
        cls = "Ca_carbonate_phosphate"
    elif "Ca" in elements and co3:
        cls = "Ca_carbonate"
    elif "Ca" in elements and po4:
        cls = "Ca_phosphate"
    elif elements & {"P", "C"}:
        cls = "other_P_or_C"
    else:
        cls = "Ca_O_H_only"
    return cls, co3, po4, len(elements)


def rmsd_and_volume_change(gen, rel):
    disp, dvol = [], []
    for g, r in zip(gen, rel):
        if g.get_chemical_symbols() != r.get_chemical_symbols():
            continue
        frac = r.get_scaled_positions(wrap=False) - g.get_scaled_positions(wrap=False)
        frac -= np.round(frac)
        disp.append(np.sqrt(np.mean(np.sum((frac @ r.cell.array) ** 2, axis=1))))
        dvol.append(r.get_volume() / g.get_volume() - 1)
    return float(np.median(disp)), float(np.median(dvol))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Summarize a MatterGen guidance-factor scan")
    parser.add_argument("--scan-dir", type=Path, required=True, help="Directory holding one subdirectory per run")
    parser.add_argument("--out-dir", type=Path, required=True, help="Directory for tracked summary files")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    runs = sorted((p for p in args.scan_dir.iterdir() if p.is_dir()), key=lambda p: guidance_from_name(p.name))
    summary, rows, frames, manifest = {}, [], {}, {}
    for run in runs:
        gen = read(run / "generated_crystals.extxyz", ":")
        rel = read(run / "relaxed_structures.extxyz", ":")
        stable = read(run / "output" / "stable" / "stable_structures.extxyz", ":")
        sun = read(run / "output" / "sun" / "sun_structures.extxyz", ":")
        n_atoms = np.array([len(a) for a in rel])
        n_elements = Counter(len(set(a.get_chemical_symbols())) for a in rel)
        rmsd, dvol = rmsd_and_volume_change(gen, rel)

        sun_classes = Counter()
        for atoms in sun:
            cls, co3, po4, nel = classify(atoms)
            sun_classes[cls] += 1
            comp = Composition(atoms.get_chemical_formula(mode="metal"))
            idx = int(atoms.info["structure_index"])
            frames[(run.name, idx)] = atoms
            rows.append(dict(
                run=run.name, gf=guidance_from_name(run.name), structure_index=idx,
                formula=comp.formula.replace(" ", ""), reduced=comp.reduced_formula,
                natoms=len(atoms), n_elements=nel, n_CO3=co3, n_PO4=po4, cls=cls,
                e_per_atom=float(atoms.info["total_energy"]) / len(atoms),
                volume_per_atom=atoms.get_volume() / len(atoms),
                charge_sum_fixed_ox=int(round(sum(FORMAL_OX[str(e)] * n for e, n in comp.items()))),
                all_C_in_CO3=int(comp["C"]) == co3,
                all_P_in_PO4=int(comp["P"]) == po4,
                motif_O_not_shared=int(comp["O"]) >= 3 * co3 + 4 * po4,
            ))

        metrics = {k: v["value"] for k, v in json.loads((run / "metrics.json").read_text()).items()}
        summary[run.name] = dict(
            guidance_factor=guidance_from_name(run.name),
            n_generated=len(gen), n_relaxed=len(rel), n_stable=len(stable), n_sun=len(sun),
            natoms_min=int(n_atoms.min()), natoms_median=float(np.median(n_atoms)), natoms_max=int(n_atoms.max()),
            n_all_5_elements=int(n_elements.get(5, 0)),
            n_with_P=sum("P" in a.get_chemical_symbols() for a in rel),
            n_with_C=sum("C" in a.get_chemical_symbols() for a in rel),
            recomputed_gen_to_relaxed_rmsd_median_A=rmsd,
            recomputed_volume_change_median=dvol,
            sun_classes=dict(sun_classes),
            mattergen_metrics=metrics,
        )
        manifest[run.name] = {
            f.name: {"bytes": f.stat().st_size, **({"sha256": sha256(f)} if f.stat().st_size < 50_000_000 else {})}
            for f in sorted(run.iterdir()) if f.is_file() and f.name != "desktop.ini"
        }
        print(run.name, dict(sun_classes), flush=True)

    for row in rows:
        row["priority"] = (row["cls"] == "Ca_carbonate_phosphate" and row["charge_sum_fixed_ox"] == 0
                           and row["all_C_in_CO3"] and row["all_P_in_PO4"] and row["motif_O_not_shared"])

    # Cross-run duplicates: StructureMatcher within each reduced formula
    matcher = StructureMatcher()
    by_formula = defaultdict(list)
    for k, row in enumerate(rows):
        by_formula[row["reduced"]].append(k)
    group_id = 0
    for ks in by_formula.values():
        structs = [AseAtomsAdaptor.get_structure(frames[(rows[k]["run"], rows[k]["structure_index"])]) for k in ks]
        for group in matcher.group_structures(structs):
            for s in group:
                rows[ks[next(n for n, t in enumerate(structs) if t is s)]]["xrun_group"] = group_id
            group_id += 1

    groups = defaultdict(list)
    for k, row in enumerate(rows):
        groups[row["xrun_group"]].append(k)
    reps = []
    for ks in groups.values():
        rep = dict(rows[min(ks, key=lambda k: rows[k]["e_per_atom"])])
        rep["n_runs_found"] = len({rows[k]["run"] for k in ks})
        rep["n_copies"] = len(ks)
        reps.append(rep)
    reps.sort(key=lambda r: (r["cls"], r["reduced"], r["e_per_atom"]))
    priority = sorted((r for r in reps if r["priority"]), key=lambda r: r["e_per_atom"])

    summary["_cross_run"] = dict(
        n_sun_total=len(rows), n_unique_after_xrun_dedup=len(reps),
        unique_classes=dict(Counter(r["cls"] for r in reps)),
        n_priority_Ca_CO3_PO4=len(priority),
        priority_by_gf=dict(Counter(str(r["gf"]) for r in priority)),
        n_runs_found=dict(Counter(str(r["n_runs_found"]) for r in reps)),
    )
    write_csv(args.out_dir / "sun_all_runs.csv", rows)
    write_csv(args.out_dir / "sun_unique_cross_run.csv", reps)
    write_csv(args.out_dir / "priority_Ca_CO3_PO4_shortlist.csv", priority)
    (args.out_dir / "scan_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (args.out_dir / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary["_cross_run"], indent=2))


if __name__ == "__main__":
    main()
