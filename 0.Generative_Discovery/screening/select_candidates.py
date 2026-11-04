#!/usr/bin/env python3
"""
S1 chemistry gate for imported MatterGen S.U.N. structures.
Input:  sun_unique_cross_run.csv from summarize_gf_scan.py + the raw scan directory
Output: s1_candidates.extxyz (atoms.info carries the CSV metadata) and s1_gate.csv

Keeps Ca carbonates, Ca phosphates and Ca carbonate-phosphates that are formal-charge
neutral, have every C in an isolated CO3 and every P in an isolated PO4 with no shared O,
and have every H bonded to O (nearest neighbour O within 1.15 A).
"""

import argparse
import csv
from collections import Counter
from pathlib import Path

import numpy as np
from ase.io import read, write
from ase.neighborlist import neighbor_list

KEEP_CLASSES = {"Ca_carbonate_phosphate", "Ca_carbonate", "Ca_phosphate"}


def h_bonded_to_o(atoms):
    sym = np.array(atoms.get_chemical_symbols())
    h_idx = np.where(sym == "H")[0]
    if len(h_idx) == 0:
        return True
    i, j, d = neighbor_list("ijd", atoms, 2.0)
    for h in h_idx:
        mask = i == h
        if not mask.any():
            return False
        nearest = j[mask][np.argmin(d[mask])]
        if sym[nearest] != "O" or d[mask].min() > 1.15:
            return False
    return True


def main():
    parser = argparse.ArgumentParser(description="S1 chemistry gate")
    parser.add_argument("--summary-csv", type=Path, required=True)
    parser.add_argument("--scan-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(open(args.summary_csv)))
    frames_by_run = {}
    kept, gate_rows = [], []
    for r in rows:
        reasons = []
        if r["cls"] not in KEEP_CLASSES:
            reasons.append("class")
        if int(r["charge_sum_fixed_ox"]) != 0:
            reasons.append("charge")
        if r["all_C_in_CO3"] != "True" or r["all_P_in_PO4"] != "True" or r["motif_O_not_shared"] != "True":
            reasons.append("motif")
        atoms = None
        if not reasons:
            if r["run"] not in frames_by_run:
                frames = read(args.scan_dir / r["run"] / "output" / "sun" / "sun_structures.extxyz", ":")
                frames_by_run[r["run"]] = {int(a.info["structure_index"]): a for a in frames}
            atoms = frames_by_run[r["run"]][int(r["structure_index"])].copy()
            if not h_bonded_to_o(atoms):
                reasons.append("H_not_on_O")
        gate_rows.append({**r, "s1_pass": not reasons, "s1_reject": ";".join(reasons)})
        if not reasons:
            atoms.calc = None
            atoms.info = {"cand_id": f"{r['run'].split('gf')[-1]}_{r['structure_index']}",
                          "run": r["run"], "structure_index": int(r["structure_index"]),
                          "reduced": r["reduced"], "cls": r["cls"],
                          "mattersim_e_per_atom": float(r["e_per_atom"])}
            kept.append(atoms)

    write(args.out_dir / "s1_candidates.extxyz", kept, format="extxyz")
    with open(args.out_dir / "s1_gate.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(gate_rows[0].keys()))
        writer.writeheader()
        writer.writerows(gate_rows)
    print(f"S1: {len(kept)}/{len(rows)} pass;", dict(Counter(a.info["cls"] for a in kept)))
    print("rejections:", dict(Counter(g["s1_reject"] for g in gate_rows if not g["s1_pass"])))


if __name__ == "__main__":
    main()
