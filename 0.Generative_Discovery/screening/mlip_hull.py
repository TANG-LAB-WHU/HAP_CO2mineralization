#!/usr/bin/env python3
"""
S2 gate: per-model convex hulls from mlip_relax.py outputs, then a committee decision.
Input:  one <model>/energies.csv per model, each covering s1_candidates and reference_phases
Output: s2_hull.csv (E_hull per model per candidate) and s2_pass.extxyz (relaxed by --geometry-model)

A candidate passes when at least --min-models models give E_hull <= --ehull-max (eV/atom)
against references relaxed with the same model, with a converged relaxation and
a volume change within --max-volume-change.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

from ase.io import read, write
from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram
from pymatgen.core import Composition


def load(model_dir):
    return {r["key"]: r for r in csv.DictReader(open(model_dir / "energies.csv"))}


def main():
    parser = argparse.ArgumentParser(description="Committee hull gate")
    parser.add_argument("--model-dirs", nargs="+", type=Path, required=True)
    parser.add_argument("--candidate-source", default="s1_candidates")
    parser.add_argument("--reference-source", default="reference_phases")
    parser.add_argument("--ehull-max", type=float, default=0.05)
    parser.add_argument("--min-models", type=int, default=2)
    parser.add_argument("--max-volume-change", type=float, default=0.05)
    parser.add_argument("--geometry-model", default=None, help="model dir name whose relaxed geometry is exported")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    per_model = {}
    for d in args.model_dirs:
        rows = load(d)
        # Unconverged references are kept: their energy is an upper bound, and dropping a hull
        # vertex would make candidates look more stable than they are.
        ref_rows = {k: r for k, r in rows.items() if r["source"] == args.reference_source}
        refs = [PDEntry(Composition(r["formula"]), float(r["e_per_atom"]) * int(r["natoms"]), name=k)
                for k, r in ref_rows.items()]
        n_unconverged = sum(r["converged"] != "True" for r in ref_rows.values())
        pd = PhaseDiagram(refs)
        res = {}
        for k, r in rows.items():
            if r["source"] != args.candidate_source:
                continue
            entry = PDEntry(Composition(r["formula"]), float(r["e_per_atom"]) * int(r["natoms"]))
            _, ehull = pd.get_decomp_and_e_above_hull(entry, allow_negative=True)
            ok = (r["converged"] == "True" and abs(float(r["volume_ratio"]) - 1) <= args.max_volume_change)
            res[k] = dict(ehull=ehull, ok=ok, formula=r["formula"])
        per_model[d.name] = res
        print(f"{d.name}: {len(refs)} references ({n_unconverged} unconverged), {len(res)} candidates")

    keys = sorted(set().union(*[set(v) for v in per_model.values()]), key=lambda s: int(s.split(":")[1]))
    out_rows, passed = [], []
    for k in keys:
        row = {"key": k}
        votes = 0
        for m, res in per_model.items():
            if k in res:
                row[f"ehull_{m}"] = round(res[k]["ehull"], 4)
                row[f"ok_{m}"] = res[k]["ok"]
                votes += int(res[k]["ok"] and res[k]["ehull"] <= args.ehull_max)
        row["formula"] = next(res[k]["formula"] for res in per_model.values() if k in res)
        row["votes"] = votes
        row["s2_pass"] = votes >= args.min_models
        out_rows.append(row)
        if row["s2_pass"]:
            passed.append(k)

    fields = ["key", "formula"] + sorted({f for r in out_rows for f in r} - {"key", "formula", "votes", "s2_pass"}) + ["votes", "s2_pass"]
    with open(args.out_dir / "s2_hull.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fields)
        writer.writeheader()
        writer.writerows(out_rows)

    geo = next(d for d in args.model_dirs if d.name == (args.geometry_model or args.model_dirs[0].name))
    relaxed = {a.info["key"]: a for a in read(geo / "relaxed.extxyz", ":")}
    write(args.out_dir / "s2_pass.extxyz", [relaxed[k] for k in passed if k in relaxed], format="extxyz")
    by_votes = defaultdict(int)
    for r in out_rows:
        by_votes[r["votes"]] += 1
    print(f"S2: {len(passed)}/{len(out_rows)} pass (>= {args.min_models} models, E_hull <= {args.ehull_max} eV/atom); votes {dict(by_votes)}")


if __name__ == "__main__":
    main()
