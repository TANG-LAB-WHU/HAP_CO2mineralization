#!/usr/bin/env python3
"""
Extract Alex-MP reference phases for every subsystem of a chemical system.
Input:  MatterGen reference_*.gz (gzip-compressed LMDB of ComputedStructureEntry dicts)
Output: reference_phases.extxyz with entry_id, DFT energy and DFT E_hull in atoms.info,
        and reference_phases.csv with the same metadata.
"""

import argparse
import csv
import gzip
import itertools
import pickle
import shutil
import tempfile
from pathlib import Path

import lmdb
from ase.io import write
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.computed_entries import ComputedStructureEntry
from pymatgen.io.ase import AseAtomsAdaptor


def main():
    parser = argparse.ArgumentParser(description="Extract reference phases for a chemical system")
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--elements", nargs="+", default=["Ca", "P", "O", "H", "C"])
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    tmp = Path(tempfile.mkdtemp())
    lmdb_path = tmp / "reference.lmdb"
    with gzip.open(args.reference, "rb") as fin, open(lmdb_path, "wb") as fout:
        shutil.copyfileobj(fin, fout, 1 << 24)

    entries = []
    env = lmdb.open(str(lmdb_path), subdir=False, readonly=True, lock=False, readahead=False)
    with env.begin() as txn:
        name = pickle.loads(txn.get(b"name"))
        systems = set(pickle.loads(txn.get(b"chemical_systems")))
        for n in range(1, len(args.elements) + 1):
            for combo in itertools.combinations(sorted(args.elements), n):
                chemsys = "-".join(combo)
                if chemsys not in systems:
                    continue
                for formula in pickle.loads(txn.get(f"{chemsys}.reduced_formulas".encode())):
                    length = pickle.loads(txn.get(f"{chemsys}.{formula}.length".encode()))
                    for i in range(length):
                        d = pickle.loads(txn.get(f"{chemsys}.{formula}.{i}".encode()))
                        entries.append(ComputedStructureEntry.from_dict(d))
    env.close()
    shutil.rmtree(tmp)

    pd = PhaseDiagram(entries)
    frames, rows = [], []
    for k, e in enumerate(entries):
        atoms = AseAtomsAdaptor.get_atoms(e.structure)
        meta = dict(
            ref_index=k, entry_id=str(e.entry_id), formula=e.composition.reduced_formula,
            natoms=len(atoms), dft_energy_per_atom=e.energy_per_atom,
            dft_ehull=pd.get_e_above_hull(e),
        )
        # "natoms" is read back by ASE as a calculator property, so keep it out of atoms.info
        atoms.info.update({k: v for k, v in meta.items() if k != "natoms"})
        frames.append(atoms)
        rows.append(meta)
    write(args.out_dir / "reference_phases.extxyz", frames, format="extxyz")
    with open(args.out_dir / "reference_phases.csv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    on_hull = sum(r["dft_ehull"] < 1e-6 for r in rows)
    print(f"{name}: {len(rows)} entries in {'-'.join(args.elements)} subsystems, {on_hull} on the DFT hull")


if __name__ == "__main__":
    main()
