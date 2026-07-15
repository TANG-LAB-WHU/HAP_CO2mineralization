# clean_structure.py
from ase.io import read, write
from ase.geometry import get_duplicate_atoms
from ase.cell import Cell
import os

# Define cell parameters
cell_params = {
    'a': 35.426316,
    'b': 22.686947,
    'c': 25.890946,
    'alpha': 90.0,
    'beta': 90.0,
    'gamma': 88.75317
}

"""
This script derives the facet-aware project base from the parent folder name.
Expected folder names look like: HAP_112_Perfect, HAP_211_THICKNESS, etc.
We use the second underscore-separated token as the hkl (e.g., 112) and
construct PROJECT_BASE = f"HAP_{hkl}". If detection fails, fallback to HAP_hkl.
"""

script_dir = os.path.dirname(os.path.abspath(__file__))
folder_name = os.path.basename(script_dir)

project_base = "HAP_hkl"
parts = folder_name.split("_")
if len(parts) >= 2 and parts[0].upper() == "HAP" and parts[1]:
    project_base = f"HAP_{parts[1]}"

original_candidates = [
    os.path.join(script_dir, f"{project_base}_Original.xyz"),
    os.path.join(script_dir, "HAP_hkl_Original.xyz"),
]

source_xyz = None
for candidate in original_candidates:
    if os.path.exists(candidate):
        source_xyz = candidate
        break

if source_xyz is None:
    raise FileNotFoundError(
        f"Could not find input XYZ. Tried: {original_candidates}"
    )

# Read the original structure
atoms = read(source_xyz)

# Set up the cell using Cell object
cell = Cell.fromcellpar([cell_params['a'], cell_params['b'], cell_params['c'],
                        cell_params['alpha'], cell_params['beta'], cell_params['gamma']])
atoms.set_cell(cell)
atoms.set_pbc([True, True, True])
atoms.wrap()

# Remove duplicate atoms
dups = get_duplicate_atoms(atoms, cutoff=0.05)
if len(dups) > 0:
    del atoms[dups]

# Write the cleaned structure
write(f"{project_base}_clean.xyz", atoms)

# Print summary
print(f"PROJECT_BASE: {project_base}")
print(f"Saved clean structure, removed {len(dups)} duplicates.")
print(f"Cell parameters: a={cell_params['a']:.3f}, b={cell_params['b']:.3f}, c={cell_params['c']:.3f}")
print(f"Cell angles: alpha={cell_params['alpha']:.1f}, beta={cell_params['beta']:.1f}, gamma={cell_params['gamma']:.1f}")
print(f"Cell volume: {atoms.get_volume():.3f} cubic Angstrom")