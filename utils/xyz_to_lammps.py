#!/usr/bin/env python3
"""
Convert XYZ file (ASE extended format) to LAMMPS data format.
Handles the ASE extended XYZ format with Lattice and Properties in comment line.

Usage:
    python xyz_to_lammps.py <input.xyz> <output.lmpdat>

Example:
    python xyz_to_lammps.py HAP_002_Perfect.xyz HAP_002.data
"""

import os
import sys
import re
import numpy as np


def parse_ase_xyz(xyz_file):
    """
    Parse ASE extended XYZ file and extract atoms and box information.
    """
    with open(xyz_file, 'r') as f:
        lines = f.readlines()
    
    if len(lines) < 2:
        raise ValueError(f"Invalid XYZ file: {xyz_file} - file too short")
    
    try:
        num_atoms = int(lines[0].strip())
    except ValueError:
        raise ValueError(f"Invalid number of atoms in first line: {lines[0].strip()}")
    
    comment = lines[1].strip()
    
    cell = np.eye(3) * 10.0  # Default 10 A cubic cell
    lattice_match = re.search(r'Lattice="([^"]+)"', comment)
    if lattice_match:
        lattice_str = lattice_match.group(1)
        lattice_values = [float(x) for x in lattice_str.split()]
        if len(lattice_values) == 9:
            cell = np.array(lattice_values).reshape(3, 3)
        elif len(lattice_values) == 3:
            cell = np.diag(lattice_values)
    
    pbc = (True, True, True)  
    pbc_match = re.search(r'pbc="([^"]+)"', comment)
    if pbc_match:
        pbc_str = pbc_match.group(1)
        pbc_parts = pbc_str.split()
        if len(pbc_parts) == 3:
            pbc = tuple(p.upper() in ('T', 'TRUE', '1') for p in pbc_parts)
    
    atoms = []
    for i in range(2, 2 + num_atoms):
        if i >= len(lines):
            raise ValueError(f"Not enough atoms in file: expected {num_atoms}, found {len(atoms)}")
        
        parts = lines[i].strip().split()
        if len(parts) < 4:
            continue
        
        symbol = parts[0]
        try:
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            atoms.append((symbol, x, y, z))
        except ValueError:
            print(f"Warning: Skipping invalid line {i+1}: {lines[i].strip()}", file=sys.stderr)
            continue
    
    if len(atoms) != num_atoms:
        print(f"Warning: Expected {num_atoms} atoms, found {len(atoms)}", file=sys.stderr)
    
    return atoms, cell, pbc


def write_lammps_data(output_file, atoms, cell, title="LAMMPS data file"):
    """
    Write LAMMPS data file for Apatite family minerals.
    """
    elements = sorted(set([atom[0] for atom in atoms]))
    element_map = {elem: i+1 for i, elem in enumerate(elements)}
    
    masses = {
        'C': 12.0107,
        'F': 18.9984,
        'O': 15.9994,
        'H': 1.00794,
        'P': 30.9738,
        'Ca': 40.078,
        'Sr': 87.62,
        'Mg': 24.305,
        'Cl': 35.45,
        'Si': 28.085,
    }
    
    off_diag = np.array([cell[0, 1], cell[0, 2], cell[1, 0], cell[1, 2], cell[2, 0], cell[2, 1]])
    is_orthogonal = np.allclose(off_diag, 0.0, atol=1e-6)
    
    positions = np.array([[a[1], a[2], a[3]] for a in atoms])
    
    if is_orthogonal:
        xlo, ylo, zlo = 0.0, 0.0, 0.0
        xhi, yhi, zhi = cell[0, 0], cell[1, 1], cell[2, 2]
        
        min_pos = positions.min(axis=0)
        if np.any(min_pos < 0):
            shift = np.minimum(min_pos, 0)
            positions -= shift
            xlo -= shift[0]
            ylo -= shift[1]
            zlo -= shift[2]
            xhi -= shift[0]
            yhi -= shift[1]
            zhi -= shift[2]
    else:
        # For triclinic boxes
        xlo, ylo, zlo = 0.0, 0.0, 0.0
        xhi = cell[0, 0]
        yhi = cell[1, 1]
        zhi = cell[2, 2]
        xy = cell[0, 1]
        xz = cell[0, 2]
        yz = cell[1, 2]
        
        min_pos = positions.min(axis=0)
        if np.any(min_pos < 0):
            shift = np.minimum(min_pos, 0)
            positions -= shift
    
    with open(output_file, 'w') as f:
        f.write(f"{title}\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{len(elements)} atom types\n\n")
        
        f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
        f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
        f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n")
        
        if not is_orthogonal:
            f.write(f"{xy:.6f} {xz:.6f} {yz:.6f} xy xz yz\n")
        f.write("\n")
        
        f.write("Masses\n\n")
        for elem in elements:
            mass = masses.get(elem, 1.0)
            f.write(f"{element_map[elem]} {mass:.4f}  # {elem}\n")
        
        f.write("\nAtoms\n\n")
        for i, (symbol, _, _, _) in enumerate(atoms):
            atom_type = element_map[symbol]
            x, y, z = positions[i]
            f.write(f"{i+1} {atom_type} {x:.6f} {y:.6f} {z:.6f}\n")
    
    return elements, element_map


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: xyz_to_lammps.py <input.xyz> <output.lmpdat>")
        sys.exit(1)
    
    xyz_file = sys.argv[1]
    lmpdat_file = sys.argv[2]
    
    if not os.path.exists(xyz_file):
        print(f"Error: Input file not found: {xyz_file}", file=sys.stderr)
        sys.exit(1)
    
    try:
        atoms, cell, pbc = parse_ase_xyz(xyz_file)
        elements, _ = write_lammps_data(lmpdat_file, atoms, cell, 
                                         title="LAMMPS data file for Apatite material")
        
        print(f"Successfully converted {xyz_file} to {lmpdat_file}")
        print(f"Elements found: {elements}")
        print(f"Total atoms: {len(atoms)}")
        
    except Exception as e:
        print(f"Error during conversion: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)
