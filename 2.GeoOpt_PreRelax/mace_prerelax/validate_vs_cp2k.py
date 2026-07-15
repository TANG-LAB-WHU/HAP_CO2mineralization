"""
Validate MACE-MH-1 pre-relaxation against CP2K GeoOpt results.
"""
import os
import sys
import numpy as np

def calculate_rmsd(xyz_file_1, xyz_file_2):
    try:
        from ase.io import read
        atoms1 = read(xyz_file_1)
        atoms2 = read(xyz_file_2)
        
        if len(atoms1) != len(atoms2):
            return "Atom count mismatch"
            
        pos1 = atoms1.get_positions()
        pos2 = atoms2.get_positions()
        
        rmsd = np.sqrt(np.mean((pos1 - pos2)**2))
        return rmsd
    except ImportError:
        return "ASE not installed"
    except Exception as e:
        return f"Error: {e}"

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python validate_vs_cp2k.py <mace_opt.xyz> <cp2k_opt.xyz>")
        sys.exit(1)
        
    rmsd = calculate_rmsd(sys.argv[1], sys.argv[2])
    print(f"RMSD: {rmsd}")
