#!/usr/bin/env python3
"""
HAP slab generator from CIF file using ASE
Creates slab models for specified Miller indices with supercell expansion
FIXED: Added CIF file parsing error handling for malformed citation sections
"""

import numpy as np
import sys
import os
from pathlib import Path
import tempfile
import re

try:
    from ase.io import read, write
    from ase.build import surface, make_supercell
    from ase import Atoms
    from ase.geometry import cell_to_cellpar, cellpar_to_cell, wrap_positions
    HAS_ASE = True
except ImportError:
    HAS_ASE = False
    print("ASE not available, please install: pip install ase")

try:
    from pymatgen.io.ase import AseAtomsAdaptor
    from pymatgen.core.surface import SlabGenerator
    from pymatgen.io.cif import CifParser
    HAS_PYMATGEN = True
except ImportError:
    HAS_PYMATGEN = False
    print("pymatgen not available, using ASE only")

def clean_cif_file(cif_file):
    """
    Clean CIF file to fix common parsing issues, especially multiline citation entries
    Returns path to cleaned temporary CIF file
    """
    print(f"Cleaning CIF file: {cif_file}")
    
    with open(cif_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Fix the specific issue with multiline journal citation
    # Pattern: primary followed by multiline text in semicolons
    multiline_pattern = r'primary\s*\n\s*;\s*\n([^;]+)\n\s*;\s*(\d+\s+\d+\s+\d+\s+\d+\s+\w+)'
    
    def replace_multiline_citation(match):
        journal_text = match.group(1).strip()
        # Clean up the journal text - remove newlines and extra spaces
        journal_text = re.sub(r'\s+', ' ', journal_text)
        # Remove commas and colons that might cause issues
        journal_text = journal_text.replace(',', '').replace(':', '')
        citation_data = match.group(2).strip()
        return f'primary "{journal_text}" {citation_data}'
    
    # Apply the fix
    cleaned_content = re.sub(multiline_pattern, replace_multiline_citation, content, flags=re.MULTILINE | re.DOTALL)
    
    # Additional cleanup: fix any remaining problematic multiline entries
    lines = cleaned_content.split('\n')
    cleaned_lines = []
    i = 0
    in_loop = False
    loop_headers = []
    
    while i < len(lines):
        line = lines[i].strip()
        
        # Track if we're in a loop section
        if line.startswith('loop_'):
            in_loop = True
            loop_headers = []
            cleaned_lines.append(lines[i])
            i += 1
            continue
        elif line.startswith('_') and in_loop:
            loop_headers.append(line)
            cleaned_lines.append(lines[i])
            i += 1
            continue
        elif in_loop and (line.startswith('#') or line == '' or line.startswith('_')):
            # End of loop data
            in_loop = False
            loop_headers = []
            cleaned_lines.append(lines[i])
            i += 1
            continue
        
        # Handle problematic citation entries in loops
        if in_loop and line.startswith('primary') and 'citation' in ' '.join(loop_headers).lower():
            # Look ahead for multiline text
            if i + 1 < len(lines) and lines[i + 1].strip().startswith(';'):
                # Found multiline entry, collect all lines until closing semicolon
                multiline_parts = [line]  # Start with 'primary'
                j = i + 1
                collecting_text = False
                text_parts = []
                
                while j < len(lines):
                    current_line = lines[j].strip()
                    
                    if current_line.startswith(';') and not collecting_text:
                        # Start of multiline text
                        collecting_text = True
                        if len(current_line) > 1:
                            text_parts.append(current_line[1:].strip())
                    elif current_line == ';' and collecting_text:
                        # End of multiline text
                        collecting_text = False
                        j += 1
                        break
                    elif collecting_text:
                        text_parts.append(current_line)
                    elif not collecting_text and current_line:
                        # This should be the rest of the citation data
                        multiline_parts.append(current_line)
                        j += 1
                        break
                    j += 1
                
                # Reconstruct the citation line
                if text_parts:
                    journal_text = ' '.join(text_parts).replace(',', '').replace(':', '').strip()
                    if len(multiline_parts) > 1:
                        cleaned_lines.append(f'primary "{journal_text}" {multiline_parts[1]}')
                    else:
                        cleaned_lines.append(f'primary "{journal_text}"')
                else:
                    cleaned_lines.append(lines[i])
                
                i = j
                continue
        
        cleaned_lines.append(lines[i])
        i += 1
    
    # Write cleaned CIF to temporary file
    with tempfile.NamedTemporaryFile(mode='w', suffix='_cleaned.cif', delete=False, encoding='utf-8') as tmp_f:
        tmp_f.write('\n'.join(cleaned_lines))
        temp_path = tmp_f.name
        
    print(f"Created cleaned CIF file: {temp_path}")
    return temp_path

def write_lammps_data(atoms, filename, title="HAP slab"):
    """Write LAMMPS data file preserving the original triclinic cell exactly
    (Lx, Ly, Lz and xy, xz, yz) without normalizing or shearing coordinates."""

    # Element type mapping
    element_types = {"Ca": 1, "P": 2, "O": 3, "H": 4, "C": 5}
    masses = {"Ca": 40.078, "P": 30.974, "O": 15.999, "H": 1.008, "C": 12.011}

    # Ensure PBC; keep the original cell as-is (do not reparameterize)
    atoms.set_pbc(True)

    # Wrap positions inside the cell
    try:
        wrapped_pos = wrap_positions(atoms.get_positions(), cell=atoms.cell, pbc=[True, True, True], eps=1e-12)
    except Exception:
        wrapped_pos = atoms.get_positions().copy()

    # Unique element types present
    elements = list(set(atoms.get_chemical_symbols()))
    used_types = {el: element_types[el] for el in elements if el in element_types}

    # Extract triclinic components from ASE cell
    cell = atoms.get_cell()
    ax, ay, az = cell[0]
    bx, by, bz = cell[1]
    cx, cy, cz = cell[2]

    # Triclinic parameters expected by LAMMPS (taken directly from the cell)
    lx, ly, lz = ax, by, cz
    xy, xz, yz = bx, cx, cy

    # Use wrapped positions as-is. Do not normalize tilt or shear coordinates.
    pos = wrapped_pos.copy()

    # Final wrap to keep positions in box
    try:
        pos = wrap_positions(pos, cell=atoms.cell, pbc=[True, True, True], eps=1e-12)
    except Exception:
        pass

    # Compute orthogonal bounding box for the triclinic cell (LAMMPS formula)
    xlo = ylo = zlo = 0.0
    xhi, yhi, zhi = lx, ly, lz
    xlo_bound = xlo + min(0.0, xy, xz, xy + xz)
    xhi_bound = xhi + max(0.0, xy, xz, xy + xz)
    ylo_bound = ylo + min(0.0, yz)
    yhi_bound = yhi + max(0.0, yz)
    zlo_bound = zlo
    zhi_bound = zhi

    with open(filename, 'w') as f:
        f.write(f"# {title}\n")
        f.write(f"# Generated from HAP CIF file\n")

        # Composition info for debugging
        composition = {}
        for symbol in atoms.get_chemical_symbols():
            composition[symbol] = composition.get(symbol, 0) + 1
        f.write(f"# Composition: {composition}\n\n")

        f.write(f"{len(atoms)} atoms\n")
        f.write("0 bonds\n")
        f.write("0 angles\n")
        f.write("0 dihedrals\n")
        f.write("0 impropers\n\n")

        f.write(f"{len(used_types)} atom types\n\n")

        # Triclinic box bounds and tilt factors (always write tilt line)
        f.write(f"{xlo_bound:.6f} {xhi_bound:.6f} xlo xhi\n")
        f.write(f"{ylo_bound:.6f} {yhi_bound:.6f} ylo yhi\n")
        f.write(f"{zlo_bound:.6f} {zhi_bound:.6f} zlo zhi\n")
        f.write(f"{xy:.6f} {xz:.6f} {yz:.6f} xy xz yz\n")
        f.write("\n")

        # Masses
        f.write("Masses\n\n")
        for element, type_id in used_types.items():
            f.write(f"{type_id} {masses[element]:.3f}  # {element}\n")
        f.write("\n")

        # Atoms section (type and Cartesian coordinates)
        f.write("Atoms\n\n")
        symbols = atoms.get_chemical_symbols()
        for i, (symbol, p) in enumerate(zip(symbols, pos), 1):
            type_id = used_types[symbol]
            f.write(f"{i} {type_id} {p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")

    # Debug info
    print(f"DEBUG: Wrote {filename} with {len(atoms)} atoms and {len(used_types)} atom types")
    print(f"DEBUG: Atom types: {list(used_types.keys())}")
    print(f"DEBUG: Composition: {composition}")

def validate_structure(atoms, min_dist=0.5, max_dist=20.0):
    """Validate structure stability"""
    print("=== STRUCTURE VALIDATION ===")
    
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    
    # Check if coordinates are reasonable
    if np.any(np.isnan(positions)) or np.any(np.isinf(positions)):
        print("ERROR: Invalid coordinates detected")
        return False
    
    # Calculate interatomic distances
    try:
        from scipy.spatial.distance import pdist, squareform
        distances = pdist(positions)
        
        min_distance = distances.min()
        max_distance = distances.max()
        
        print(f"Minimum interatomic distance: {min_distance:.3f} Å")
        print(f"Maximum interatomic distance: {max_distance:.3f} Å")
        
        if min_distance < min_dist:
            print(f"WARNING: Atoms too close! Minimum distance {min_distance:.3f} Å < {min_dist} Å")
            return False
        
        if max_distance > max_dist:
            print(f"WARNING: Atoms too far! Maximum distance {max_distance:.3f} Å > {max_dist} Å")
    except ImportError:
        print("WARNING: scipy not available, skipping distance validation")
    
    # Check CO2 molecule position
    co2_indices = [i for i, symbol in enumerate(symbols) if symbol == 'C']
    if co2_indices:
        co2_pos = positions[co2_indices[0]]
        cell = atoms.get_cell()
        center_x = cell[0, 0] / 2
        center_y = cell[1, 1] / 2
        
        distance_from_center = np.sqrt((co2_pos[0] - center_x)**2 + (co2_pos[1] - center_y)**2)
        print(f"CO2 distance from surface center: {distance_from_center:.3f} Å")
        
        if distance_from_center > 2.0:
            print("WARNING: CO2 molecule is far from surface center")
    
    print("✓ Structure validation passed")
    return True

def create_slab_ase(cif_file, miller_index, layers=4, vacuum=15.0, supercell=(2, 2, 1)):
    """Create slab using ASE with error handling for CIF parsing"""
    
    if not HAS_ASE:
        raise ImportError("ASE is required for slab generation")
    
    # Try to read structure from CIF
    try:
        print(f"Attempting to read CIF file with ASE: {cif_file}")
        atoms = read(cif_file)
        print("Successfully read CIF file with ASE")
    except Exception as e:
        print(f"Error reading CIF with ASE: {e}")
        print("Attempting to clean CIF file and retry...")
        
        # Try to clean the CIF file by fixing common formatting issues
        try:
            cleaned_cif = clean_cif_file(cif_file)
            atoms = read(cleaned_cif)
            print("Successfully read cleaned CIF file with ASE")
            # Clean up temporary file
            os.unlink(cleaned_cif)
        except Exception as e2:
            print(f"Failed to read even cleaned CIF file: {e2}")
            raise RuntimeError(f"Could not parse CIF file even after cleaning. Original error: {e}, Cleaned error: {e2}")
    
    # Create surface slab
    print(f"Creating surface slab with Miller index {miller_index}")
    slab = surface(atoms, miller_index, layers, vacuum=vacuum, periodic=True)
    
    # Create supercell for larger system
    if supercell != (1, 1, 1):
        print(f"Creating supercell: {supercell}")
        supercell_matrix = np.diag(supercell)
        slab = make_supercell(slab, supercell_matrix)
    
    # Center slab in z-direction
    slab.center(vacuum=vacuum, axis=2)
    
    return slab

def create_slab_pymatgen(cif_file, miller_index, min_slab_size=15.0, min_vacuum_size=15.0, supercell=(2, 2, 1)):
    """Create slab using pymatgen with error handling for CIF parsing"""
    
    if not HAS_PYMATGEN:
        raise ImportError("pymatgen is required for advanced slab generation")
    
    # Try to parse structure
    try:
        print(f"Attempting to read CIF file with pymatgen: {cif_file}")
        parser = CifParser(cif_file)
        structure = parser.get_structures()[0]
        print("Successfully read CIF file with pymatgen")
    except Exception as e:
        print(f"Error reading CIF with pymatgen: {e}")
        print("Attempting to clean CIF file and retry...")
        
        # Try to clean the CIF file by fixing common formatting issues
        try:
            cleaned_cif = clean_cif_file(cif_file)
            parser = CifParser(cleaned_cif)
            structure = parser.get_structures()[0]
            print("Successfully read cleaned CIF file with pymatgen")
            # Clean up temporary file
            os.unlink(cleaned_cif)
        except Exception as e2:
            print(f"Failed to read even cleaned CIF file: {e2}")
            raise RuntimeError(f"Could not parse CIF file even after cleaning. Original error: {e}, Cleaned error: {e2}")
    
    # Generate slab
    print(f"Generating slab with pymatgen SlabGenerator")
    slabgen = SlabGenerator(
        initial_structure=structure,
        miller_index=miller_index,
        min_slab_size=min_slab_size,
        min_vacuum_size=min_vacuum_size,
        center_slab=True
    )
    
    slabs = slabgen.get_slabs()
    if not slabs:
        raise ValueError(f"No valid slabs generated for {miller_index}")
    
    slab = slabs[0]  # Take first termination
    
    # Convert to ASE for supercell creation
    adaptor = AseAtomsAdaptor()
    ase_slab = adaptor.get_atoms(slab)
    
    # Create supercell
    if supercell != (1, 1, 1):
        print(f"Creating supercell: {supercell}")
        supercell_matrix = np.diag(supercell)
        ase_slab = make_supercell(ase_slab, supercell_matrix)
    
    return ase_slab

def add_co2_molecule(slab, height=3.0, position=None):
    """Add CO2 molecule above the surface with improved positioning"""
    
    positions = slab.get_positions()
    cell = slab.get_cell()
    
    # Find surface center and top with improved positioning
    if position is None:
        # Use surface geometric center instead of simple cell center
        surface_threshold = positions[:, 2].max() - 5.0  # Select atoms near surface
        surface_atoms = positions[:, 2] > surface_threshold
        
        if np.any(surface_atoms):
            surface_positions = positions[surface_atoms]
            center_x = np.mean(surface_positions[:, 0])
            center_y = np.mean(surface_positions[:, 1])
        else:
            # If surface atoms cannot be identified, use cell center
            center_x = cell[0, 0] / 2
            center_y = cell[1, 1] / 2
    else:
        center_x, center_y = position
    
    top_z = positions[:, 2].max()
    co2_z = top_z + height
    
    # CO2 geometry (linear molecule, C=O bond length ~1.16 Å)
    co2_positions = [
        [center_x, center_y, co2_z],                    # C
        [center_x - 1.16, center_y, co2_z],            # O
        [center_x + 1.16, center_y, co2_z]             # O
    ]
    
    # Add CO2 atoms
    symbols = list(slab.get_chemical_symbols()) + ['C', 'O', 'O']
    all_positions = np.vstack([positions, co2_positions])
    
    new_slab = Atoms(symbols=symbols, positions=all_positions, cell=slab.get_cell(), pbc=slab.get_pbc())
    
    # Verify CO2 was added correctly
    original_symbols = set(slab.get_chemical_symbols())
    new_symbols = set(new_slab.get_chemical_symbols())
    added_symbols = new_symbols - original_symbols
    
    print(f"Original atoms: {len(slab)}")
    print(f"New atoms: {len(new_slab)}")
    print(f"Added symbols: {added_symbols}")
    print(f"CO2 position: center_x={center_x:.3f}, center_y={center_y:.3f}, height={height:.3f}")
    
    if 'C' not in added_symbols:
        raise RuntimeError("CO2 molecule was not added correctly - C atom missing")
    
    return new_slab

def add_multiple_co2_molecules(slab, n_co2=3, height=3.0):
    """Add multiple CO2 molecules above the surface with uniform distribution"""
    
    positions = slab.get_positions()
    cell = slab.get_cell()
    top_z = positions[:, 2].max()
    
    # Calculate grid positions
    x_range = cell[0, 0]
    y_range = cell[1, 1]
    
    co2_positions = []
    co2_symbols = []
    
    for i in range(n_co2):
        for j in range(n_co2):
            x_pos = (i + 0.5) * x_range / n_co2
            y_pos = (j + 0.5) * y_range / n_co2
            z_pos = top_z + height
            
            # CO2 molecule positions
            co2_positions.extend([
                [x_pos, y_pos, z_pos],                    # C
                [x_pos - 1.16, y_pos, z_pos],            # O
                [x_pos + 1.16, y_pos, z_pos]             # O
            ])
            co2_symbols.extend(['C', 'O', 'O'])
    
    # Merge structures
    all_symbols = list(slab.get_chemical_symbols()) + co2_symbols
    all_positions = np.vstack([positions, co2_positions])
    
    new_slab = Atoms(symbols=all_symbols, positions=all_positions, 
                    cell=slab.get_cell(), pbc=slab.get_pbc())
    
    print(f"Added {n_co2**2} CO2 molecules ({len(co2_symbols)} atoms)")
    print(f"CO2 distribution: {n_co2}x{n_co2} grid")
    return new_slab

def main():
    """Main function"""
    
    print("HAP Slab Generator")
    print("=" * 40)
    
    # Parse command line arguments
    if len(sys.argv) < 2:
        print("Usage: python create_hap_slabs.py <miller_h> <miller_k> <miller_l> [options]")
        print("\nExamples:")
        print("  python create_hap_slabs.py 0 0 1          # (001) surface")
        print("  python create_hap_slabs.py 0 1 0          # (010) surface") 
        print("  python create_hap_slabs.py 1 0 0          # (100) surface")
        print("  python create_hap_slabs.py 1 1 0          # (110) surface")
        print("\nOptions:")
        print("  --layers N       Number of atomic layers (default: 4)")
        print("  --vacuum V       Vacuum thickness in Å (default: 15.0)")
        print("  --supercell A B  Supercell expansion (default: 2 2)")
        print("  --add-co2        Add CO2 molecule above surface")
        print("  --n-co2 N        Number of CO2 molecules per side (default: 1, creates NxN grid)")
        print("  --co2-height H   CO2 height above top surface in Å (default: 10.0)")
        print("  --method METHOD  Use 'ase' or 'pymatgen' (default: ase)")
        print("  --output NAME    Output filename prefix (default: auto)")
        sys.exit(1)
    
    # Parse Miller indices
    try:
        h, k, l = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
        miller_index = (h, k, l)
    except (ValueError, IndexError):
        print("Error: Please provide three integers for Miller indices")
        sys.exit(1)
    
    # Parse options
    layers = 4
    vacuum = 15.0
    supercell = (2, 2, 1)
    add_co2 = False
    n_co2 = 1  # Default number of CO2 molecules per side
    co2_height = 10.0
    method = 'ase'
    output_prefix = None
    
    i = 4
    while i < len(sys.argv):
        if sys.argv[i] == '--layers' and i + 1 < len(sys.argv):
            layers = int(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == '--vacuum' and i + 1 < len(sys.argv):
            vacuum = float(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == '--supercell' and i + 2 < len(sys.argv):
            supercell = (int(sys.argv[i + 1]), int(sys.argv[i + 2]), 1)
            i += 3
        elif sys.argv[i] == '--add-co2':
            add_co2 = True
            i += 1
        elif sys.argv[i] == '--n-co2' and i + 1 < len(sys.argv):
            n_co2 = int(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == '--co2-height' and i + 1 < len(sys.argv):
            co2_height = float(sys.argv[i + 1])
            i += 2
        elif sys.argv[i] == '--method' and i + 1 < len(sys.argv):
            method = sys.argv[i + 1]
            i += 2
        elif sys.argv[i] == '--output' and i + 1 < len(sys.argv):
            output_prefix = sys.argv[i + 1]
            i += 2
        else:
            print(f"Unknown option: {sys.argv[i]}")
            i += 1
    
    # Check for CIF file
    cif_file = "HAP_ConventionalCell.cif"
    if not os.path.exists(cif_file):
        print(f"Error: {cif_file} not found in current directory")
        sys.exit(1)
    
    print(f"Miller index: {miller_index}")
    print(f"Layers: {layers}")
    print(f"Vacuum: {vacuum} Å")
    print(f"Supercell: {supercell}")
    print(f"Method: {method}")
    print(f"Add CO2: {add_co2}")
    if add_co2:
        print(f"Number of CO2 molecules per side: {n_co2} (total: {n_co2**2})")
        print(f"CO2 height: {co2_height} Å")
    
    # Create slab
    try:
        if method == 'pymatgen' and HAS_PYMATGEN:
            print("Using pymatgen method for slab generation")
            slab = create_slab_pymatgen(cif_file, miller_index, 
                                     min_slab_size=layers*3, 
                                     min_vacuum_size=vacuum,
                                     supercell=supercell)
        elif method == 'ase' and HAS_ASE:
            print("Using ASE method for slab generation")
            slab = create_slab_ase(cif_file, miller_index, 
                                 layers=layers, 
                                 vacuum=vacuum,
                                 supercell=supercell)
        else:
            print(f"Error: {method} method not available")
            sys.exit(1)
            
        print(f"Created slab: {len(slab)} atoms")
        print(f"Cell dimensions: {slab.get_cell().lengths()}")
        
        # Generate output filenames - FIX: Use clean base name
        if output_prefix:
            base_name = output_prefix.strip()  # Remove any whitespace/newlines
        else:
            base_name = f"hap_{h}{k}{l}".strip()  # Remove any whitespace/newlines
        
        # Clean surface
        data_file = f"{base_name}.data"
        xyz_file = f"{base_name}.xyz"
        
        write_lammps_data(slab, data_file, f"HAP {miller_index} surface")
        write(xyz_file, slab)
        
        print(f"✓ Generated: {data_file}")
        print(f"✓ Generated: {xyz_file}")
        
        # Variables to track which file to use for workflow
        workflow_source_file = data_file
        workflow_description = "without CO2"
        
        # Add CO2 if requested
        if add_co2:
            print("Adding CO2 molecule above surface...")
            
            # Validate original structure
            if not validate_structure(slab):
                print("WARNING: Original structure has issues, attempting to continue...")
            
            if n_co2 == 1:
                slab_co2 = add_co2_molecule(slab, height=co2_height)  # Increased height to 10.0Å
            else:
                slab_co2 = add_multiple_co2_molecules(slab, n_co2=n_co2, height=co2_height)
            
            # Validate structure with CO2
            if not validate_structure(slab_co2):
                print("WARNING: Structure with CO2 has issues")
            
            data_co2_file = f"{base_name}_co2.data"
            xyz_co2_file = f"{base_name}_co2.xyz"
            
            write_lammps_data(slab_co2, data_co2_file, f"HAP {miller_index} surface with CO2")
            write(xyz_co2_file, slab_co2)
            
            print(f"✓ Generated: {data_co2_file}")
            print(f"✓ Generated: {xyz_co2_file}")
            
            # Update workflow source to CO2 version
            workflow_source_file = data_co2_file
            workflow_description = "with CO2"
        
        # Create standard workflow compatible file
        # CRITICAL FIX: Always create workflow file based on add_co2 flag
        workflow_file = "hap_hkl.data"
        import shutil
        
        # Copy the appropriate file based on CO2 setting
        shutil.copy2(workflow_source_file, workflow_file)
        print(f"✓ Created workflow file: {workflow_file} ({workflow_description})")
        
        # Verify the workflow file content
        with open(workflow_file, 'r') as f:
            first_few_lines = [f.readline().strip() for _ in range(10)]
        
        # Count atom types in the workflow file
        atom_types_line = None
        for line in first_few_lines:
            if 'atom types' in line:
                atom_types_line = line
                break
        
        if atom_types_line:
            atom_types_count = int(atom_types_line.split()[0])
            print(f"✓ Workflow file contains {atom_types_count} atom types")
            if add_co2 and atom_types_count != 5:
                print(f"WARNING: Expected 5 atom types with CO2, but found {atom_types_count}")
            elif not add_co2 and atom_types_count != 4:
                print(f"WARNING: Expected 4 atom types without CO2, but found {atom_types_count}")
        
        # Print summary - FIX: Use the correct slab for composition
        summary_slab = slab_co2 if add_co2 else slab
        composition = {}
        for symbol in summary_slab.get_chemical_symbols():
            composition[symbol] = composition.get(symbol, 0) + 1
        
        print(f"\nSlab composition: {composition}")
        print(f"Cell parameters: a={summary_slab.get_cell()[0,0]:.2f}, b={summary_slab.get_cell()[1,1]:.2f}, c={summary_slab.get_cell()[2,2]:.2f} Å")
        print(f"\nTo use in workflow:")
        print(f"  Workflow file: {workflow_file} ({workflow_description})")
        
    except Exception as e:
        print(f"Error creating slab: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()