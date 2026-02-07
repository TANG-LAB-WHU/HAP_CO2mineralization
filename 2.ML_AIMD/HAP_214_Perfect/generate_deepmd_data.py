import dpdata
import pathlib
import os
import shutil
import sys
import traceback
import numpy as np

# -----------------------------------------------------------------------------
# Work-around for cp2kdata <0.8 TypeError bug when parsing forces
# -----------------------------------------------------------------------------

try:
    # Only patch if the module is available – keep runtime flexibility
    from cp2kdata.block_parser import forces as _cp2k_forces_mod  # type: ignore

    if not hasattr(_cp2k_forces_mod, "_dp_force_patch_applied"):
        _orig_parse_atomic_forces_list = _cp2k_forces_mod.parse_atomic_forces_list

        def _patched_parse_atomic_forces_list(output_file, *args, **kwargs):  # type: ignore
            """Shim that converts *output_file* to a raw string when a file-like
            object is mistakenly provided by older cp2kdata versions.
            """

            # 1) If it's a file handle, read() its entire content
            if hasattr(output_file, "read"):
                try:
                    output_file = output_file.read()
                except Exception:
                    output_file = str(output_file)

            # 2) If it's a pathlib.Path / os.PathLike, read its text content
            elif isinstance(output_file, (pathlib.Path, os.PathLike)):
                try:
                    output_file = pathlib.Path(output_file).read_text()
                except Exception:
                    output_file = str(output_file)

            # 3) Fallback – ensure it's at least str so regex doesn't fail
            if not isinstance(output_file, str):
                output_file = str(output_file)

            return _orig_parse_atomic_forces_list(output_file, *args, **kwargs)

        _cp2k_forces_mod.parse_atomic_forces_list = _patched_parse_atomic_forces_list  # type: ignore
        _cp2k_forces_mod._dp_force_patch_applied = True  # type: ignore

        # ALSO replace any stale aliases already imported (e.g. in cp2kdata.output)
        import importlib, sys
        for mod_name in ("cp2kdata.output", "cp2kdata.dpdata_plugin"):
            try:
                mod = sys.modules.get(mod_name)
                if mod is None:
                    mod = importlib.import_module(mod_name)
                if hasattr(mod, "parse_atomic_forces_list"):
                    setattr(mod, "parse_atomic_forces_list", _patched_parse_atomic_forces_list)  # type: ignore
            except ImportError:
                pass

        # Optional: notify user once for easier debugging
        print("Applied compatibility patch for cp2kdata.parse_atomic_forces_list (file-handle bug)")
except ImportError:
    # Fallback: try the pluralized module path introduced in newer cp2kdata versions
    try:
        from cp2kdata.block_parsers import forces as _cp2k_forces_mod  # type: ignore

        if not hasattr(_cp2k_forces_mod, "_dp_force_patch_applied"):
            _orig_parse_atomic_forces_list = _cp2k_forces_mod.parse_atomic_forces_list

            def _patched_parse_atomic_forces_list(output_file, *args, **kwargs):  # type: ignore
                try:
                    if hasattr(output_file, "read"):
                        output_file = output_file.read()
                except Exception:
                    output_file = str(output_file)

                return _orig_parse_atomic_forces_list(output_file, *args, **kwargs)

            _cp2k_forces_mod.parse_atomic_forces_list = _patched_parse_atomic_forces_list  # type: ignore
            _cp2k_forces_mod._dp_force_patch_applied = True  # type: ignore
            # Also patch any module-level aliases (e.g., within cp2kdata.output)
            try:
                import importlib
                _cp2k_output_mod = importlib.import_module("cp2kdata.output")
                if hasattr(_cp2k_output_mod, "parse_atomic_forces_list"):
                    _cp2k_output_mod.parse_atomic_forces_list = _patched_parse_atomic_forces_list  # type: ignore
            except (ModuleNotFoundError, ImportError):
                pass
            print("Applied compatibility patch for cp2kdata.block_parsers.forces (file-handle bug)")
    except ImportError:
        # cp2kdata not installed – nothing to patch
        pass


# -----------------------------------------------------------------------------
# Helper utilities
# -----------------------------------------------------------------------------

# Physical conversion factors
HARTREE_TO_EV = 27.211386245988  # eV / Hartree
BOHR_TO_ANGSTROM = 0.529177210903  # Å / Bohr
FORCE_HARTREE_PER_BOHR_TO_EV_PER_ANG = HARTREE_TO_EV / BOHR_TO_ANGSTROM

def detect_ensemble_type(cp2k_dir: pathlib.Path, project: str) -> str:
    """Determine the MD ensemble type to pass to dpdata."""
    # 1. Explicit override via environment variable.
    env_val = os.getenv("MD_ENSEMBLE_TYPE")
    if env_val:
        return env_val.strip().upper()

    # 2. Try to read <project>.inp and look for the ENSEMBLE line.
    inp_file = cp2k_dir / f"{project}.inp"
    if inp_file.is_file():
        try:
            with inp_file.open("r", encoding="utf-8") as fh:
                for raw_line in fh:
                    # Remove inline comments starting with ! or #
                    line = raw_line.split("!")[0].split("#")[0].strip()
                    if line.upper().startswith("ENSEMBLE"):
                        parts = line.split()
                        if len(parts) >= 2:
                            return parts[1].upper()
        except Exception as exc:
            # Non-fatal – fall back to default, but inform the user.
            print(f"WARNING: Failed to parse ensemble type from {inp_file}: {exc}")

    # 3. Default fallback.
    return "NVT"


def manually_parse_xyz_trajectory(xyz_file):
    """Manually parse multi-frame XYZ trajectory when dpdata fails."""
    if not xyz_file.exists():
        return None, None, None
    
    try:
        print(f"[DEBUG] Manually parsing XYZ trajectory: {xyz_file}")
        
        with open(xyz_file, 'r') as f:
            lines = f.readlines()
        
        if len(lines) < 2:
            return None, None, None
        
        # First line should contain number of atoms
        natoms = int(lines[0].strip())
        print(f"[DEBUG] Atoms per frame: {natoms}")
        
        # Calculate expected lines per frame (natoms + 2 header lines)
        lines_per_frame = natoms + 2
        total_frames = len(lines) // lines_per_frame
        
        print(f"[DEBUG] Total lines: {len(lines)}")
        print(f"[DEBUG] Lines per frame: {lines_per_frame}")
        print(f"[DEBUG] Calculated frames: {total_frames}")
        print(f"[DEBUG] Remainder lines: {len(lines) % lines_per_frame}")
        
        if total_frames == 0:
            return None, None, None
        
        # Show sample frame headers for verification
        print(f"[DEBUG] Sample frame headers:")
        for frame in range(min(3, total_frames)):
            frame_start = frame * lines_per_frame
            comment_line = lines[frame_start + 1].strip()
            print(f"  Frame {frame}: {comment_line[:60]}...")
        
        # Parse all frames
        all_coords = []
        all_elements = []
        
        for frame in range(total_frames):
            frame_start = frame * lines_per_frame
            
            # Skip header lines for this frame (natoms line + comment line)
            coord_start = frame_start + 2
            coord_end = coord_start + natoms
            
            if coord_end > len(lines):
                print(f"[WARNING] Frame {frame} is incomplete, stopping at frame {frame}")
                break
            
            frame_coords = []
            frame_elements = []
            
            for i in range(coord_start, coord_end):
                line = lines[i].strip()
                if line:
                    parts = line.split()
                    if len(parts) >= 4:
                        element = parts[0]
                        x, y, z = map(float, parts[1:4])
                        
                        frame_coords.append([x, y, z])
                        if frame == 0:  # Only collect elements from first frame
                            frame_elements.append(element)
            
            if len(frame_coords) == natoms:
                all_coords.append(frame_coords)
                if frame == 0:
                    all_elements = frame_elements
                    print(f"[DEBUG] First frame elements: {set(frame_elements)}")
            else:
                print(f"[WARNING] Frame {frame} has {len(frame_coords)} atoms, expected {natoms}")
                break
        
        if all_coords:
            coords_array = np.array(all_coords, dtype=np.float64)
            print(f"[DEBUG] Successfully parsed {len(all_coords)} frames")
            print(f"[DEBUG] Coordinates shape: {coords_array.shape}")
            
            # Verify coordinate ranges are reasonable
            x_range = coords_array[:, :, 0].min(), coords_array[:, :, 0].max()
            y_range = coords_array[:, :, 1].min(), coords_array[:, :, 1].max()
            z_range = coords_array[:, :, 2].min(), coords_array[:, :, 2].max()
            print(f"[DEBUG] Coordinate ranges: X({x_range[0]:.1f}, {x_range[1]:.1f}), Y({y_range[0]:.1f}, {y_range[1]:.1f}), Z({z_range[0]:.1f}, {z_range[1]:.1f})")
            
            return coords_array, all_elements, natoms
        else:
            return None, None, None
            
    except Exception as e:
        print(f"[ERROR] Failed to manually parse XYZ: {e}")
        traceback.print_exc()
        return None, None, None


def parse_cp2k_forces_file(forces_file, natoms, expected_frames):
    """Parse CP2K forces file with proper handling of the specific CP2K format."""
    if not forces_file.exists():
        return None
    
    try:
        print(f"[DEBUG] Parsing CP2K forces file: {forces_file}")
        
        with open(forces_file, 'r') as f:
            lines = f.readlines()
        
        print(f"[DEBUG] Forces file has {len(lines)} total lines")
        
        force_data = []
        current_frame_forces = []
        frame_count = 0
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # Look for the start of a new force section
            if line.startswith("FORCES|") and "Atomic forces" in line:
                if frame_count < 5:  # Only debug first few frames
                    print(f"[DEBUG] Found force section start at line {i+1}")
                
                # If we have a previous frame, save it
                if len(current_frame_forces) == natoms:
                    force_data.append(np.array(current_frame_forces))
                    frame_count += 1
                    if frame_count <= 5:
                        print(f"[DEBUG] Completed frame {frame_count} with {len(current_frame_forces)} atoms")
                
                # Start a new frame
                current_frame_forces = []
                
                # Skip the header lines
                i += 1
                # Skip the "Atom x y z |f|" line if it exists
                if i < len(lines) and lines[i].strip().startswith("FORCES|") and "Atom" in lines[i]:
                    i += 1
                
                continue
                
            # Parse force data lines - but skip summary lines
            elif line.startswith("FORCES|") and len(line.split()) >= 5:
                parts = line.split()
                
                # Skip summary lines that contain "Sum" or "Total"
                if len(parts) > 1 and ("Sum" in parts[1] or "Total" in parts[1]):
                    i += 1
                    continue
                
                try:
                    # Format: FORCES| atom_id fx fy fz |f|
                    atom_id = int(parts[1])
                    fx = float(parts[2])
                    fy = float(parts[3])
                    fz = float(parts[4])
                    
                    current_frame_forces.append([fx, fy, fz])
                    
                    # Debug first few entries of first frame
                    if frame_count == 0 and len(current_frame_forces) <= 5:
                        print(f"[DEBUG] Atom {atom_id}: [{fx:.6e}, {fy:.6e}, {fz:.6e}]")
                        
                except (ValueError, IndexError) as e:
                    if frame_count < 5:  # Only debug first few frames
                        print(f"[DEBUG] Failed to parse force line {i+1}: {line[:50]}... Error: {e}")
                    
            i += 1
        
        # Add the last frame if complete
        if len(current_frame_forces) == natoms:
            force_data.append(np.array(current_frame_forces))
            frame_count += 1
            print(f"[DEBUG] Completed final frame {frame_count} with {len(current_frame_forces)} atoms")
        elif len(current_frame_forces) > 0:
            print(f"[DEBUG] Incomplete final frame: got {len(current_frame_forces)} atoms, expected {natoms}")
        
        if force_data:
            forces_array = np.array(force_data)
            print(f"[DEBUG] Successfully parsed {len(force_data)} frames of forces")
            print(f"[DEBUG] Forces shape: {forces_array.shape}")
            print(f"[DEBUG] Expected shape: ({expected_frames}, {natoms}, 3)")
            
            # Convert from Hartree/Bohr to eV/Å for DeepMD compatibility
            forces_array *= FORCE_HARTREE_PER_BOHR_TO_EV_PER_ANG
            print(f"[DEBUG] Converted forces to eV/Å (factor {FORCE_HARTREE_PER_BOHR_TO_EV_PER_ANG:.6f})")
            
            # Check force magnitudes
            force_magnitudes = np.linalg.norm(forces_array, axis=2)
            print(f"[DEBUG] Force magnitude range (eV/Å): {force_magnitudes.min():.6e} to {force_magnitudes.max():.6e}")
            
            return forces_array
        else:
            print(f"[DEBUG] No valid force frames found")
            return None
            
    except Exception as e:
        print(f"[ERROR] Failed to parse forces file: {e}")
        traceback.print_exc()
        return None


def parse_cp2k_atomic_forces_file(forces_file, natoms, expected_frames):
    """Parse CP2K forces file with 'ATOMIC FORCES in [a.u.]' header format.
    
    Format:
    ATOMIC FORCES in [a.u.]
    
    # Atom   Kind   Element          X              Y              Z
         1      1      O           0.00386741     0.00342249    -0.02093797
         ...
       387      1      O           0.00024085     0.00014984     0.00193903
    SUM OF ATOMIC FORCES          -0.00537890     0.10385758    -0.01592790     0.10520945
    """
    if not forces_file.exists():
        return None
    
    try:
        print(f"[DEBUG] Parsing CP2K atomic forces file: {forces_file}")
        
        with open(forces_file, 'r') as f:
            lines = f.readlines()
        
        print(f"[DEBUG] Forces file has {len(lines)} total lines")
        
        force_data = []
        current_frame_forces = []
        frame_count = 0
        i = 0
        
        while i < len(lines):
            line = lines[i].strip()
            
            # Detect start of new force section
            if "ATOMIC FORCES" in line.upper() and "in [a.u.]" in line:
                if frame_count < 5:
                    print(f"[DEBUG] Found ATOMIC FORCES section at line {i+1}")
                
                # If we have a previous frame, save it
                if len(current_frame_forces) == natoms:
                    force_data.append(np.array(current_frame_forces))
                    frame_count += 1
                    if frame_count <= 5:
                        print(f"[DEBUG] Completed frame {frame_count} with {len(current_frame_forces)} atoms")
                
                # Start a new frame
                current_frame_forces = []
                
                # Skip header lines: empty line and comment line
                i += 1  # Skip "ATOMIC FORCES in [a.u.]" line
                # Skip empty line
                if i < len(lines) and not lines[i].strip():
                    i += 1
                # Skip comment line "# Atom   Kind   Element          X              Y              Z"
                if i < len(lines) and lines[i].strip().startswith("#"):
                    i += 1
                
                continue
            
            # Skip SUM OF ATOMIC FORCES line
            elif "SUM OF ATOMIC FORCES" in line.upper():
                i += 1
                continue
            
            # Parse force data lines (skip empty lines and comments)
            elif line and not line.startswith("#"):
                parts = line.split()
                
                # Expected format: atom_id kind element fx fy fz
                if len(parts) >= 6:
                    try:
                        atom_id = int(parts[0])
                        # Skip kind (parts[1]) and element (parts[2]), get forces (parts[3:6])
                        fx = float(parts[3])
                        fy = float(parts[4])
                        fz = float(parts[5])
                        
                        current_frame_forces.append([fx, fy, fz])
                        
                        # Debug first few entries of first frame
                        if frame_count == 0 and len(current_frame_forces) <= 5:
                            print(f"[DEBUG] Atom {atom_id}: [{fx:.6e}, {fy:.6e}, {fz:.6e}]")
                            
                    except (ValueError, IndexError) as e:
                        if frame_count < 5:
                            print(f"[DEBUG] Failed to parse force line {i+1}: {line[:50]}... Error: {e}")
            
            i += 1
        
        # Add the last frame if complete
        if len(current_frame_forces) == natoms:
            force_data.append(np.array(current_frame_forces))
            frame_count += 1
            print(f"[DEBUG] Completed final frame {frame_count} with {len(current_frame_forces)} atoms")
        elif len(current_frame_forces) > 0:
            print(f"[DEBUG] Incomplete final frame: got {len(current_frame_forces)} atoms, expected {natoms}")
        
        if force_data:
            forces_array = np.array(force_data)
            print(f"[DEBUG] Successfully parsed {len(force_data)} frames of forces")
            print(f"[DEBUG] Forces shape: {forces_array.shape}")
            print(f"[DEBUG] Expected shape: ({expected_frames}, {natoms}, 3)")
            
            # Convert from Hartree/Bohr to eV/Å for DeepMD compatibility
            forces_array *= FORCE_HARTREE_PER_BOHR_TO_EV_PER_ANG
            print(f"[DEBUG] Converted forces to eV/Å (factor {FORCE_HARTREE_PER_BOHR_TO_EV_PER_ANG:.6f})")
            
            # Check force magnitudes
            force_magnitudes = np.linalg.norm(forces_array, axis=2)
            print(f"[DEBUG] Force magnitude range (eV/Å): {force_magnitudes.min():.6e} to {force_magnitudes.max():.6e}")
            
            return forces_array
        else:
            print(f"[DEBUG] No valid force frames found")
            return None
            
    except Exception as e:
        print(f"[ERROR] Failed to parse atomic forces file: {e}")
        traceback.print_exc()
        return None


def parse_standard_forces_file(forces_file, natoms, expected_frames):
    """Parse standard forces file format (simple numerical format)."""
    if not forces_file.exists():
        return None
    
    try:
        print(f"[DEBUG] Parsing standard forces file: {forces_file}")
        
        # Try to load as simple numerical data
        try:
            data = np.loadtxt(forces_file)
            print(f"[DEBUG] Loaded data shape: {data.shape}")
            
            if len(data.shape) == 2 and data.shape[1] >= 3:
                # Reshape to frames
                total_atoms = data.shape[0]
                if total_atoms % natoms == 0:
                    nframes = total_atoms // natoms
                    forces = data[:, :3].reshape(nframes, natoms, 3)
                    forces *= FORCE_HARTREE_PER_BOHR_TO_EV_PER_ANG
                    print(f"[DEBUG] Reshaped to {nframes} frames and converted to eV/Å")
                    return forces
                else:
                    print(f"[DEBUG] Cannot reshape: {total_atoms} atoms not divisible by {natoms}")
                    return None
            else:
                print(f"[DEBUG] Unexpected data format")
                return None
                
        except Exception as e:
            print(f"[DEBUG] Failed to load as numerical data: {e}")
            return None
            
    except Exception as e:
        print(f"[ERROR] Failed to parse standard forces file: {e}")
        return None


def parse_cp2k_virial_file(virial_file, expected_frames):
    """Parse CP2K virial tensor file to extract stress information.
    
    CP2K format example:
    STRESS| Analytical stress tensor [bar]
    STRESS|                        x                   y                   z
    STRESS|      x       -2.39627363376E+03  -6.96139579662E+02  -1.65851624361E+01
    STRESS|      y       -6.96139579662E+02  -3.36301965788E+03  -1.04709895422E+02
    STRESS|      z       -1.65851624361E+01  -1.04709895422E+02   1.57031582075E+02
    """
    if not virial_file.exists():
        return None
    
    try:
        print(f"[DEBUG] Parsing CP2K virial file: {virial_file}")
        
        with open(virial_file, 'r') as f:
            lines = f.readlines()
        
        print(f"[DEBUG] Virial file has {len(lines)} total lines")
        
        virial_data = []
        current_stress_matrix = []
        in_stress_block = False
        frame_count = 0
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # Detect start of new stress tensor block
            if line.startswith("STRESS|") and "Analytical stress tensor" in line:
                if frame_count < 5:  # Only debug first few frames
                    print(f"[DEBUG] Found stress tensor block at line {i+1}")
                
                # If we have a complete previous stress matrix, save it
                if len(current_stress_matrix) == 3:
                    stress_tensor = np.array(current_stress_matrix)
                    virial_data.append(stress_tensor)
                    frame_count += 1
                    
                    # Debug first few stress tensors
                    if frame_count <= 3:
                        print(f"[DEBUG] Frame {frame_count}: Stress tensor (bar)")
                        print(f"        [[{stress_tensor[0,0]:.6e}, {stress_tensor[0,1]:.6e}, {stress_tensor[0,2]:.6e}],")
                        print(f"         [{stress_tensor[1,0]:.6e}, {stress_tensor[1,1]:.6e}, {stress_tensor[1,2]:.6e}],")
                        print(f"         [{stress_tensor[2,0]:.6e}, {stress_tensor[2,1]:.6e}, {stress_tensor[2,2]:.6e}]]")
                
                # Reset for new stress tensor
                current_stress_matrix = []
                in_stress_block = True
                
                # Skip the header line with "x y z"
                i += 1
                if i < len(lines) and "x" in lines[i] and "y" in lines[i] and "z" in lines[i]:
                    i += 1
                
                continue
            
            # Parse stress tensor components
            elif in_stress_block and line.startswith("STRESS|"):
                parts = line.split()
                
                # Check if this is a stress tensor component line
                # Format: STRESS|      x/y/z    val1    val2    val3
                if len(parts) >= 5 and parts[1] in ['x', 'y', 'z']:
                    try:
                        axis = parts[1]  # x, y, or z
                        val1 = float(parts[2])
                        val2 = float(parts[3])
                        val3 = float(parts[4])
                        
                        stress_row = [val1, val2, val3]
                        current_stress_matrix.append(stress_row)
                        
                        # Debug parsing
                        if frame_count == 0 and len(current_stress_matrix) <= 3:
                            print(f"[DEBUG] Stress row {axis}: [{val1:.6e}, {val2:.6e}, {val3:.6e}]")
                        
                    except (ValueError, IndexError) as e:
                        if frame_count < 3:
                            print(f"[DEBUG] Failed to parse stress line {i+1}: {line[:50]}... Error: {e}")
                
                # Check for end of stress block
                elif ("Trace" in line or "Determinant" in line or "Eigenvectors" in line):
                    in_stress_block = False
                    
            i += 1
        
        # Add the last stress tensor if complete
        if len(current_stress_matrix) == 3:
            stress_tensor = np.array(current_stress_matrix)
            virial_data.append(stress_tensor)
            frame_count += 1
            print(f"[DEBUG] Added final stress tensor frame {frame_count}")
        
        if virial_data:
            virial_array = np.array(virial_data)
            print(f"[DEBUG] Successfully parsed {len(virial_data)} frames of virial data")
            print(f"[DEBUG] Virial shape: {virial_array.shape}")
            print(f"[DEBUG] Expected shape: ({expected_frames}, 3, 3)")
            
            # Check stress magnitudes and units
            stress_traces = np.trace(virial_array, axis1=1, axis2=2)  # Trace of stress tensor
            print(f"[DEBUG] Stress trace range: {stress_traces.min():.6e} to {stress_traces.max():.6e} bar")
            
            # Convert from bar to eV/Å³ for DeepMD compatibility
            # 1 bar = 1e5 Pa = 6.241506e-7 eV/Å³
            bar_to_ev_per_ang3 = 6.241506e-7
            virial_array_ev = virial_array * bar_to_ev_per_ang3
            
            print(f"[DEBUG] Converted stress range: {(stress_traces * bar_to_ev_per_ang3).min():.6e} to {(stress_traces * bar_to_ev_per_ang3).max():.6e} eV/Å³")
            
            # Convert to DeepMD format: flatten to 9 components per frame
            # DeepMD expects virial in the format: [xx, xy, xz, yx, yy, yz, zx, zy, zz]
            virial_flat = virial_array_ev.reshape(len(virial_data), 9)
            
            return virial_flat
        else:
            print(f"[DEBUG] No valid virial data found")
            return None
            
    except Exception as e:
        print(f"[ERROR] Failed to parse virial file: {e}")
        traceback.print_exc()
        return None


def read_energies_from_file_advanced(energy_file, expected_frames):
    """Advanced energy file reader for CP2K format."""
    if not energy_file.exists():
        return None
    
    try:
        print(f"[DEBUG] Reading energy file: {energy_file}")
        
        energies = []
        with open(energy_file, 'r') as f:
            for line_num, line in enumerate(f):
                line = line.strip()
                
                # Skip comment lines and empty lines
                if line and not line.startswith('#') and not line.startswith('!'):
                    parts = line.split()
                    
                    # CP2K energy file format (based on your actual data):
                    # # Step Nr. Time[fs] Kin.[a.u.] Temp[K] Pot.[a.u.] Cons Qty[a.u.] UsedTime[s]
                    # Column indices: 0=Step, 1=Time, 2=Kinetic, 3=Temp, 4=Potential, 5=Conserved, 6=UsedTime
                    
                    if len(parts) >= 6:  # Ensure we have enough columns (need at least 6 for potential energy)
                        try:
                            step = int(parts[0])
                            time = float(parts[1])
                            potential_energy = float(parts[4])  # 5th column (0-indexed: 4) is Pot.[a.u.]
                            energies.append(potential_energy)
                            
                            # Debug first few entries
                            if len(energies) <= 5:
                                print(f"[DEBUG] Step {step:3d}: Time={time:8.3f}fs, Potential={potential_energy:.6f} a.u.")
                                
                        except ValueError:
                            if line_num < 10:  # Only show first few parsing failures
                                print(f"[DEBUG] Skipping non-numeric energy line {line_num}: {line[:50]}...")
                            continue
                    else:
                        if line_num < 10:
                            print(f"[DEBUG] Skipping short line {line_num}: {line[:50]}...")
        
        if energies:
            energies_array = np.array(energies)
            print(f"[DEBUG] Successfully parsed {len(energies)} energy values")
            print(f"[DEBUG] Expected {expected_frames} values")
            print(f"[DEBUG] Energy range: {energies_array.min():.6f} to {energies_array.max():.6f} a.u.")
            print(f"[DEBUG] Mean energy: {energies_array.mean():.6f} a.u.")
            
            # Verify energy values are reasonable
            if all(e < 0 for e in energies):
                print(f"[DEBUG] ✅ All energies are negative (expected for this system)")
            else:
                print(f"[DEBUG] ⚠️  Some energies are positive (unexpected)")
            
            energies_array *= HARTREE_TO_EV
            print(f"[DEBUG] Converted energies to eV (factor {HARTREE_TO_EV:.6f})")
                
            return energies_array
        else:
            print(f"[DEBUG] No valid energy data found")
            return None
            
    except Exception as e:
        print(f"[WARNING] Failed to read energies: {e}")
        return None


def read_cell_matrix_from_file(cell_file):
    """Manually read the cell matrix from CP2K cell trajectory file."""
    if not cell_file.exists():
        return None
    
    try:
        print(f"[DEBUG] Reading CP2K cell trajectory file: {cell_file}")
        
        with open(cell_file, 'r') as f:
            lines = f.readlines()
        
        # Parse the CP2K cell trajectory format
        # Format: Step Time Ax Ay Az Bx By Bz Cx Cy Cz Volume
        # Columns: 0=Step, 1=Time, 2-4=A vector, 5-7=B vector, 8-10=C vector, 11=Volume
        
        cell_matrices = []
        
        for line in lines:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 11:  # Ensure we have all required columns
                    try:
                        # Extract the three lattice vectors
                        ax, ay, az = float(parts[2]), float(parts[3]), float(parts[4])
                        bx, by, bz = float(parts[5]), float(parts[6]), float(parts[7])
                        cx, cy, cz = float(parts[8]), float(parts[9]), float(parts[10])
                        
                        # Create 3x3 cell matrix
                        cell_matrix = np.array([
                            [ax, ay, az],
                            [bx, by, bz],
                            [cx, cy, cz]
                        ])
                        
                        cell_matrices.append(cell_matrix)
                        
                    except ValueError as e:
                        print(f"[WARNING] Failed to parse cell line: {line[:50]}... Error: {e}")
                        continue
        
        if cell_matrices:
            cells_array = np.array(cell_matrices)
            print(f"[DEBUG] Successfully parsed {len(cell_matrices)} cell frames")
            print(f"[DEBUG] Cell array shape: {cells_array.shape}")
            print(f"[DEBUG] First cell matrix:\n{cells_array[0]}")
            
            # Check if all cells are the same (rigid cell)
            if len(cell_matrices) > 1:
                if np.allclose(cells_array[0], cells_array[-1], atol=1e-10):
                    print(f"[DEBUG] Cell appears to be rigid (constant throughout trajectory)")
                else:
                    print(f"[DEBUG] Cell varies during trajectory")
            
            return cells_array
        else:
            print(f"[DEBUG] No valid cell data found")
            return None
        
    except Exception as e:
        print(f"[WARNING] Failed to read cell trajectory: {e}")
        return None


def create_element_type_mapping(elements):
    """Create element to type mapping for DeepMD."""
    # Preferred order for HAP system
    preferred_order = ["Ca", "P", "O", "H", "C"]
    unique_elements = list(dict.fromkeys(elements))  # Preserve order, remove duplicates
    
    # Reorder according to preferred order
    ordered_elements = []
    for el in preferred_order:
        if el in unique_elements:
            ordered_elements.append(el)
    
    # Add any remaining elements not in preferred order
    for el in unique_elements:
        if el not in ordered_elements:
            ordered_elements.append(el)
    
    # Create type mapping (DeepMD uses 0-based indexing)
    type_map = {el: idx for idx, el in enumerate(ordered_elements)}
    
    print(f"[DEBUG] Element type mapping: {type_map}")
    return type_map, ordered_elements


def create_deepmd_dataset_v3(output_dir, coords, elements, cells=None, forces=None, energies=None, virial=None):
    """Create DeepMD dataset compatible with DeepMD-kit v3.x format."""
    
    try:
        print(f"[DEBUG] Creating DeepMD v3 dataset...")
        
        # Create set.000 directory
        set_dir = output_dir / "set.000"
        set_dir.mkdir(exist_ok=True)
        
        nframes, natoms, _ = coords.shape
        print(f"[DEBUG] Dataset dimensions: {nframes} frames, {natoms} atoms")
        
        # Create element type mapping
        type_map, ordered_elements = create_element_type_mapping(elements)
        
        # Create atom types array
        atom_types = np.array([type_map[el] for el in elements], dtype=np.int32)
        
        # Save coordinates in both .npy and .raw formats for compatibility
        coord_flat = coords.reshape(nframes, natoms * 3).astype(np.float64)
        
        # Save as .npy files (newer DeepMD versions prefer this)
        np.save(set_dir / "coord.npy", coord_flat)
        
        # Also save as .raw files for older compatibility
        np.savetxt(set_dir / "coord.raw", coord_flat, fmt="%.12f")
        
        # Save cell information
        if cells is not None:
            cell_flat = cells.reshape(nframes, 9).astype(np.float64)
            np.save(set_dir / "box.npy", cell_flat)
            np.savetxt(set_dir / "box.raw", cell_flat, fmt="%.12f")
        
        # Save forces if available
        if forces is not None:
            force_flat = forces.reshape(nframes, natoms * 3).astype(np.float64)
            np.save(set_dir / "force.npy", force_flat)
            np.savetxt(set_dir / "force.raw", force_flat, fmt="%.12f")
        
        # Save energies if available
        if energies is not None:
            # Ensure we have the right number of energies
            if len(energies) >= nframes:
                energies_subset = energies[:nframes].astype(np.float64)
                np.save(set_dir / "energy.npy", energies_subset)
                np.savetxt(set_dir / "energy.raw", energies_subset, fmt="%.12f")
            else:
                print(f"[WARNING] Energy count mismatch: {len(energies)} < {nframes}")
        
        # Save virial tensor if available
        if virial is not None:
            # Ensure we have the right number of virial frames
            if virial.shape[0] >= nframes:
                virial_subset = virial[:nframes].astype(np.float64)
                np.save(set_dir / "virial.npy", virial_subset)
                np.savetxt(set_dir / "virial.raw", virial_subset, fmt="%.12f")
                print(f"[DEBUG] Saved virial data with shape: {virial_subset.shape}")
            else:
                print(f"[WARNING] Virial count mismatch: {virial.shape[0]} < {nframes}")
        
        # Create type.raw (atom types for the system)
        np.savetxt(output_dir / "type.raw", atom_types, fmt="%d")
        
        # Create type_map.raw (mapping from type index to element symbol)
        with open(output_dir / "type_map.raw", 'w') as f:
            for element in ordered_elements:
                f.write(f"{element}\n")
        
        print(f"[DEBUG] Successfully created DeepMD v3 dataset")
        print(f"  - Frames: {nframes}")
        print(f"  - Atoms: {natoms}")
        print(f"  - Elements: {ordered_elements}")
        print(f"  - Coordinates: YES (.npy + .raw)")
        print(f"  - Cell: {'YES (.npy + .raw)' if cells is not None else 'NO'}")
        print(f"  - Forces: {'YES (.npy + .raw)' if forces is not None else 'NO'}")
        print(f"  - Energies: {'YES (.npy + .raw)' if energies is not None else 'NO'}")
        print(f"  - Virial: {'YES (.npy + .raw)' if virial is not None else 'NO'}")
        
        return True
        
    except Exception as e:
        print(f"[ERROR] Failed to create DeepMD dataset: {e}")
        traceback.print_exc()
        return False


def detect_project_name(cp2k_dir: pathlib.Path) -> str:
    """Auto-detect project name from CP2K output files by finding *_md-pos.xyz files."""
    # Method 1: Look for *_md-pos.xyz files (primary method)
    pos_files = list(cp2k_dir.glob("*_md-pos.xyz"))
    if pos_files:
        # Extract project name from filename (e.g., "HAP_300_md-pos.xyz" -> "HAP_300")
        filename = pos_files[0].stem  # "HAP_300_md-pos"
        project_name = filename.replace("_md-pos", "")
        print(f"[DEBUG] Detected project name from trajectory file: {project_name}")
        return project_name
    
    # Method 2: Fallback to *_md-1.xyz files
    xyz_files = list(cp2k_dir.glob("*_md-1.xyz"))
    if xyz_files:
        filename = xyz_files[0].stem  # "HAP_300_md-1"
        project_name = filename.replace("_md-1", "")
        print(f"[DEBUG] Detected project name from XYZ file: {project_name}")
        return project_name
    
    # Method 3: Look for HAP_*_md.inp file
    inp_files = list(cp2k_dir.glob("HAP_*_md.inp"))
    if inp_files:
        filename = inp_files[0].stem  # "HAP_300_md"
        project_name = filename.replace("_md", "")
        print(f"[DEBUG] Detected project name from input file: {project_name}")
        return project_name
    
    # Fallback to default
    print(f"[WARNING] Could not auto-detect project name, using default 'HAP_hkl'")
    return "HAP_hkl"


def main():
    """
    Generates a DeepMD dataset from CP2K AIMD output files, including virial tensor support.
    """
    cp2k_source_dir = pathlib.Path("cp2k_run")
    output_dir = pathlib.Path("deepmd_data")

    print(f"Output directory: {output_dir}")
    print(f"CP2K source directory: {cp2k_source_dir}")

    # Auto-detect project name from files in cp2k_run directory
    project_name = detect_project_name(cp2k_source_dir)
    print(f"Using project name: {project_name}")

    # Clean and recreate the output directory
    if output_dir.exists():
        print(f"Removing existing output directory: {output_dir}")
        shutil.rmtree(output_dir)
    output_dir.mkdir(exist_ok=True)

    try:
        print("Generating DeepMD dataset from CP2K AIMD output...")

        # First, analyze all available files
        print("=== Available Files Analysis ===")
        all_files = list(cp2k_source_dir.glob("*"))
        relevant_files = [f for f in all_files if f.is_file() and project_name in f.name]
        
        print(f"Found {len(relevant_files)} relevant files:")
        for f in sorted(relevant_files):
            size = f.stat().st_size
            print(f"  - {f.name}: {size:,} bytes")

        # Find the main trajectory file
        print(f"\n=== Coordinate File Analysis ===")
        
        # Use the original CP2K trajectory output file
        coord_file = cp2k_source_dir / f"{project_name}_md-pos.xyz"
        
        if not coord_file.exists():
            # Fallback to alternative naming
            coord_file = cp2k_source_dir / f"{project_name}_md-1.xyz"
            if not coord_file.exists():
                raise FileNotFoundError(f"No coordinate file found. Expected: {project_name}_md-pos.xyz")
        
        file_size = coord_file.stat().st_size
        print(f"Using original CP2K trajectory file: {coord_file.name} (size: {file_size:,} bytes)")
        
        # Manually parse the XYZ trajectory to get all frames
        coords, elements, natoms = manually_parse_xyz_trajectory(coord_file)
        
        if coords is None:
            raise ValueError("Failed to parse trajectory coordinates")
        
        nframes = coords.shape[0]
        print(f"Successfully loaded: {nframes} frames, {natoms} atoms per frame")
        print(f"Elements found: {set(elements)}")
        
        # Try to load cell information
        print(f"\n=== Cell Information ===")
        
        # Use the original CP2K cell output file
        cell_file = cp2k_source_dir / f"{project_name}_md-1.cell"
        
        if not cell_file.exists():
            print(f"❌ CP2K cell file not found: {cell_file}")
            print("[WARNING] No cell data will be loaded")
            cells = None
        else:
            file_size = cell_file.stat().st_size
            print(f"Found original CP2K cell file: {cell_file.name} (size: {file_size:,} bytes)")
            
            cell_matrices = read_cell_matrix_from_file(cell_file)
            if cell_matrices is not None:
                # Check if we have the right number of cell frames
                if len(cell_matrices) >= nframes:
                    cells = cell_matrices[:nframes]  # Take only the frames we need
                    print(f"✅ Added cell information for {nframes} frames")
                    print(f"Cell array shape: {cells.shape}")
                elif len(cell_matrices) == 1:
                    # Single cell matrix - replicate for all frames
                    cells = np.tile(cell_matrices[0], (nframes, 1, 1))
                    print(f"✅ Replicated single cell matrix for {nframes} frames")
                    print(f"Cell array shape: {cells.shape}")
                else:
                    print(f"[WARNING] Cell frame count mismatch: {len(cell_matrices)} cells vs {nframes} coord frames")
                    # Use what we have and pad/truncate as needed
                    if len(cell_matrices) < nframes:
                        # Pad with the last cell
                        last_cell = cell_matrices[-1]
                        padding = [last_cell] * (nframes - len(cell_matrices))
                        cell_matrices_padded = np.concatenate([cell_matrices, padding])
                        cells = cell_matrices_padded
                    else:
                        cells = cell_matrices[:nframes]
                    print(f"✅ Adjusted cell data to {nframes} frames")
                    print(f"Cell array shape: {cells.shape}")
            else:
                print("❌ Could not read cell matrix")
                cells = None
        
        # Try to load forces
        print(f"\n=== Forces Information ===")
        
        # Use the original CP2K output file, not the script-generated copy
        forces_file = cp2k_source_dir / f"{project_name}_forces.dat"
        
        if not forces_file.exists():
            print(f"❌ Original CP2K forces file not found: {forces_file}")
            print("[WARNING] No forces data will be loaded")
            forces = None
        else:
            file_size = forces_file.stat().st_size
            print(f"Found original CP2K forces file: {forces_file.name} (size: {file_size:,} bytes)")
            
            print(f"Attempting to parse forces from: {forces_file}")
            
            # Check the file format
            with open(forces_file, 'r') as f:
                first_lines = [f.readline().strip() for _ in range(5)]
            
            print(f"[DEBUG] First few lines of {forces_file.name}:")
            for i, line in enumerate(first_lines):
                print(f"  Line {i+1}: {line[:60]}...")
            
            # Check the file format and use appropriate parser
            if any("FORCES|" in line for line in first_lines):
                print(f"[DEBUG] Confirmed CP2K FORCES| format in {forces_file.name}")
                forces = parse_cp2k_forces_file(forces_file, natoms, nframes)
            elif any("ATOMIC FORCES" in line.upper() for line in first_lines):
                print(f"[DEBUG] Confirmed CP2K ATOMIC FORCES format in {forces_file.name}")
                forces = parse_cp2k_atomic_forces_file(forces_file, natoms, nframes)
            else:
                print(f"[WARNING] Unexpected format in {forces_file.name}, trying standard parser")
                forces = parse_standard_forces_file(forces_file, natoms, nframes)
            
            if forces is not None:
                print(f"✅ Successfully loaded forces: {forces.shape}")
                # Match the number of frames
                if forces.shape[0] >= nframes:
                    forces = forces[:nframes]
                    print(f"Trimmed forces to {nframes} frames")
                elif forces.shape[0] < nframes:
                    print(f"[WARNING] Forces have fewer frames ({forces.shape[0]}) than coordinates ({nframes})")
            else:
                print(f"❌ Failed to parse forces from {forces_file}")
                forces = None

        if forces is None:
            print("[WARNING] No forces data loaded")

        # Try to load virial tensor (stress tensor)
        print(f"\n=== Virial Tensor Information ===")
        
        # Use the original CP2K virial output file
        virial_file = cp2k_source_dir / f"{project_name}_virial.dat"
        
        if not virial_file.exists():
            print(f"❌ CP2K virial file not found: {virial_file}")
            print("[WARNING] No virial tensor data will be loaded")
            virial = None
        else:
            file_size = virial_file.stat().st_size
            print(f"Found CP2K virial file: {virial_file.name} (size: {file_size:,} bytes)")
            
            # Show sample virial file content for debugging
            print(f"[DEBUG] First few lines of {virial_file.name}:")
            with open(virial_file, 'r') as f:
                for i, line in enumerate(f):
                    if i >= 10:
                        break
                    print(f"  Line {i+1}: {line.strip()[:70]}...")
            
            print(f"Attempting to parse virial tensor from: {virial_file}")
            virial = parse_cp2k_virial_file(virial_file, nframes)
            
            if virial is not None:
                print(f"✅ Successfully loaded virial tensor: {virial.shape}")
                if virial.shape[0] >= nframes:
                    print(f"Virial frame count matches or exceeds coordinate frame count ({nframes})")
                else:
                    print(f"[WARNING] Virial frame count ({virial.shape[0]}) less than coordinate frame count ({nframes})")
            else:
                print(f"❌ Failed to parse virial tensor from {virial_file}")
                virial = None

        if virial is None:
            print("[WARNING] No virial tensor data loaded - this will limit model accuracy for stress/pressure predictions")

        # Try to load energies
        print(f"\n=== Energy Information ===")
        
        # Use the original CP2K energy output file
        energy_file = cp2k_source_dir / f"{project_name}_md-1.ener"
        
        if not energy_file.exists():
            print(f"❌ CP2K energy file not found: {energy_file}")
            print("[WARNING] No energy data will be loaded")
            energies = None
        else:
            file_size = energy_file.stat().st_size
            print(f"Found CP2K energy file: {energy_file.name} (size: {file_size:,} bytes)")
            
            print(f"Attempting to parse energies from: {energy_file}")
            energies = read_energies_from_file_advanced(energy_file, nframes)
            
            if energies is not None:
                print(f"✅ Successfully loaded energies: {len(energies)} values")
                if len(energies) >= nframes:
                    print(f"Energy count matches or exceeds frame count ({nframes})")
                else:
                    print(f"[WARNING] Energy count ({len(energies)}) less than frame count ({nframes})")
            else:
                print(f"❌ Failed to parse energies from {energy_file}")
                energies = None
        
        # Create the DeepMD dataset
        print(f"\n=== Creating DeepMD Dataset ===")
        success = create_deepmd_dataset_v3(
            output_dir, coords, elements, cells, forces, energies, virial
        )
        
        if success:
            print(f"\n=== SUCCESS ===")
            print(f"DeepMD dataset generated successfully!")
            print(f"Dataset saved to: {output_dir}")
            
            # Verify the generated files
            set_dir = output_dir / "set.000"
            if set_dir.exists():
                npy_files = list(set_dir.glob("*.npy"))
                raw_files = list(set_dir.glob("*.raw"))
                print(f"Generated .npy files: {[f.name for f in npy_files]}")
                print(f"Generated .raw files: {[f.name for f in raw_files]}")
                
                # Check coord.npy specifically
                coord_npy = set_dir / "coord.npy"
                if coord_npy.exists():
                    coord_data = np.load(coord_npy)
                    print(f"coord.npy verification: shape = {coord_data.shape}")
                else:
                    print("[ERROR] coord.npy not found!")
                    
                # Check other key files
                for key_file in ["box.npy", "force.npy", "energy.npy", "virial.npy"]:
                    file_path = set_dir / key_file
                    if file_path.exists():
                        data = np.load(file_path)
                        print(f"{key_file} verification: shape = {data.shape}")
                    else:
                        print(f"[INFO] {key_file} not present (optional)")
            
            # Check root files
            root_files = ["type.raw", "type_map.raw"]
            for root_file in root_files:
                file_path = output_dir / root_file
                if file_path.exists():
                    print(f"✅ {root_file} created")
                    # Show content for verification
                    if root_file == "type_map.raw":
                        with open(file_path, 'r') as f:
                            content = f.read().strip()
                            print(f"   Content: {content.replace(chr(10), ', ')}")
                else:
                    print(f"❌ {root_file} missing")
            
            # Summary of what was included
            print(f"\n=== Dataset Summary ===")
            print(f"📊 Total frames: {nframes}")
            print(f"🔬 Total atoms: {natoms}")
            print(f"⚛️  Elements: {set(elements)}")
            print(f"📍 Coordinates: ✅ Included")
            print(f"📦 Cell vectors: {'✅ Included' if cells is not None else '❌ Missing'}")
            print(f"⚡ Forces: {'✅ Included' if forces is not None else '❌ Missing'}")
            print(f"🔋 Energies: {'✅ Included' if energies is not None else '❌ Missing'}")
            print(f"📈 Virial tensor: {'✅ Included' if virial is not None else '❌ Missing'}")
            
            if virial is not None:
                print(f"\n📈 Virial Tensor Benefits:")
                print(f"   • Improved accuracy for pressure/stress predictions")
                print(f"   • Better performance in NPT ensemble simulations")
                print(f"   • Enhanced mechanical property predictions")
                print(f"   • More robust model extrapolation to different pressure conditions")
                print(f"   • Support for accurate bulk modulus and elastic constant calculations")
            else:
                print(f"\n⚠️  Missing Virial Tensor:")
                print(f"   • Model may have reduced accuracy for stress/pressure")
                print(f"   • NPT simulations may be less reliable")
                print(f"   • Consider ensuring CP2K outputs STRESS information")
                print(f"   • Check if virial.dat file contains stress tensor data")
            
            # Final recommendations
            print(f"\n💡 Recommendations:")
            if forces is not None and energies is not None and cells is not None:
                if virial is not None:
                    print(f"   🎯 Excellent! All required data present for high-quality training")
                    print(f"   🚀 Your model will support both NVT and NPT simulations")
                else:
                    print(f"   ✅ Good dataset for NVT simulations")
                    print(f"   📝 Consider adding virial data for NPT support")
            else:
                missing = []
                if forces is None:
                    missing.append("forces")
                if energies is None:
                    missing.append("energies")
                if cells is None:
                    missing.append("cell vectors")
                print(f"   ⚠️  Missing critical data: {', '.join(missing)}")
                print(f"   📝 Ensure CP2K is configured to output all required properties")
            
        else:
            raise RuntimeError("Failed to create DeepMD dataset")

    except Exception as exc:
        print(f"\n=== ERROR ===")
        print("Failed to generate DeepMD dataset.")
        print(f"Error: {exc}")
        traceback.print_exc()
        
        print(f"\nTroubleshooting suggestions:")
        print("1. Check that MD simulation completed successfully")
        print("2. Verify that trajectory files are not corrupted")
        print("3. Check that forces, energy, and virial files are properly formatted")
        print("4. Consider running with a smaller subset of frames for testing")
        print("5. Check the debug output above for specific parsing errors")
        print("6. Ensure virial.dat file contains stress tensor data in the expected format")
        print("7. Verify file permissions and accessibility")
        
        sys.exit(1)


if __name__ == "__main__":
    main()