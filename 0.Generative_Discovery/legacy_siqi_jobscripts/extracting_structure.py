#!/usr/bin/env python3

import argparse
import sys
import re
from datetime import datetime
from pathlib import Path
import ase.io
from mattergen.common.utils.eval_utils import load_structures, save_structures
from mattergen.evaluation.metrics.evaluator import MetricsEvaluator
from mattergen.evaluation.reference.reference_dataset_serializer import LMDBGZSerializer
from pymatgen.io.ase import AseAtomsAdaptor

def sanitize_filename(formula: str) -> str:
    """Sanitize chemical formula for use as filename."""
    # Remove spaces and replace special characters with underscores
    # Keep alphanumeric characters, parentheses, and common formula characters
    sanitized = re.sub(r'[^\w\-\(\)]', '_', formula)
    # Remove multiple consecutive underscores
    sanitized = re.sub(r'_+', '_', sanitized)
    # Remove leading/trailing underscores
    sanitized = sanitized.strip('_')
    return sanitized

def save_structures_as_individual_cifs(output_dir: Path, structures: list) -> int:
    """Save each structure as an individual CIF file with chemical formula as filename.
    
    Args:
        output_dir: Directory to save CIF files
        structures: List of pymatgen Structure objects
        
    Returns:
        Number of structures saved
    """
    if len(structures) == 0:
        return 0
    
    formula_count = {}
    saved_count = 0
    
    for structure in structures:
        # Get chemical formula from structure
        formula = structure.composition.reduced_formula
        
        # Sanitize formula for filename
        base_filename = sanitize_filename(formula)
        
        # Handle duplicate formulas by adding index
        if base_filename in formula_count:
            formula_count[base_filename] += 1
            filename = f"{base_filename}_{formula_count[base_filename]}.cif"
        else:
            formula_count[base_filename] = 0
            filename = f"{base_filename}.cif"
        
        # Convert Structure to ASE Atoms
        ase_atoms = AseAtomsAdaptor.get_atoms(structure)
        
        # Save as CIF file
        cif_path = output_dir / filename
        ase.io.write(str(cif_path), ase_atoms, format="cif")
        saved_count += 1
    
    return saved_count

def save_structures_as_extxyz(output_dir: Path, structures: list, filename: str = "structures.extxyz") -> bool:
    """Save all structures as a single extxyz file.
    
    Args:
        output_dir: Directory to save extxyz file
        structures: List of pymatgen Structure objects
        filename: Name of the output extxyz file
        
    Returns:
        True if saved successfully, False otherwise
    """
    if len(structures) == 0:
        return False
    
    # Convert all structures to ASE Atoms
    ase_atoms_list = [AseAtomsAdaptor.get_atoms(s) for s in structures]
    
    # Save as extxyz file
    extxyz_path = output_dir / filename
    ase.io.write(str(extxyz_path), ase_atoms_list, format="extxyz")
    
    return True

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Extract and classify structures based on stability, uniqueness, and novelty'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Output directory for extracted structures (default: script_dir/output)'
    )
    parser.add_argument(
        '--structures-file',
        type=str,
        default=None,
        help='Path to relaxed structures file (default: script_dir/relaxed_structures.extxyz)'
    )
    args = parser.parse_args()

    # Get script directory for relative path resolution
    script_dir = Path(__file__).parent.resolve()

    # Set output directory
    if args.output_dir:
        output_base = Path(args.output_dir).resolve()
    else:
        output_base = script_dir / "output"
    output_base.mkdir(parents=True, exist_ok=True)

    # Setup logging: create log file with timestamp in output directory
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = output_base / f"extraction_log_{timestamp}.txt"

    # Custom class to write to both console and log file
    class TeeOutput:
        def __init__(self, log_file_path):
            self.log_file = open(log_file_path, 'w', encoding='utf-8')
            self.terminal = sys.stdout
            
        def write(self, message):
            self.terminal.write(message)
            self.log_file.write(message)
            self.log_file.flush()
            
        def flush(self):
            self.terminal.flush()
            self.log_file.flush()
            
        def close(self):
            if self.log_file:
                self.log_file.close()

    # Redirect stdout to both console and file
    tee = TeeOutput(log_file)
    original_stdout = sys.stdout
    sys.stdout = tee

    # Print header with timestamp
    print("=" * 60)
    print(f"Structure Extraction Log - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    print(f"Log file: {log_file}")
    print(f"Output directory: {output_base}")
    print()

    # Set structures file path
    if args.structures_file:
        structures_path = Path(args.structures_file).resolve()
    else:
        structures_path = script_dir / "relaxed_structures.extxyz"

    # 1. Load relaxed structure file
    try:
        if not structures_path.exists():
            print(f"ERROR: Structure file not found: {structures_path}")
            sys.stdout = original_stdout
            tee.close()
            raise FileNotFoundError(f"Structure file not found: {structures_path}")

        structures = load_structures(structures_path)

        if len(structures) == 0:
            print("ERROR: No structures loaded from file")
            sys.stdout = original_stdout
            tee.close()
            raise ValueError("No structures loaded from file")

        print(f"Loaded {len(structures)} structures from {structures_path}")

        # 2. Extract energy information from extxyz file
        # ASE stores energy in atoms.info['total_energy']
        print("Extracting energy information...")
        ase_atoms = ase.io.read(str(structures_path), ":")
        energies = []
        for i, atoms in enumerate(ase_atoms):
            if 'total_energy' not in atoms.info:
                print(f"ERROR: Missing 'total_energy' in structure {i}")
                sys.stdout = original_stdout
                tee.close()
                raise KeyError(f"Missing 'total_energy' in structure {i}")
            energies.append(atoms.info['total_energy'])

        if len(energies) != len(structures):
            print(f"ERROR: Energy count ({len(energies)}) doesn't match structure count ({len(structures)})")
            sys.stdout = original_stdout
            tee.close()
            raise ValueError(f"Energy count ({len(energies)}) doesn't match structure count ({len(structures)})")
        
        print(f"Extracted {len(energies)} energy values")

        # 3. Load reference dataset
        print("\nLoading reference dataset...")
        mattergen_root = script_dir.parent.parent.parent / "mattergen"
        reference_path = mattergen_root / "data-release" / "alex-mp" / "reference_MP2020correction.gz"

        if not reference_path.exists():
            print(f"ERROR: Reference dataset not found: {reference_path}")
            sys.stdout = original_stdout
            tee.close()
            raise FileNotFoundError(f"Reference dataset not found: {reference_path}")

        print(f"Reference dataset path: {reference_path}")
        reference = LMDBGZSerializer().deserialize(str(reference_path))
        print("Reference dataset loaded successfully")

        # 4. Create evaluator
        print("\nCreating evaluator...")
        evaluator = MetricsEvaluator.from_structures_and_energies(
            structures=structures,
            energies=energies,
            reference=reference,
        )
        print("Evaluator created successfully")

        # 5. Get classification masks
        print("\nComputing classification masks...")
        is_stable = evaluator.is_stable
        is_unique = evaluator.is_unique
        is_novel = evaluator.is_novel
        print("Classification masks computed")

        # 6. Extract classified structures
        print("\nExtracting classified structures...")
        stable_structures = [s for i, s in enumerate(structures) if is_stable[i]]
        unique_structures = [s for i, s in enumerate(structures) if is_unique[i]]
        novel_structures = [s for i, s in enumerate(structures) if is_novel[i]]

        # S.U.N. structures (satisfying all three conditions: stable, unique, novel)
        sun_mask = is_stable & is_unique & is_novel
        sun_structures = [s for i, s in enumerate(structures) if sun_mask[i]]

        # Novel & Unique structures (satisfying novel and unique, but not necessarily stable)
        novel_unique_mask = is_novel & is_unique
        novel_unique_structures = [s for i, s in enumerate(structures) if novel_unique_mask[i]]
        
        print(f"Extracted: {len(stable_structures)} stable, {len(unique_structures)} unique, "
              f"{len(novel_structures)} novel, {len(novel_unique_structures)} novel_unique, "
              f"{len(sun_structures)} S.U.N. structures")

        # 7. Save results - create directories first
        structure_categories = {
            "stable": stable_structures,
            "unique": unique_structures,
            "novel": novel_structures,
            "novel_unique": novel_unique_structures,
            "sun": sun_structures,
        }

        print(f"\nSaving results to: {output_base}")
        for category, structs in structure_categories.items():
            output_dir = output_base / category
            output_dir.mkdir(parents=True, exist_ok=True)
            if len(structs) > 0:
                # Save individual CIF files
                saved_count = save_structures_as_individual_cifs(output_dir, structs)
                print(f"  Saved {saved_count} {category} structures as individual CIF files to {output_dir}")
                
                # Save all structures as a single extxyz file
                if save_structures_as_extxyz(output_dir, structs, f"{category}_structures.extxyz"):
                    print(f"  Saved {len(structs)} {category} structures as {category}_structures.extxyz to {output_dir}")
            else:
                print(f"  No {category} structures to save (skipping)")

        # 8. Print summary
        total = len(structures)
        print("\n" + "=" * 60)
        print("SUMMARY")
        print("=" * 60)
        print(f"Total structures: {total}")
        
        # Calculate percentages
        stable_pct = (len(stable_structures) / total * 100) if total > 0 else 0.0
        unique_pct = (len(unique_structures) / total * 100) if total > 0 else 0.0
        novel_pct = (len(novel_structures) / total * 100) if total > 0 else 0.0
        novel_unique_pct = (len(novel_unique_structures) / total * 100) if total > 0 else 0.0
        sun_pct = (len(sun_structures) / total * 100) if total > 0 else 0.0
        
        print(f"Stable: {len(stable_structures)}/{total} ({stable_pct:.2f}%)")
        print(f"Unique: {len(unique_structures)}/{total} ({unique_pct:.2f}%)")
        print(f"Novel: {len(novel_structures)}/{total} ({novel_pct:.2f}%)")
        print(f"Novel & Unique: {len(novel_unique_structures)}/{total} ({novel_unique_pct:.2f}%)")
        print(f"S.U.N. (Stable & Unique & Novel): {len(sun_structures)}/{total} ({sun_pct:.2f}%)")
        print("=" * 60)
        print(f"\nLog file saved to: {log_file}")
        print(f"Completed at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
    except Exception as e:
        # Restore stdout before raising exception
        sys.stdout = original_stdout
        tee.close()
        raise
    
    # Restore stdout and close log file
    sys.stdout = original_stdout
    tee.close()

if __name__ == '__main__':
    main()