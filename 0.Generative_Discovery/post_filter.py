#!/usr/bin/env python3
"""
Filter MatterGen candidates by stability and novelty.
Input:  detailed_metrics.json + generated_crystals.extxyz
Output: filtered_candidates/ directory with CIF files
"""

import os
import math
import json
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Filter generated candidates")
    parser.add_argument("--results-dir", type=str, required=True, help="Directory containing MatterGen results")
    parser.add_argument("--ehull-threshold", type=float, default=0.1, help="Max energy above hull threshold (eV/atom)")
    parser.add_argument("--require-novelty", action=argparse.BooleanOptionalAction, default=True, help="Require structures to be novel")
    parser.add_argument("--required-elements", nargs="+", default=["Ca", "P", "O", "H"], help="Necessary chemistry gate; default HAP elements, not a structure-type proof")
    parser.add_argument("--structures", type=Path, help="Exact evaluated structures; defaults to relaxed_structures.extxyz")
    args = parser.parse_args()

    if not math.isfinite(args.ehull_threshold) or args.ehull_threshold < 0:
        parser.error("Hull threshold must be finite and nonnegative")
    results_dir = Path(args.results_dir)
    detailed_metrics_file = results_dir / "detailed_metrics.json"
    structures_file = args.structures or results_dir / "relaxed_structures.extxyz"
    output_dir = results_dir / "filtered_candidates"
    
    if not detailed_metrics_file.exists():
        print(f"Error: {detailed_metrics_file} not found.")
        raise SystemExit(1)
        
    if not structures_file.exists():
        print(f"Error: {structures_file} not found.")
        raise SystemExit(1)

    print(f"Loading metrics from {detailed_metrics_file}...")
    with open(detailed_metrics_file, "r") as f:
        metrics = json.load(f)
        
    try:
        from ase.io import read, write
        frames = read(structures_file, index=':')
    except ImportError:
        print("Error: 'ase' is required to process structures. Please install it.")
        raise SystemExit(1)
        
    ehull_list = metrics.get("energy_above_hull", [])
    novelty_list = metrics.get("is_novel", [])
    
    num_metrics = len(ehull_list) if ehull_list else 0
    if not frames or len(frames) != num_metrics:
        raise ValueError("Nonempty structures and metrics must have identical lengths")
    if args.require_novelty and (len(novelty_list) != len(frames) or not all(type(x) is bool for x in novelty_list)):
        raise ValueError("One boolean novelty label per structure is required")
    if not all(isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x) for x in ehull_list):
        raise ValueError("Nonfinite or invalid hull energy")
        
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError("Preserving existing filtered candidates; use a fresh result directory")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    saved_count = 0
    chemistry_rejected = 0
    for idx, structure in enumerate(frames):
        if not set(args.required_elements).issubset(structure.get_chemical_symbols()):
            chemistry_rejected += 1
            continue
        # Extract metrics for this specific structure
        ehull = ehull_list[idx] if idx < len(ehull_list) else 100.0
        is_novel = novelty_list[idx] if idx < len(novelty_list) else False
        
        # Check criteria
        if ehull <= args.ehull_threshold:
            if not args.require_novelty or is_novel:
                formula = structure.get_chemical_formula()
                filename = f"candidate_{idx:04d}_{formula}_ehull{ehull:.3f}.cif"
                filepath = output_dir / filename
                
                # Save as CIF
                write(filepath, structure)
                saved_count += 1
                
    import hashlib
    summary = {"input_count":len(frames), "saved_count":saved_count, "chemistry_rejected":chemistry_rejected,
        "required_elements":args.required_elements, "manuscript_eligible":False,
        "structures_sha256":hashlib.sha256(structures_file.read_bytes()).hexdigest(),
        "metrics_sha256":hashlib.sha256(detailed_metrics_file.read_bytes()).hexdigest(),
        "note":"Necessary chemistry gate only; apatite structure identification and scientific validation remain required"}
    (results_dir/'filter_summary.json').write_text(json.dumps(summary,indent=2))
    print(f"Filtering complete.")
    print(f"Saved {saved_count} candidates to {output_dir}")

if __name__ == "__main__":
    main()
