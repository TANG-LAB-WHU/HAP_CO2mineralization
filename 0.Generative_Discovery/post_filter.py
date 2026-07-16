#!/usr/bin/env python3
"""
Filter MatterGen candidates by stability and novelty.
Input:  detailed_metrics.json + generated_crystals.extxyz
Output: filtered_candidates/ directory with CIF files
"""

import os
import json
import argparse
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Filter generated candidates")
    parser.add_argument("--results-dir", type=str, required=True, help="Directory containing MatterGen results")
    parser.add_argument("--ehull-threshold", type=float, default=0.1, help="Max energy above hull threshold (eV/atom)")
    parser.add_argument("--require-novelty", action="store_true", default=True, help="Require structures to be novel")
    args = parser.parse_args()

    results_dir = Path(args.results_dir)
    detailed_metrics_file = results_dir / "detailed_metrics.json"
    structures_file = results_dir / "generated_crystals.extxyz"
    output_dir = results_dir / "filtered_candidates"
    
    if not detailed_metrics_file.exists():
        print(f"Error: {detailed_metrics_file} not found.")
        return
        
    if not structures_file.exists():
        print(f"Error: {structures_file} not found.")
        return

    print(f"Loading metrics from {detailed_metrics_file}...")
    with open(detailed_metrics_file, "r") as f:
        metrics = json.load(f)
        
    try:
        from ase.io import read, write
        frames = read(structures_file, index=':')
    except ImportError:
        print("Error: 'ase' is required to process structures. Please install it.")
        return
        
    ehull_list = metrics.get("energy_above_hull", [])
    novelty_list = metrics.get("novelty", [])
    
    num_metrics = len(ehull_list) if ehull_list else 0
    if len(frames) != num_metrics:
        print(f"Warning: Number of structures ({len(frames)}) does not match metrics length ({num_metrics})")
        
    output_dir.mkdir(parents=True, exist_ok=True)
    
    saved_count = 0
    for idx, structure in enumerate(frames):
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
                
    print(f"Filtering complete.")
    print(f"Saved {saved_count} candidates to {output_dir}")

if __name__ == "__main__":
    main()
