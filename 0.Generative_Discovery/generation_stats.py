#!/usr/bin/env python3
"""
Generate summary statistics for MatterGen outputs.
"""

import sys
from pathlib import Path
import argparse

def print_stats(structures_path):
    try:
        from ase.io import read
        from collections import Counter
        
        frames = read(structures_path, index=':')
        print(f"Total structures generated: {len(frames)}")
        
        if len(frames) == 0:
            return
            
        elements = Counter()
        for f in frames:
            elements.update(f.get_chemical_symbols())
        print(f"\nElement distribution across all structures:")
        for elem, count in elements.most_common():
            print(f"  {elem}: {count}")
            
        sizes = [len(f) for f in frames]
        print(f"\nAtoms per structure:")
        print(f"  Min: {min(sizes)}")
        print(f"  Max: {max(sizes)}")
        print(f"  Avg: {sum(sizes)/len(sizes):.1f}")
        
    except ImportError:
        print("Error: 'ase' is required to run this script. Please install it.")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(
        description='Generate summary statistics for MatterGen outputs'
    )
    parser.add_argument(
        'structures_path',
        type=str,
        help='Path to the .extxyz file containing generated structures'
    )
    args = parser.parse_args()

    path = Path(args.structures_path)
    if not path.exists():
        print(f"Error: {path} not found.")
        sys.exit(1)

    print_stats(path)

if __name__ == "__main__":
    main()
