#!/usr/bin/env python3
"""Extract the last frame from a multi-frame XYZ file.

Usage:
    python extract_last_frame.py input.xyz output.xyz
This writes the final configuration (number of atoms + comment + coordinates)
from *input.xyz* into *output.xyz*.
"""
import sys
from pathlib import Path

def extract_last_frame(src: Path, dst: Path):
    with src.open('r') as fh:
        lines = fh.readlines()

    idx = 0
    last_frame = None
    total_lines = len(lines)
    while idx < total_lines:
        try:
            natoms = int(lines[idx].strip())
        except ValueError:
            break  # Malformed file
        frame_lines = lines[idx: idx + natoms + 2]
        if len(frame_lines) < natoms + 2:
            break  # Incomplete frame at EOF
        last_frame = frame_lines
        idx += natoms + 2

    if last_frame is None:
        raise RuntimeError(f"No valid frame detected in {src}")

    # CP2K expects plain coordinate lines without the XYZ header (natoms + comment)
    # Therefore we drop the first two lines when writing the output file.
    # If you need the original XYZ frame with header for other purposes,
    # simply keep a copy before this transformation.

    # Remove natoms and comment line
    coordinate_lines = last_frame[2:]

    # Safety check to ensure the coordinate count matches `natoms`
    if len(coordinate_lines) != natoms:
        raise RuntimeError(
            f"Frame seems malformed: expected {natoms} coordinate lines, got {len(coordinate_lines)}"
        )

    dst.write_text(''.join(coordinate_lines))
    print(f"Wrote last frame (no XYZ header) with {natoms} atoms to {dst}")


def main():
    if len(sys.argv) != 3:
        print("Usage: python extract_last_frame.py <input.xyz> <output.xyz>")
        sys.exit(1)

    src = Path(sys.argv[1])
    dst = Path(sys.argv[2])
    if not src.is_file():
        print(f"Input file {src} not found")
        sys.exit(1)

    extract_last_frame(src, dst)

if __name__ == "__main__":
    main() 