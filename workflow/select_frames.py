"""Pick time-decorrelated frames from exploration trajectories for CP2K labelling.

Frames are taken evenly in time within every (system, temperature) block so each block
contributes the same number. The sealed test set must come from a separate exploration
run with a different seed; this script never mixes the two.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import numpy as np
from ase.io import read, write


def main():
    parser = argparse.ArgumentParser(description="Even-in-time frame selection")
    parser.add_argument("--frames", nargs="+", type=Path, required=True)
    parser.add_argument("--per-block", type=int, default=20)
    parser.add_argument("--skip-first", type=int, default=5, help="frames dropped at the start of each block")
    parser.add_argument("--strain", type=float, default=0.04, help="max diagonal strain of the bulk copies")
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    total = 0
    for path in args.frames:
        blocks = defaultdict(list)
        for f in read(path, ":"):
            blocks[f.info["temperature_K"]].append(f)
        picked = []
        rng = np.random.default_rng(args.seed)
        for t, frames in sorted(blocks.items()):
            frames = frames[args.skip_first:]
            idx = np.unique(np.linspace(0, len(frames) - 1, min(args.per_block, len(frames))).round().astype(int))
            picked += [frames[i] for i in idx]
            # Exploration MD is NVT; bulk cells also get strained copies so NPT has labels to learn from
            if not frames[0].info.get("system", "").startswith("hap001"):
                for i in idx[::3]:
                    s = frames[i].copy()
                    eps = np.diag(rng.uniform(-args.strain, args.strain, 3))
                    eps[np.triu_indices(3, 1)] = eps[np.tril_indices(3, -1)] = rng.uniform(-args.strain / 2, args.strain / 2, 3)
                    s.set_cell(s.cell.array @ (np.eye(3) + eps), scale_atoms=True)
                    s.info = {**s.info, "strain_max": float(np.abs(eps).max())}
                    picked.append(s)
        write(args.out / path.name.replace("_frames", "_selected"), picked)
        total += len(picked)
        print(path.name, len(picked))
    print("total", total)


if __name__ == "__main__":
    main()
