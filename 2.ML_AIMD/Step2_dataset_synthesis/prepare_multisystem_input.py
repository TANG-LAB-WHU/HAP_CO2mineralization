#!/usr/bin/env python3
"""Build DeepMD multi-system input JSON from pooled datasets."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def discover_systems(pool_dir: Path) -> list[Path]:
    systems: list[Path] = []
    if not pool_dir.exists():
        return systems
    for child in sorted(pool_dir.iterdir()):
        if not child.is_dir():
            continue
        if (child / "set.000").is_dir() and (child / "type.raw").is_file():
            systems.append(child)
    return systems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pool-dir", required=True, help="Local deepmd_pool directory")
    parser.add_argument("--output", required=True, help="Output JSON file path")
    parser.add_argument(
        "--template",
        default="Step3_mlip_deepmd/input.json",
        help="Template input JSON to copy model/training hyperparameters from",
    )
    parser.add_argument(
        "--container-prefix",
        default="/mnt/shared/workspace/deepmd_pool",
        help="Container-visible pool path prefix",
    )
    parser.add_argument(
        "--min-systems",
        type=int,
        default=1,
        help="Minimum required pooled systems for training",
    )
    args = parser.parse_args()

    pool_dir = Path(args.pool_dir)
    template_path = Path(args.template)
    output_path = Path(args.output)

    if not template_path.is_file():
        raise FileNotFoundError(f"Template file not found: {template_path}")

    systems = discover_systems(pool_dir)
    if len(systems) < args.min_systems:
        raise RuntimeError(
            f"Need at least {args.min_systems} systems in {pool_dir}, found {len(systems)}"
        )

    with template_path.open("r", encoding="utf-8") as f:
        conf = json.load(f)

    container_systems = [f"{args.container_prefix}/{sys_dir.name}" for sys_dir in systems]
    if len(container_systems) == 1:
        train_systems = container_systems
        valid_systems = container_systems
    else:
        train_systems = container_systems[:-1]
        valid_systems = [container_systems[-1]]

    training = conf.setdefault("training", {})
    training_data = training.setdefault("training_data", {})
    validation_data = training.setdefault("validation_data", {})

    training_data["systems"] = train_systems
    training_data["batch_size"] = training_data.get("batch_size", "auto")
    training_data["auto_prob"] = training_data.get("auto_prob", "prob_sys_size")

    validation_data["systems"] = valid_systems
    validation_data["batch_size"] = validation_data.get("batch_size", "auto")
    validation_data["numb_btch"] = validation_data.get("numb_btch", 3)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        json.dump(conf, f, indent=4)
        f.write("\n")

    print(f"Discovered pooled systems: {len(container_systems)}")
    print(f"Training systems: {len(train_systems)}")
    print(f"Validation systems: {len(valid_systems)}")
    print(f"Wrote multi-system input: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
