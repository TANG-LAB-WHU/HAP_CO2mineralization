#!/usr/bin/env python3
"""Build DeepMD multi-system input JSON from pooled datasets."""

from __future__ import annotations

import argparse
import json
import math
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


def split_systems_by_ratio(
    systems: list[str], train_ratio: float, valid_ratio: float, test_ratio: float
) -> tuple[list[str], list[str], list[str]]:
    total = len(systems)
    if total == 0:
        return [], [], []
    if total == 1:
        return systems[:], systems[:], systems[:]
    if total == 2:
        return [systems[0]], [systems[1]], [systems[1]]

    train_count = int(math.floor(total * train_ratio))
    valid_count = int(math.floor(total * valid_ratio))
    test_count = int(math.floor(total * test_ratio))

    # Ensure non-empty split for >=3 systems.
    train_count = max(1, train_count)
    valid_count = max(1, valid_count)
    test_count = max(1, test_count)

    assigned = train_count + valid_count + test_count
    if assigned > total:
        overflow = assigned - total
        reduce_order = ["test", "valid", "train"]
        for key in reduce_order:
            if overflow <= 0:
                break
            if key == "test" and test_count > 1:
                delta = min(overflow, test_count - 1)
                test_count -= delta
                overflow -= delta
            elif key == "valid" and valid_count > 1:
                delta = min(overflow, valid_count - 1)
                valid_count -= delta
                overflow -= delta
            elif key == "train" and train_count > 1:
                delta = min(overflow, train_count - 1)
                train_count -= delta
                overflow -= delta
    elif assigned < total:
        train_count += total - assigned

    train_end = train_count
    valid_end = train_end + valid_count
    train_systems = systems[:train_end]
    valid_systems = systems[train_end:valid_end]
    test_systems = systems[valid_end : valid_end + test_count]

    if not test_systems:
        test_systems = valid_systems[-1:]
    return train_systems, valid_systems, test_systems


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pool-dir", required=True, help="Local deepmd_pool directory")
    parser.add_argument("--output", required=True, help="Output JSON file path")
    parser.add_argument(
        "--template",
        default="Step3_mlip_deepmd/input_template.json",
        help="Template input JSON to copy model/training hyperparameters from",
    )
    parser.add_argument(
        "--container-prefix",
        default="/mnt/deepmd_pool",
        help="Container-visible pool path prefix",
    )
    parser.add_argument(
        "--min-systems",
        type=int,
        default=1,
        help="Minimum required pooled systems for training",
    )
    parser.add_argument(
        "--train-ratio",
        type=float,
        default=0.8,
        help="System-level ratio allocated to training_data",
    )
    parser.add_argument(
        "--valid-ratio",
        type=float,
        default=0.1,
        help="System-level ratio allocated to validation_data",
    )
    parser.add_argument(
        "--test-ratio",
        type=float,
        default=0.1,
        help="System-level ratio allocated to test split summary",
    )
    parser.add_argument(
        "--test-output",
        default="",
        help="Optional path to write test systems JSON list",
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

    ratio_sum = args.train_ratio + args.valid_ratio + args.test_ratio
    if ratio_sum <= 0:
        raise ValueError("Sum of --train-ratio/--valid-ratio/--test-ratio must be > 0")
    train_ratio = args.train_ratio / ratio_sum
    valid_ratio = args.valid_ratio / ratio_sum
    test_ratio = args.test_ratio / ratio_sum

    container_systems = [f"{args.container_prefix}/{sys_dir.name}" for sys_dir in systems]
    train_systems, valid_systems, test_systems = split_systems_by_ratio(
        container_systems, train_ratio, valid_ratio, test_ratio
    )

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

    if args.test_output:
        test_output_path = Path(args.test_output)
        test_output_path.parent.mkdir(parents=True, exist_ok=True)
        with test_output_path.open("w", encoding="utf-8") as f:
            json.dump({"test_systems": test_systems}, f, indent=4)
            f.write("\n")

    print(f"Discovered pooled systems: {len(container_systems)}")
    print(f"Training systems: {len(train_systems)}")
    print(f"Validation systems: {len(valid_systems)}")
    print(f"Test systems: {len(test_systems)}")
    print(f"Wrote multi-system input: {output_path}")
    if args.test_output:
        print(f"Wrote test split list: {args.test_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
