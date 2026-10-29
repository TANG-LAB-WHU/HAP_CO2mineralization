"""Pipeline v2 DeePMD ensemble: write inputs, train/freeze/compress 4 seeds, test on a sealed set.

Architecture follows the smooth Deep Potential used for the SCAN-ML CaCO3 model
(se_e2_a, rcut 6 A, embedding 25/50/100, fitting 3x120 there; 3x240 here for the wider
Ca-P-O-H-C chemistry). Run inside hap-mlip (DeePMD-kit 3.1.3, PyTorch backend, CUDA).

  write:  python dp_ensemble.py write --train D1 D2 --valid V1 --out iter-000
  train:  python dp_ensemble.py train --out iter-000          (sequential seeds on one GPU)
  test:   python dp_ensemble.py test --out iter-000 --test T1 T2
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path

TYPE_MAP = ["Ca", "P", "O", "H", "C"]
SEEDS = [11, 23, 37, 53]


def config(train, valid, seed, steps):
    return {
        "model": {
            "type_map": TYPE_MAP,
            "descriptor": {"type": "se_e2_a", "sel": "auto", "rcut_smth": 0.5, "rcut": 6.0,
                           "neuron": [25, 50, 100], "axis_neuron": 16, "resnet_dt": False,
                           "type_one_side": False, "seed": seed},
            "fitting_net": {"neuron": [240, 240, 240], "resnet_dt": True, "seed": seed},
        },
        "learning_rate": {"type": "exp", "start_lr": 1e-3, "stop_lr": 3e-8, "decay_steps": 5000},
        "loss": {"type": "ener", "start_pref_e": 0.02, "limit_pref_e": 1.0, "start_pref_f": 1000.0,
                 "limit_pref_f": 1.0, "start_pref_v": 0.0, "limit_pref_v": 0.02},
        "training": {
            "training_data": {"systems": [Path(p).as_posix() for p in train], "batch_size": "auto"},
            "validation_data": {"systems": [Path(p).as_posix() for p in valid], "batch_size": "auto", "numb_btch": 4},
            "numb_steps": steps, "seed": seed, "disp_file": "lcurve.out", "disp_freq": 1000,
            "save_freq": 20000,
        },
    }


def run(cmd, cwd):
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=cwd, check=True)


def main():
    parser = argparse.ArgumentParser(description="DeePMD ensemble driver")
    sub = parser.add_subparsers(dest="cmd", required=True)
    w = sub.add_parser("write")
    w.add_argument("--train", nargs="+", type=Path, required=True)
    w.add_argument("--valid", nargs="+", type=Path, required=True)
    w.add_argument("--steps", type=int, default=400_000)
    w.add_argument("--out", type=Path, required=True)
    t = sub.add_parser("train")
    t.add_argument("--out", type=Path, required=True)
    t.add_argument("--init-from", type=Path, default=None, help="previous iteration dir to restart weights from")
    t.add_argument("--index", type=int, nargs="*", default=None, help="train only these members (run members in parallel)")
    e = sub.add_parser("test")
    e.add_argument("--out", type=Path, required=True)
    e.add_argument("--test", nargs="+", type=Path, required=True)
    args = parser.parse_args()

    if args.cmd == "write":
        for k, seed in enumerate(SEEDS):
            d = args.out / f"model-{k}"
            d.mkdir(parents=True, exist_ok=True)
            (d / "input.json").write_text(json.dumps(config(args.train, args.valid, seed, args.steps), indent=2) + "\n")
        print(f"wrote {len(SEEDS)} inputs under {args.out}")
    elif args.cmd == "train":
        for k in (args.index if args.index is not None else range(len(SEEDS))):
            d = (args.out / f"model-{k}").resolve()
            if (d / "graph-compress.pth").exists():
                continue
            cmd = ["dp", "--pt", "train", "input.json"]
            if args.init_from:
                cmd += ["--init-model", str((args.init_from / f"model-{k}" / "model.ckpt.pt").resolve())]
            run(cmd, d)
            run(["dp", "--pt", "freeze", "-o", "graph.pth"], d)
            try:
                run(["dp", "--pt", "compress", "-i", "graph.pth", "-o", "graph-compress.pth"], d)
            except subprocess.CalledProcessError:
                # Compression only speeds up MD; record the fallback instead of failing the iteration
                (d / "graph-compress.pth").write_bytes((d / "graph.pth").read_bytes())
                (d / "COMPRESSION_FAILED").write_text("dp --pt compress failed; graph-compress.pth is the uncompressed model\n")
    else:
        report = {}
        for k in range(len(SEEDS)):
            d = (args.out / f"model-{k}").resolve()
            if not (d / "graph.pth").exists():
                report[f"model-{k}"] = "not trained"
                continue
            res = {}
            for sysdir in args.test:
                # dp test reports through logging (stderr), so read both streams
                proc = subprocess.run(["dp", "--pt", "test", "-m", "graph.pth", "-s", sysdir.resolve().as_posix(), "-n", "100000"],
                                      cwd=d, capture_output=True, text=True)
                text = proc.stdout + proc.stderr
                if proc.returncode != 0:
                    res[sysdir.name] = {"error": text.strip().splitlines()[-1] if text.strip() else "dp test failed"}
                    continue
                vals = dict(re.findall(r"(Energy RMSE/Natoms|Force  RMSE|Force RMSE)\s*:\s*([0-9.eE+-]+)", text))
                res[sysdir.name] = {"energy_rmse_eV_atom": float(vals.get("Energy RMSE/Natoms", "nan")),
                                    "force_rmse_eV_A": float(vals.get("Force  RMSE", vals.get("Force RMSE", "nan"))),
                                    "frames": int(re.search(r"number of test data\s*:\s*(\d+)", text).group(1))
                                    if re.search(r"number of test data\s*:\s*(\d+)", text) else None}
            report[f"model-{k}"] = res
        (args.out / "sealed_test.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
