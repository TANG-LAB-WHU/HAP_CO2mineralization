"""Pipeline v2 orchestrator on Shawn-PC: CP2K labels on the CPU, DeePMD + LAMMPS on the GPU.

Each iteration: label frames -> collect -> train 4 models in parallel -> sealed test ->
LAMMPS model-deviation exploration -> next frames. Progress is kept in <run>/pipeline_state.json
so the script resumes after an interruption. It stops when exploration is converged
(accurate fraction >= --converged) or after --iters iterations; gates are reported,
never relaxed.

  python workflow/pipeline_v2.py --run .research-runs/pipeline-v2 --cutoff 600 --rel-cutoff 60
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CP2K_IMAGE = "cp2k/cp2k:2026.2_mpich_generic_psmp"
MLIP_IMAGE = "hap-mlip:2026.10"


def sh(cmd, **kw):
    print("+", " ".join(map(str, cmd)), flush=True)
    return subprocess.run(list(map(str, cmd)), check=True, **kw)


def docker_mlip(args_list, name=None, detach=False):
    cmd = ["docker", "run", "--gpus", "all", "-e", "OMP_NUM_THREADS=2", "-v", f"{ROOT.as_posix()}:/work", "-w", "/work"]
    cmd += ["-d", "--name", name] if detach else ["--rm"]
    return sh(cmd + [MLIP_IMAGE] + args_list, env={"MSYS_NO_PATHCONV": "1", **os.environ})


def in_container(path: Path) -> str:
    return "/work/" + path.resolve().relative_to(ROOT).as_posix()


def label(run: Path, label_dir: Path, frames, data_out: Path, args, valid_fraction=0.1):
    py = sys.executable
    if not (label_dir / "index.csv").exists():
        sh([py, ROOT / "workflow/cp2k_label.py", "prepare", "--frames", *frames, "--label-dir", label_dir,
            "--cutoff", args.cutoff, "--rel-cutoff", args.rel_cutoff])
    sh(["docker", "run", "--rm", "--shm-size=4g", "-v", f"{label_dir.resolve().as_posix()}:/data", CP2K_IMAGE,
        "bash", "/data/run_labels.sh", "/data", args.ranks, 1], env={"MSYS_NO_PATHCONV": "1", **os.environ})
    sh([py, ROOT / "workflow/cp2k_label.py", "collect", "--label-dir", label_dir, "--out", data_out,
        "--valid-fraction", valid_fraction])


def main():
    parser = argparse.ArgumentParser(description="Pipeline v2 orchestrator")
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--cutoff", type=int, required=True)
    parser.add_argument("--rel-cutoff", type=int, required=True)
    parser.add_argument("--ranks", type=int, default=24)
    parser.add_argument("--iters", type=int, default=5)
    parser.add_argument("--steps", type=int, default=300_000, help="DeePMD training steps per iteration")
    parser.add_argument("--converged", type=float, default=0.99)
    args = parser.parse_args()
    run = args.run.resolve()
    state_path = run / "pipeline_state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else {"done": []}

    def mark(step, **info):
        state["done"].append(step)
        state.setdefault("info", {})[step] = info
        state_path.write_text(json.dumps(state, indent=2) + "\n")

    dp = ROOT / "3.ML_AIMD/Step3_mlip_deepmd"

    def sealed_test():
        # Labelled on the CPU while the first ensemble trains on the GPU
        if "sealed_test" not in state["done"]:
            label(run, run / "sealed-test/labels", sorted((run / "exploration/w8-test/selected").glob("*_selected.extxyz")),
                  run / "sealed-test/data", args, valid_fraction=0.0)
            mark("sealed_test")
        return sorted((run / "sealed-test/data/train").glob("*__*"))

    for it in range(args.iters):
        d = run / f"iter-{it:03d}"
        frames = (sorted((run / "exploration/w8/selected").glob("*_selected.extxyz")) if it == 0
                  else [run / f"iter-{it - 1:03d}/explore/candidates.extxyz"])
        if f"label-{it}" not in state["done"]:
            label(run, d / "labels", frames, d / "data", args)
            mark(f"label-{it}")
        if f"train-{it}" not in state["done"]:
            train = sorted(p for k in range(it + 1) for p in (run / f"iter-{k:03d}/data/train").glob("*__*"))
            valid = sorted(p for k in range(it + 1) for p in (run / f"iter-{k:03d}/data/valid").glob("*__*"))
            sh([sys.executable, dp / "dp_ensemble.py", "write", "--train", *map(in_container, train),
                "--valid", *map(in_container, valid), "--steps", args.steps, "--out", d / "train"])
            names = [f"dp-train-{it}-{k}" for k in range(4)]
            for k, name in enumerate(names):
                subprocess.run(["docker", "rm", "-f", name], capture_output=True)
                docker_mlip(["python", in_container(dp / "dp_ensemble.py"), "train", "--out", in_container(d / "train"),
                             "--index", str(k)], name=name, detach=True)
            if it == 0:
                sealed_test()
            for name in names:
                code = int(subprocess.run(["docker", "wait", name], capture_output=True, text=True).stdout.strip() or 1)
                if code != 0:
                    raise SystemExit(f"{name} exited with {code}; see docker logs {name}")
            mark(f"train-{it}", n_train_systems=len(train))
        test_systems = sealed_test()
        if f"test-{it}" not in state["done"] and test_systems:
            docker_mlip(["python", in_container(dp / "dp_ensemble.py"), "test", "--out", in_container(d / "train"),
                         "--test", *map(in_container, test_systems)])
            mark(f"test-{it}")
        if f"explore-{it}" not in state["done"]:
            models = [in_container(d / f"train/model-{k}/graph-compress.pth") for k in range(4)]
            inits = sorted((run / "exploration/w8").glob("*_initial.extxyz"))
            docker_mlip(["python", in_container(dp / "dp_explore.py"), "--models", *models,
                         "--structures", *map(in_container, inits), "--temperatures", "300", "500",
                         "--steps", str(20000 * (it + 1)), "--out", in_container(d / "explore")])
            summary = json.loads((d / "explore/explore_summary.json").read_text())
            mark(f"explore-{it}", mean_accurate=summary["mean_accurate"], n_candidates=summary["n_candidates"])
            if summary["mean_accurate"] >= args.converged:
                print(f"converged at iteration {it}")
                break


if __name__ == "__main__":
    main()
