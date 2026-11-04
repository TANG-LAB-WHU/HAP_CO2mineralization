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
    # Without expandable segments two DeePMD trainings on iteration-2 data filled the 24 GB GPU (13.8 + 9.0 GB) and
    # slowed to 0.8 s/batch or stalled; with them they used 5.8 + 4.1 GB and ran at 0.25 s/batch (2026-11-03)
    cmd = ["docker", "run", "--gpus", "all", "-e", "OMP_NUM_THREADS=2", "-e", "PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True",
           "-v", f"{ROOT.as_posix()}:/work", "-w", "/work"]
    cmd += ["-d", "--name", name] if detach else ["--rm"]
    return sh(cmd + [MLIP_IMAGE] + args_list, env={"MSYS_NO_PATHCONV": "1", **os.environ})


def in_container(path: Path) -> str:
    return "/work/" + path.resolve().relative_to(ROOT).as_posix()


def seed_containers(*statuses, name="seed-"):
    """Names of seed-label containers (started by hand, e.g. seed-s3) in the given states."""
    cmd = ["docker", "ps", "-a", "--format", "{{.Names}}", "--filter", f"name=^{name}"]
    for s in statuses:
        cmd += ["--filter", f"status={s}"]
    return subprocess.run(cmd, capture_output=True, text=True).stdout.split()


def label(run: Path, label_dir: Path, frames, data_out: Path, args, valid_fraction=0.1):
    py = sys.executable
    if not (label_dir / "index.csv").exists():
        sh([py, ROOT / "workflow/cp2k_label.py", "prepare", "--frames", *frames, "--label-dir", label_dir,
            "--cutoff", args.cutoff, "--rel-cutoff", args.rel_cutoff])
    # Seed-label containers (seed-<n>) pause while the pipeline labels, so the two never share the CPU. A seed
    # already paused was left so by an orchestrator stopped during a label step; it resumes here afterwards too
    paused = seed_containers("paused")
    seeds = seed_containers("running")
    for s in seeds:
        sh(["docker", "pause", s])
    seeds += paused
    try:
        sh(["docker", "run", "--rm", "--shm-size=4g", "-v", f"{label_dir.resolve().as_posix()}:/data", CP2K_IMAGE,
            "bash", "/data/run_labels.sh", "/data", args.ranks, 1], env={"MSYS_NO_PATHCONV": "1", **os.environ})
    finally:
        for s in seeds:
            sh(["docker", "unpause", s])
    sh([py, ROOT / "workflow/cp2k_label.py", "collect", "--label-dir", label_dir, "--out", data_out,
        "--valid-fraction", valid_fraction])


def main():
    parser = argparse.ArgumentParser(description="Pipeline v2 orchestrator")
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--cutoff", type=int, required=True)
    parser.add_argument("--rel-cutoff", type=int, required=True)
    parser.add_argument("--ranks", type=int, default=20)
    parser.add_argument("--gpu-jobs", type=int, default=1, help="DeePMD models trained at once on the GPU")
    parser.add_argument("--explore-jobs", type=int, default=1, help="LAMMPS exploration runs at once on the GPU")
    parser.add_argument("--iters", type=int, default=5)
    parser.add_argument("--steps", type=int, default=300_000, help="DeePMD steps from scratch (iteration 0 and final)")
    parser.add_argument("--warm-steps", type=int, default=100_000,
                        help="steps for iterations >= 1, warm-started from the previous ensemble; 0 trains from scratch")
    parser.add_argument("--converged", type=float, default=0.99)
    parser.add_argument("--extra-systems", nargs="*", default=[],
                        help="exploration sets added mid-run (e.g. s3): <run>/exploration/<name>/selected is labelled "
                             "for training and <name>-test/selected for the sealed test before the next training step; "
                             "<name>/*_initial.extxyz then joins the LAMMPS exploration")
    args = parser.parse_args()
    run = args.run.resolve()
    state_path = run / "pipeline_state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else {"done": []}

    def mark(step, **info):
        state["done"].append(step)
        state.setdefault("info", {})[step] = info
        state_path.write_text(json.dumps(state, indent=2) + "\n")

    dp = ROOT / "3.ML_AIMD/Step3_mlip_deepmd"

    def extras():
        return [n for n in args.extra_systems if f"label-extra-{n}" in state["done"]]

    def label_extras(wait=False):
        for n in args.extra_systems:
            if f"label-extra-{n}" in state["done"]:
                continue
            # A seed-<n> container labels these frames on the CPU while the GPU trains; join the next training
            # step instead of blocking this one (the final ensemble waits for it). Outside label() nothing should
            # hold it paused
            for s in seed_containers("paused", name=f"seed-{n}$"):
                sh(["docker", "unpause", s])
            if seed_containers("running", name=f"seed-{n}$"):
                if not wait:
                    print(f"extra {n}: seed-{n} still labelling; deferred to the next training step", flush=True)
                    continue
                print(f"extra {n}: waiting for seed-{n}", flush=True)
                subprocess.run(["docker", "wait", f"seed-{n}"], capture_output=True)
            label(run, run / f"extra/{n}/labels", sorted((run / f"exploration/{n}/selected").glob("*_selected.extxyz")),
                  run / f"extra/{n}/data", args)
            label(run, run / f"extra/{n}-test/labels", sorted((run / f"exploration/{n}-test/selected").glob("*_selected.extxyz")),
                  run / f"extra/{n}-test/data", args, valid_fraction=0.0)
            mark(f"label-extra-{n}")

    def sealed_test():
        # Labelled on the CPU while the first ensemble trains on the GPU
        if "sealed_test" not in state["done"]:
            label(run, run / "sealed-test/labels", sorted((run / "exploration/w8-test/selected").glob("*_selected.extxyz")),
                  run / "sealed-test/data", args, valid_fraction=0.0)
            mark("sealed_test")
        return sorted((run / "sealed-test/data/train").glob("*__*")) + sorted(
            p for n in extras() for p in (run / f"extra/{n}-test/data/train").glob("*__*"))

    def train_ensemble(tag, out, last_iter, steps, init_from=None, during=None):
        dirs = [run / f"iter-{k:03d}/data" for k in range(last_iter + 1)] + [run / f"extra/{n}/data" for n in extras()]
        train = sorted(p for d in dirs for p in (d / "train").glob("*__*"))
        valid = sorted(p for d in dirs for p in (d / "valid").glob("*__*"))
        sh([sys.executable, dp / "dp_ensemble.py", "write", "--train", *map(in_container, train),
            "--valid", *map(in_container, valid), "--steps", steps, "--out", out,
            *(["--sel-from", init_from] if init_from else [])])
        extra = ["--init-from", in_container(init_from)] if init_from else []
        # Concurrent GPU containers coincided with NVIDIA FECS errors and Docker VM crashes on
        # 2026-10-29, so models train --gpu-jobs at a time
        for start in range(0, 4, args.gpu_jobs):
            names = [f"dp-{tag}-{k}" for k in range(start, min(4, start + args.gpu_jobs))]
            for name in names:
                # After an orchestrator restart a model still training keeps going; it is waited for below
                st = subprocess.run(["docker", "inspect", "-f", "{{.State.Running}}", name], capture_output=True, text=True)
                if st.stdout.strip() == "true":
                    print(f"{name} still running; waiting for it", flush=True)
                    continue
                subprocess.run(["docker", "rm", "-f", name], capture_output=True)
                docker_mlip(["python", in_container(dp / "dp_ensemble.py"), "train", "--out", in_container(out),
                             "--index", name.rsplit("-", 1)[1], *extra], name=name, detach=True)
            if during and start == 0:
                during()
            for name in names:
                code = int(subprocess.run(["docker", "wait", name], capture_output=True, text=True).stdout.strip() or 1)
                if code != 0:
                    raise SystemExit(f"{name} exited with {code}; see docker logs {name}")
        mark(tag, n_train_systems=len(train), steps=steps, init_from=str(init_from) if init_from else None)

    def test_ensemble(out):
        docker_mlip(["python", in_container(dp / "dp_ensemble.py"), "test", "--out", in_container(out),
                     "--test", *map(in_container, sealed_test())])

    for it in range(args.iters):
        d = run / f"iter-{it:03d}"
        frames = (sorted((run / "exploration/w8/selected").glob("*_selected.extxyz")) if it == 0
                  else [run / f"iter-{it - 1:03d}/explore/candidates.extxyz"])
        if f"label-{it}" not in state["done"]:
            label(run, d / "labels", frames, d / "data", args)
            mark(f"label-{it}")
        if f"train-{it}" not in state["done"]:
            label_extras()
            # Intermediate iterations warm-start from the previous ensemble; the final ensemble is
            # retrained from scratch on all labels after the loop
            warm = it > 0 and args.warm_steps > 0
            train_ensemble(f"train-{it}", d / "train", it, args.warm_steps if warm else args.steps,
                           init_from=run / f"iter-{it - 1:03d}/train" if warm else None,
                           during=sealed_test if it == 0 else None)
        if f"test-{it}" not in state["done"]:
            test_ensemble(d / "train")
            mark(f"test-{it}")
        if f"explore-{it}" not in state["done"]:
            models = [in_container(d / f"train/model-{k}/graph-compress.pth") for k in range(4)]
            inits = sorted((run / "exploration/w8").glob("*_initial.extxyz")) + sorted(
                p for n in extras() for p in (run / f"exploration/{n}").glob("*_initial.extxyz"))
            # Converged only once every extra system was part of the exploration
            pending = [n for n in args.extra_systems if n not in extras()]
            docker_mlip(["python", in_container(dp / "dp_explore.py"), "--models", *models,
                         "--structures", *map(in_container, inits), "--temperatures", "300", "500",
                         "--steps", str(20000 * (it + 1)), "--parallel", str(args.explore_jobs),
                         "--out", in_container(d / "explore")])
            summary = json.loads((d / "explore/explore_summary.json").read_text())
            mark(f"explore-{it}", mean_accurate=summary["mean_accurate"], n_candidates=summary["n_candidates"])
            if summary["mean_accurate"] >= args.converged:
                if not pending:
                    print(f"converged at iteration {it}")
                    break
                print(f"iteration {it}: accurate, but {pending} not explored yet; continuing")
            if summary["n_candidates"] == 0:
                raise SystemExit(f"iteration {it}: not converged and no candidates; inspect {d / 'explore'}")

    # Final ensemble: from scratch, full step count, every label collected so far
    label_extras(wait=True)
    last = max(k for k in range(args.iters) if f"label-{k}" in state["done"])
    if "train-final" not in state["done"]:
        train_ensemble("train-final", run / "final/train", last, args.steps)
    if "test-final" not in state["done"]:
        test_ensemble(run / "final/train")
        mark("test-final")

    # Stage 7: lattice, stability and NVE-drift gates of the final ensemble (dp_validate.py). Lattice references
    # are the CP2K cell optimisations of the matching DFT-hull structures. The w8 systems carry the formal gates
    # (500 ps); the extra systems are checked for 200 ps
    groups = [("w8", sorted((run / "exploration/w8").glob("*_initial.extxyz")), 500)] + [
        (n, sorted((run / f"exploration/{n}").glob("*_initial.extxyz")), 200) for n in extras()]
    models = [in_container(run / f"final/train/model-{k}/graph-compress.pth") for k in range(4)]
    for tag, structures, ps in groups:
        if f"validate-{tag}" in state["done"]:
            continue
        v = run / f"final/validate/{tag}"
        v.mkdir(parents=True, exist_ok=True)
        ref = []
        if (run / "dft-hull/index.json").exists():
            sh([sys.executable, ROOT / "0.Generative_Discovery/screening/cp2k_hull.py", "cells", "--out", run / "dft-hull",
                "--match", *structures, "--write", v / "reference_cells.extxyz"])
            ref = ["--reference", in_container(v / "reference_cells.extxyz")]
        docker_mlip(["python", in_container(dp / "dp_validate.py"), "--models", *models,
                     "--structures", *map(in_container, structures), *ref, "--md-ps", str(ps),
                     "--parallel", str(args.explore_jobs), "--out", in_container(v)])
        mark(f"validate-{tag}", **json.loads((v / "validation_summary.json").read_text())["formal"])


if __name__ == "__main__":
    main()
