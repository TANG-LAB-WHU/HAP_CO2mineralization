"""Pipeline v2 orchestrator on Shawn-PC: CP2K labels on the CPU, DeePMD + LAMMPS on the GPU.

Each iteration: label frames -> collect -> train 4 models in parallel -> sealed test ->
LAMMPS model-deviation exploration -> next frames. Progress is kept in <run>/pipeline_state.json
so the script resumes after an interruption. It stops when exploration is converged
(accurate fraction >= --converged) or after --iters iterations; gates are reported,
never relaxed.

  python workflow/pipeline_v2.py --run .research-runs/pipeline-v2 --cutoff 600 --rel-cutoff 60

On a Slurm cluster (--runtime apptainer) the same steps run in Apptainer images of the two Docker images, and
--phase cpu|gpu stops at the first step of the other phase, so CP2K labels and GPU steps go to separate jobs;
<run>/next_phase then names the phase to submit next (or "done"). See workflow/hpc/PIPELINE_V2_HPC.md.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[1]
CP2K_IMAGE = "cp2k/cp2k:2026.2_mpich_generic_psmp"
MLIP_IMAGE = "hap-mlip:2026.10"
# Set by main(): container runtime and images (Docker tags or Apptainer .sif files)
RUNTIME = {"runtime": "docker", "cp2k": CP2K_IMAGE, "mlip": MLIP_IMAGE}
# Host directory -> path inside the MLIP container; main() adds the run directory when it lies outside ROOT
MOUNTS = [(ROOT, "/work")]
# The CP2K image's entrypoint puts cp2k.psmp and mpirun on PATH; Docker runs it, Apptainer exec does not
CP2K_LAUNCH = "/opt/cp2k/install/bin/launch"
# Inside a Slurm job MPICH's mpirun starts its ranks with srun, which the image lacks ("execvp error on file srun",
# checked with the 2026.2 image under a Slurm environment on 2026-11-08); fork keeps them on the job's node
HYDRA_ENV = "HYDRA_LAUNCHER=fork,HYDRA_RMK=user"
# Without expandable segments two DeePMD trainings on iteration-2 data filled the 24 GB GPU (13.8 + 9.0 GB) and
# slowed to 0.8 s/batch or stalled; with them they used 5.8 + 4.1 GB and ran at 0.25 s/batch (2026-11-03)
MLIP_ENV = {"OMP_NUM_THREADS": "2", "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True"}


def sh(cmd, **kw):
    print("+", " ".join(map(str, cmd)), flush=True)
    return subprocess.run(list(map(str, cmd)), check=True, **kw)


def apptainer() -> bool:
    return RUNTIME["runtime"] == "apptainer"


def mlip(args_list, name=None, detach=False, gpu=None):
    """Run args_list in the MLIP image. detach starts it in the background: Docker as container `name`, Apptainer
    as a child process on GPU `gpu` (returned)."""
    if apptainer():
        cmd = ["apptainer", "exec", "--nv", "--pwd", "/work", "--env", ",".join(f"{k}={v}" for k, v in MLIP_ENV.items())]
        for host, mount in MOUNTS:
            cmd += ["--bind", f"{host.as_posix()}:{mount}"]
        cmd += [RUNTIME["mlip"], *args_list]
        if not detach:
            return sh(cmd)
        print("+", " ".join(map(str, cmd)), flush=True)
        return subprocess.Popen(list(map(str, cmd)), env={**os.environ, **({"CUDA_VISIBLE_DEVICES": gpu} if gpu else {})})
    cmd = ["docker", "run", "--gpus", "all"]
    for k, v in MLIP_ENV.items():
        cmd += ["-e", f"{k}={v}"]
    for host, mount in MOUNTS:
        cmd += ["-v", f"{host.as_posix()}:{mount}"]
    cmd += ["-w", "/work"] + (["-d", "--name", name] if detach else ["--rm"])
    return sh(cmd + [RUNTIME["mlip"]] + args_list, env={"MSYS_NO_PATHCONV": "1", **os.environ})


def mlip_start(args_list, name, gpu=None):
    """Start a background MLIP job; returns what mlip_wait() takes. After an orchestrator restart a Docker model
    still training keeps going and is waited for; Apptainer jobs end with the Slurm job that started them."""
    if apptainer():
        return mlip(args_list, detach=True, gpu=gpu)
    st = subprocess.run(["docker", "inspect", "-f", "{{.State.Running}}", name], capture_output=True, text=True)
    if st.stdout.strip() == "true":
        print(f"{name} still running; waiting for it", flush=True)
        return name
    subprocess.run(["docker", "rm", "-f", name], capture_output=True)
    mlip(args_list, name=name, detach=True)
    return name


def mlip_wait(job) -> int:
    if isinstance(job, subprocess.Popen):
        return job.wait()
    return int(subprocess.run(["docker", "wait", job], capture_output=True, text=True).stdout.strip() or 1)


def in_container(path: Path) -> str:
    p = path.resolve()
    for host, mount in MOUNTS:
        try:
            return str(PurePosixPath(mount) / p.relative_to(host).as_posix())
        except ValueError:
            continue
    raise SystemExit(f"{path} is outside the directories mounted in the MLIP container")


def seed_containers(*statuses, name="seed-"):
    """Names of seed-label containers (started by hand, e.g. seed-s3) in the given states; none under Apptainer."""
    if apptainer():
        return []
    cmd = ["docker", "ps", "-a", "--format", "{{.Names}}", "--filter", f"name=^{name}"]
    for s in statuses:
        cmd += ["--filter", f"status={s}"]
    return subprocess.run(cmd, capture_output=True, text=True).stdout.split()


@contextmanager
def seeds_paused():
    """Seed-label containers (seed-<n>) pause while the pipeline labels, so the two never share the CPU. A seed
    already paused was left so by an orchestrator stopped during a label step; it resumes here afterwards too."""
    paused = seed_containers("paused")
    seeds = seed_containers("running")
    for s in seeds:
        sh(["docker", "pause", s])
    try:
        yield
    finally:
        for s in seeds + paused:
            sh(["docker", "unpause", s])


def cp2k(label_dir: Path, args, script, *extra):
    run = ["bash", f"/data/{script}", "/data", args.ranks, 1, *extra]
    if apptainer():
        sh(["apptainer", "exec", "--env", HYDRA_ENV, "--bind", f"{label_dir.resolve().as_posix()}:/data",
            RUNTIME["cp2k"], CP2K_LAUNCH, *run])
        return
    sh(["docker", "run", "--rm", "--shm-size=4g", "-v", f"{label_dir.resolve().as_posix()}:/data",
        RUNTIME["cp2k"], *run], env={"MSYS_NO_PATHCONV": "1", **os.environ})


def screen(label_dir: Path, models, out="energy_screen.json"):
    mlip(["python", in_container(ROOT / "workflow/screen_labels.py"), "--label-dir", in_container(label_dir),
          "--models", *models, "--out", out])
    return json.loads((label_dir / out).read_text())["flagged"]


def label(run: Path, label_dir: Path, frames, data_out: Path, args, valid_fraction=0.1, screen_models=None):
    py = sys.executable
    if not (label_dir / "index.csv").exists():
        sh([py, ROOT / "workflow/cp2k_label.py", "prepare", "--frames", *frames, "--label-dir", label_dir,
            "--cutoff", args.cutoff, "--rel-cutoff", args.rel_cutoff])
    # A directory prepared earlier gets the current run_labels.sh (e.g. the shard arguments of 2026-11-08)
    sh([py, ROOT / "workflow/cp2k_label.py", "scripts", "--label-dir", label_dir])
    with seeds_paused():
        cp2k(label_dir, args, "run_labels.sh")
        if screen_models:
            # A converged SCF can still sit in a wrong electronic state (see screen_labels.py): flagged frames are
            # relabelled once from the atomic guess with OT CG, then screened again; collect leaves out the rest
            if screen(label_dir, screen_models):
                sh([py, ROOT / "workflow/cp2k_label.py", "relabel", "--label-dir", label_dir])
                cp2k(label_dir, args, "relabel.sh")
                screen(label_dir, screen_models)
    sh([py, ROOT / "workflow/cp2k_label.py", "collect", "--label-dir", label_dir, "--out", data_out,
        "--valid-fraction", valid_fraction, *(["--exclude-screen"] if screen_models else [])])


def recheck(label_dir: Path, data_out: Path, models, args, tag, valid_fraction=0.1):
    """Screen a collected training set again with newer models; relabel what they flag and collect it again.

    label() screens with the ensemble that explored the frames. Where that ensemble has not learned a system its
    models disagree by more than --max-spread and a wrong SCF state passes: iter-004 candidates_00304 (s3_055),
    converged by OT DIIS 0.325 Ha above the state OT CG reached, was found only by the iteration-4 models that had
    trained on it (2026-11-07). Acts only on frames the earlier screens did not already leave out; a frame is still
    relabelled once. The previous data is kept as <data>.pre-<tag>. Returns whether the set was collected again."""
    out = f"energy_screen.{tag}.json"
    earlier = {n for s in label_dir.glob("energy_screen.json") for n in json.loads(s.read_text())["flagged"]}
    if not set(screen(label_dir, models, out)) - earlier:
        return False
    py = sys.executable
    with seeds_paused():
        sh([py, ROOT / "workflow/cp2k_label.py", "relabel", "--label-dir", label_dir, "--screen", out,
            "--list", f"relabel.{tag}.txt"])
        cp2k(label_dir, args, "relabel.sh", f"relabel.{tag}.txt")
        screen(label_dir, models, out)
    keep = data_out.with_name(f"{data_out.name}.pre-{tag}")
    if keep.exists():  # an earlier attempt was interrupted after moving the data; its new collection is redone
        shutil.rmtree(data_out, ignore_errors=True)
    else:
        data_out.rename(keep)
    sh([py, ROOT / "workflow/cp2k_label.py", "collect", "--label-dir", label_dir, "--out", data_out,
        "--valid-fraction", valid_fraction, "--exclude-screen", "--screen", "energy_screen.json", out])
    return True


def step_phase(step: str) -> str:
    """CP2K labels, with the energy screens and relabels around them, form the cpu phase; DeePMD training and
    tests, LAMMPS exploration and validation the gpu phase."""
    return "cpu" if step.startswith(("label-", "recheck-", "sealed_test")) else "gpu"


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
    parser.add_argument("--npt-temperatures", nargs="*", default=[],
                        help="exploration temperatures at which bulk systems run NPT instead of NVT (dp_explore.py)")
    parser.add_argument("--extra-systems", nargs="*", default=[],
                        help="exploration sets added mid-run (e.g. s3): <run>/exploration/<name>/selected is labelled "
                             "for training and <name>-test/selected for the sealed test before the next training step; "
                             "<name>/*_initial.extxyz, if any, then joins the LAMMPS exploration")
    parser.add_argument("--runtime", choices=["docker", "apptainer"], default="docker")
    parser.add_argument("--cp2k-image", help=f"Docker tag or Apptainer .sif of {CP2K_IMAGE}")
    parser.add_argument("--mlip-image", help=f"Docker tag or Apptainer .sif of {MLIP_IMAGE}")
    parser.add_argument("--phase", choices=["all", "cpu", "gpu"], default="all",
                        help="cpu: CP2K labels with their screens and rechecks; gpu: training, tests, exploration and "
                             "validation. Stops before the first step of the other phase and writes it to <run>/next_phase")
    args = parser.parse_args()
    run = args.run.resolve()
    RUNTIME.update(runtime=args.runtime, cp2k=args.cp2k_image or CP2K_IMAGE, mlip=args.mlip_image or MLIP_IMAGE)
    if apptainer() and not (args.cp2k_image and args.mlip_image):
        raise SystemExit("--runtime apptainer needs --cp2k-image and --mlip-image (.sif files)")
    try:
        run.relative_to(ROOT)
    except ValueError:
        MOUNTS.insert(0, (run, "/hap-run"))
    state_path = run / "pipeline_state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else {"done": []}

    def mark(step, **info):
        state["done"].append(step)
        state.setdefault("info", {})[step] = info
        state_path.write_text(json.dumps(state, indent=2) + "\n")

    def begin(step):
        """Called before a step runs; with --phase, the first step of the other phase ends this run."""
        phase = step_phase(step)
        if args.phase not in ("all", phase):
            (run / "next_phase").write_text(phase + "\n", newline="\n")
            print(f"--phase {args.phase}: stopping before {step} ({phase} phase)", flush=True)
            raise SystemExit(0)

    dp = ROOT / "3.ML_AIMD/Step3_mlip_deepmd"

    def extras():
        return [n for n in args.extra_systems if f"label-extra-{n}" in state["done"]]

    def screen_models(before):
        """TF models of the newest ensemble explored before iteration `before`, for the label energy screen."""
        for k in range(before - 1, -1, -1):
            models = [run / f"iter-{k:03d}/explore/models/model-{i}.pb" for i in range(4)]
            if all(m.exists() for m in models):
                return [in_container(m) for m in models]
        return None

    def label_extras(before, wait=False):
        for n in args.extra_systems:
            if f"label-extra-{n}" in state["done"]:
                continue
            # A seed-<n> container labels these frames on the CPU while the GPU trains; join the next training
            # step instead of blocking this one (the final ensemble waits for it). Outside label() nothing should
            # hold it paused
            for s in seed_containers("paused", name=f"seed-{n}$"):
                sh(["docker", "unpause", s])
            # A set still being explored (the biased runs write their selection at the end) waits for its frames
            if not wait and not (run / f"exploration/{n}/selected").is_dir():
                print(f"extra {n}: no selected frames yet; deferred to the next training step", flush=True)
                continue
            if seed_containers("running", name=f"seed-{n}$"):
                # The last iteration waits, so that its exploration tests the ensemble that has seen these frames
                if not wait and before < args.iters - 1:
                    print(f"extra {n}: seed-{n} still labelling; deferred to the next training step", flush=True)
                    continue
                print(f"extra {n}: waiting for seed-{n}", flush=True)
                subprocess.run(["docker", "wait", f"seed-{n}"], capture_output=True)
            begin(f"label-extra-{n}")
            # Either part of a set may be empty
            for sub, vf in ((n, 0.1), (f"{n}-test", 0.0)):
                frames = sorted((run / f"exploration/{sub}/selected").glob("*_selected.extxyz"))
                if frames:
                    label(run, run / f"extra/{sub}/labels", frames, run / f"extra/{sub}/data", args, valid_fraction=vf,
                          screen_models=screen_models(before))
            mark(f"label-extra-{n}")

    def recheck_all(before):
        """Before a training step: every training set collected so far, screened with the models that screened
        label-`before`. Test sets are left alone: frames left out on the models' word would bias the test."""
        models = screen_models(before)
        if not models or f"recheck-{before}" in state["done"]:
            return
        begin(f"recheck-{before}")
        sets =[run / f"iter-{k:03d}" for k in range(before) if f"label-{k}" in state["done"]] + [
            run / f"extra/{n}" for n in extras()]
        changed = [s.name for s in sets if recheck(s / "labels", s / "data", models, args, f"recheck-{before}")]
        mark(f"recheck-{before}", collected_again=changed)

    def sealed_test():
        # Labelled on the CPU while the first ensemble trains on the GPU
        if "sealed_test" not in state["done"]:
            begin("sealed_test")
            label(run, run / "sealed-test/labels", sorted((run / "exploration/w8-test/selected").glob("*_selected.extxyz")),
                  run / "sealed-test/data", args, valid_fraction=0.0)
            mark("sealed_test")
        return sorted((run / "sealed-test/data/train").glob("*__*")) + sorted(
            p for n in extras() for p in (run / f"extra/{n}-test/data/train").glob("*__*"))

    def train_ensemble(tag, out, last_iter, steps, init_from=None, during=None):
        begin(tag)
        dirs =[run / f"iter-{k:03d}/data" for k in range(last_iter + 1)] + [run / f"extra/{n}/data" for n in extras()]
        train = sorted(p for d in dirs for p in (d / "train").glob("*__*"))
        valid = sorted(p for d in dirs for p in (d / "valid").glob("*__*"))
        sh([sys.executable, dp / "dp_ensemble.py", "write", "--train", *map(in_container, train),
            "--valid", *map(in_container, valid), "--steps", steps, "--out", out,
            *(["--sel-from", init_from] if init_from else [])])
        extra = ["--init-from", in_container(init_from)] if init_from else []
        # Under Apptainer each model gets one of the job's GPUs (CUDA_VISIBLE_DEVICES set by Slurm), in turn
        gpus = [g for g in os.environ.get("CUDA_VISIBLE_DEVICES", "").split(",") if g] if apptainer() else []
        # Concurrent GPU containers coincided with NVIDIA FECS errors and Docker VM crashes on
        # 2026-10-29, so models train --gpu-jobs at a time
        for start in range(0, 4, args.gpu_jobs):
            jobs = {f"dp-{tag}-{k}": mlip_start(["python", in_container(dp / "dp_ensemble.py"), "train", "--out",
                                                 in_container(out), "--index", k, *extra], f"dp-{tag}-{k}",
                                                gpu=gpus[k % len(gpus)] if gpus else None)
                    for k in range(start, min(4, start + args.gpu_jobs))}
            if during and start == 0:
                during()
            for name, job in jobs.items():
                code = mlip_wait(job)
                if code != 0:
                    raise SystemExit(f"{name} exited with {code}" + ("" if apptainer() else f"; see docker logs {name}"))
        mark(tag, n_train_systems=len(train), steps=steps, init_from=str(init_from) if init_from else None)

    def test_ensemble(out):
        mlip(["python", in_container(dp / "dp_ensemble.py"), "test", "--out", in_container(out),
              "--test", *map(in_container, sealed_test())])

    for it in range(args.iters):
        d = run / f"iter-{it:03d}"
        frames = (sorted((run / "exploration/w8/selected").glob("*_selected.extxyz")) if it == 0
                  else [run / f"iter-{it - 1:03d}/explore/candidates.extxyz"])
        if f"label-{it}" not in state["done"]:
            begin(f"label-{it}")
            label(run, d / "labels", frames, d / "data", args, screen_models=screen_models(it))
            mark(f"label-{it}")
        if f"train-{it}" not in state["done"]:
            label_extras(it)
            recheck_all(it)
            # Intermediate iterations warm-start from the previous ensemble; the final ensemble is
            # retrained from scratch on all labels after the loop
            warm = it > 0 and args.warm_steps > 0
            train_ensemble(f"train-{it}", d / "train", it, args.warm_steps if warm else args.steps,
                           init_from=run / f"iter-{it - 1:03d}/train" if warm else None,
                           during=sealed_test if it == 0 else None)
        if f"test-{it}" not in state["done"]:
            begin(f"test-{it}")
            test_ensemble(d / "train")
            mark(f"test-{it}")
        if f"explore-{it}" not in state["done"]:
            begin(f"explore-{it}")
            models = [in_container(d / f"train/model-{k}/graph-compress.pth") for k in range(4)]
            inits = sorted((run / "exploration/w8").glob("*_initial.extxyz")) + sorted(
                p for n in extras() for p in (run / f"exploration/{n}").glob("*_initial.extxyz"))
            # Converged only once every extra system was part of the exploration
            pending = [n for n in args.extra_systems if n not in extras()]
            mlip(["python", in_container(dp / "dp_explore.py"), "--models", *models,
                         "--structures", *map(in_container, inits), "--temperatures", "300", "500",
                         "--steps", str(20000 * (it + 1)), "--parallel", str(args.explore_jobs),
                         "--npt-temperatures", *args.npt_temperatures, "--out", in_container(d / "explore")])
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
    label_extras(args.iters, wait=True)
    last = max(k for k in range(args.iters) if f"label-{k}" in state["done"])
    if "train-final" not in state["done"]:
        recheck_all(args.iters)
        train_ensemble("train-final", run / "final/train", last, args.steps)
    if "test-final" not in state["done"]:
        begin("test-final")
        test_ensemble(run / "final/train")
        mark("test-final")

    # Stage 7: lattice, stability and NVE-drift gates of the final ensemble (dp_validate.py). Lattice references
    # are the CP2K cell optimisations of the matching DFT-hull structures. The w8 systems carry the formal gates
    # (500 ps); the extra systems are checked for 200 ps
    groups = [("w8", sorted((run / "exploration/w8").glob("*_initial.extxyz")), 500)] + [
        (n, sorted((run / f"exploration/{n}").glob("*_initial.extxyz")), 200) for n in extras()]
    models = [in_container(run / f"final/train/model-{k}/graph-compress.pth") for k in range(4)]
    for tag, structures, ps in groups:
        # An extra set without starting structures (the biased metad frames) has nothing to run MD from
        if f"validate-{tag}" in state["done"] or not structures:
            continue
        begin(f"validate-{tag}")
        v = run / f"final/validate/{tag}"
        v.mkdir(parents=True, exist_ok=True)
        ref = []
        if (run / "dft-hull/index.json").exists():
            sh([sys.executable, ROOT / "0.Generative_Discovery/screening/cp2k_hull.py", "cells", "--out", run / "dft-hull",
                "--match", *structures, "--write", v / "reference_cells.extxyz"])
            ref = ["--reference", in_container(v / "reference_cells.extxyz")]
        mlip(["python", in_container(dp / "dp_validate.py"), "--models", *models,
              "--structures", *map(in_container, structures), *ref, "--md-ps", str(ps),
              "--parallel", str(args.explore_jobs), "--out", in_container(v)])
        mark(f"validate-{tag}", **json.loads((v / "validation_summary.json").read_text())["formal"])
    (run / "next_phase").write_text("done\n", newline="\n")


if __name__ == "__main__":
    main()
