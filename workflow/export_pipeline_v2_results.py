"""Export a compact, hash-checked snapshot of a pipeline-v2 run directory for Git.

Copies the small records of every stage (state, protocol report, label audits, training inputs and learning
curves, sealed-test and exploration summaries, exploration start structures), rewrites host and container path
prefixes to repository-relative paths, and records the SHA-256 of every source file. Training and sealed-test
labels, models, checkpoints and trajectories stay in the run directory; only their identities are recorded.
Label steps that are still running are counted from their directories. Never starts a calculation.

  python workflow/export_pipeline_v2_results.py --run .research-runs/pipeline-v2 \
      --finite-size finite_size_pairs.csv --out results/pipeline-v2-20261104
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
IMAGES = ["hap-mlip:2026.10", "hap-cp2k:2026.2-cuda-sm89", "cp2k/cp2k:2026.2_mpich_generic_psmp"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def tree_identity(d):
    """SHA-256 over the sorted (relative path, file SHA-256) list of a directory, plus its size."""
    files = sorted(p for p in d.rglob("*") if p.is_file())
    h = hashlib.sha256()
    for p in files:
        h.update(f"{p.relative_to(d).as_posix()}\t{sha256(p)}\n".encode())
    return dict(files=len(files), bytes=sum(p.stat().st_size for p in files), tree_sha256=h.hexdigest())


def path_rewriter(repo):
    """Replace absolute host paths under the repository, and the container mount /work, by relative paths."""
    win, posix = str(repo).replace("/", "\\"), repo.as_posix()
    msys = "/" + posix[0].lower() + posix[2:] if posix[1:2] == ":" else posix
    # Windows paths appear plain in CSV files and with doubled backslashes inside JSON strings
    prefixes = [re.escape(p) for p in sorted({win + "\\", win.replace("\\", "\\\\") + "\\\\", posix + "/", msys + "/",
                                              "/work/"}, key=len, reverse=True)]
    pattern = re.compile("(?:" + "|".join(prefixes) + r")([^\s,\"'\]]*)")
    return lambda text: pattern.sub(lambda m: m.group(1).replace("\\\\", "/").replace("\\", "/"), text)


def label_progress(d):
    """Counts for a label directory whose step has not been collected yet."""
    frames = sorted(p for p in d.iterdir() if p.is_dir())
    failed = [p.name for p in frames if (p / "output.out").exists()
              and "SCF run NOT converged" in (p / "output.out").read_text(errors="replace")]
    done = [p.name for p in frames if (p / "done").exists()]
    return dict(frames=len(frames), done_marker=len(done), scf_not_converged=failed,
                converged_done=len([n for n in done if n not in failed]))


def docker_images():
    out = {}
    for tag in IMAGES:
        r = subprocess.run(["docker", "image", "inspect", "--format",
                            "{{.Id}}|{{json .RepoDigests}}|{{.Created}}", tag], capture_output=True, text=True)
        if r.returncode == 0:
            iid, digests, created = r.stdout.strip().split("|")
            out[tag] = dict(image_id=iid, repo_digests=json.loads(digests), created=created)
        else:
            out[tag] = "not inspected (image or docker CLI unavailable)"
    return out


PACKAGES = """echo "## base (DeePMD-kit) python: $(python --version 2>&1)"; python -m pip freeze 2>/dev/null
for v in mace mattersim sevenn; do echo "## venv-$v python: $(/opt/venv-$v/bin/python --version 2>&1)"
/opt/venv-$v/bin/pip freeze 2>/dev/null; done
echo "## lammps: $(lmp -h | sed -n 2p)"; echo "## plumed: $(plumed --no-mpi info --long-version 2>/dev/null | tail -1)"
"""


def image_packages(out):
    """Installed versions inside the images; Dockerfile.mlip does not pin its pip packages, so a rebuild may differ."""
    written = []
    r = subprocess.run(["docker", "run", "--rm", "--entrypoint", "bash", IMAGES[0], "-lc", PACKAGES],
                       capture_output=True, text=True)
    if r.returncode == 0:
        (out / "environment/hap-mlip-2026.10-packages.txt").write_text(r.stdout, encoding="utf-8", newline="")
        written.append("environment/hap-mlip-2026.10-packages.txt")
    r = subprocess.run(["docker", "run", "--rm", IMAGES[2], "cp2k", "--version"], capture_output=True, text=True)
    if r.returncode == 0:
        (out / "environment/cp2k-2026.2-version.txt").write_text(r.stdout, encoding="utf-8", newline="")
        written.append("environment/cp2k-2026.2-version.txt")
    return written


def main():
    parser = argparse.ArgumentParser(description="Compact snapshot of a pipeline-v2 run")
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--finite-size", type=Path, help="CSV written by finite_size_compare.py")
    args = parser.parse_args()
    run, out = args.run.resolve(), args.out.resolve()
    repo = run.parents[1] if run.parent.name == ".research-runs" else ROOT
    rewrite = path_rewriter(repo)
    # Rewrite the exported records; keep the hand-written README, REPRODUCE, scripts/ and tables/
    if out.exists():
        for p in out.iterdir():
            if p.name not in ("README.md", "REPRODUCE.md", "scripts", "tables"):
                shutil.rmtree(p) if p.is_dir() else p.unlink()
    files, excluded = [], {}

    def copy(src, dst):
        src, dst = run / src, out / dst
        dst.parent.mkdir(parents=True, exist_ok=True)
        data = src.read_bytes()
        if src.suffix in (".json", ".csv", ".md", ".out", ".txt"):
            text = data.decode("utf-8")
            new = rewrite(text)
            dst.write_text(new, encoding="utf-8", newline="")
            transformed = new != text
        else:
            dst.write_bytes(data)
            transformed = False
        files.append(dict(path=dst.relative_to(out).as_posix(), source=src.relative_to(repo).as_posix(),
                          source_sha256=hashlib.sha256(data).hexdigest(), host_paths_rewritten=transformed))

    state = json.loads((run / "pipeline_state.json").read_text())
    copy("pipeline_state.json", "state/pipeline_state.json")
    copy("cp2k-protocol/protocol_report.json", "protocol/protocol_report.json")
    copy("exploration/w8/summary.json", "exploration/w8_mace_summary.json")
    for group in ("w8", "w8-test", "s3", "s3-test"):
        for p in sorted((run / "exploration" / group).glob("*_initial.extxyz")):
            copy(p.relative_to(run), f"exploration/{group}/{p.name}")
    for name in ("label_summary.json", "label_audit.csv"):
        copy(f"sealed-test/data/{name}", f"sealed-test/{name}")
    excluded["sealed-test/data"] = tree_identity(run / "sealed-test/data")
    if args.finite_size:
        dst = out / "finite-size/finite_size_pairs.csv"
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(args.finite_size, dst)
        files.append(dict(path="finite-size/finite_size_pairs.csv", source="workflow/finite_size_compare.py output",
                          source_sha256=sha256(args.finite_size), host_paths_rewritten=False))
    progress = {}
    for it in range(100):
        d = run / f"iter-{it:03d}"
        if not d.exists():
            break
        tag = f"iter-{it:03d}"
        if f"label-{it}" in state["done"]:
            for name in ("label_summary.json", "label_audit.csv"):
                copy(f"{tag}/data/{name}", f"{tag}/{name}")
            for part in ("train", "valid"):
                if (d / "data" / part).exists():
                    excluded[f"{tag}/data/{part}"] = tree_identity(d / "data" / part)
        elif (d / "labels").exists():
            progress[f"label-{it}"] = label_progress(d / "labels")
        if f"train-{it}" in state["done"]:
            for k in range(4):
                for name in ("input.json", "lcurve.out"):
                    copy(f"{tag}/train/model-{k}/{name}", f"{tag}/train/model-{k}/{name}")
                for name in ("graph.pth", "graph-compress.pth", "model.ckpt.pt"):
                    p = d / f"train/model-{k}/{name}"
                    excluded[p.relative_to(run).as_posix()] = dict(bytes=p.stat().st_size, sha256=sha256(p))
        if f"test-{it}" in state["done"]:
            copy(f"{tag}/train/sealed_test.json", f"{tag}/sealed_test.json")
        if f"explore-{it}" in state["done"]:
            copy(f"{tag}/explore/explore_summary.json", f"{tag}/explore_summary.json")
            p = d / "explore/candidates.extxyz"
            excluded[p.relative_to(run).as_posix()] = dict(bytes=p.stat().st_size, sha256=sha256(p))
        for inval in sorted(d.glob("explore_INVALID_*")):
            copy(f"{tag}/{inval.name}/explore_summary.json", f"{tag}/{inval.name}_summary.json")
    for name in ("s3", "s3-test"):
        if (run / f"extra/{name}/labels").exists():
            progress[f"extra-{name}"] = label_progress(run / f"extra/{name}/labels")
    snapshot = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat()
    (out / "state").mkdir(parents=True, exist_ok=True)
    (out / "state/progress.json").write_text(json.dumps(dict(snapshot_utc=snapshot, stages_done=state["done"],
                                                             uncollected_label_steps=progress), indent=2) + "\n")
    commit = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain", "--", "workflow", "3.ML_AIMD",
                            "0.Generative_Discovery"], capture_output=True, text=True).stdout.split("\n")
    environment = dict(docker_images=docker_images(),
                       gpu=subprocess.run(["nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader"],
                                          capture_output=True, text=True).stdout.strip() or "not queried")
    (out / "environment.json").write_text(json.dumps(environment, indent=2) + "\n")
    (out / "environment").mkdir(exist_ok=True)
    generated = ["state/progress.json", "environment.json", ".gitattributes"] + image_packages(out)
    (out / ".gitattributes").write_text("* -text\n")
    manifest = dict(schema="hap-pipeline-v2-snapshot/v1", snapshot_utc=snapshot,
                    run_directory=run.relative_to(repo).as_posix(),
                    run_checkout_commit=commit, run_checkout_uncommitted_code=[l.strip() for l in dirty if l.strip()],
                    scientific_calculations_run_for_export=False, manuscript_eligible=False,
                    production_qualified=False, files=[], excluded_assets=excluded)
    for f in files + [dict(path=p) for p in generated]:
        p = out / f["path"]
        f.update(bytes=p.stat().st_size, sha256=sha256(p))
        manifest["files"].append(f)
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"{len(manifest['files'])} files written to {out}; {len(excluded)} excluded assets identified")


if __name__ == "__main__":
    main()
