"""Run manifest (workflow/schemas/run_manifest.example.json) for one pipeline-v2 Slurm job.

  start:  before the orchestrator: source commit, arguments file, images and the steps already done
  finish: after it: return code, the steps this job completed (with their pipeline_state.json info) and the phase
          the orchestrator wrote to <run>/next_phase

Manifests go to <run>/manifests/job-<SLURM_JOB_ID>.manifest.json and are never overwritten. A normal exit never
sets validation.passed: the gates recorded in pipeline_state.json and reviewed in PIPELINE_V2.md decide that.
Standard library only, so it runs with the host Python outside the images.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
STAGE = {"cpu": "aimd", "gpu": "mlip_training", "all": "mlip_training"}  # CP2K labels are DFT single points
UNITS = {"mean_accurate": "fraction of explored frames", "n_candidates": "frames", "n_train_systems": "systems",
         "steps": "training steps"}


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def image_ref(p: Path) -> dict:
    """SIF images are tens of GB; their checksum is kept beside them (<image>.sha256) and reused while newer."""
    p = p.resolve()
    side = p.with_name(p.name + ".sha256")
    if side.exists() and side.stat().st_mtime >= p.stat().st_mtime:
        digest = side.read_text().split()[0]
    else:
        digest = sha(p)
        try:
            side.write_text(f"{digest}  {p.name}\n")
        except OSError:
            pass
    return {"path": str(p), "sha256": digest}


def git(*args) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("action", choices=["start", "finish"])
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--phase", choices=sorted(STAGE), help="start: the --phase given to the orchestrator")
    p.add_argument("--args-file", type=Path, help="start: the sourced arguments file")
    p.add_argument("--command", help="start: the orchestrator command line")
    p.add_argument("--image", type=Path, action="append", default=[], help="start: SIF images used")
    p.add_argument("--returncode", type=int, help="finish: the orchestrator's exit code")
    a = p.parse_args()
    run = a.run.resolve()
    job = os.environ.get("SLURM_JOB_ID", "manual")
    if not job.replace("_", "").replace("-", "").isalnum():
        raise SystemExit(f"unexpected job id {job!r}")
    path = run / "manifests" / f"job-{job}.manifest.json"
    state_path = run / "pipeline_state.json"
    state = json.loads(state_path.read_text()) if state_path.exists() else {"done": []}
    if a.action == "start":
        if path.exists():
            raise SystemExit(f"{path} exists; a manifest is never overwritten")
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "schema_version": 1,
            "run_id": f"{run.name}-{job}",
            "status": "running",
            "backend": "slurm",
            "stage": STAGE[a.phase],
            "git": {"repository": "TANG-LAB-WHU/HAP_CO2mineralization", "branch": git("branch", "--show-current"),
                    "commit": git("rev-parse", "HEAD"), "dirty": bool(git("status", "--porcelain"))},
            "configuration": {"path": str(a.args_file.resolve()), "sha256": sha(a.args_file)},
            "inputs": [{"path": str(state_path), "sha256": sha(state_path)}] if state_path.exists() else [],
            "environment": {
                "software": {"host": platform.node(), "python": platform.python_version()},
                "container": [image_ref(i) for i in a.image],
                "scheduler_job_id": job,
                "allocation": {k: os.environ.get(k) for k in ("SLURM_JOB_PARTITION", "SLURM_NTASKS",
                                                               "SLURM_CPUS_PER_TASK", "SLURM_MEM_PER_NODE",
                                                               "SLURM_JOB_GPUS", "CUDA_VISIBLE_DEVICES")},
            },
            "outputs": [],
            "metrics": {},
            "validation": {"criteria": ["Gates of PIPELINE_V2.md, checked on the recorded outputs after the job"],
                           "passed": False, "reviewer": None, "reviewed_at": None},
            "started_at": now(),
            "finished_at": None,
            "pipeline_v2": {"phase": a.phase, "command": a.command, "done_at_start": state["done"]},
            "manuscript_eligible": False,
            "notes": "Process status only; scientific validity is reviewed separately.",
        }
    else:
        data = json.loads(path.read_text())
        before = set(data["pipeline_v2"]["done_at_start"])
        new = [s for s in state["done"] if s not in before]
        nxt = run / "next_phase"
        data.update(status="completed" if a.returncode == 0 else "failed", returncode=a.returncode, finished_at=now())
        data["pipeline_v2"].update(completed_steps=new, next_phase=nxt.read_text().strip() if nxt.exists() else None)
        data["outputs"] = [{"path": str(f), "sha256": sha(f)} for f in (state_path, nxt) if f.exists()]
        data["metrics"] = {f"{s}.{k}": {"value": v, "unit": UNITS.get(k, "as recorded in pipeline_state.json")}
                           for s in new for k, v in state.get("info", {}).get(s, {}).items()
                           if isinstance(v, (int, float)) and not isinstance(v, bool)}
    path.write_text(json.dumps(data, indent=2) + "\n")


if __name__ == "__main__":
    main()
