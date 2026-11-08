"""2026-11-07: the iteration-4 sealed test improved S3 forces (0.082-0.084 -> 0.066-0.068 eV/A) but not S3 energies
(1.5-2.6 -> 3.0-3.7 meV/atom). Per frame and per system: is the energy error a constant offset per system (mean) or
scatter (std)? Iteration-3 and -4 models (graph.pth, CPU) on the S3 sealed-test frames, and the iteration-4 model 0
on the S3 training frames (MACE-explored in extra/s3, DP-explored in iter-004). Run in hap-mlip:
  CUDA_VISIBLE_DEVICES= python bias.py > bias.json"""
import json
from pathlib import Path

import numpy as np
from deepmd.infer import DeepPot

run = Path("/work/.research-runs/pipeline-v2")


def errors(pot, sysdir):
    types = np.loadtxt(sysdir / "type.raw", dtype=int).reshape(-1)
    s = sysdir / "set.000"
    c, b, e = np.load(s / "coord.npy"), np.load(s / "box.npy"), np.load(s / "energy.npy").reshape(-1)
    pred = pot.eval(c, b, types)[0].reshape(-1)
    return 1000 * (pred - e) / len(types)  # meV/atom, model - DFT


def stats(d):
    return dict(frames=len(d), mean=float(d.mean()), std=float(d.std()), rmse=float(np.sqrt((d ** 2).mean())))


out = {"note": "meV/atom, model - DFT"}
test = sorted((run / "extra/s3-test/data/train").glob("s3_*"))
for it in (3, 4):
    for k in range(4):
        pot = DeepPot(str(run / f"iter-{it:03d}/train/model-{k}/graph.pth"))
        out[f"test it{it} model-{k}"] = {p.name.split("__")[0]: stats(errors(pot, p)) for p in test}
pot = DeepPot(str(run / "iter-004/train/model-0/graph.pth"))
for tag, root in (("train MACE (extra/s3)", run / "extra/s3/data"), ("train DP (iter-004)", run / "iter-004/data")):
    out[f"{tag} it4 model-0"] = {p.name.split("__")[0]: stats(np.concatenate([errors(pot, q) for q in (root / sub / p.name
                                                                               for sub in ("train", "valid"))
                                                                               if q.exists()]))
                                 for p in sorted((root / "train").glob("s3_*"))}
print(json.dumps(out, indent=1))
