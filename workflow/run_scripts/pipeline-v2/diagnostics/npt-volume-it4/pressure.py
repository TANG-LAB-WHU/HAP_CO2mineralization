"""2026-11-08: explore-4 ran the bulk systems at 300 K in NPT (1 atm) for the first time. On the old bulk systems the
mean volume of the second half agrees within about 1 % with where a linear fit of the DFT pressure of the labels
crosses -P_kin; on several S3 systems the cells shrank 3-7 % further than that fit (s3_055 0.909 vs 0.953, s3_302
0.932 vs 1.000), and no S3 label lies below 0.984 V0. Two questions for the iteration-4 ensemble (graph.pth):
  1. labels: model - DFT potential pressure on every labelled bulk frame (train + valid), per system;
  2. extrapolation: model pressure along an isotropic scan 0.88-1.06 V0 of the first S3 label (MACE-explored, near V0).
P_pot = trace(virial) / (3 V) in GPa; the sign was checked on hap_bulk (DFT +2.55 GPa at V0, NPT expanded 2.0 %).
Run in hap-mlip: python pressure.py > pressure.json"""
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from deepmd.infer import DeepPot

run = Path("/work/.research-runs/pipeline-v2")
GPa = 160.21766208
pots = [DeepPot(str(run / f"iter-004/train/model-{k}/graph.pth")) for k in range(4)]


def p_pot(virial, box):
    v = np.abs(np.linalg.det(box.reshape(-1, 3, 3)))
    return np.trace(virial.reshape(-1, 3, 3), axis1=1, axis2=2) / 3 / v * GPa


def model_p(c, b, types):
    return np.array([p_pot(pot.eval(c, b, types)[2], b) for pot in pots])  # (models, frames)


labels = defaultdict(lambda: defaultdict(list))
for d in sorted([*run.glob("iter-*/data/*/*"), *run.glob("extra/*/data/*/*")]):
    s = d / "set.000"
    if not (s / "virial.npy").exists():
        continue
    types = np.loadtxt(d / "type.raw", dtype=int).reshape(-1)
    c, b, w = np.load(s / "coord.npy"), np.load(s / "box.npy"), np.load(s / "virial.npy")
    g = labels[d.name.split("__")[0]]
    g["dft"].append(p_pot(w, b))
    g["model"].append(model_p(c, b, types))
    g["volume"].append(np.abs(np.linalg.det(b.reshape(-1, 3, 3))))

out = {"note": "GPa; potential pressure only; model = mean of the 4 iteration-4 models", "labels": {}, "scan": {}}
for name, g in sorted(labels.items()):
    dft, mod = np.concatenate(g["dft"]), np.concatenate(g["model"], axis=1)
    diff = mod.mean(0) - dft
    out["labels"][name] = dict(frames=len(dft), dft_mean=float(dft.mean()), model_mean=float(mod.mean()),
                               mean_diff=float(diff.mean()), rmse=float(np.sqrt((diff ** 2).mean())),
                               model_spread=float(mod.std(0).mean()))

for d in sorted((run / "extra/s3/data/train").glob("s3_*")):
    types = np.loadtxt(d / "type.raw", dtype=int).reshape(-1)
    c0 = np.load(d / "set.000/coord.npy")[:1]
    b0 = np.load(d / "set.000/box.npy")[:1]
    rows = []
    for x in np.arange(0.88, 1.0601, 0.02):
        f = x ** (1 / 3)
        p = model_p(c0 * f, b0 * f, types)[:, 0]
        rows.append(dict(V_over_V0=round(float(x), 2), model_mean=float(p.mean()), model_spread=float(p.std())))
    out["scan"][d.name.split("__")[0]] = rows
print(json.dumps(out, indent=1))
