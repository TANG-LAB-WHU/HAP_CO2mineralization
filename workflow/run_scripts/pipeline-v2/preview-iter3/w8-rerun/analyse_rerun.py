"""2026-11-06: the gates of dp_validate.py on the standalone rerun of hap001_water_co2 (slab, so no lattice).
Same functions, same models (preview-iter3/w8/models), same trust_lo, 100 frames for the deviation. Run in hap-mlip:
  python analyse_rerun.py > rerun_summary.json"""
import json
import sys
from pathlib import Path

import numpy as np
from ase.io import read

sys.path.insert(0, "/work/3.ML_AIMD/Step3_mlip_deepmd")
import dp_validate  # noqa: E402
from dp_explore import TYPE_MAP, python_devi  # noqa: E402

here = Path(__file__).resolve().parent / "hap001_water_co2"
models = sorted((here.parents[1] / "w8/models").glob("model-*.pb"))
frames = read(here / "traj.dump", ":", format="lammps-dump-text")
for f in frames:
    f.set_chemical_symbols([TYPE_MAP[z - 1] for z in f.get_atomic_numbers()])
natoms = len(frames[0])
r = dict(models=[str(m) for m in models], frames=len(frames), md_ps=len(frames) - 1)
for lim in (1.2, dp_validate.H_O_MAX):
    dp_validate.H_O_MAX = lim
    fails = [(k, fl) for k, f in enumerate(frames) if (fl := dp_validate.geometry_gate(f))]
    r[f"geometry_H_O_max_{lim}"] = dict(failed_frames=len(fails), first_failures=[dict(time_ps=k, first=fl[:3])
                                                                                   for k, fl in fails[:5]])
blocks = dp_validate.thermo_blocks(here / "log.lammps")
r["thermo_blocks"] = len(blocks)
nve = dp_validate.thermo_blocks(here / "nve.log")
if nve and len(nve[0]["Step"]) > 1:
    t = (nve[0]["Step"] - nve[0]["Step"][0]) * 0.5 / 1000
    r["nve"] = dict(ps=float(t[-1]), drift_meV_atom_ps=dp_validate.drift(t, nve[0]["TotEng"], natoms),
                    etotal_std_meV_atom=float(nve[0]["TotEng"].std() / natoms * 1000))
    r["nve"]["passed"] = abs(r["nve"]["drift_meV_atom_ps"]) <= dp_validate.GATES["nve_drift_meV_atom_ps"]
from deepmd.infer import DeepPot  # noqa: E402
pots = [DeepPot(str(m)) for m in models]
pick = np.unique(np.linspace(0, len(frames) - 1, min(100, len(frames))).round().astype(int))
dev = python_devi(pots, [frames[k] for k in pick])
r["model_devi"] = dict(frames=len(pick), max_devi_f_max=float(dev.max()), frac_ge_trust_lo=float(np.mean(dev >= 0.10)),
                       frac_ge_trust_lo_first_154ps=float(np.mean(dev[pick <= 154] >= 0.10)),
                       frac_ge_trust_lo_after_154ps=float(np.mean(dev[pick > 154] >= 0.10)))
print(json.dumps(r, indent=1))
