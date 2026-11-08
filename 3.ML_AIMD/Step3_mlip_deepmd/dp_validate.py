"""Pipeline v2 validation gates with LAMMPS + DeePMD (stage 7 of workflow/PIPELINE_V2.md).

Per system, with model 0 of the ensemble (TensorFlow backend, converted and checked as in dp_explore):
  bulk  - 0 K cell relaxation (box/relax tri, 1 atm) of the input cell, then NPT 300 K, 1 atm for --md-ps on the
          input cell (the >= 8 A label cell, the same cell as the CP2K reference) or a supercell with every width
          >= --min-width;
  slab  - NVT 300 K for --md-ps (systems named hap001*, as in cp2k_label);
  both  - an NVE segment of --nve-ps from the end of that run.
Gates (thresholds from PIPELINE_V2.md, never relaxed here):
  lattice   - NPT averages after --equil-ps, per formula cell, within 1 % of the CP2K cell optimisation
              (--reference); the formal gate is HAP a and c, the other phases are reported the same way;
  stability - every dumped frame passes the geometry gate: minimum distance > 0.7 A, every P has 4 O within
              1.9 A, every H has its nearest O at 0.8-1.45 A (H_O_MAX) and every C has 2 or 3 O within 1.6 A (CO2 or
              CO3); the run must also finish;
  NVE drift - |slope of the total energy| <= 1 meV/atom/ps.
The 4-model deviation on the dumped frames is reported, not gated. Run inside hap-mlip.

  python dp_validate.py --models iter-005/train/model-*/graph-compress.pth --structures init/*.extxyz \
      --reference refs.extxyz --out final/validate
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dp_explore import TYPE_MAP, is_slab, python_devi, tf_models  # noqa: E402

GATES = {"lattice_rel": 0.01, "nve_drift_meV_atom_ps": 1.0}
# Was 1.2 A (the v1 gate of refine_pipeline.py). Shawn set 1.45 A on 2026-11-06, after the iteration-3 preview and
# before the formal run: monetite H are shared in short hydrogen bonds (O...O 2.47/2.58 A) and sit at 1.21-1.27 A
# from their nearest O (workflow/PIPELINE_V2.md, "Stage-7 preview"). Same limit as the metad selection gate.
H_O_MAX = 1.45

HEAD = """units metal
boundary p p p
atom_style atomic
read_data {data}
{masses}
pair_style deepmd {model}
pair_coeff * *
# Skin 0, as in dp_explore (the PyTorch plugin miscounts neighbours with a skin; kept for the TF backend too)
neighbor 0.0 bin
neigh_modify every 1 delay 0 check yes
"""

RELAX = HEAD + """thermo_style custom step pe press vol cella cellb cellc cellalpha cellbeta cellgamma
thermo 100
fix 1 all box/relax tri 1.01325 vmax 0.001
min_style cg
minimize 0 1e-4 5000 50000
minimize 0 1e-4 5000 50000
minimize 0 1e-4 5000 50000
print "CELL $(cella) $(cellb) $(cellc) $(cellalpha) $(cellbeta) $(cellgamma) $(press) $(pe)" file cell.txt
"""

MD = HEAD + """velocity all create {temp} {seed} mom yes rot yes dist gaussian
timestep {dt}
thermo_style custom step temp pe ke etotal press vol cella cellb cellc cellalpha cellbeta cellgamma
# The default 8 significant digits round the total energy of a 640-atom cell to 0.01 eV
thermo_modify format float %.12g
thermo {thermo}
fix 1 all {ensemble}
dump 1 all custom {dump} traj.dump id type x y z
dump_modify 1 sort id
run {steps}
unfix 1
undump 1
log nve.log
thermo_style custom step temp pe ke etotal
thermo_modify format float %.12g
thermo 20
fix 2 all nve
run {nve_steps}
"""


def repeats(atoms, min_width):
    cell = atoms.cell.array
    vol = abs(np.linalg.det(cell))
    widths = [vol / np.linalg.norm(np.cross(cell[(k + 1) % 3], cell[(k + 2) % 3])) for k in range(3)]
    return [max(1, math.ceil(min_width / w - 1e-9)) for w in widths]


def thermo_blocks(log):
    """Thermo tables of a LAMMPS log, one dict of column arrays per run."""
    blocks, cols, rows = [], None, []
    for line in Path(log).read_text().splitlines():
        parts = line.split()
        if parts and parts[0] == "Step":
            cols, rows = parts, []
        elif cols and line.startswith("Loop time"):
            blocks.append({c: np.array([r[k] for r in rows]) for k, c in enumerate(cols)})
            cols = None
        elif cols:
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if len(vals) == len(cols):
                rows.append(vals)
    return blocks


def drift(t_ps, etotal, natoms):
    """Linear drift of the total energy in meV/atom/ps."""
    return float(np.polyfit(t_ps, etotal, 1)[0] / natoms * 1000)


def geometry_gate(a):
    from ase.neighborlist import neighbor_list
    sym = np.array(a.get_chemical_symbols())
    i, j, d = neighbor_list("ijd", a, 2.0)
    fails = []
    if len(d) and d.min() <= 0.7:
        fails.append(f"min distance {d.min():.2f} A")
    for k in np.where(sym == "P")[0]:
        n = int(np.sum((i == k) & (sym[j] == "O") & (d < 1.9)))
        if n != 4:
            fails.append(f"P{k} has {n} O within 1.9 A")
    for k in np.where(sym == "H")[0]:
        m = d[(i == k) & (sym[j] == "O")]
        if not len(m) or not 0.8 < m.min() < H_O_MAX:
            fails.append(f"H{k} nearest O {m.min() if len(m) else float('inf'):.2f} A")
    for k in np.where(sym == "C")[0]:
        n = int(np.sum((i == k) & (sym[j] == "O") & (d < 1.6)))
        if n not in (2, 3):
            fails.append(f"C{k} has {n} O within 1.6 A")
    return fails


def lattice_compare(npt, ref, rep):
    """Relative deviation of the NPT lattice (per input cell) from the reference cell."""
    par = np.array(npt) / np.r_[rep, 1, 1, 1]
    refp = np.array(ref)
    vol = lambda c: float(c[0] * c[1] * c[2] * np.sqrt(1 - np.sum(np.cos(np.radians(c[3:])) ** 2)
                                                      + 2 * np.prod(np.cos(np.radians(c[3:])))))
    return dict(md=par.round(4).tolist(), reference=refp.round(4).tolist(),
                rel_volume_dev=round(vol(par) / vol(refp) - 1, 5),
                rel_length_dev=((par[:3] - refp[:3]) / refp[:3]).round(5).tolist(),
                angle_dev_deg=(par[3:] - refp[3:]).round(3).tolist())


def main():
    parser = argparse.ArgumentParser(description="DeePMD validation gates (lattice, stability, NVE drift)")
    parser.add_argument("--models", nargs=4, type=Path, required=True)
    parser.add_argument("--structures", nargs="+", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=None,
                        help="extxyz of CP2K-relaxed cells, info['system'] = structure name (cp2k_hull.py cells)")
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--md-ps", type=float, default=500.0)
    parser.add_argument("--equil-ps", type=float, default=100.0, help="dropped before lattice averages")
    parser.add_argument("--nve-ps", type=float, default=50.0)
    parser.add_argument("--dt-fs", type=float, default=0.5)
    parser.add_argument("--dump-ps", type=float, default=1.0)
    # A 14 A supercell (640-720 atoms) costs 50-100 ms per step on the RTX 4090, about 22 h per 500 ps run
    parser.add_argument("--min-width", type=float, default=0.0, help="A, supercell for NPT (0: the input cell)")
    parser.add_argument("--trust-lo", type=float, default=0.10)
    parser.add_argument("--seed", type=int, default=20261101)
    parser.add_argument("--parallel", type=int, default=3)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    structures = {p.stem.replace("_initial", ""): read(p) for p in args.structures}
    # cp2k_hull.py cells writes an empty file when no structure has a finished reference yet
    refs = ({a.info["system"]: a for a in read(args.reference, ":")}
            if args.reference and args.reference.exists() and args.reference.stat().st_size else {})
    pbs, convert = tf_models([m.resolve() for m in args.models], out / "models", list(structures.values()), 1e-6)
    masses = "\n".join(f"mass {i + 1} {atomic_masses[atomic_numbers[s]]:.4f}" for i, s in enumerate(TYPE_MAP))
    steps = round(args.md_ps * 1000 / args.dt_fs)
    jobs = []
    for name, a in structures.items():
        d = out / name
        rep = [1, 1, 1] if is_slab(name) else repeats(a, args.min_width)
        md = d / "md"
        md.mkdir(parents=True, exist_ok=True)
        write(md / "conf.lmp", a.repeat(rep), format="lammps-data", specorder=TYPE_MAP, masses=False)
        p = 1.01325
        ensemble = (f"nvt temp {args.temperature} {args.temperature} 0.05" if is_slab(name)
                    else f"npt temp {args.temperature} {args.temperature} 0.05 tri {p} {p} 0.5 flip no")
        (md / "in.lmp").write_text(MD.format(data="conf.lmp", masses=masses, model=pbs[0], temp=args.temperature,
                                             seed=args.seed, dt=args.dt_fs / 1000, thermo=200, ensemble=ensemble,
                                             dump=round(args.dump_ps * 1000 / args.dt_fs), steps=steps,
                                             nve_steps=round(args.nve_ps * 1000 / args.dt_fs)))
        jobs.append(md)
        if not is_slab(name):
            rl = d / "relax"
            rl.mkdir(parents=True, exist_ok=True)
            write(rl / "conf.lmp", a, format="lammps-data", specorder=TYPE_MAP, masses=False)
            (rl / "in.lmp").write_text(RELAX.format(data="conf.lmp", masses=masses, model=pbs[0]))
            jobs.append(rl)
        (d / "repeat.json").write_text(json.dumps(rep) + "\n")
    env = {"DP_INTRA_OP_PARALLELISM_THREADS": "2", "DP_INTER_OP_PARALLELISM_THREADS": "1", **os.environ}

    # Output and exit code are kept: the preview run of hap001_water_co2 died at 154.6 ps (2026-11-05) with its output
    # discarded; the kernel log showed the WSL GPU driver failing an allocation (dxg ioctl -12) and lmp aborting.
    def run_lmp(w):
        with open(w / "lmp.out", "w") as log:
            rc = subprocess.run(["lmp", "-in", "in.lmp", "-log", "log.lammps"], cwd=w, env=env, stdout=log,
                                stderr=subprocess.STDOUT).returncode
        (w / "exit_code.txt").write_text(f"{rc}\n")

    with ThreadPoolExecutor(args.parallel) as pool:
        list(pool.map(run_lmp, [w for w in jobs if not (w / ("cell.txt" if w.name == "relax" else "nve.log")).exists()]))

    from deepmd.infer import DeepPot
    pots = [DeepPot(str(m)) for m in pbs]
    report = {}
    for name, a in structures.items():
        d = out / name
        rep = json.loads((d / "repeat.json").read_text())
        natoms = len(a) * int(np.prod(rep))
        r = dict(kind="slab" if is_slab(name) else "bulk", natoms_md=natoms, repeat=rep)
        blocks = thermo_blocks(d / "md/log.lammps") if (d / "md/log.lammps").exists() else []
        nve = thermo_blocks(d / "md/nve.log") if (d / "md/nve.log").exists() else []
        finished = len(blocks) >= 1 and len(nve) >= 1 and len(nve[0]["Step"]) > 1
        r["finished"] = finished
        if (d / "md/exit_code.txt").exists():
            r["lmp_exit_code"] = int((d / "md/exit_code.txt").read_text())
        frames = read(d / "md/traj.dump", ":", format="lammps-dump-text") if (d / "md/traj.dump").exists() else []
        symbols = a.repeat(rep).get_chemical_symbols()
        fails = []
        for k, f in enumerate(frames):
            f.set_chemical_symbols(symbols)
            fl = geometry_gate(f)
            if fl:
                fails.append(dict(time_ps=k * args.dump_ps, first=fl[:3], n=len(fl)))
        r["stability"] = dict(frames=len(frames), md_ps=(len(frames) - 1) * args.dump_ps if frames else 0.0,
                              failed_frames=len(fails), first_failures=fails[:5],
                              passed=bool(finished and frames and not fails and (len(frames) - 1) * args.dump_ps >= args.md_ps - 1e-9))
        if frames:
            pick = np.unique(np.linspace(0, len(frames) - 1, min(100, len(frames))).round().astype(int))
            dev = python_devi(pots, [frames[k] for k in pick])
            r["model_devi"] = dict(frames=len(pick), max_devi_f_max=float(dev.max()),
                                   frac_ge_trust_lo=float(np.mean(dev >= args.trust_lo)))
        if nve and len(nve[0]["Step"]) > 1:
            t = (nve[0]["Step"] - nve[0]["Step"][0]) * args.dt_fs / 1000
            r["nve"] = dict(ps=float(t[-1]), drift_meV_atom_ps=drift(t, nve[0]["TotEng"], natoms),
                            etotal_std_meV_atom=float(nve[0]["TotEng"].std() / natoms * 1000))
            r["nve"]["passed"] = abs(r["nve"]["drift_meV_atom_ps"]) <= GATES["nve_drift_meV_atom_ps"]
        if not is_slab(name):
            if (d / "relax/cell.txt").exists():
                c = [float(x) for x in (d / "relax/cell.txt").read_text().split()[1:]]
                r["relax_0K"] = dict(cellpar=c[:6], press_bar=c[6])
            if blocks:
                b = blocks[0]
                t = (b["Step"] - b["Step"][0]) * args.dt_fs / 1000
                keep = t >= args.equil_ps
                cols = ["Cella", "Cellb", "Cellc", "CellAlpha", "CellBeta", "CellGamma"]
                if keep.sum() >= 10:
                    vals = np.array([b[c][keep] for c in cols])
                    # 5 contiguous blocks give the standard error of the mean
                    sem = np.array([np.std([x.mean() for x in np.array_split(v, 5)], ddof=1) / math.sqrt(5) for v in vals])
                    r["npt"] = dict(ps_averaged=float(t[keep][-1] - t[keep][0]), cellpar=vals.mean(1).tolist(),
                                    sem=sem.tolist())
            if name in refs and "npt" in r:
                ref = refs[name].cell.cellpar()
                r["lattice"] = lattice_compare(r["npt"]["cellpar"], ref, rep)
                r["lattice"]["reference_source"] = refs[name].info.get("source", "")
                r["lattice"]["passed"] = bool(np.all(np.abs(r["lattice"]["rel_length_dev"]) <= GATES["lattice_rel"]))
                if "relax_0K" in r:
                    r["lattice_0K"] = lattice_compare(r["relax_0K"]["cellpar"], ref, [1, 1, 1])
        report[name] = r
        print(name, json.dumps({k: v for k, v in r.items() if k in ("finished", "stability", "nve", "lattice")})[:400],
              flush=True)
    hap = report.get("hap_bulk", {}).get("lattice")
    summary = dict(gates=GATES, md_ps=args.md_ps, nve_ps=args.nve_ps, temperature_K=args.temperature, dt_fs=args.dt_fs,
                   model=str(pbs[0]), conversion_check=convert, systems=report,
                   formal=dict(hap_lattice_a_c=None if hap is None else bool(abs(hap["rel_length_dev"][0]) <= 0.01
                                                                              and abs(hap["rel_length_dev"][2]) <= 0.01),
                               stability_all=all(r["stability"]["passed"] for r in report.values()),
                               nve_all=all(r.get("nve", {}).get("passed", False) for r in report.values())))
    (out / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("formal gates", summary["formal"])


if __name__ == "__main__":
    main()
