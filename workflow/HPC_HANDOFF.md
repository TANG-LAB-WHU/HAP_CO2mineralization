# HPC handoff: prepared, disabled, and not submitted

The local run is a bounded implementation and small-model test. Its force labels,
datasets and learned models do not become production-ready by moving them to an HPC.
No connection, transfer, Slurm submission or cluster execution has been performed.

The existing `hpc/hap-stage.sbatch` covers MACE and CP2K with explicit reviewed restart
inputs. `hpc/hap-mlip.sbatch` adds DeepMD training/freeze and LAMMPS. It requires a new
attempt directory, a clean source commit, a pinned local Apptainer image, and a validated
predecessor manifest. It emits `schemas/run_manifest.example.json`-compatible provenance.
The predecessor must list output paths/checksums, and those actual files are rehashed before use. Restart files and standard model/trajectory outputs are also hashed. Its process-success result deliberately leaves scientific validation false. An operator
must run/review the relevant geometry, labeling, split, model or dynamics audit before
the next stage accepts that manifest. `hpc/dependencies.sh` prints a dependency plan;
it never submits anything. Slurm `afterok` alone does not perform scientific validation.

The scripts are templates and are not cluster-tested. DeepMD `--restart` syntax was
checked against the pinned 3.1.0 local executable. That does not establish HPC checkpoint
compatibility across versions/hardware. Supply both the actual checkpoint file/hash and
its CLI path prefix; never use `--force-load` to bypass missing state. CP2K continuation
uses matching coordinates, cell, wavefunction and Hessian with unchanged electronic and
convergence settings. LAMMPS restart inputs must keep the exact model and type mapping.

| Parameter | Measured local pilot / handoff meaning |
|---|---|
| Structure | Dry one-repeat (001), explicit antiparallel OH, all atoms mobile; no claim of stable termination |
| Numerical labels | PBE, DZVP-MOLOPT-SR-GTH, GTH Ca10/P5/O6/H1/C4; 1100/100 Ry and SCF 1e-8 selected by one-geometry force-grid sensitivity |
| Reference comparison | 300/40, 500/60, 700/60 Ry PBE; separate 500/60 Ry PBE-D3 method comparison; matched reference geometries |
| Geometry gate | Maximum vector force <0.05 eV/A plus all CP2K optimizer criteria; retained failed attempts are not accepted minima |
| AIMD pilot | Eight separately seeded local starts, 16 x0.25 fs =4 fs each; neither equilibrium nor decorrelation established |
| Training pilot | One-frame sanity fit and two compact architectures; complete trajectories assigned to separate splits before fitting |
| Dynamics pilot | Same-state engine comparison, 250 fs NVE, binary restart and half-timestep comparison; no long-time stability claim |
| Production candidates | Converge basis, cutoff/relative cutoff, k points, slab thickness/lateral size/vacuum and dipole treatment; justify hydration, protonation, sites and charge |
| Production sampling | Select physical temperature/ensemble; establish equilibration, correlation times and independent coverage; reactive labels and long-range validation need separate work |

Use the actual run's `verification.json` for outcomes. The table defines the pilot and
remaining studies, not a universal standard input or adequate production duration.
The immutable Docker digests are recorded in the local runner; preserve their resolved
image metadata and record the SHA256 of the corresponding reviewed Apptainer SIF.
Module-based use requires equivalent verified software and explicit CP2K data paths;
do not silently assume that a host module matches the tested container.

Only these site details remain necessary before adapting/enabling the handoff:

- Scheduler and allocation/account/partition (including required QoS).
- Available modules or Apptainer and the approved image/runtime route.
- Allocated CPU/GPU hardware and RAM, walltime, storage and job-limit constraints.

Repository/run paths are operator-selected parameters. Keep credentials outside the
repository. The example resource requests are placeholders, not measured cluster sizing.
