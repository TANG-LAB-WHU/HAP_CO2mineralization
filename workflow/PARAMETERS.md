# Local smoke and proposed HPC parameters

All values below are **functional smoke settings or unvalidated production candidates**.
There is no production-converged adsorption ranking, AIMD ensemble, MLIP, or mechanism.
All new and legacy computational artifacts remain manuscript-ineligible.

## Structure and matched references

The repository ICSD 26205 CIF is retained unchanged. Its half-occupied OH sites cannot
be used as simultaneous full atoms. The smoke model explicitly selects O fractional z
0.1978/0.6978 and H z 0.0608/0.5608, giving Ca10(PO4)6(OH)2 (44 atoms).
This is a modeling assumption, **not human-reviewed ground truth**. Other OH orderings,
surface stoichiometries/terminations, polarity, protonation and hydration must be examined.
The (001) one-cell slab has 12 A vacuum on each side. Bottom-half atoms are fixed for
MACE relaxation; the three-reference cells and periodicity are identical. One CO2 is
placed horizontally over a top Ca at 2.8 A. No claim of exhaustive site selection is made.
The AIMD plumbing examples allow all atoms to move and do not preserve the MACE bottom
constraint; they intentionally establish execution only, not a surface-physics ensemble.

Adsorption convention: E(slab+CO2)-E(slab)-E(CO2), in eV per CO2. Negative is favorable
under this **electronic-energy convention**, not a finite-temperature free energy.
Use matching functional, pseudopotentials, basis, dispersion, cell, charge and settings.
Never mix MACE absolute energies with CP2K absolute energies. The observed sign difference
between MACE and CP2K on the first smoke geometries is a failed validation of interchangeability.

## Parameter table

| Quantity | Executed smoke | Proposed production investigation (not enabled) |
|---|---|---|
| MACE | MP-0 small, float64, CUDA, fmax 0.05 eV/A; 12 then max 200 FIRE steps | Evaluate newer suitable models and committees against DFT; verify elemental/domain coverage; fmax 0.02 eV/A trial |
| Slab | (001), one unit-cell repeat, 24 A total vacuum | Thickness series, lateral 1x1/2x2/3x3 and vacuum series; enumerate terminations, sites, orientations, coverages |
| CP2K | 2025.1, PBE, DZVP-MOLOPT-SR-GTH, Ca q10/P q5/O q6/H q1/C q4 | Consistent dispersion choice and basis comparison; choose only after adsorption-energy convergence |
| Cutoff / relative cutoff | 300/40 Ry | Test 400/600/800 Ry and 40/60/80 Ry; seek <0.02 eV change in adsorption difference (proposed gate) |
| SCF | 1e-6, OT DIIS, outer SCF; 1e-8 for stress test | 1e-7 to 1e-8 trials; all frames converged; never ignore failures |
| Boundary / k points | XYZ / Gamma | Test in-plane k sampling; slab dipole treatment and image interactions; isolation convergence for molecular reference |
| AIMD | NVE, 100/200/300 K initial velocities, seeds 710/711/712, 5 x 0.25 fs = 1.25 fs each | Trial NVT at scientifically chosen temperature, CSVR tau 100 fs, dt 0.25 fs; compare smaller dt, equilibrate, run independent replicas until observables/autocorrelation converge |
| MLIP | 20 GPU updates, tiny se_e2_a, 5 frames per run after excluding shared initial geometries | Split by independent trajectory/composition; test holdout force RMSE <0.1 eV/A and energy <10 meV/atom as provisional gates; assess tails, reactions and extrapolation |
| LAMMPS | 10 x 0.25 fs, uninterrupted vs 5+restart+5 | Validate forces/energies against the training engine, NVE drift, timestep, system-size and sampling convergence before scale-up |

200 steps x 0.25 fs is only 50 fs. Neither that duration nor normal termination proves
equilibration or adequate sampling. The executed five-step runs are even shorter.
No pressure/virial validation can be inferred from DeepMD's printed zero virial errors:
the tiny training sets intentionally contain no virial labels.

## Stress convention

CP2K's printed stress is compressive-positive; ASE uses the opposite sign.
DeepMD virial (eV) = V(A^3) x CP2K stress(bar) x 1e5/1.602176634e11.
Equivalently virial = -V x ASE stress(eV/A^3). The tensor order is xx,xy,xz,yx,yy,yz,zx,zy,zz.
The legacy conversion now rejects missing/misaligned volumes and unexpected units.
Counts alone do not establish temporal alignment for legacy files: review step IDs before
scientific use. New dataset export verifies coordinate, force and cell step IDs explicitly.
The actual +/-0.1% isotropic strain test compares dE/dstrain with -trace(virial).

## HPC handoff

`workflow/hpc/hap-stage.sbatch` and `dependencies.sh` are templates, **not submitted or
cluster-tested**. HPC remains disabled. Fill only the actual site's scheduler, account,
partition/QoS, modules or Apptainer availability, image paths, CPU/GPU/RAM/walltime limits,
and shared/scratch paths. No credentials belong in this repository.

Pin the container digest, model hash, Git commit, inputs and configuration. Preserve failed
and interrupted runs. Slurm `afterok` proves process success only; independently validate
the predecessor's manifest, SCF, completeness and scientific gates before consuming data.
CP2K restarts require matching input/cell/coordinates, WFN and restart files. MACE checkpoints
here support a saved-geometry energy repeatability check; optimizer-state continuation is
not implemented. Do not claim optimizer resume. LAMMPS binary restart was actually tested.

Useful primary implementation references: [CP2K stress output](https://manual.cp2k.org/cp2k-2025_1-branch/CP2K_INPUT/FORCE_EVAL/PRINT/STRESS_TENSOR.html),
[MACE](https://github.com/ACEsuit/mace), [MatterGen](https://github.com/microsoft/mattergen).
