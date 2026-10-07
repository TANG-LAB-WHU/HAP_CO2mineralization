#!/bin/bash
# Print a dependency plan only. This script never submits jobs or connects remotely.
set -euo pipefail
cat <<'PLAN'
After filling account/partition/modules or pinned Apptainer images:
  preparation -> MACE matched references -> independent relaxation/geometry gates
  -> CP2K matched references + cutoff/basis/slab/vacuum checks
  -> independent AIMD replicas -> frame-aligned dataset export
  -> disjoint trajectory/system splits -> training -> held-out validation
  -> LAMMPS scale-up -> analysis

Slurm pattern for an operator-authorized submission:
  first=$(sbatch --parsable --account=ACCOUNT --partition=PARTITION hap-stage.sbatch)
  next=$(sbatch --parsable --dependency=afterok:$first hap-stage.sbatch)

afterok checks process status only. Each dependent stage must also validate the
upstream manifest and its scientific gates. Failed/partial outputs are retained.
No inference/training/MD result becomes manuscript-eligible automatically.
PLAN
