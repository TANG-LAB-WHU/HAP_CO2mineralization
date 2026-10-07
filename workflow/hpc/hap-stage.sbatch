#!/bin/bash
#SBATCH --job-name=hap-stage
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
# Set partition/account/GPU requests at submission for the actual cluster.
set -euo pipefail
: "${HAP_REPO:?absolute shared repository path required}"
: "${HAP_RUN:?new absolute run directory required}"
: "${HAP_STAGE:?mace or cp2k required}"
: "${HAP_IMAGE:?pinned Apptainer image required}"
cd "$HAP_REPO"
test -z "$(git status --porcelain)" || { echo 'Clean committed source required'; exit 2; }
mkdir -p "$HAP_RUN"
exec > >(tee -a "$HAP_RUN/job-${SLURM_JOB_ID:-manual}.log") 2>&1
git rev-parse HEAD > "$HAP_RUN/source-${SLURM_JOB_ID:-manual}.txt"
sha256sum "$HAP_IMAGE" > "$HAP_RUN/image-${SLURM_JOB_ID:-manual}.sha256"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
trap 'echo "Interrupted: retain all output and restart files; do not mark complete"; exit 143' TERM INT
if [[ "$HAP_STAGE" == mace ]]; then
  : "${HAP_MODEL:?pinned model checkpoint required}"
  test ! -e "$HAP_RUN/mace.manifest.json" || { echo 'Use new run directory; preserve prior manifest'; exit 2; }
  sha256sum "$HAP_MODEL" > "$HAP_RUN/model.sha256"
  srun apptainer exec --nv --bind "$HAP_REPO,$HAP_RUN" "$HAP_IMAGE" python workflow/local_compute.py mace --out "$HAP_RUN" --model "$HAP_MODEL" --steps "${HAP_STEPS:-12}"
elif [[ "$HAP_STAGE" == cp2k ]]; then
  : "${HAP_INPUT:?explicit CP2K input or reviewed restart input required}"
  # Restart is explicit: pass reviewed *.restart with matching WFN/cell/coordinates.
  # Never infer completion from a trajectory or restart-file presence.
  sha256sum "$HAP_INPUT" > "$HAP_RUN/input-${SLURM_JOB_ID:-manual}.sha256"
  cd "$HAP_RUN"
  srun apptainer exec --nv --bind "$HAP_REPO,$HAP_RUN" "$HAP_IMAGE" cp2k -i "$HAP_INPUT" -o "cp2k-${SLURM_JOB_ID:-manual}.log"
else
  echo "Unsupported stage $HAP_STAGE"; exit 2
fi
echo 'Execution finished; independent convergence and evidence review remain required.'
