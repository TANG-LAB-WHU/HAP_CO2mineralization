#!/bin/bash
#SBATCH --job-name=hap-mlip
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
# EXAMPLE ONLY. Set the real account, partition and GPU request at submission.
set -euo pipefail
: "${HAP_REPO:?absolute committed repository path required}"
: "${HAP_RUN:?new attempt directory with reviewed staged inputs required}"
: "${HAP_INPUT:?explicit training JSON or LAMMPS input required}"
: "${HAP_IMAGE:?existing pinned Apptainer image required; no automatic download}"
: "${HAP_UPSTREAM_MANIFEST:?completed and scientifically validated predecessor required}"
: "${HAP_STAGE:?train or lammps required}"
test "$HAP_STAGE" = train || test "$HAP_STAGE" = lammps || exit 2
cd "$HAP_RUN"
record="$HAP_REPO/workflow/hpc/job_record.py"
stage=mlip_training
test "$HAP_STAGE" = train || stage=scaleup
args=(start --run "$HAP_RUN" --stage "$stage" --input "$HAP_INPUT" --image "$HAP_IMAGE" --upstream "$HAP_UPSTREAM_MANIFEST")
if [[ "$HAP_STAGE" == lammps ]]; then
  : "${HAP_MODEL:?exact trained model used by the LAMMPS input required}"
  args+=(--model "$HAP_MODEL")
fi
if [[ -n "${HAP_RESTART_FILE:-}" ]]; then args+=(--restart-file "$HAP_RESTART_FILE"); fi
python3 "$record" "${args[@]}"
log="$HAP_RUN/job-${SLURM_JOB_ID:-manual}.log"
exec > >(tee "$log") 2>&1
finish() {
  code=$?
  trap - EXIT
  outputs=(--output "$log")
  if [[ "$HAP_STAGE" == train ]]; then
    outputs+=(--output "$HAP_RUN/model.pth" --output "$HAP_RUN/lcurve.out")
  else
    outputs+=(--output "$HAP_RUN/lammps.log")
  fi
  # Preserve actual standard checkpoint/trajectory outputs when produced by the reviewed input.
  shopt -s nullglob
  for artifact in "$HAP_RUN"/*.dump "$HAP_RUN"/*.restart "$HAP_RUN"/model.ckpt*.pt; do
    outputs+=(--output "$artifact")
  done
  python3 "$record" finish --run "$HAP_RUN" --returncode "$code" "${outputs[@]}"
  exit "$code"
}
trap finish EXIT
trap 'exit 143' TERM INT
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-4}"
export DP_INTRA_OP_PARALLELISM_THREADS="$OMP_NUM_THREADS"
export DP_INTER_OP_PARALLELISM_THREADS=1
engine=(apptainer exec --nv --bind "$HAP_REPO,$HAP_RUN" "$HAP_IMAGE")
"${engine[@]}" python -c 'import torch; assert torch.cuda.is_available(), "Expected allocated CUDA GPU"; print(torch.__version__, torch.cuda.get_device_name(0))'
if [[ "$HAP_STAGE" == train ]]; then
  restart=()
  if [[ -n "${HAP_RESTART:-}" ]]; then
    : "${HAP_RESTART_FILE:?actual checkpoint file for checksum verification required}"
    : "${HAP_RESTART_SHA256:?reviewed checkpoint checksum required}"
    actual=$(sha256sum "$HAP_RESTART_FILE" | cut -d' ' -f1)
    test "$actual" = "$HAP_RESTART_SHA256" || { echo 'Checkpoint checksum mismatch'; exit 2; }
    restart=(--restart "$HAP_RESTART")
  fi
  srun "${engine[@]}" dp --pt train "$HAP_INPUT" "${restart[@]}"
  srun "${engine[@]}" dp --pt freeze -o model.pth
else
  # A restart is a separate reviewed input containing read_restart, matching model and types.
  # Never infer completion from a trajectory/restart file's existence.
  srun "${engine[@]}" lmp -in "$HAP_INPUT" -log lammps.log
fi
echo 'Engine finished. Complete independent validation before enabling the dependent stage.'
