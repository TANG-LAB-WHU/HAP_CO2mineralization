#!/bin/bash
#SBATCH --job-name=hap-pv2-cp2k
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --array=0-7
# EXAMPLE ONLY, not cluster-tested. Optional: spreads the CP2K runs of one prepared label directory over an array;
# task k labels every N-th frame directory from k (run_labels.sh <shard> <shards>), skipping done and failed frames.
# Afterwards a HAP_PHASE=cpu pipeline-v2.sbatch job finds the frames labelled and only screens, relabels and
# collects. Never run it beside a cpu-phase job on the same directory: both would start the same frames.
# The directory must be prepared and carry the sharding run_labels.sh (workflow/hpc/PIPELINE_V2_HPC.md):
#   python3 workflow/cp2k_label.py scripts --label-dir "$HAP_LABEL_DIR"
set -euo pipefail
: "${HAP_REPO:?absolute path of the committed repository checkout}"
: "${HAP_RUN:?absolute path of the run directory (for the job manifest)}"
: "${HAP_LABEL_DIR:?prepared label directory, e.g. \$HAP_RUN/iter-005/labels}"
: "${HAP_CP2K_SIF:?Apptainer image of cp2k/cp2k:2026.2_mpich_generic_psmp}"
test -f "$HAP_LABEL_DIR/index.csv" || { echo "$HAP_LABEL_DIR is not prepared"; exit 2; }
grep -q SHARDS "$HAP_LABEL_DIR/run_labels.sh" || { echo "run_labels.sh predates sharding; run cp2k_label.py scripts"; exit 2; }
# module load apptainer   # site specific
cd "$HAP_REPO"
test -z "$(git status --porcelain)" || { echo 'Clean committed source required'; exit 2; }
shards=${HAP_SHARDS:-${SLURM_ARRAY_TASK_COUNT:-1}}
shard=$(( ${SLURM_ARRAY_TASK_ID:-0} - ${SLURM_ARRAY_TASK_MIN:-0} ))
ranks=${HAP_RANKS:-${SLURM_CPUS_PER_TASK:-32}}
# MPICH in the image would start its ranks with srun, which the image lacks; fork keeps them on this node
cmd=(apptainer exec --env HYDRA_LAUNCHER=fork,HYDRA_RMK=user --bind "$HAP_LABEL_DIR:/data" "$HAP_CP2K_SIF"
     /opt/cp2k/install/bin/launch bash /data/run_labels.sh /data "$ranks" 1 "$shard" "$shards")
mkdir -p "$HAP_RUN/logs"
exec > >(tee -a "$HAP_RUN/logs/cp2k-labels-${SLURM_JOB_ID:-manual}.log") 2>&1
record=(python3 workflow/hpc/pipeline_v2_record.py)
"${record[@]}" start --run "$HAP_RUN" --phase cpu --args-file workflow/hpc/pipeline_v2_args.sh \
  --command "${cmd[*]}" --image "$HAP_CP2K_SIF"
code=0
"${cmd[@]}" || code=$?
"${record[@]}" finish --run "$HAP_RUN" --returncode "$code"
echo "shard $shard/$shards exit $code"
exit "$code"
