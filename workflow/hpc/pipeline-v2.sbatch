#!/bin/bash
#SBATCH --job-name=hap-pv2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=48:00:00
# EXAMPLE ONLY, not cluster-tested. One phase of the pipeline v2 orchestrator per job (workflow/hpc/PIPELINE_V2_HPC.md):
#   HAP_PHASE=cpu  CP2K labels, energy screens, relabels and rechecks; CPU node, --cpus-per-task = MPI ranks
#   HAP_PHASE=gpu  DeePMD training/tests, LAMMPS exploration, stage-7 validation; GPU node, e.g. --gpus=4
# The orchestrator stops before the first step of the other phase and writes it to $HAP_RUN/next_phase
# (cpu, gpu or done). Set account, partition, QoS, GPUs and walltime at submission for the actual cluster.
set -euo pipefail
: "${HAP_REPO:?absolute path of the committed repository checkout}"
: "${HAP_RUN:?absolute path of the transferred run directory (pipeline_state.json inside)}"
: "${HAP_PHASE:?cpu or gpu}"
: "${HAP_CP2K_SIF:?Apptainer image of cp2k/cp2k:2026.2_mpich_generic_psmp}"
: "${HAP_MLIP_SIF:?Apptainer image of hap-mlip:2026.10}"
test "$HAP_PHASE" = cpu || test "$HAP_PHASE" = gpu || { echo "HAP_PHASE must be cpu or gpu"; exit 2; }
test -f "$HAP_RUN/pipeline_state.json" || { echo "No pipeline_state.json in $HAP_RUN"; exit 2; }
# module load apptainer   # site specific
cd "$HAP_REPO"
test -z "$(git status --porcelain)" || { echo 'Clean committed source required'; exit 2; }
python3 -c 'import numpy, ase' || { echo 'Host python3 needs numpy and ase (pymatgen and scipy for stage 7)'; exit 2; }
source workflow/hpc/pipeline_v2_args.sh
hw=()
if [[ "$HAP_PHASE" == cpu ]]; then
  hw=(--ranks "${HAP_RANKS:-${SLURM_CPUS_PER_TASK:-32}}")
else
  # One model per allocated GPU (at most 4); explore-jobs LAMMPS runs share a GPU (3 on the local RTX 4090)
  ngpu=$(echo "${CUDA_VISIBLE_DEVICES:-0}" | tr ',' '\n' | grep -c .)
  hw=(--gpu-jobs "${HAP_GPU_JOBS:-$(( ngpu < 4 ? ngpu : 4 ))}" --explore-jobs "${HAP_EXPLORE_JOBS:-3}")
fi
cmd=(python3 workflow/pipeline_v2.py --run "$HAP_RUN" "${PIPELINE_V2_ARGS[@]}" "${hw[@]}"
     --runtime apptainer --cp2k-image "$HAP_CP2K_SIF" --mlip-image "$HAP_MLIP_SIF" --phase "$HAP_PHASE")
mkdir -p "$HAP_RUN/logs"
exec > >(tee -a "$HAP_RUN/logs/pipeline_v2-${SLURM_JOB_ID:-manual}.log") 2>&1
record=(python3 workflow/hpc/pipeline_v2_record.py)
"${record[@]}" start --run "$HAP_RUN" --phase "$HAP_PHASE" --args-file workflow/hpc/pipeline_v2_args.sh \
  --command "${cmd[*]}" --image "$HAP_CP2K_SIF" --image "$HAP_MLIP_SIF"
code=0
"${cmd[@]}" || code=$?
"${record[@]}" finish --run "$HAP_RUN" --returncode "$code"
echo "orchestrator exit $code; next phase: $(cat "$HAP_RUN/next_phase" 2>/dev/null || echo unknown)"
echo 'Process status only: review the step results against the PIPELINE_V2.md gates before relying on them.'
exit "$code"
