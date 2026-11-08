# Pipeline v2 science arguments for every orchestrator job (sourced by pipeline-v2.sbatch). They match the local run
# (logs of 2026-11-08 01:24), except --iters: 6 -> 8 (decided 2026-11-08), so explore-5 and explore-6 candidates are
# labelled and trained on before the final ensemble. Hardware arguments (--ranks, --gpu-jobs, --explore-jobs) are
# set per job from the Slurm allocation. Changing any value here changes the protocol: record it in PIPELINE_V2.md.
PIPELINE_V2_ARGS=(
  --cutoff 800 --rel-cutoff 60
  --iters 8
  --warm-steps 100000
  --npt-temperatures 300
  --extra-systems s3 metad
)
