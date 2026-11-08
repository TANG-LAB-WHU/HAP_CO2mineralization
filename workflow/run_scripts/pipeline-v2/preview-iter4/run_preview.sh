#!/bin/bash
# 2026-11-08: stage-7 gates (dp_validate.py, unchanged) on the iteration-4 ensemble while the GPU waits for label-5,
# the same early look as preview-iter3 (200 ps for every group; the formal w8 run is 500 ps). Not the stage-7 result.
# Groups as arguments, S3 first: explore-4 NPT moved several S3 cells 3-9 % in volume (s3_055 from 6- to about
# 9-coordinated Ca within 25 ps); this repeats that over 200 ps with the 0 K relaxation and NVE gates. w8 follows only
# if label-5 leaves the GPU free long enough.
set -e
ROOT=/d/project/HAP_CO2mineralization
cd "$ROOT"
RUN=.research-runs/pipeline-v2
P=$RUN/preview-iter4
MODELS=$(for k in 0 1 2 3; do echo /work/$RUN/iter-004/train/model-$k/graph-compress.pth; done)
for g in "$@"; do
  mkdir -p $P/$g
  python -I 0.Generative_Discovery/screening/cp2k_hull.py cells --out $RUN/dft-hull --match $RUN/exploration/$g/*_initial.extxyz --write $P/$g/reference_cells.extxyz
  S=$(for f in $RUN/exploration/$g/*_initial.extxyz; do echo /work/$f; done)
  MSYS_NO_PATHCONV=1 docker run --rm --name preview-iter4-$g --gpus all -e OMP_NUM_THREADS=2 \
    -e PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True -v "$ROOT:/work" -w /work hap-mlip:2026.10 \
    python /work/3.ML_AIMD/Step3_mlip_deepmd/dp_validate.py --models $MODELS --structures $S \
      --reference /work/$P/$g/reference_cells.extxyz --md-ps 200 --parallel 3 --out /work/$P/$g
done
echo PREVIEW_FINISHED
