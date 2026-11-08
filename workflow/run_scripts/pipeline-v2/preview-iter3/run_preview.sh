#!/bin/bash
# 2026-11-05: stage-7 gates (dp_validate.py) on the iteration-3 ensemble while the GPU waits for label-4, as an early
# look before the final ensemble; 200 ps for every group (the formal w8 run is 500 ps). Not the stage-7 result.
set -e
ROOT=/d/project/HAP_CO2mineralization
cd "$ROOT"
RUN=.research-runs/pipeline-v2
P=$RUN/preview-iter3
MODELS=$(for k in 0 1 2 3; do echo /work/$RUN/iter-003/train/model-$k/graph-compress.pth; done)
for g in w8 s3; do
  mkdir -p $P/$g
  python -I 0.Generative_Discovery/screening/cp2k_hull.py cells --out $RUN/dft-hull --match $RUN/exploration/$g/*_initial.extxyz --write $P/$g/reference_cells.extxyz
  S=$(for f in $RUN/exploration/$g/*_initial.extxyz; do echo /work/$f; done)
  MSYS_NO_PATHCONV=1 docker run --rm --name preview-iter3-$g --gpus all -e OMP_NUM_THREADS=2 \
    -e PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True -v "$ROOT:/work" -w /work hap-mlip:2026.10 \
    python /work/3.ML_AIMD/Step3_mlip_deepmd/dp_validate.py --models $MODELS --structures $S \
      --reference /work/$P/$g/reference_cells.extxyz --md-ps 200 --parallel 3 --out /work/$P/$g
done
echo PREVIEW_FINISHED
