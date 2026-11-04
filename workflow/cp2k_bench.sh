#!/bin/bash
# CP2K GPU vs CPU timing on the protocol-test inputs (same hap-cp2k binary for every case).
# Usage (Git Bash on Shawn-PC): bash workflow/cp2k_bench.sh D:/path/to/cp2k-protocol <case>...
# Each case runs in a fresh container; results go to <protocol_dir>/bench/<layout>__<case>/.
set -u
IMAGE=${IMAGE:-hap-cp2k:2026.2-cuda-sm89}
DIR=$1; shift
LAYOUTS=${LAYOUTS:-"cpu_24x1 cpu_16x2 gpu_4x8 gpu_8x4"}
for name in "$@"; do
  for layout in $LAYOUTS; do
    kind=${layout%%_*}; ranks=${layout#*_}; np=${ranks%x*}; omp=${ranks#*x}
    out="bench/${TAG:-}${layout}__${name}"
    gpu=""; [[ $kind == gpu ]] && gpu="--gpus all"
    t0=$(date +%s)
    MSYS_NO_PATHCONV=1 docker run --rm $gpu --shm-size=4g -v "$DIR:/data" -w /data $IMAGE \
      bash -c "mkdir -p $out && cp $name/input.inp $out/ && cd $out && \
               OMP_NUM_THREADS=$omp mpirun -np $np cp2k.psmp -i input.inp -o output.out > run.log 2>&1"
    t1=$(date +%s)
    e=$(grep -m1 "Total FORCE_EVAL" "$DIR/$out/output.out" | awk '{print $NF}')
    echo "$layout $name wall_s=$((t1 - t0)) energy=$e"
  done
done
