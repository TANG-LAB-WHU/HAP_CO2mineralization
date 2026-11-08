#!/bin/bash
# Protocol tests on the official CPU image, 24 MPI ranks x 1 thread (fastest layout on Shawn-PC).
cd /data
for d in cut400_rel60 cut500_rel60 cut600_rel40 cut600_rel60 cut600_rel80 cut800_rel60 cut1000_rel60 width12_2x2x2 width16_2x2x3; do
  [ -f $d/run/done ] && continue
  mkdir -p $d/run && cp $d/input.inp $d/run/ && cd $d/run
  OMP_NUM_THREADS=1 mpirun -np 24 cp2k.psmp -i input.inp -o output.out > run.log 2>&1 && grep -q "PROGRAM ENDED AT" output.out && touch done
  echo "$d $(grep -m1 'Total FORCE_EVAL' output.out | awk '{print $NF}')"
  cd /data
done
mkdir -p width8_1x1x2/run && cp cut600_rel60/run/output.out width8_1x1x2/run/ 2>/dev/null
echo PROTOCOL_DONE
