#!/bin/bash
# Usage inside the CP2K image: run_job.sh <job_dir> <mpi_ranks>
cd "$1" || exit 1
inp=input.inp
# A stopped cell optimisation resumes from its last step
# (a restart file cut off by the stop falls back to the backup of the step before)
for f in opt-1.restart opt-1.restart.bak-1; do
  [ -f $f ] && grep -q "&END FORCE_EVAL" $f && cp $f resume.inp && inp=resume.inp && break
done
n=$(ls output_*.out 2>/dev/null | wc -l)
OMP_NUM_THREADS=1 mpirun -np $2 cp2k.psmp -i $inp -o output_$n.out && grep -q "PROGRAM ENDED AT" output_$n.out && touch done
