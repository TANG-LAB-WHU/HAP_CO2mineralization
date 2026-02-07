#!/bin/sh
#SBATCH -p cp4
#SBATCH -N 6
#SBATCH -n 324

export OMP_NUM_THREADS=1

source /fs2/software/cp2k/pkg/cp2k/tools/toolchain/install/setup

module load GCC/12.2.0 
module load MPI/mpich/4.0.2-mpi-x-gcc12.2 
module load scalapack/2.2.0-gcc12.2-mpi-x
module load openblas/0.3.28-gcc12.2 
module load dbcsr/2.8.0-gcc12.2-mpi-x

module load cp2k/2025.2-all-gcc12.2-mpich

yhrun -c $OMP_NUM_THREADS cp2k.psmp -i HAP_431_opt.inp