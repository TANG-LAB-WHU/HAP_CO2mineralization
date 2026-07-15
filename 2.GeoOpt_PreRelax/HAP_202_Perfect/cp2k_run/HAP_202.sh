#!/bin/sh
#SBATCH -p cp4
#SBATCH -N 6
#SBATCH -n 324

# Set number of OpenMP threads to 1
export OMP_NUM_THREADS=1

# Load CP2K environment dependencies
source /fs2/software/cp2k/pkg/cp2k-gcc12.2/tools/toolchain/install/setup
source /fs2/software/cp2k/pkg/cp2k-gcc12.2/tools/toolchain/env.sh

# Add CP2K executable path to PATH
export PATH="/fs2/software/cp2k/pkg/cp2k-gcc12.2/exe/local:$PATH"

# Disable UCX TLS to avoid memory allocation issues
export UCX_TLS=^glex

# Run CP2K
yhrun -c $OMP_NUM_THREADS cp2k.psmp -i HAP_202_opt.inp -o HAP_202_opt.log