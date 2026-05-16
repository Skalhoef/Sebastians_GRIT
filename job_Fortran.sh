#!/bin/bash -l

# Set the allocation to be charged for this job
# not required if you have set a default allocation
#SBATCH -A naiss2025-1-28

# The name of the script is myjob
#SBATCH -J GRIT

# partition
#SBATCH -p shared

# 10 hours wall-clock time will be given to this job
#SBATCH -t 23:59:59

# Number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=8

# load modules
ml PDC/24.11
ml cray-fftw

### set number of OMP threads
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=false
export OMP_STACKSIZE=256M
ulimit -Ss unlimited

echo "Script initiated at `date` on `hostname`"

srun -n 32 grit.x > out.log

echo "Script finished at `date` on `hostname`"

