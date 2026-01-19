#!/bin/sh
#SBATCH  --job-name=mayuan
#SBATCH  --output=log.out
#SBATCH  --error=log.err
#SBATCH  --partition=liuhanyu
#SBATCH  --nodes=1
#SBATCH  --ntasks=48
#SBATCH  --ntasks-per-node=48
#SBATCH  --cpus-per-task=1
#SBATCH  --exclude=node49,node98,node57

source /work/home/mayuan/intel/oneapi/setvars.sh --force
ulimit -s unlimited
export I_MPI_ADJUST_REDUCE=3
export MPIR_CVAR_COLL_ALIAS_CHECK=0

echo 'Job started at' `date` 
mpirun -np 48 vasp_std
echo 'Job finished at' `date` 
