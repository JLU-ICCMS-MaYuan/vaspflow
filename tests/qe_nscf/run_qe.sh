#!/bin/sh
#SBATCH  --job-name=qe_job
#SBATCH  --output=log.out
#SBATCH  --error=log.err
#SBATCH  --partition=liuhanyu
#SBATCH  --nodes=1
#SBATCH  --ntasks=48
#SBATCH  --ntasks-per-node=48
#SBATCH  --cpus-per-task=1

source /work/home/mayuan/intel/oneapi/setvars.sh --force
ulimit -s unlimited

echo 'Job started at' `date` 
mpirun -np 48 pw.x < nscf.in > nscf.out
echo 'Job finished at' `date` 
