#!/bin/bash
#SBATCH -p skx-dev
#SBATCH -N 1 
#SBATCH -n 1
#SBATCH -t 0-2:00 
#SBATCH -o slurm.%N.%j.out
#SBATCH -e slurm.%N.%j.err


echo "Working" 
