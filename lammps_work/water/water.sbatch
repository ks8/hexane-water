#!/bin/bash
#SBATCH --job-name=water
#SBATCH --time=24:00:00
#SBATCH --partition=sandyb
#SBATCH --nodes=8
#SBATCH --output=output
#SBATCH --error=errors

module load gcc
module load openmpi
mpirun -np 8 /home/swansonk1/lammps-31Mar17/src/lmp_mpi -in in.water



