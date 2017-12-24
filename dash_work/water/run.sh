#!/bin/bash
#SBATCH --job-name=water
#SBATCH --partition=depablo-gpu
#SBATCH --gres=gpu:1
#SBATCH --output=output
#SBATCH --error=errors
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time 24:00:00

# LOAD NECESSARY MODULES
module load boost/1.62.0+openmpi-1.6+gcc-4.7
module load python
module load cuda/8.0

# SET REQUIRED PATHS
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/swansonk1/OLD_DASH/md_engine/build/


python in.tip4pF.py

