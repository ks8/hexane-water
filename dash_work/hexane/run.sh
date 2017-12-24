#!/bin/bash

#SBATCH --job-name=hexane
#SBATCH --output=output
#SBATCH --error=errors
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --account=pi-depablo
#SBATCH --partition=depablo-gpu
#SBATCH --gres=gpu:1
#SBATCH --mail-user=swansonk1@uchicago.edu
PATH=$PATH:/home/swansonk1/OLD_DASH/md_engine
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/swansonk1/OLD_DASH/md_engine/build

module load cuda/8.0
module load boost/1.62.0+openmpi-1.6+gcc-4.7


python in.webb_hexane.py

