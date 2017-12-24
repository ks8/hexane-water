#!/bin/bash
#SBATCH --job-name=interface
#SBATCH --partition=depablo-gpu
#SBATCH --gres=gpu:1
#SBATCH --output=output
#SBATCH --error=errors
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time 24:00:00
# SBATCH --exclude=midway[230-232]

# LOAD NECESSARY MODULES
module load boost/1.62.0+openmpi-1.6+gcc-4.7
module load python
module load cuda/8.0
#module load cuda/7.5

# SET REQUIRED PATHS
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/swansonk1/OLD_DASH/md_engine/build/


# CHANGE INTO WORKING DIRECTORIES
echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"

#echo "Hello WORLD!"
python test.py -water_restart_file=kirk_tip4p_restart1950000.xml -hexane_lammps_file=hexane_restart.txt

