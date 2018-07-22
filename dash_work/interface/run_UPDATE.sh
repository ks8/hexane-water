#!/bin/bash
#SBATCH --job-name=full-255-32bead
#SBATCH --partition=depablo-gpu
#SBATCH --gres=gpu:1
#SBATCH --output=output-full-255-32bead
#SBATCH --error=errors-full-255-32bead
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time 48:00:00
# SBATCH --exclude=midway[230-232]

# LOAD NECESSARY MODULES
module load boost/1.62.0+openmpi-1.6+gcc-4.7
module load python
module load cuda/8.0
#module load cuda/7.5

# SET REQUIRED PATHS
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/swansonk1/NEW_DASH/md_engine/build/


# CHANGE INTO WORKING DIRECTORIES
echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"

#echo "Hello WORLD!"
python interface_UPDATE.py -water_restart_file=kirk_tip4p_restart1950000.xml -hexane_restart_file=hexane_restart.txt -hexane_settings_file=hexane.in.settings -equil_ensemble=NVT -nSteps_equilibration=1000000 -prod_ensemble=NVT -nSteps_production=6000000 -T=255.0 -PI=True -nBeads=32 -ptensorFreq=20 -mix=waldman -record_restart=True -restartFreq=10000 -record_traj=True -trajFreq=10000 -zlo_change=-20.0 -zhi_change=20.0 -filename=full-255-32bead -restart=False -restart_file=
