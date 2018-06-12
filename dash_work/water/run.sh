#!/bin/bash
#SBATCH --job-name=final-restart-test
#SBATCH --partition=depablo-gpu
#SBATCH --gres=gpu:1
#SBATCH --output=output-final-restart-test
#SBATCH --error=errors-final-restart-test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time 24:00:00

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
python tip4pF_UPDATE.py -equil_ensemble=NVT -nSteps_equilibration=3000 -prod_ensemble=NVT -nSteps_production=0 -numMolecules=216 -T=298.0 -P=1.0 -PI=True -nBeads=16 -dataFreq=100 -record_restart=False -restartFreq=1000 -record_traj=True -trajFreq1=100 -trajFreq2=50 -zlo_change=0.0 -zhi_change=0.0 -filename=final-restart-test -restart=True -restart_file=final-restart-test2000.xml
