#!/bin/bash
#SBATCH --account=rrg-awachs
#SBATCH --job-name=Basilisk 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1

#SBATCH --time=00-00:20:00
#SBATCH --output log.%j
#SBATCH --error err.%j

#SBATCH --partition=debug


#SBATCH --array=1-1%1   # Run a job array, one job at a time.
# ---------------------------------------------------------------------
echo "Current working directory: `pwd`"
echo "Starting run at: `date`"
# ---------------------------------------------------------------------
echo ""
echo "Job Array ID / Job ID: $SLURM_ARRAY_JOB_ID / $SLURM_JOB_ID"
echo "This is job $SLURM_ARRAY_TASK_ID out of $SLURM_ARRAY_TASK_COUNT jobs."
echo ""
# ---------------------------------------------------------------------


# Number of cores
# ---------------
export nprocs=$SLURM_NTASKS


# USER DEFINED Working directory
# ------------------------------
local_dir__=$SLURM_SUBMIT_DIR
cd $local_dir__


# PacIFiC variables
# -----------------
source ${HOME}/NonCubicOctree2026/PacIFiC/env/trillium.env.sh


source /scinet/vast/etc/vastpreload-openmpi.bash # important for MPI-IO


# Job info
# --------
echo " "
echo "Job ID:" $SLURM_JOBID
echo "-----------------"
echo " "
echo "Node numbers" 
echo "------------"
echo $SLURM_JOB_NODELIST
echo " "
echo "Numbers of tasks per node" 
echo "-------------------------"
echo $SLURM_TASKS_PER_NODE
echo " "


# Run job
# -------
${PACIFIC_MPIRUN} -n $nprocs ./$__Basilisk_app__
