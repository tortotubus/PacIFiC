#!/bin/bash
#SBATCH --account=rrg-awachs
#SBATCH -J L2r1R10C0.2
#SBATCH --nodes=16
##SBATCH --ntasks=512
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=3900mb
#SBATCH --time=3-00:00:00
#SBATCH --output basilisk.log.%j
#SBATCH -e basilisk.err.%j
#SBATCH --signal=B:INT@1800


# Number of cores
# ---------------
export nprocs=$SLURM_NTASKS


# USER DEFINED Working directory
# ------------------------------
local_dir__=$SLURM_SUBMIT_DIR

termination_handler(){
  echo ""
  echo "Function termination_handler called. Cleaning up, resubmitting and exiting."
  echo ""

# Log the termination event
  echo "$(date): Job was interrupted due to time limit. Resubmitting job." >> termination.log

  #ssh login01 "cd $SLURM_SUBMIT_DIR; sbatch --export=ALL basilisk-fir.sh"
  if [ $? -ne 0 ]; then
    echo "Failed to resubmit the job"
  fi

  # Exit the script
  exit 0
}

trap termination_handler INT


# PacIFiC variables
# -----------------
#source template.env.sh
## PacIFiC environment variables


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
num=$NUM
#num=$(($num+1))
if [ "$NUM" -eq 1 ]; then
   echo "Running a fresh simulation"
   echo "--------------------------"
   ####chmod +x ./bsuspension
   srun ./bsuspension &
   wait
else
   echo "Running a simulation using the last restart file"
   echo "------------------------------------------------"
   srun ./bsuspension &
   wait
fi

# End of compile and submit
