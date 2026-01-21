#!/bin/bash
#SBATCH --job-name=pacific-mpi
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --partition=compute
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

module purge
module load mpi/openmpi

IMAGE=/path/to/pacific-openmpi-hybrid-el9.sif
APP=grains

cd "$SLURM_SUBMIT_DIR"

# Make sure the file exists in the submit directory
test -f Grains/Init/insert.xml

srun --mpi=pmix \
  apptainer exec \
    --bind "$SLURM_SUBMIT_DIR:/work" \
    --pwd /work \
    "$IMAGE" \
    "$APP" Grains/Init/insert.xml
