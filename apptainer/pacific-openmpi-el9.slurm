#!/bin/bash
#SBATCH --job-name=grains-container-mpi
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --time=00:10:00
#SBATCH --partition=compute
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

set -euo pipefail

IMAGE=/path/to/pacific-openmpi-el9.sif
APP=grains

cd "$SLURM_SUBMIT_DIR"
test -f Grains/Init/insert.xml

apptainer exec \
  --bind "$SLURM_SUBMIT_DIR:/work" \
  --pwd /work \
  "$IMAGE" \
  mpirun -np $SLURM_NTASKS \
    "$APP" Grains/Init/insert.xml
