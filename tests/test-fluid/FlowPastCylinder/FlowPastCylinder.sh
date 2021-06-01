#!/bin/bash

# Run test case on four cores
./pacific ./InputFiles/time.mac 4

# Assess success: compare to reference file (this command's exit status returns 0 if the two files are identical)
cmp -s DS_results/particle_forces.csv DS_results/particle_forces_reference.csv
SUCCESS=$?
if [ SUCCESS != "0" ]
then
  echo "Error: output file and reference file differ"
  echo "reference file:"
  cat DS_results/particle_forces_reference.csv
  echo "output file:"
  cat DS_results/particle_forces.csv
fi
echo "$SUCCESS" | grep -q "0"
