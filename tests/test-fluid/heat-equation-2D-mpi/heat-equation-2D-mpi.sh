#!/bin/bash

# Run test case on four cores
./pacific ./InputFiles/time.mac 4

# Assess success: compare to reference file (this command's exit status returns 0 if the two files are identical)
cd Res
for time in {1..2}
do
  for proc in {1..4}
  do
    cmp save${time}_${proc}.vtr ref_save${time}_${proc}.vtr >> heat-equation-2D-mpi-success.txt
  done
done

$(! grep -q "1" heat-equation-2D-mpi-success.txt)
