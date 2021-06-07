#!/bin/bash

# 1) Run test case on four cores
./pacific ./InputFiles/time.mac 4

# 2) Assess sucess
echo "Success assessment output:"

# 2.i) Compare to reference file (this command's exit status returns 0 if the two files are identical)
for time in {1..2}
do
  for proc in {0..3}
  do
    cmp Res/saveT${time}_${proc}.vtr Reference_results/ref_saveT${time}_${proc}.vtr ; echo $? >> heat-equation-2D-mpi-success.txt
    output="ref_saveT${time}_${proc}.vtr: "
    output+=$(tail -n 1 heat-equation-2D-mpi-success.txt)
    echo $output
  done
done

# 2.ii) We run a command that fails if any of the above success assessments failed 
$(! grep -q "1" heat-equation-2D-mpi-success.txt)
