# source ~/pacific-caps/Env/PacIFiC-OpenMPI-4.1.2-GNU-11.4.0_Ubuntu22.env.sh

PACIFIC_THIRDPARTY=${PACIFIC_BUILDDIR_ABS}/third_party
PACIFIC_OCTREE=${PACIFIC_SRCDIR_ABS}/Octree

CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -DINITIAL_FLOW_FIELD=0 \
        -DDISPLAY=1 \
        -DMTRACE=2 -DTRACE=2 \
        $1 -o "${1%.*}" \
	-I${PACIFIC_OCTREE} \
	-L${PACIFIC_THIRDPARTY}/lib64 \
	-lglutils -lfb_tiny -lwsserver -lm \
        -grid=multigrid3D 



chmod +x "${1%.*}"
mkdir -p Res
mkdir -p Savings
