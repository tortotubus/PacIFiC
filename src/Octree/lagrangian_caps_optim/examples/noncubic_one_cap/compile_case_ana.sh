
CC99='mpicc -std=c99' qcc -events -Wall -O2 -D_MPI=1 \
	-D_DEFAULT_SOURCE \
	-DEBUG=1 \
	-DDISPLAY=1 \
        -DMTRACE=2 -DTRACE=2 \
        $1 -o "${1%.*}" \
        -I${BASILISK} -I${OCTREE_HOME} \
        -L${BASILISK}/gl -L${BASILISK}/wsServer \
        -grid=multigrid3D \
        -lglutils -lfb_tiny -lws -lm -Wno-unused-result

chmod +x "${1%.*}"
mkdir -p Res
mkdir -p Savings
