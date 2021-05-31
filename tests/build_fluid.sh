source ../Env/PacIFiC-CI-RUNNER-OpenMPI-2.1.1-GNU-8.2.1.env.sh

# Compilation of FLUID
cd ${PACIFIC_HOME}/Cartesian/FLUID
./compil
find ${PACIFIC_HOME}/Cartesian/FLUID/lib/Linux-${MAC_FULL_EXT}/ -name "exe0"
find ${PACIFIC_HOME}/Cartesian/FLUID/lib/Linux-${MAC_FULL_EXT}/ -name "exe2"
