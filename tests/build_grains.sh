source ../Env/PacIFiC-CI-RUNNER-OpenMPI-2.1.1-GNU-8.2.1.env.sh

# Compilation of Grains3D
cd $GRAINS_HOME
./makeARCH create ; make update ; make dtd

# Check if the compilation was successful
find ${GRAINS_HOME}/Main/bin${GRAINS_FULL_EXT}/ -name "grains" >> build_success.txt
