source ../Env/PacIFiC-CI-RUNNER-OpenMPI-2.1.1-GNU-8.2.1.env.sh

# Compilation of MacWorld
cd $MACWORLD_ROOT/MAC
./install-mac.sh

# Check if compilation was successful
find ${MAC_HOME}/lib/Linux-${MAC_FULL_EXT}/ -name "libmac0.so"
find ${MAC_HOME}/lib/Linux-${MAC_FULL_EXT}/ -name "libmac2.so"
