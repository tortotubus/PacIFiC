cd ${MACWORLD_ROOT}/MAC
./install-mac.sh
find "${MAC_HOME}/lib/Linux-${MAC_FULL_EXT}/" -name "libmac0.so"
find "${MAC_HOME}/lib/Linux-${MAC_FULL_EXT}/" -name "libmac2.so"
