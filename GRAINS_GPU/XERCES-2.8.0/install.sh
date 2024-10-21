#!/bin/bash
echo "AA"
echo ${GRAINS_CPP_COMPILER_VERSION}
if [ -d lib${GRAINS_BITS_DEFAULT}-${GRAINS_SERCOMPIL_ENV}-${GRAINS_SERCOMPIL_VERSION} ]
  then
    echo 'Xerces is built for' ${GRAINS_SERCOMPIL_ENV}-${GRAINS_SERCOMPIL_VERSION}
  else
    echo 'Building Xerces for' ${GRAINS_SERCOMPIL_ENV}-${GRAINS_SERCOMPIL_VERSION}
    cd src/xercesc
    make clean
    ./runConfigure -plinux -c${GRAINS_C} -x${GRAINS_CPP} -b${GRAINS_BITS_DEFAULT} -minmem -nsocket -tnative -rpthread
    make
    cd ../../
    cp -r lib lib${GRAINS_BITS_DEFAULT}-${GRAINS_SERCOMPIL_ENV}-${GRAINS_SERCOMPIL_VERSION}
    cd lib${GRAINS_BITS_DEFAULT}-${GRAINS_SERCOMPIL_ENV}-${GRAINS_SERCOMPIL_VERSION}
    rm libxerces-c.so libxerces-c.so.28 libxerces-depdom.so libxerces-depdom.so.28
    ln -s libxerces-c.so.28.0 libxerces-c.so.28
    ln -s libxerces-c.so.28 libxerces-c.so
    ln -s libxerces-depdom.so.28.0 libxerces-depdom.so.28          
    ln -s libxerces-depdom.so.28 libxerces-depdom.so
    cd ../
    rm -rf lib
fi