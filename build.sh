sed -i "s|TO_BE_OVERWRITTEN|$(echo $PWD)|" Env/PacIFiC-CI-RUNNER-OpenMPI-2.1.1-GNU-8.2.1.env.sh
sed -i 's|${GRAINS_HOME}/XERCES-2.8.0|/home/gitlab-runner/dependencies/XERCES-2.8.0|' GRAINS/Env/grains_default.env.sh
sed -i 's|${GRAINS_HOME}/XERCES-2.8.0|/home/gitlab-runner/dependencies/XERCES-2.8.0|' GRAINS/Env/grains_default.env.csh
cd Env/
source PacIFiC-CI-RUNNER-OpenMPI-2.1.1-GNU-8.2.1.env.sh
cp ${PACIFIC_HOME}/GRAINS/Env/grains_env_template.env.sh ${PACIFIC_HOME}/GRAINS/Env/grains-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
cp ${PACIFIC_HOME}/Cartesian/MacWorld/Env/macworld_env_template.env.sh ${PACIFIC_HOME}/Cartesian/MacWorld/Env/macworld-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
cp ${PACIFIC_HOME}/Cartesian/FLUID/Env/fluid_env_template.env.sh ${PACIFIC_HOME}/Cartesian/FLUID/Env/fluid-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh
source PacIFiC-CI-RUNNER-${PACIFIC_MPI_DISTRIB}-${PACIFIC_MPI_VERSION}-${PACIFIC_SERCOMPIL_ENV}-${PACIFIC_SERCOMPIL_VERSION}.env.sh

cd $GRAINS_HOME
./makeARCH create ; make update ; make dtd
