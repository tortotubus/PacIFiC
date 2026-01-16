
#
# Declare all dependencies for our Octree codebase
#

set(HDF5_PREFER_PARALLEL TRUE) 
find_package(HDF5 COMPONENTS C HL REQUIRED)
find_package(MPI REQUIRED)

#
# Pull in Basilisk as an ExternalProject and 'install' it to our build folder 
# so that it can be used for transpiling Basilisk C to ISO99 C as part of our
# toolchain
#

include(FetchContent)

if(DEFINED BASILISK_SRC)
  FetchContent_Declare(Basilisk
    SOURCE_DIR ${BASILISK_SRC}
  )
else()
  FetchContent_Declare(Basilisk
    GIT_REPOSITORY https://github.com/tortotubus/beam-basilisk.git
  )
endif()

FetchContent_MakeAvailable(Basilisk)


#
#
#
#
#

set(OCTREE_DLMFD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/DLMFD" CACHE PATH "Octree DLMFD directory" FORCE)
set(OCTREE_VTKHYPERTREE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/VTKHyperTree" CACHE PATH "Octree VTKHyperTree directory" FORCE)

#
#
#
#
#

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/BasiliskCTargets.cmake)


#
#
#
#
#
include(GNUInstallDirs)

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/DLMFD"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/Octree"
  FILES_MATCHING
  PATTERN "*.h"
)

install(
  DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/VTKHyperTree"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/Octree"
  FILES_MATCHING
  PATTERN "*.h"
)

set(QCC_INSTALL_LOCATION "qcc")

# Example extra flags to always add
set(QCC_EXTRA_OPTS "-lm -I${HDF5_INCLUDE_DIRS} -I${CMAKE_INSTALL_FULL_INCLUDEDIR}/Octree/DLMFD -I${CMAKE_INSTALL_FULL_INCLUDEDIR}/Octree/VTKHyperTree -L${CMAKE_INSTALL_FULL_LIBDIR} -Wl,-rpath,${CMAKE_INSTALL_FULL_LIBDIR} -lGrains ${HDF5_C_LIBRARIES} ${HDF5_C_HL_LIBRARIES}")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/qcc-pacific.in
  ${CMAKE_CURRENT_BINARY_DIR}/qcc-pacific
  @ONLY
)

install(
  PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/qcc-pacific"
  DESTINATION "${CMAKE_INSTALL_BINDIR}"
)

 