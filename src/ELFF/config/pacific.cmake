# ELFF settings when ELFF is built as a PacIFiC subproject.
#
# PacIFiC owns global build policy and third-party dependency setup. This file
# maps those parent-level choices onto the ELFF variables consumed by
# src/ELFF/CMakeLists.txt, while keeping ELFF-specific features configurable.

set(ELFF_PRECISION "double" CACHE STRING
  "Floating-point precision to use for ELFF: single or double")
set_property(CACHE ELFF_PRECISION PROPERTY STRINGS single double)

option(ELFF_USE_EXCEPTIONS "Enable exceptions in ELFF" OFF)
option(ELFF_USE_OPENMP "Enable OpenMP in ELFF" OFF)
option(ELFF_USE_EIGEN3 "Enable Eigen3 support in ELFF" ON)
option(ELFF_BUILD_TESTS "Build ELFF tests as part of PacIFiC" OFF)
option(ELFF_BUILD_DOC "Build ELFF documentation as part of PacIFiC" OFF)
option(ELFF_USE_COVERAGE "Enable compiler coverage instrumentation for ELFF" OFF)

# Install directory mappings. ELFF/CMakeLists.txt already prefers the
# CMAKE_MPI_INSTALL_* variables when they are set.
if(DEFINED PACIFIC_MPI_INSTALL_LIBDIR)
  set(CMAKE_MPI_INSTALL_LIBDIR "${PACIFIC_MPI_INSTALL_LIBDIR}")
endif()

if(DEFINED PACIFIC_MPI_INSTALL_BINDIR)
  set(CMAKE_MPI_INSTALL_BINDIR "${PACIFIC_MPI_INSTALL_BINDIR}")
endif()

if(DEFINED PACIFIC_MPI_INSTALL_INCLUDEDIR)
  set(CMAKE_MPI_INSTALL_INCLUDEDIR "${PACIFIC_MPI_INSTALL_INCLUDEDIR}")
endif()

if(DEFINED PACIFIC_INSTALL_LIBDIR)
  set(CMAKE_INSTALL_LIBDIR "${PACIFIC_INSTALL_LIBDIR}")
endif()

if(DEFINED PACIFIC_INSTALL_BINDIR)
  set(CMAKE_INSTALL_BINDIR "${PACIFIC_INSTALL_BINDIR}")
endif()

if(DEFINED PACIFIC_INSTALL_INCLUDEDIR)
  set(CMAKE_INSTALL_INCLUDEDIR "${PACIFIC_INSTALL_INCLUDEDIR}")
endif()

# PacIFiC-level policy mappings.
if(DEFINED PACIFIC_ENABLE_MPI)
  set(ELFF_USE_MPI "${PACIFIC_ENABLE_MPI}")
else()
  set(ELFF_USE_MPI ON)
endif()

if(DEFINED BUILD_THIRD_PARTY_EIGEN3)
  set(ELFF_BUILD_THIRD_PARTY_EIGEN3 "${BUILD_THIRD_PARTY_EIGEN3}")
else()
  set(ELFF_BUILD_THIRD_PARTY_EIGEN3 OFF)
endif()

if(DEFINED BUILD_THIRD_PARTY_HDF5)
  set(ELFF_BUILD_THIRD_PARTY_HDF5 "${BUILD_THIRD_PARTY_HDF5}")
else()
  set(ELFF_BUILD_THIRD_PARTY_HDF5 OFF)
endif()


if(DEFINED USE_SUBMODULES)
  set(ELFF_USE_SUBMODULES "${USE_SUBMODULES}")
else()
  set(ELFF_USE_SUBMODULES OFF)
endif()

# HDF5 is a PacIFiC dependency in the current tree; there is no top-level
# PACIFIC_ENABLE_HDF5 switch to mirror.
set(ELFF_USE_HDF5 ON)

# PacIFiC owns third-party dependency setup. ELFF should consume dependencies
# from the parent build or from find_package(), not create its own dependency
# subtree while embedded. 
set(ELFF_BUILD_THIRD_PARTY_GTEST OFF)
set(ELFF_BUILD_THIRD_PARTY_BASILISK OFF)

# Avoid ELFF's Basilisk examples/subproject by default inside PacIFiC. PacIFiC's
# Basilisk/qcc setup is handled elsewhere, primarily by Octree.
option(ELFF_USE_BASILISK "Enable ELFF Basilisk examples/subproject in PacIFiC" OFF)
