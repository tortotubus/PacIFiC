# ================================================================================
#
# Top level project settings
#
# ================================================================================

cmake_minimum_required(VERSION 3.12.0...4.0.0)
message(STATUS "CMake version: ${CMAKE_VERSION}")

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force the use of the chosen C++ standard.")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Enable C++ standard extensions.")

project(Pacific 
  VERSION 1.0.0
  LANGUAGES C CXX Fortran
)

# set(CMAKE_C_VISIBILITY_PRESET default
#     CACHE STRING "" FORCE)
# set(CMAKE_CXX_VISIBILITY_PRESET default
#     CACHE STRING "" FORCE)
# set(CMAKE_VISIBILITY_INLINES_HIDDEN OFF
#     CACHE BOOL "" FORCE)



# ================================================================================
#
# USER PROVIDED DEPENDENCIES
#
# ================================================================================

#
# Require MPI
#

find_package(MPI REQUIRED COMPONENTS C CXX Fortran)

#
# Require ZLib
#

find_package(ZLIB REQUIRED)

#
# Require HDF5
#

find_package(HDF5 REQUIRED)

#
# Require BLAS/LAPACK
#

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

set(BLAS_LAPACK_LIBS
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
)

if (DEFINED BLAS_INCLUDE_DIR)      # some FindBLAS modules define this
  set(BLAS_LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIR})
else()
  find_path(BLAS_LAPACK_INCLUDE_DIRS
    NAMES   cblas.h
    HINTS   ${BLAS_VENDOR_ROOT}/include            # e.g. /opt/openblas
            ${LAPACK_VENDOR_ROOT}/include          # if different
            /usr/include                           # fallback
  )
  set(BLAS_LAPACK_INCLUDE_DIRS "${BLAS_LAPACK_INCLUDE_DIRS}" CACHE INTERNAL "BLAS+LAPACK include directories")
  if (NOT BLAS_LAPACK_INCLUDE_DIRS)
    message(WARNING
      "Could not locate cblas.h; you may need to set BLAS_LAPACK_INCLUDE_DIRS by hand"
    )
  endif()
  
endif()

message(STATUS "Found BLAS Includes: ${BLAS_LAPACK_INCLUDE_DIRS}")
# ================================================================================
#
# External project declarations
#
# ================================================================================

include(ExternalProject)

# ================================================================================
#
# Add Xerces-C 2.8.0 as an external project; 
# build this in ${CMAKE_BINARY_DIR}/_ext/xerces-c/
#
# ================================================================================

set(XERCESCROOT ${CMAKE_BINARY_DIR}/_ext/xerces-c/src/xerces-c)

# force CMake to create these paths on project configuration:
file(MAKE_DIRECTORY "${XERCESCROOT}/include")
file(MAKE_DIRECTORY "${XERCESCROOT}/lib")

# External project declaration
ExternalProject_Add(
  xerces-c
  PREFIX ${CMAKE_BINARY_DIR}/_ext/xerces-c
  URL https://archive.apache.org/dist/xerces/c/2/sources/xerces-c-src_2_8_0.tar.gz

  # 1. No CMake “configure” step – we do our own:
  CONFIGURE_COMMAND  
    ${CMAKE_COMMAND} -E env XERCESCROOT=${XERCESCROOT}
    ${CMAKE_COMMAND} -E chdir ${XERCESCROOT}/src/xercesc
    ./runConfigure -p linux -x g++ -b 64

  # 2. The build step is just 'make' in that same directory:
  BUILD_COMMAND
    ${CMAKE_COMMAND} -E env XERCESCROOT=${XERCESCROOT}
    ${CMAKE_COMMAND} -E chdir ${XERCESCROOT}/src/xercesc
    make

  INSTALL_COMMAND ""

  BUILD_BYPRODUCTS
    ${XERCESCROOT}/lib/libxerces-c.so
    ${XERCESCROOT}/lib/libxerces-depdom.so
)

add_library(xerces-c_lib SHARED IMPORTED GLOBAL)

set_target_properties(xerces-c_lib PROPERTIES
  IMPORTED_LOCATION              "${CMAKE_BINARY_DIR}/_ext/xerces-c/src/xerces-c/lib/libxerces-c.so"
  INTERFACE_SONAME               "libxerces-c.so"
  INTERFACE_INCLUDE_DIRECTORIES  # only active during build
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_ext/xerces-c/src/xerces-c/include>
)

add_dependencies(xerces-c_lib xerces-c)

add_library(xerces-depdom_lib SHARED IMPORTED GLOBAL)

set_target_properties(xerces-depdom_lib PROPERTIES
  IMPORTED_LOCATION              "${CMAKE_BINARY_DIR}/_ext/xerces-c/src/xerces-c/lib/libxerces-depdom.so"
  INTERFACE_SONAME               "libxerces-depdom.so"
  INTERFACE_INCLUDE_DIRECTORIES
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_ext/xerces-c/src/xerces-c/include>
)

# ================================================================================
#
# Add Hypre 2.32.0 as an external project; 
# build this in ${CMAKE_BINARY_DIR}/_ext/hypre
#
# ================================================================================

set(HYPRE_PREFIX "${CMAKE_BINARY_DIR}/_ext/hypre")
set(HYPRE_INSTALL_DIR "${HYPRE_PREFIX}/install")
set(HYPRE_LIBRARY_FILES "${HYPRE_INSTALL_DIR}/lib64/libHYPRE${CMAKE_SHARED_LIBRARY_SUFFIX}")
set(HYPRE_INCLUDE_DIR "${HYPRE_INSTALL_DIR}/include/")

# Define the external project
ExternalProject_Add(
  hypre
  PREFIX ${HYPRE_PREFIX}

  INSTALL_DIR ${HYPRE_PREFIX}/install

  SOURCE_SUBDIR "src"
  URL    https://github.com/hypre-space/hypre/archive/refs/tags/v2.32.0.zip

  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DHYPRE_ENABLE_SHARED=ON
)

# ================================================================================
#
# Add PETSc 3.23.4 as an external project; 
# build this in ${CMAKE_BINARY_DIR}/_ext/petsc
#
# ================================================================================

find_package(Python3 REQUIRED COMPONENTS Interpreter)
if (NOT Python3_Interpreter_FOUND)
  message(FATAL_ERROR "Python 3 interpreter not found; PETSc configure.py needs Python")
endif()

set(PETSC_PREFIX "${CMAKE_BINARY_DIR}/_ext/petsc")
set(PETSC_ARCH "arch-cmake-generated")

# Define external project for PETSc
ExternalProject_Add(
  petsc
  DEPENDS hypre
  PREFIX ${PETSC_PREFIX}
  INSTALL_DIR "${PETSC_PREFIX}/install"

  GIT_REPOSITORY  https://gitlab.com/petsc/petsc.git
  GIT_TAG         v3.23.4

  CONFIGURE_COMMAND
    ${CMAKE_COMMAND} -E env
      PETSC_DIR=<SOURCE_DIR>
      PETSC_ARCH=${PETSC_ARCH}
    ${CMAKE_COMMAND} -E chdir <SOURCE_DIR>
    ${Python3_EXECUTABLE} ./configure
      --prefix=<INSTALL_DIR>
      --with-cc=${MPI_C_COMPILER}
      --with-cxx=${MPI_CXX_COMPILER}
      --with-fc=${MPI_Fortran_COMPILER}
      --CFLAGS="-pthread"
      --CXXFLAGS="-pthread"
      --FFLAGS="-pthread"
      --with-debugging=0
      --with-mpi=1
      --with-mpi-compilers=1
      --with-mpi-include=${MPI_C_INCLUDE_PATH}
      --with-mpi-lib=${MPI_LIBS}
      --with-mpiexec=${MPIEXEC_EXECUTABLE}
      --with-blaslapack=1
      --with-blaslapack-lib=${BLAS_LAPACK_LIBS}
      --with-blaslapack-include=${BLAS_LAPACK_INCLUDE_DIRS}
      --with-hypre=1
      --with-hypre-lib=${HYPRE_LIBRARY_FILES}
      --with-hypre-include=${HYPRE_INCLUDE_DIR}

  BUILD_COMMAND
    make
      -C <SOURCE_DIR>
      PETSC_DIR=<SOURCE_DIR>
      PETSC_ARCH=${PETSC_ARCH}
      all 

  INSTALL_COMMAND 
    make 
      -C <SOURCE_DIR>
      PETSC_DIR=<SOURCE_DIR>
      PETSC_ARCH=${PETSC_ARCH}
      install

  UPDATE_COMMAND ""
  PATCH_COMMAND ""

  BUILD_BYPRODUCTS
    <INSTALL_DIR>/lib/libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX}

)

# force CMake to create these paths now:
file(MAKE_DIRECTORY "${PETSC_PREFIX}/install/include")
file(MAKE_DIRECTORY "${PETSC_PREFIX}/install/lib")


# Define shared library target that PETSc generates
add_library(libpetsc SHARED IMPORTED GLOBAL)

set_target_properties(libpetsc PROPERTIES
  IMPORTED_LOCATION              "${CMAKE_BINARY_DIR}/_ext/petsc/install/lib/libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX}"
  INTERFACE_SONAME               "libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX}"
  INTERFACE_INCLUDE_DIRECTORIES  # only active during build
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_ext/petsc/install/include>
)

# Define alias for the library
add_library(PETSC::PETSC ALIAS libpetsc)

# Make sure that when targets link to libpetsc, the external project is built first
add_dependencies(libpetsc petsc)


# ================================================================================
#
# Source tree targets/projects
#
# ================================================================================

#
# Add Grains 
#

add_subdirectory(Grains3D)

#
# Add Cartesian 
#

add_subdirectory(Cartesian)
