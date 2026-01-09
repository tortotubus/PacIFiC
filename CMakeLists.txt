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

message(STATUS "Found BLAS: ${BLAS_LIBRARIES}")
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")

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

# set(HYPRE_PREFIX "${CMAKE_BINARY_DIR}/_ext/hypre")
# set(HYPRE_INSTALL_DIR "${HYPRE_PREFIX}/install")
# set(HYPRE_LIBRARY_FILES "${HYPRE_INSTALL_DIR}/lib64/libHYPRE${CMAKE_SHARED_LIBRARY_SUFFIX}")
# set(HYPRE_INCLUDE_DIR "${HYPRE_INSTALL_DIR}/include/")

# # Define the external project
# ExternalProject_Add(
#   hypre
#   PREFIX ${HYPRE_PREFIX}

#   INSTALL_DIR ${HYPRE_PREFIX}/install

#   SOURCE_SUBDIR "src"
#   URL    https://github.com/hypre-space/hypre/archive/refs/tags/v2.32.0.zip

#   CMAKE_ARGS
#     -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
#     -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
#     -DHYPRE_ENABLE_SHARED=ON
# )

# ================================================================================
#
# Add PETSc 3.23.4 as an external project; 
# build this in ${CMAKE_BINARY_DIR}/_ext/petsc
#
# ================================================================================

# find_package(Python3 REQUIRED COMPONENTS Interpreter)
# if (NOT Python3_Interpreter_FOUND)
#   message(FATAL_ERROR "Python 3 interpreter not found; PETSc configure.py needs Python")
# endif()

# set(PETSC_PREFIX "${CMAKE_BINARY_DIR}/_ext/petsc")
# set(PETSC_ARCH "arch-cmake-generated")

# # Define external project for PETSc
# ExternalProject_Add(
#   petsc
#   DEPENDS hypre
#   PREFIX ${PETSC_PREFIX}
#   INSTALL_DIR "${PETSC_PREFIX}/install"

#   GIT_REPOSITORY  https://gitlab.com/petsc/petsc.git
#   GIT_TAG         v3.23.0

#   CONFIGURE_COMMAND
#     ${CMAKE_COMMAND} -E env
#       #--unset=CUDA_HOME
#       #--unset=CUDA_ROOT
#       #--unset=COMPILER_PATH
#       #--unset=GCC_ROOT
#       #--unset=GMP_ROOT
#       #--unset=HDF5_ROOT
#       #--unset=LIBSZIP_ROOT
#       #--unset=LIBXML2_ROOT
#       #--unset=MPIF77
#       #--unset=MPIF90
#       #--unset=MPICC
#       #--unset=MPICXX
#       #--unset=OPENBLAS_ROOT
#       #--unset=OPENMPI_ROOT
#       #--unset=ZLIB_ROOT
#       #--unset=ZLIB_NG_ROOT
#       PETSC_DIR=<SOURCE_DIR>
#       PETSC_ARCH=${PETSC_ARCH}
#     ${CMAKE_COMMAND} -E chdir <SOURCE_DIR>
#     ${Python3_EXECUTABLE} ./configure
#       --prefix=<INSTALL_DIR>

#       --with-cc=${MPI_C_COMPILER}
#       --with-cxx=${MPI_CXX_COMPILER}
#       --with-fc=${MPI_Fortran_COMPILER}

#       --with-debugging=0

#       --CFLAGS="-pthread"
#       --CXXFLAGS="-pthread"
#       --FFLAGS="-pthread"

#       --with-blaslapack=1
#       --with-blas-lib=${BLAS_LIBRARIES}
#       --with-lapack-lib=${LAPACK_LIBRARIES}

#       --with-mpi=1
#       --with-mpi-compilers=1
#       --with-mpi-include=${MPI_C_INCLUDE_PATH}
#       --with-mpi-lib=${MPI_LIBS}
#       --with-mpiexec=${MPIEXEC_EXECUTABLE}

#       --with-hypre=1
#       --with-hypre-lib=${HYPRE_LIBRARY_FILES}
#       --with-hypre-include=${HYPRE_INCLUDE_DIR}
#   BUILD_COMMAND
#     make
#       -C <SOURCE_DIR>
#       PETSC_DIR=<SOURCE_DIR>
#       PETSC_ARCH=${PETSC_ARCH}
#       all 

#   INSTALL_COMMAND 
#     make 
#       -C <SOURCE_DIR>
#       PETSC_DIR=<SOURCE_DIR>
#       PETSC_ARCH=${PETSC_ARCH}
#       install

#   UPDATE_COMMAND ""
#   PATCH_COMMAND ""

#   BUILD_BYPRODUCTS
#     <INSTALL_DIR>/lib/libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX}

# )

# # force CMake to create these paths now:
# file(MAKE_DIRECTORY "${PETSC_PREFIX}/install/include")
# file(MAKE_DIRECTORY "${PETSC_PREFIX}/install/lib")


# # Define shared library target that PETSc generates
# add_library(libpetsc SHARED IMPORTED GLOBAL)

# set_target_properties(libpetsc PROPERTIES
#   IMPORTED_LOCATION              "${CMAKE_BINARY_DIR}/_ext/petsc/install/lib/libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX}"
#   INTERFACE_SONAME               "libpetsc${CMAKE_SHARED_LIBRARY_SUFFIX}"
#   INTERFACE_INCLUDE_DIRECTORIES  # only active during build
#     $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/_ext/petsc/install/include>
# )

# # Define alias for the library
# add_library(PETSC::PETSC ALIAS libpetsc)

# # Make sure that when targets link to libpetsc, the external project is built first
# add_dependencies(libpetsc petsc)

# ================================================================================
#
# Add Basilisk as an external project; 
# build this in ${CMAKE_BINARY_DIR}/_ext/basilisk
#
# ================================================================================

# ExternalProject_Add(basilisk
# 	GIT_REPOSITORY        https://github.com/tortotubus/beam-basilisk
#   UPDATE_DISCONNECTED   1
#   PREFIX ${CMAKE_BINARY_DIR}/_ext/basilisk/
# 	# SOURCE_DIR            ${CMAKE_BINARY_DIR}/_ext/basilisk/src/
# 	# BINARY_DIR            ${CMAKE_BINARY_DIR}/_ext/basilisk/bin/
# 	INSTALL_DIR           ${CMAKE_BINARY_DIR}/_ext/basilisk/install/

# 	CMAKE_ARGS
# 	-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
# 	-DCMAKE_BUILD_TYPE=Debug

# 	BUILD_COMMAND 
# 	${CMAKE_COMMAND}
# 	--build <BINARY_DIR>
# 	--config Debug
# 	--target qcc

# 	INSTALL_COMMAND
# 	${CMAKE_COMMAND}
# 	--install <BINARY_DIR>
# 	--prefix <INSTALL_DIR>
# 	--config ${CMAKE_BUILD_TYPE}
# 	--component QCC
# )

# ExternalProject_Get_Property(basilisk install_dir)
# ExternalProject_Get_Property(basilisk binary_dir)

# set(QCC_EXECUTABLE "${install_dir}/bin/qcc" CACHE FILEPATH "Path to built qcc")


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

# add_subdirectory(Cartesian)

#
# Add Octree
#

# add_subdirectory(Octree)

#
# Add tests
#

enable_testing()
add_subdirectory(gtest)