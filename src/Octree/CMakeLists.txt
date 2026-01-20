cmake_minimum_required(VERSION 3.20)
project(Octree LANGUAGES C)

# ------------------------------------------------------------------------------
# Core deps
# ------------------------------------------------------------------------------
set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(MPI REQUIRED)

include(GNUInstallDirs)
include(FetchContent)

# ------------------------------------------------------------------------------
# Basilisk / qcc provider selection
#
#   SYSTEM : use qcc found on PATH
#   VENDORED : use user-provided Basilisk source dir via -DBASILISK_DIR=/path/to/basilisk
#   GIT    : fetch Basilisk from git
# ------------------------------------------------------------------------------
set(OCTREE_BASILISK_PROVIDER "GIT" CACHE STRING "How to provide Basilisk/qcc: SYSTEM, VENDORED, or GIT")
set_property(CACHE OCTREE_BASILISK_PROVIDER PROPERTY STRINGS SYSTEM VENDORED GIT)

# User-provided Basilisk source directory for VENDORED mode
set(BASILISK_DIR "" CACHE PATH "Path to a Basilisk source tree (must contain CMakeLists.txt) used when OCTREE_BASILISK_PROVIDER=VENDORED")

# What the wrapper should execute. Default is plain 'qcc' (works if we install qcc to bindir).
set(QCC_PROGRAM "qcc")

if(OCTREE_BASILISK_PROVIDER STREQUAL "SYSTEM")
  find_program(QCC_PROGRAM qcc REQUIRED)
  message(STATUS "Using system qcc: ${QCC_PROGRAM}")

elseif(OCTREE_BASILISK_PROVIDER STREQUAL "VENDORED")
  if(BASILISK_DIR STREQUAL "")
    message(FATAL_ERROR
      "OCTREE_BASILISK_PROVIDER=VENDORED requires -DBASILISK_DIR=/path/to/basilisk")
  endif()
  if(NOT EXISTS "${BASILISK_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
      "BASILISK_DIR='${BASILISK_DIR}' does not look like a CMake Basilisk tree (missing CMakeLists.txt)")
  endif()

  # Tell FetchContent to use the user-provided source dir for Basilisk (no network).
  # NOTE: The variable name must match the declared content name: Basilisk -> BASILISK
  set(FETCHCONTENT_SOURCE_DIR_BASILISK "${BASILISK_DIR}" CACHE PATH "" FORCE)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_BASILISK ON CACHE BOOL "" FORCE)

  FetchContent_Declare(Basilisk)
  FetchContent_MakeAvailable(Basilisk)

  # If Basilisk provides a target called qcc, install it with us.
  if(TARGET qcc-toolchain)
    install(TARGETS qcc-toolchain RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
    message(STATUS "Vendored Basilisk provides target 'qcc-toolchain' (will be installed).")
  else()
    message(WARNING "Vendored Basilisk did not define a 'qcc-toolchain' target. "
                    "Your wrapper will call '${QCC_PROGRAM}' so qcc must be available at runtime.")
  endif()

elseif(OCTREE_BASILISK_PROVIDER STREQUAL "GIT")
  FetchContent_Declare(
    Basilisk
    GIT_REPOSITORY https://github.com/tortotubus/pacific-basilisk.git
    # Optional pin:
    # GIT_TAG v0.0.1
  )
  FetchContent_MakeAvailable(Basilisk)

  if(TARGET qcc-toolchain)
    install(TARGETS qcc-toolchain RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
    message(STATUS "Git Basilisk provides target 'qcc-toolchain' (will be installed).")
  else()
    message(WARNING "Git Basilisk did not define a 'qcc-toolchain' target. "
                    "Your wrapper will call '${QCC_PROGRAM}' so qcc must be available at runtime.")
  endif()

else()
  message(FATAL_ERROR
    "Unknown OCTREE_BASILISK_PROVIDER='${OCTREE_BASILISK_PROVIDER}'. Use SYSTEM, VENDORED, or GIT.")
endif()

# ------------------------------------------------------------------------------
# Octree headers
# ------------------------------------------------------------------------------
set(OCTREE_DLMFD_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/DLMFD"        CACHE PATH "Octree DLMFD directory")
set(OCTREE_VTKHYPERTREE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/VTKHyperTree" CACHE PATH "Octree VTKHyperTree directory")

install(
  DIRECTORY "${OCTREE_DLMFD_DIR}/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/Octree/DLMFD"
  COMPONENT Octree_DLMFD_Devel
  FILES_MATCHING PATTERN "*.h"
)

install(
  DIRECTORY "${OCTREE_VTKHYPERTREE_DIR}/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/Octree/VTKHyperTree"
  COMPONENT Octree_VTKHypertree_Devel
  FILES_MATCHING PATTERN "*.h"
)


include("${CMAKE_CURRENT_SOURCE_DIR}/cmake/BasiliskCTargets.cmake")

# ------------------------------------------------------------------------------
# qcc wrapper script install
# ------------------------------------------------------------------------------
# For SYSTEM, QCC_PROGRAM is an absolute path (from find_program).
# For VENDORED/GIT, we install 'qcc' into our bindir, so calling "qcc" is fine.
# set(QCC_INSTALL_LOCATION "${QCC_PROGRAM}")

# set(QCC_EXTRA_OPTS
#   "-lm "
#   "-I${CMAKE_INSTALL_FULL_INCLUDEDIR}/Octree/DLMFD "
#   "-I${CMAKE_INSTALL_FULL_INCLUDEDIR}/Octree/VTKHyperTree "
#   "-L${CMAKE_INSTALL_FULL_LIBDIR} "
#   "-Wl,-rpath,${CMAKE_INSTALL_FULL_LIBDIR}"
#   "-lGrains ${HDF5_C_LIBRARIES} ${HDF5_C_HL_LIBRARIES}"
# )

# configure_file(
#   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/qcc-pacific.in"
#   "${CMAKE_CURRENT_BINARY_DIR}/qcc-pacific"
#   @ONLY
# )

# install(
#   PROGRAMS "${CMAKE_CURRENT_BINARY_DIR}/qcc-pacific"
#   DESTINATION "${CMAKE_INSTALL_BINDIR}"
# )
