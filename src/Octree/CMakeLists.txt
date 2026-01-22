cmake_minimum_required(VERSION 3.20)
project(Octree LANGUAGES C)

# ------------------------------------------------------------------------------
# Core deps
# ------------------------------------------------------------------------------

if (PACIFIC_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

if (BUILD_THIRD_PARTY_HDF5)
else()
  set(HDF5_PREFER_PARALLEL TRUE)
  find_package(HDF5 REQUIRED COMPONENTS C HL)
endif()

include(GNUInstallDirs)
include(FetchContent)

# ------------------------------------------------------------------------------
# Basilisk / qcc provider selection
#
#   SYSTEM : use qcc found on PATH
#   VENDORED : use user-provided Basilisk source dir via -DBASILISK_DIR=/path/to/basilisk
#   GIT    : fetch Basilisk from git
# ------------------------------------------------------------------------------
 
set(OCTREE_BASILISK_PROVIDER "GIT" CACHE STRING "How to provide Basilisk/qcc: SYSTEM, VENDORED, GIT")
set_property(CACHE OCTREE_BASILISK_PROVIDER PROPERTY STRINGS SYSTEM VENDORED GIT)


# What the wrapper should execute. Default is plain 'qcc' (works if we install qcc to bindir).
set(QCC_PROGRAM "qcc")

if(OCTREE_BASILISK_PROVIDER STREQUAL "SYSTEM")
  find_program(QCC_PROGRAM qcc REQUIRED)
  message(STATUS "Using system qcc: ${QCC_PROGRAM}")

elseif(OCTREE_BASILISK_PROVIDER STREQUAL "VENDORED")
elseif(OCTREE_BASILISK_PROVIDER STREQUAL "VENDORED")

  # Define the cache variable ONCE
  set(BASILISK_DIR "" CACHE PATH
      "Path to a Basilisk source tree (must contain CMakeLists.txt)")

  # Compute effective directory
  set(_basilisk_dir "${BASILISK_DIR}")

  if(USE_SUBMODULES AND _basilisk_dir STREQUAL "")
    set(_basilisk_dir
        "${CMAKE_SOURCE_DIR}/third_party/basilisk/basilisk_submodule")
    message(STATUS "USE_SUBMODULES=ON: using Basilisk submodule at ${_basilisk_dir}")
  endif()

  if(_basilisk_dir STREQUAL "")
    message(FATAL_ERROR
      "OCTREE_BASILISK_PROVIDER=VENDORED requires:\n"
      "  -DBASILISK_DIR=/path/to/basilisk\n"
      "or\n"
      "  -DUSE_SUBMODULES=ON (expects submodule checkout)")
  endif()

  if(NOT EXISTS "${_basilisk_dir}/CMakeLists.txt")
    message(FATAL_ERROR
      "Basilisk dir '${_basilisk_dir}' is missing CMakeLists.txt.\n"
      "If using submodules, run:\n"
      "  git submodule update --init --recursive")
  endif()

  # Correct FetchContent override
  set(FETCHCONTENT_SOURCE_DIR_BASILISK "${_basilisk_dir}" CACHE PATH "" FORCE)
  set(FETCHCONTENT_UPDATES_DISCONNECTED_BASILISK ON CACHE BOOL "" FORCE)
  mark_as_advanced(
    FETCHCONTENT_SOURCE_DIR_BASILISK
    FETCHCONTENT_UPDATES_DISCONNECTED_BASILISK
  )

  FetchContent_Declare(Basilisk)
  FetchContent_MakeAvailable(Basilisk)

  if(TARGET qcc-toolchain)
    install(TARGETS qcc-toolchain RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
  endif()

elseif(OCTREE_BASILISK_PROVIDER STREQUAL "GIT")
  FetchContent_Declare(
    Basilisk
    GIT_REPOSITORY https://github.com/tortotubus/pacific-basilisk.git
  )
  FetchContent_MakeAvailable(Basilisk)

  mark_as_advanced(
    FETCHCONTENT_SOURCE_DIR_BASILISK
    FETCHCONTENT_UPDATES_DISCONNECTED_BASILISK
  )

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
mark_as_advanced(OCTREE_DLMFD_DIR OCTREE_VTKHYPERTREE_DIR)

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
