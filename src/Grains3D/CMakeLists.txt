cmake_minimum_required(VERSION 3.12.0...4.0.0)
message(STATUS "CMake version: ${CMAKE_VERSION}")

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force the use of the chosen C++ standard.")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Enable C++ standard extensions.")

project(Grains
  VERSION 1.0.0
  LANGUAGES CXX
)

# ================================================================================
#
# USER PROVIDED DEPENDENCIES
#
# ================================================================================

find_package(MPI REQUIRED COMPONENTS C)
find_package(ZLIB REQUIRED)
find_package(XercesC REQUIRED)


# ================================================================================
#
# Source tree targets/projects
#
# ================================================================================

add_subdirectory(Grains)
add_subdirectory(Main)
# add_subdirectory(Tools)

option(GRAINS_BUILD_TESTS "Build Grains tests" OFF)

if(GRAINS_BUILD_TESTS)
  set(GRAINS_SOURCE_ROOT_DIR "${PROJECT_SOURCE_DIR}")
  set(GRAINS_BINARY_ROOT_DIR "${PROJECT_BINARY_DIR}")
  enable_testing()
  add_subdirectory(Tests)
endif()

# ================================================================================
#
# Install 
#
# ================================================================================

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ---- 1. Install the shared library ------------------------------------------------
install(TARGETS Grains Grains_main
  EXPORT GrainsTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}      # libgrains.so
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}      # static lib (if any)
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}      # Windows DLLs / executables
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} # for consumers
)

# ---- 2. Install all public headers -------------------------------------------------
# If your project has multiple include roots, add them all.
install(DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/Base/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/Component/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/Geometry/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/Fluid/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/PostProcessing/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/TimeIntegration/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/XML/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/InterfaceCBasilisk/include/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/Grains
  FILES_MATCHING PATTERN "*.h*"      # match .h, .hh, .hpp
)

# ---- 3. Install extra data/config files (optional) --------------------------------
if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/Grains/Dtd")
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Grains/Dtd/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/Grains/Dtd
    FILES_MATCHING PATTERN "*.dtd"
  )
endif()

if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/Grains/Tools")
  install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Grains/Tools/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/Grains/Tools
    FILES_MATCHING PATTERN "*.exec"
  )
endif()

# ---- 4. Generate CMake package configuration --------------------------------------
# This makes your project discoverable via find_package(Grains CONFIG)
set(GRAINS_INSTALL_CMAKEDIR ${CMAKE_INSTALL_LIBDIR}/cmake/Grains)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/GrainsConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/GrainsConfig.cmake
  INSTALL_DESTINATION ${GRAINS_INSTALL_CMAKEDIR}
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/GrainsConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# ---- 5. Install the package + targets ---------------------------------------------
install(EXPORT GrainsTargets
  FILE GrainsTargets.cmake
  NAMESPACE PacIFiC::
  DESTINATION ${GRAINS_INSTALL_CMAKEDIR}
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/GrainsConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/GrainsConfigVersion.cmake
  DESTINATION ${GRAINS_INSTALL_CMAKEDIR}
)

# ---- 6. Export from the build tree for development -------------------------------
export(EXPORT GrainsTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/GrainsTargets.cmake
  NAMESPACE PacIFiC::
)



