cmake_minimum_required(VERSION 3.12.0...4.0.0)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force the use of the chosen C++ standard.")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Enable C++ standard extensions.")

project(Grains
  VERSION ${PACIFIC_VERSION_SEMVER}
  LANGUAGES CXX
)

# ================================================================================
#
# USER PROVIDED DEPENDENCIES
#
# ================================================================================

if (PACIFIC_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS C)
endif()

if (BUILD_THIRD_PARTY_ZLIB)
  # find_package(ZLIB REQUIRED)
else()
  find_package(ZLIB REQUIRED)
endif()

if (BUILD_THIRD_PARTY_XERCESC)
  # find_package(XercesC REQUIRED)
else()
  find_package(XercesC REQUIRED)
endif()


# ================================================================================
#
# Source tree targets/projects
#
# ================================================================================

add_subdirectory(Grains)
add_subdirectory(Main) 

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
  LIBRARY DESTINATION ${PACIFIC_MPI_INSTALL_LIBDIR}      # libgrains.so
  COMPONENT Grains_Runtime
  RUNTIME DESTINATION ${PACIFIC_MPI_INSTALL_BINDIR}
  COMPONENT Grains_Runtime
  INCLUDES DESTINATION ${PACIFIC_MPI_INSTALL_INCLUDEDIR} # for consumers
)

# ---- 2. Install all public headers -------------------------------------------------
install(DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/Base/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/Component/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/Geometry/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/Fluid/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/PostProcessing/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/TimeIntegration/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/XML/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Grains/InterfaceCBasilisk/include/
  DESTINATION ${PACIFIC_MPI_INSTALL_INCLUDEDIR}/Grains
  COMPONENT Grains_Devel
  FILES_MATCHING PATTERN "*.h*" # match .h, .hh, .hpp
)

# ---- 3. Install data/config files --------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Grains/Dtd/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/Grains/Dtd
  COMPONENT Grains_Data
  FILES_MATCHING PATTERN "*.dtd"
)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Grains/Tools/
  DESTINATION ${CMAKE_INSTALL_DATADIR}/Grains/Tools
  COMPONENT Grains_Data
  FILES_MATCHING PATTERN "*.exec"
)

# ---- 4. Generate CMake package configuration --------------------------------------
# This makes the project discoverable via find_package(Grains CONFIG)
set(GRAINS_INSTALL_CMAKEDIR ${PACIFIC_MPI_INSTALL_LIBDIR}/cmake/Grains)

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
  COMPONENT Grains_Devel
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/GrainsConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/GrainsConfigVersion.cmake
  DESTINATION ${GRAINS_INSTALL_CMAKEDIR}
  COMPONENT Grains_Devel
)

# ---- 6. Export from the build tree for development -------------------------------
export(EXPORT GrainsTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/GrainsTargets.cmake
  NAMESPACE PacIFiC::
)



