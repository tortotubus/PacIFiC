cmake_minimum_required(VERSION 3.12.0...4.0.0)

include(CheckLanguage)
check_language(CUDA)

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force the use of the chosen C++ standard.")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Enable C++ standard extensions.")

project(Grains3D_GPU
  VERSION ${PACIFIC_VERSION_SEMVER}
  LANGUAGES CUDA CXX 
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
