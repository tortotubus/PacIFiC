include(ExternalProject)
include(GNUInstallDirs)

if (USE_SUBMODULES)
  ExternalProject_Add(eigen_ep
    UPDATE_COMMAND ""
    SOURCE_DIR "${CMAKE_SOURCE_DIR}/third_party/eigen/eigen_submodule" 
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/eigen-build"
    INSTALL_DIR "${EIGEN3_PREFIX}"
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
  )
else()
  ExternalProject_Add(eigen_ep
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG        3.4.1

    SOURCE_DIR  "${CMAKE_BINARY_DIR}/_deps/eigen-src"
    BINARY_DIR  "${CMAKE_BINARY_DIR}/_deps/eigen-build"
    INSTALL_DIR "${EIGEN3_PREFIX}"
  
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      
    UPDATE_DISCONNECTED 1
    BUILD_ALWAYS OFF
  )
endif()

file(MAKE_DIRECTORY "${EIGEN3_PREFIX}/include/eigen3")
add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
set_target_properties(Eigen3::Eigen PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_PREFIX}/include/eigen3"
)
add_dependencies(Eigen3::Eigen eigen_ep)

list(PREPEND CMAKE_PREFIX_PATH "${EIGEN3_PREFIX}")
