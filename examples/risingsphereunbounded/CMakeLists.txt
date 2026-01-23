cmake_minimum_required(VERSION 3.14...4.0)
project(RisingSphereUnbounded)
enable_language(C)

# Needed for octree_add_executable
find_package(Octree CONFIG REQUIRED) 

set(GRAINS_DATA_SRC_DIR "${CMAKE_SOURCE_DIR}/Grains")
set(GRAINS_DATA_BIN_DIR "${CMAKE_BINARY_DIR}/Grains")

add_custom_target(risingsphereunbounded_data ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory "${GRAINS_DATA_BIN_DIR}"
  COMMAND ${CMAKE_COMMAND} -E make_directory "Res"
  COMMAND ${CMAKE_COMMAND} -E make_directory "Savings"
  COMMAND ${CMAKE_COMMAND} -E copy_directory "${GRAINS_DATA_SRC_DIR}" "${GRAINS_DATA_BIN_DIR}"
  COMMENT "Copying Grains/ into build tree"
  VERBATIM
) 
 
octree_add_executable(${CMAKE_SOURCE_DIR}/risingsphereunbounded.c)

add_dependencies(risingsphereunbounded risingsphereunbounded_data)
