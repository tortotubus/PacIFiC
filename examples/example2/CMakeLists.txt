cmake_minimum_required(VERSION 3.14...4.0)
project(GrainsExample2)

#
#
#
#
#

include(FetchContent)


FetchContent_Declare(PacIFiC
  SOURCE_DIR ../../../
)

# FetchContent_Declare(PacIFiC
#   GIT_REPOSITORY https://github.com/anthonywachs/PacIFiC.git
#)

FetchContent_MakeAvailable(PacIFiC)

set_target_properties(FLUID  PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(Grains PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(Grains_main PROPERTIES EXCLUDE_FROM_ALL FALSE)
set_target_properties(Grains_position PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(Grains_rotationmatrix PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(Grains_shapefile PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(MAC PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(qcc-grammar PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(qcc-install PROPERTIES EXCLUDE_FROM_ALL TRUE)
set_target_properties(qcc-toolchain PROPERTIES EXCLUDE_FROM_ALL TRUE)

message(STATUS "PacIFiC source: ${PacIFiC_SOURCE_DIR}")

#
#
#
#
#
  
