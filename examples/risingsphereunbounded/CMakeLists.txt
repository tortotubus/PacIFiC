cmake_minimum_required(VERSION 3.14...4.0)
project(RisingSphereUnbounded)

include(FetchContent)

FetchContent_Declare(PacIFiC
  SOURCE_DIR ../../../
)

FetchContent_MakeAvailable(PacIFiC)

message(STATUS "PacIFiC source: ${PacIFiC_SOURCE_DIR}")

# Enable C language
enable_language(C)

find_package(HDF5 REQUIRED COMPONENTS C HL)
find_package(MPI REQUIRED)
#
#
#
#
#

set(GRAINS_DATA_SRC_DIR "${CMAKE_SOURCE_DIR}/Grains")
set(GRAINS_DATA_BIN_DIR "${CMAKE_BINARY_DIR}/Grains")

add_custom_target(Grains_ensure_runtime_data ALL
  COMMAND ${CMAKE_COMMAND} -E make_directory "${GRAINS_DATA_BIN_DIR}"
  COMMAND ${CMAKE_COMMAND} -E make_directory "Res"
  COMMAND ${CMAKE_COMMAND} -E make_directory "Savings"
  COMMAND ${CMAKE_COMMAND} -E copy_directory "${GRAINS_DATA_SRC_DIR}" "${GRAINS_DATA_BIN_DIR}"
  COMMENT "Copying Grains/ into build tree"
  VERBATIM
)

#
#
#
#
# 

set (BASILISK_C_SRCS
  "${CMAKE_SOURCE_DIR}/risingsphereunbounded.c"
)

foreach(basilisk_source ${BASILISK_C_SRCS})
  get_filename_component(source_name ${basilisk_source} NAME_WE)

  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/_${source_name}.c"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    # 1) ensure the build dir exists
    COMMAND ${CMAKE_COMMAND} -E make_directory
            "${CMAKE_BINARY_DIR}"
    # 2) copy the .c file over
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_SOURCE_DIR}/${source_name}.c"
            "${CMAKE_BINARY_DIR}/${source_name}.c"
    # 3) transpile with qcc
    COMMAND ${QCC_EXECUTABLE}
      "${source_name}.c" 
      -I"${CMAKE_SOURCE_DIR}"
      -I"${OCTREE_DLMFD_DIR}"
      -I"${OCTREE_VTKHYPERTREE_DIR}"
      -DTRACE=2
      -source
    DEPENDS
      basilisk #source
      # ${BASILISK_LIBRARY_HDRS}
      ${basilisk_source}
      Grains_ensure_runtime_data
  )

  add_executable(${source_name} "_${source_name}.c")

  set_target_properties(${source_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
  )

  # 4. Include your headers
  target_include_directories(${source_name} PRIVATE
    "${CMAKE_SOURCE_DIR}"
  )

  target_compile_options(${source_name} PRIVATE
    -Wall
    -D_FORTIFY_SOURCE=2
    -g
    -pipe
  )

  target_link_libraries(${source_name}
    PRIVATE
      hdf5::hdf5
      hdf5::hdf5_hl
      MPI::MPI_C
      PacIFiC::Grains
      m
  )

  set_target_properties(${source_name} PROPERTIES
    BUILD_RPATH "${CMAKE_BINARY_DIR}"
  )

endforeach()

