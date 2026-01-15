
function(add_basilisk_executable SOURCE_FILE)
  find_package(HDF5 COMPONENTS C HL REQUIRED)
  find_package(MPI REQUIRED)
  
  get_filename_component(source_name ${SOURCE_FILE} NAME_WE)
  set(output_c "${CMAKE_BINARY_DIR}/_${source_name}.c")

  add_custom_command(
    OUTPUT "${CMAKE_BINARY_DIR}/_${source_name}.c"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${SOURCE_FILE}" "${CMAKE_BINARY_DIR}/${source_name}.c"
    COMMAND ${QCC_EXECUTABLE}
      "${source_name}.c"
      -I"${CMAKE_SOURCE_DIR}"
      -I"${OCTREE_DLMFD_DIR}"
      -I"${OCTREE_VTKHYPERTREE_DIR}"
      -DTRACE=2
      -source
    DEPENDS ${SOURCE_FILE}
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}" 
  )

  add_executable(${source_name} "_${source_name}.c")

  target_include_directories(${source_name} PRIVATE "${CMAKE_SOURCE_DIR}")
  target_compile_options(${source_name} PRIVATE -Wall -D_FORTIFY_SOURCE=2 -g -pipe)

  target_link_libraries(${source_name}
    PRIVATE
      hdf5::hdf5
      hdf5::hdf5_hl
      MPI::MPI_C
      m
  )

  set_target_properties(${source_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    BUILD_RPATH "${CMAKE_BINARY_DIR}"
  )
endfunction()

function(add_basilisk_grains_executable SOURCE_FILE)
  find_package(HDF5 COMPONENTS C HL REQUIRED)
  find_package(MPI REQUIRED)

  get_filename_component(source_name ${SOURCE_FILE} NAME_WE)
  set(output_c "${CMAKE_BINARY_DIR}/_${source_name}.c")

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

  
  add_custom_command(
    OUTPUT ${output_c}
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${SOURCE_FILE}" "${CMAKE_BINARY_DIR}/${source_name}.c"
    COMMAND $<TARGET_FILE:qcc-toolchain>
      "${source_name}.c"
      -I"${CMAKE_SOURCE_DIR}"
      -I"${OCTREE_DLMFD_DIR}"
      -I"${OCTREE_VTKHYPERTREE_DIR}"
      -DTRACE=2
      -source
    DEPENDS
      qcc-toolchain
      Grains_ensure_runtime_data
      ${SOURCE_FILE}

    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    COMMENT "Basilisk: transpiling ${SOURCE_FILE} -> ${output_c}"
  )

  add_executable(${source_name} ${output_c})

  target_include_directories(${source_name} PRIVATE "${CMAKE_SOURCE_DIR}")
  target_compile_options(${source_name} PRIVATE -Wall -D_FORTIFY_SOURCE=2 -g -pipe)

  target_link_libraries(${source_name}
    PRIVATE
      hdf5::hdf5
      hdf5::hdf5_hl
      MPI::MPI_C 
      PacIFiC::Grains
      m
  )

  set_target_properties(${source_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}"
    BUILD_RPATH "${CMAKE_BINARY_DIR}"
  )
endfunction()
