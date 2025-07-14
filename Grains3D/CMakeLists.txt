add_subdirectory(Grains)
add_subdirectory(Tools)



# Add Doxygen documentation
option(BUILD_GRAINS_DOC "Build documentation" OFF)
if (BUILD_GRAINS_DOC)
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
      add_custom_target(grains_doc ALL
          ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/DoxyGRAINS/Doxyfile
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/DoxyGRAINS/
          COMMENT "Generating API documentation with Doxygen"
          VERBATIM)
  endif()
endif()