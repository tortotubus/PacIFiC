cmake_minimum_required(VERSION 3.12)
project(grains_xml_editor LANGUAGES CXX)
include(GNUInstallDirs)

if(NOT TARGET LibXml2::LibXml2)
  find_package(LibXml2 2.14 REQUIRED)
endif()
set(CURSES_NEED_NCURSES TRUE)
find_package(Curses REQUIRED)

# Embed the standard 3D schema, expanding its local parameter entities so the
# installed executable does not depend on DTD files or source-tree paths.
set(_dtd_dir "${CMAKE_CURRENT_SOURCE_DIR}/../../Dtd")
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${_dtd_dir}/Grains3D.dtd" "${_dtd_dir}/Convex.dtd" "${_dtd_dir}/Contact.dtd")
file(READ "${_dtd_dir}/Grains3D.dtd" EDITOR_DTD_TEXT)
foreach(_part Convex Contact)
  file(READ "${_dtd_dir}/${_part}.dtd" _part_text)
  string(REPLACE "<!ENTITY % ${_part}.dtd SYSTEM \"${_part}.dtd\">" ""
    EDITOR_DTD_TEXT "${EDITOR_DTD_TEXT}")
  string(REPLACE "%${_part}.dtd;" "${_part_text}"
    EDITOR_DTD_TEXT "${EDITOR_DTD_TEXT}")
endforeach()
configure_file(schema.hpp.in "${CMAKE_CURRENT_BINARY_DIR}/generated/schema.hpp" @ONLY)

add_executable(Grains3D_xml
  main.cpp document.cpp constraints.cpp terminal.cpp editor.cpp projection.cpp cursor.cpp)
target_compile_features(Grains3D_xml PRIVATE cxx_std_17)
target_include_directories(Grains3D_xml PRIVATE
  "${CMAKE_CURRENT_BINARY_DIR}/generated" ${CURSES_INCLUDE_DIRS})
set_target_properties(Grains3D_xml PROPERTIES CXX_EXTENSIONS OFF OUTPUT_NAME grains-xml)
target_link_libraries(Grains3D_xml PRIVATE LibXml2::LibXml2 ${CURSES_LIBRARIES})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(Grains3D_xml PRIVATE -Wall -Wextra -Wpedantic -Wno-deprecated-declarations)
endif()

if(PACIFIC_MPI_INSTALL_BINDIR)
  set(_editor_bindir "${PACIFIC_MPI_INSTALL_BINDIR}")
else()
  set(_editor_bindir "${CMAKE_INSTALL_BINDIR}")
endif()
install(TARGETS Grains3D_xml RUNTIME DESTINATION "${_editor_bindir}" COMPONENT Grains3D_Runtime)
