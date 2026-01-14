find_package(GTest REQUIRED)

enable_testing()

# Create the Grains test suite
file(GLOB_RECURSE GRAINS_TEST_SOURCES *.cpp)
add_executable(grains_tests ${GRAINS_TEST_SOURCES})
target_include_directories(grains_tests PRIVATE ${CMAKE_SOURCE_DIR})
target_link_libraries(grains_tests PRIVATE grains GTest::gtest_main)

# Set compile-time definitions
target_compile_definitions(grains_tests PRIVATE
    GRAINS_HOME="${GRAINS_BINARY_ROOT_DIR}"
    GRAINS_TEST_DATA_DIR="${GRAINS_BINARY_ROOT_DIR}/Tests/data"
)

# Set compile-time options
set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
 
target_compile_options(grains_tests PRIVATE
    -O0
    -g
    -fno-inline
    -fno-omit-frame-pointer
    -fno-optimize-sibling-calls
)

gtest_discover_tests(grains_tests)

# Copy data files into the build directory
set(GRAINS_TEST_DATA_SRC "${CMAKE_CURRENT_SOURCE_DIR}/data")      # your test folder
set(GRAINS_TEST_DATA_DST "${CMAKE_CURRENT_BINARY_DIR}/data")      # where you want it in build tree

add_custom_command(
  OUTPUT "${GRAINS_TEST_DATA_DST}/.copied"
  COMMAND ${CMAKE_COMMAND} -E make_directory "${GRAINS_TEST_DATA_DST}"
  COMMAND ${CMAKE_COMMAND} -E copy_directory "${GRAINS_TEST_DATA_SRC}" "${GRAINS_TEST_DATA_DST}"
  COMMAND ${CMAKE_COMMAND} -E touch "${GRAINS_TEST_DATA_DST}/.copied"
  DEPENDS ${GRAINS_TEST_DATA_SRC}    # (note: directory-dep is coarse; see below)
  COMMENT "Copying gtest test data into build tree"
)

add_custom_target(grains_tests_copy_data ALL
  DEPENDS "${GRAINS_TEST_DATA_DST}/.copied"
)

add_dependencies(grains_tests grains_tests_copy_data)