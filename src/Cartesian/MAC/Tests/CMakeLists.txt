find_package(GTest REQUIRED)

set(MAC_TEST_SOURCES
  main_mpi_gtest.cc
  FV_ParaviewPostProcessingWriter_gtest.cc
  FV_ParaviewPostProcessingWriterVTKHDF_gtest.cc
)

add_executable(MAC_tests ${MAC_TEST_SOURCES})

set(_MAC_GTEST_LIBS "")
if(TARGET GTest::gtest)
  list(APPEND _MAC_GTEST_LIBS GTest::gtest)
elseif(TARGET GTest::GTest)
  list(APPEND _MAC_GTEST_LIBS GTest::GTest)
else()
  message(FATAL_ERROR "Unsupported GTest target names from find_package(GTest)")
endif()

target_link_libraries(MAC_tests
  PRIVATE
    PacIFiC::MAC
    MPI::MPI_C
    ${_MAC_GTEST_LIBS}
)

include(GoogleTest)

gtest_discover_tests(MAC_tests
  TEST_PREFIX "MAC."
  DISCOVERY_TIMEOUT 60
)
