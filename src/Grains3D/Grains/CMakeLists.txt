#
# Subdirectories
#

add_subdirectory(Dtd)
add_subdirectory(Tools)

#
#
#

# ---- Generate Version.hh ----
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Utility/include/Version.hh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Utility/include/Version.hh"
  @ONLY
)


#
# Collect Grains3D sources
#

set(GRAINS3D_INCLUDE_DIRS "")
set(GRAINS3D_SOURCES "")

# Base
list(APPEND GRAINS3D_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Base/include>")
file(GLOB GRAINS3D_BASE_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/Base/src/*.cpp"
)
list(SORT GRAINS3D_BASE_SOURCES)
list(APPEND GRAINS3D_SOURCES
  ${GRAINS3D_BASE_SOURCES}
)

# Component
list(APPEND GRAINS3D_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Component/include>")
file(GLOB GRAINS3D_COMPONENT_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/Component/src/*.cpp"
)
list(SORT GRAINS3D_COMPONENT_SOURCES)
list(APPEND GRAINS3D_SOURCES
  ${GRAINS3D_COMPONENT_SOURCES}
)

# Fluid
list(APPEND GRAINS3D_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Fluid/include>")
file(GLOB GRAINS3D_FLUID_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/Fluid/src/*.cpp"
)
list(SORT GRAINS3D_FLUID_SOURCES)
list(APPEND GRAINS3D_SOURCES
  ${GRAINS3D_FLUID_SOURCES}
)

# Geometry
list(APPEND GRAINS3D_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Geometry/include>")
file(GLOB GRAINS3D_GEOMETRY_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/Geometry/src/*.cpp"
)
list(SORT GRAINS3D_GEOMETRY_SOURCES)
list(APPEND GRAINS3D_SOURCES
  ${GRAINS3D_GEOMETRY_SOURCES}
)

# InterfaceCBasilisk
list(APPEND GRAINS3D_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/InterfaceCBasilisk/include>")
file(GLOB GRAINS3D_INTERFACECBASILISK_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/InterfaceCBasilisk/src/*.cpp"
)
list(SORT GRAINS3D_INTERFACECBASILISK_SOURCES)
list(APPEND GRAINS3D_SOURCES
  ${GRAINS3D_INTERFACECBASILISK_SOURCES}
)

# PostProcessing
list(APPEND GRAINS3D_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/PostProcessing/include>")
file(GLOB GRAINS3D_POSTPROCESSING_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/PostProcessing/src/*.cpp"
)
list(SORT GRAINS3D_POSTPROCESSING_SOURCES)
list(APPEND GRAINS3D_SOURCES
  ${GRAINS3D_POSTPROCESSING_SOURCES}
)

# TimeIntegration
list(APPEND GRAINS3D_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/TimeIntegration/include>")
file(GLOB GRAINS3D_TIMEINTEGRATION_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/TimeIntegration/src/*.cpp"
)
list(SORT GRAINS3D_TIMEINTEGRATION_SOURCES)
list(APPEND GRAINS3D_SOURCES
  ${GRAINS3D_TIMEINTEGRATION_SOURCES}
)

# Utility
list(APPEND GRAINS3D_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Utility/include>")
file(GLOB GRAINS3D_UTILITY_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/Utility/src/*.cpp"
)
list(SORT GRAINS3D_UTILITY_SOURCES)
list(APPEND GRAINS3D_SOURCES
  ${GRAINS3D_UTILITY_SOURCES}
)

# XML3
list(APPEND GRAINS3D_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/XML3/include>")
file(GLOB GRAINS3D_XML3_SOURCES CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/XML3/src/*.cpp"
)
list(SORT GRAINS3D_XML3_SOURCES)
list(APPEND GRAINS3D_SOURCES
  ${GRAINS3D_XML3_SOURCES}
)



#
# Build Grains3D shared object
#

# message(STATUS ">>> GRAINS3D_SOURCES = ${GRAINS3D_SOURCES}")

set(CMAKE_CXX_VISIBILITY_PRESET default)
set(CMAKE_VISIBILITY_INLINES_HIDDEN  OFF)

add_library(Grains3D SHARED ${GRAINS3D_SOURCES})
add_library(PacIFiC::Grains3D ALIAS Grains3D) 

# Pendulum's implicit structural solver uses Eigen, kept out of public headers.
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(Grains3D PRIVATE Eigen3::Eigen)

set_target_properties(Grains3D PROPERTIES
  OUTPUT_NAME "Grains3D"
)

set_target_properties(Grains3D PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

set_target_properties(Grains3D PROPERTIES
  CXX_VISIBILITY_PRESET default
)
 
target_include_directories(Grains3D
  PUBLIC
    ${GRAINS3D_INCLUDE_DIRS}
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if (PACIFIC_USE_MPI)
  target_link_libraries(Grains3D 
    PUBLIC
      MPI::MPI_C
      ZLIB::ZLIB
      XercesC::XercesC
  )
else()
  target_link_libraries(Grains3D 
    PUBLIC
      ZLIB::ZLIB
      XercesC::XercesC
  )
endif()

if(PACIFIC_BUILD_TESTS)
  add_executable(grains3d_pendulum_solver_test ../Tests/BendingRodImplicit.cpp)
  target_include_directories(grains3d_pendulum_solver_test PRIVATE Base/src)
  target_link_libraries(grains3d_pendulum_solver_test PRIVATE Grains3D Eigen3::Eigen)
  add_test(NAME grains3d.pendulum_solver COMMAND grains3d_pendulum_solver_test)
  add_executable(grains3d_pendulum_test ../Tests/Pendulum.cpp)
  target_link_libraries(grains3d_pendulum_test PRIVATE Grains3D)
  set(PENDULUM_DTD "${CMAKE_CURRENT_SOURCE_DIR}/Dtd/Grains3D-v2.dtd")
  foreach(PENDULUM_TYPE Standard MPI)
    set(PENDULUM_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/pendulum-${PENDULUM_TYPE}")
    file(MAKE_DIRECTORY "${PENDULUM_OUTPUT}")
    set(PENDULUM_DECOMPOSITION "")
    if(PENDULUM_TYPE STREQUAL "MPI")
      set(PENDULUM_DECOMPOSITION "<DomainDecomposition NX=\"2\" NY=\"1\" NZ=\"1\"/><Verbosity Level=\"0\"/>")
    endif()
    configure_file(../Tests/fixtures/pendulum.xml.in "${PENDULUM_OUTPUT}/input.xml" @ONLY)
  endforeach()
  add_test(NAME grains3d.pendulum_xml COMMAND grains3d_pendulum_test "${CMAKE_CURRENT_BINARY_DIR}/pendulum-Standard/input.xml")
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pendulum-example")
  add_test(NAME grains3d.pendulum_example COMMAND grains3d_pendulum_test
    "${CMAKE_CURRENT_SOURCE_DIR}/../Examples/Pendulum/pendulum.xml")
  set_tests_properties(grains3d.pendulum_example PROPERTIES
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pendulum-example")
  set_tests_properties(grains3d.pendulum_example grains3d.pendulum_xml PROPERTIES
    FAIL_REGULAR_EXPRESSION "ERR RigidBodyWithCrust;Handling contact error")
  foreach(geometry geometry single dynamics)
    add_test(NAME grains3d.pendulum_${geometry} COMMAND grains3d_pendulum_test
      "${CMAKE_CURRENT_BINARY_DIR}/pendulum-Standard/input.xml" ${geometry})
  endforeach()
  foreach(invalid rigidity damping radius nonfinite duplicate direction empty count fraction overflow length missing)
    add_test(NAME grains3d.pendulum_invalid_${invalid} COMMAND grains3d_pendulum_test "${CMAKE_CURRENT_BINARY_DIR}/pendulum-Standard/input.xml" ${invalid})
  endforeach()
  if(PACIFIC_USE_MPI)
    add_test(NAME grains3d.pendulum_xml_mpi
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:grains3d_pendulum_test> ${MPIEXEC_POSTFLAGS} "${CMAKE_CURRENT_BINARY_DIR}/pendulum-MPI/input.xml")
    set_tests_properties(grains3d.pendulum_xml_mpi PROPERTIES PROCESSORS 2 TIMEOUT 60)
  endif()
  add_executable(grains3d_external_model_test ../Tests/ExternalModel.cpp)
  target_link_libraries(grains3d_external_model_test PRIVATE Grains3D)
  target_compile_definitions(grains3d_external_model_test PRIVATE
    GRAINS_EXTERNAL_MODEL_LAW="${CMAKE_CURRENT_SOURCE_DIR}/../Tests/fixtures/external-model-contact.xml")
  add_test(NAME grains3d.external_model COMMAND grains3d_external_model_test)
  foreach(scheme FirstOrderExplicit SecondOrderExplicit SecondOrderAdamsBashforth)
    add_test(NAME grains3d.external_model_${scheme}
      COMMAND grains3d_external_model_test --scheme=${scheme})
  endforeach()
  foreach(failure missing duplicate foreign oversized throw)
    add_test(NAME grains3d.external_model_failure_${failure}
      COMMAND grains3d_external_model_test --bad=${failure})
  endforeach()
  if(PACIFIC_USE_MPI)
    foreach(mode mpi mpi-periodic)
      add_test(NAME grains3d.external_model_${mode}
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
          ${MPIEXEC_PREFLAGS} $<TARGET_FILE:grains3d_external_model_test>
          ${MPIEXEC_POSTFLAGS} --${mode} --bad=foreign)
      set_tests_properties(grains3d.external_model_${mode} PROPERTIES PROCESSORS 2 TIMEOUT 60)
    endforeach()
  endif()
  add_executable(grains3d_load_driver_test ../Tests/LoadsAndKinematicsDriver.cpp)
  target_link_libraries(grains3d_load_driver_test PRIVATE Grains3D)
  add_test(NAME grains3d.loads_and_kinematics_driver COMMAND grains3d_load_driver_test)
  add_executable(grains3d_external_state_test ../Tests/ExternalParticleState.cpp)
  target_link_libraries(grains3d_external_state_test PRIVATE Grains3D)
  target_compile_definitions(grains3d_external_state_test PRIVATE
    GRAINS_COMPOSITE_FIXTURE="${CMAKE_CURRENT_SOURCE_DIR}/../Tests/fixtures/external-composite.xml")
  add_test(NAME grains3d.external_state COMMAND grains3d_external_state_test)
  if(PACIFIC_USE_MPI)
    add_test(NAME grains3d.external_state_mpi
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:grains3d_external_state_test>
        ${MPIEXEC_POSTFLAGS} --mpi)
    set_tests_properties(grains3d.external_state_mpi PROPERTIES PROCESSORS 2 TIMEOUT 60)
    add_test(NAME grains3d.external_state_mpi_periodic
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:grains3d_external_state_test>
        ${MPIEXEC_POSTFLAGS} --mpi-periodic)
    set_tests_properties(grains3d.external_state_mpi_periodic
      PROPERTIES PROCESSORS 2 TIMEOUT 60)
  endif()
endif()
 
