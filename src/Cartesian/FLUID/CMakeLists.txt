
cmake_minimum_required(VERSION 3.12.0...4.0.0)
message(STATUS "CMake version: ${CMAKE_VERSION}")

set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to use")
set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Force the use of the chosen C++ standard.")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Enable C++ standard extensions.")

project(FLUID
  VERSION 0.0.1
  LANGUAGES CXX
)

#
# Collect FLUID sources
#

set(FLUID_SOURCES "")

#
# FluidSolidInterface
#

list(APPEND FLUID_SOURCES
  FluidSolidInterface/src/FS_2Dcylinder.cc
  FluidSolidInterface/src/FS_3Dbox.cc
  FluidSolidInterface/src/FS_3Dcylinder.cc
  FluidSolidInterface/src/FS_AllRigidBodies.cc
  FluidSolidInterface/src/FS_GeneralPolyhedron.cc
  FluidSolidInterface/src/FS_Grains3DPlugIn.cc
  FluidSolidInterface/src/FS_RigidBody_BuilderFactory.cc
  FluidSolidInterface/src/FS_RigidBody.cc
  FluidSolidInterface/src/FS_SolidPlugIn_BuilderFactory.cc
  FluidSolidInterface/src/FS_SolidPlugIn.cc
  FluidSolidInterface/src/FS_Sphere.cc
)

#
# Misc
#

list(APPEND FLUID_SOURCES
  Misc/src/PAC_computingtime.cc  
  Misc/src/PAC_solvercomputingtime.cc
  Misc/src/PAC_Misc.cc
)

#
# PostProcessing
#

list(APPEND FLUID_SOURCES
  PostProcessing/src/PostProcessing.cc
)

#
# Solvers
#

# list(APPEND FLUID_SOURCES
#   Solvers/DirectionSplitting/src/DS_2Dcylinder.cc
#   Solvers/DirectionSplitting/src/DS_AllRigidBodies.cc
#   Solvers/DirectionSplitting/src/DS_HeatTransfer.cc
#   Solvers/DirectionSplitting/src/DS_NavierStokesSystem.cc
#   Solvers/DirectionSplitting/src/DS_RigidBody.cc
#   Solvers/DirectionSplitting/src/DS_3Dbox.cc
#   Solvers/DirectionSplitting/src/DS_DirectionSplitting.cc
#   Solvers/DirectionSplitting/src/DS_HeatTransferSystem.cc
#   Solvers/DirectionSplitting/src/DS_PID.cc
#   Solvers/DirectionSplitting/src/DS_Sphere.cc
#   Solvers/DirectionSplitting/src/DS_3Dcylinder.cc
#   Solvers/DirectionSplitting/src/DS_GeneralPolyhedron.cc
#   Solvers/DirectionSplitting/src/DS_NavierStokes.cc
#   Solvers/DirectionSplitting/src/DS_RigidBody_BuilderFactory.cc
#   Solvers/DirectionSplitting/src/DS_STL.cc
# )

list(APPEND FLUID_SOURCES
  Solvers/DLMFD/src/DLMFD_3Dbox.cc
  Solvers/DLMFD/src/DLMFD_BoundaryMultiplierPoint.cc
  Solvers/DLMFD/src/DLMFD_InteriorMultiplierPoint.cc
  Solvers/DLMFD/src/DLMFD_RigidBody_BuilderFactory.cc
  Solvers/DLMFD/src/DLMFD_3Dcylinder.cc
  Solvers/DLMFD/src/DLMFD_FictitiousDomain.cc
  Solvers/DLMFD/src/DLMFD_ParticlePoint.cc
  Solvers/DLMFD/src/DLMFD_RigidBody.cc
  Solvers/DLMFD/src/DLMFD_AllGeomBoundaries.cc
  Solvers/DLMFD/src/DLMFD_GeneralPolyhedron.cc
  Solvers/DLMFD/src/DLMFD_ProjectionNavierStokes.cc
  Solvers/DLMFD/src/DLMFD_Sphere.cc
  Solvers/DLMFD/src/DLMFD_AllRigidBodies.cc
  Solvers/DLMFD/src/DLMFD_GeomBoundary.cc
  Solvers/DLMFD/src/DLMFD_ProjectionNavierStokesSystem.cc
  Solvers/DLMFD/src/DLMFD_XYZBoxPlan.cc
)

list(APPEND FLUID_SOURCES
  Solvers/Viscoplastic_HeatTransfer/src/VPH_HeatTransfer.cc
  Solvers/Viscoplastic_HeatTransfer/src/VPH_HeatTransferSystem.cc
  Solvers/Viscoplastic_HeatTransfer/src/VPH_Viscoplastic.cc
  Solvers/Viscoplastic_HeatTransfer/src/VPH_Viscoplastic_Misc.cc
  Solvers/Viscoplastic_HeatTransfer/src/VPH_ViscoplasticSystem.cc
)


#
# Build FLUID shared object
#

add_library(FLUID SHARED ${FLUID_SOURCES} )
add_library(PacIFiC::FLUID ALIAS FLUID)

set_target_properties(FLUID
  PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_include_directories(FLUID
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/FluidSolidInterface/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Misc/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/PostProcessing/include>
    # $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solvers/DirectionSplitting/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solvers/DLMFD/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solvers/Viscoplastic_HeatTransfer/include>

    # Install tree include dir (must NOT be in the source tree)
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/FLUID>
)

target_link_libraries(FLUID 
  PUBLIC
    PacIFiC::MAC
    PacIFiC::Grains
)


# ================================================================================
#
# Install 
#
# ================================================================================

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# ---- 1. Install the shared library ------------------------------------------------
install(TARGETS FLUID
  EXPORT FLUIDTargets
  LIBRARY DESTINATION ${PACIFIC_MPI_INSTALL_LIBDIR}      
  COMPONENT FLUID_Runtime
  RUNTIME DESTINATION ${PACIFIC_MPI_INSTALL_BINDIR}      
  COMPONENT FLUID_Runtime
  INCLUDES DESTINATION ${PACIFIC_MPI_INSTALL_INCLUDEDIR} 
)

# ---- 2. Install all public headers -------------------------------------------------
# If your project has multiple include roots, add them all.
install(DIRECTORY
  ${CMAKE_CURRENT_SOURCE_DIR}/FluidSolidInterface/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Misc/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/PostProcessing/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Solvers/DirectionSplitting/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Solvers/DLMFD/include/
  ${CMAKE_CURRENT_SOURCE_DIR}/Solvers/Viscoplastic_HeatTransfer/include/
  DESTINATION ${PACIFIC_MPI_INSTALL_INCLUDEDIR}/FLUID
  COMPONENT FLUID_Devel
  FILES_MATCHING PATTERN "*.h*"      # match .h, .hh, .hpp
)

# ---- 4. Generate CMake package configuration --------------------------------------
set(FLUID_INSTALL_CMAKEDIR ${PACIFIC_MPI_INSTALL_LIBDIR}/cmake/FLUID)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/FLUIDConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/FLUIDConfig.cmake
  INSTALL_DESTINATION ${FLUID_INSTALL_CMAKEDIR}
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/FLUIDConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

# ---- 5. Install the package + targets ---------------------------------------------
install(EXPORT FLUIDTargets
  FILE FLUIDTargets.cmake
  NAMESPACE PacIFiC::
  DESTINATION ${FLUID_INSTALL_CMAKEDIR}
  COMPONENT FLUID_Devel
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/FLUIDConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/FLUIDConfigVersion.cmake
  DESTINATION ${FLUID_INSTALL_CMAKEDIR}
  COMPONENT FLUID_Devel
)

# ---- 6. Export from the build tree for development -------------------------------
export(EXPORT FLUIDTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/FLUIDTargets.cmake
  NAMESPACE PacIFiC::
)
