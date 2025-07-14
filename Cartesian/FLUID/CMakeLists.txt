

#
# Collect FLUID sources
#

set(FLUID_HEADERS "")
set(FLUID_INCLUDE_DIRS
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
set(FLUID_SOURCES "")

#
# FluidSolidInterface
#

list(APPEND FLUID_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/FluidSolidInterface/include>")

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

list(APPEND FLUID_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Misc/include>")

list(APPEND FLUID_SOURCES
  Misc/src/PAC_computingtime.cc  
  Misc/src/PAC_solvercomputingtime.cc
  Misc/src/PAC_Misc.cc
)

#
# PostProcessing
#

list(APPEND FLUID_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/PostProcessing/include>")

list(APPEND FLUID_SOURCES
  PostProcessing/src/PostProcessing.cc
)

#
# Solvers
#

list(APPEND FLUID_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solvers/DirectionSplitting/include>")

list(APPEND FLUID_SOURCES
  Solvers/DirectionSplitting/src/DS_2Dcylinder.cc
  Solvers/DirectionSplitting/src/DS_AllRigidBodies.cc
  Solvers/DirectionSplitting/src/DS_HeatTransfer.cc
  Solvers/DirectionSplitting/src/DS_NavierStokesSystem.cc
  Solvers/DirectionSplitting/src/DS_RigidBody.cc
  Solvers/DirectionSplitting/src/DS_3Dbox.cc
  Solvers/DirectionSplitting/src/DS_DirectionSplitting.cc
  Solvers/DirectionSplitting/src/DS_HeatTransferSystem.cc
  Solvers/DirectionSplitting/src/DS_PID.cc
  Solvers/DirectionSplitting/src/DS_Sphere.cc
  Solvers/DirectionSplitting/src/DS_3Dcylinder.cc
  Solvers/DirectionSplitting/src/DS_GeneralPolyhedron.cc
  Solvers/DirectionSplitting/src/DS_NavierStokes.cc
  Solvers/DirectionSplitting/src/DS_RigidBody_BuilderFactory.cc
  Solvers/DirectionSplitting/src/DS_STL.cc
)

list(APPEND FLUID_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Solvers/DLMFD/include>")

list(APPEND FLUID_SOURCES
  Solvers/DLMFD/src/DLMFD_3Dbox.cc
  Solvers/DLMFD/src/DLMFD_GeneralPolyhedron.cc
  Solvers/DLMFD/src/DLMFD_RigidBody_BuilderFactory.cc
  Solvers/DLMFD/src/DLMFD_3Dcylinder.cc
  Solvers/DLMFD/src/DLMFD_GeomBoundary.cc
  Solvers/DLMFD/src/DLMFD_RigidBody.cc
  Solvers/DLMFD/src/DLMFD_AllGeomBoundaries.cc
  Solvers/DLMFD/src/DLMFD_InteriorMultiplierPoint.cc
  Solvers/DLMFD/src/DLMFD_Sphere.cc
  Solvers/DLMFD/src/DLMFD_AllRigidBodies.cc
  Solvers/DLMFD/src/DLMFD_ParticlePoint.cc
  Solvers/DLMFD/src/DLMFD_XYZBoxPlan.cc
  Solvers/DLMFD/src/DLMFD_BoundaryMultiplierPoint.cc
  Solvers/DLMFD/src/DLMFD_ProjectionNavierStokes.cc
  Solvers/DLMFD/src/DLMFD_FictitiousDomain.cc
  Solvers/DLMFD/src/DLMFD_ProjectionNavierStokesSystem.cc
)

#
# Build FLUID shared object
#

# message(STATUS ">>> FLUID_SOURCES = ${FLUID_SOURCES}")

add_library(fluid SHARED ${FLUID_SOURCES} )
add_library(Pacific::fluid ALIAS fluid)

target_include_directories(fluid
  PRIVATE
    ${FLUID_INCLUDE_DIRS}
    $<INSTALL_INTERFACE:include>
)

target_link_libraries(fluid 
  PUBLIC
    Pacific::mac
    Pacific::grains
)
