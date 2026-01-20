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
# Collect Grains sources
#

set(GRAINS_HEADERS "")
set(GRAINS_INCLUDE_DIRS
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
set(GRAINS_SOURCES "")

list(APPEND GRAINS_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Base/include>")

list(APPEND GRAINS_SOURCES
  Base/src/AppBrownian.cpp
  Base/src/AppCollision.cpp
  Base/src/App.cpp
  Base/src/Cell.cpp
  Base/src/computingtime.cpp
  Base/src/Error.cpp
  Base/src/GrainsBuilderFactory.cpp
  Base/src/Grains.cpp
  Base/src/GrainsCRBFeatures.cpp
  Base/src/GrainsExec.cpp
  Base/src/GrainsMPI.cpp
  Base/src/GrainsMPIWrapper.cpp
  Base/src/GrainsParameters.cpp
  Base/src/GrainsPostProcessing.cpp
  Base/src/GrainsTestDev.cpp
  Base/src/LinkedCell.cpp
  Base/src/MPINeighbors.cpp
  Base/src/solvercomputingtime.cpp
)

list(APPEND GRAINS_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Component/include>")

list(APPEND GRAINS_SOURCES
  Component/src/AllComponents.cpp
  Component/src/CompositeObstacle.cpp
  Component/src/CylindricalShell.cpp
  Component/src/HookeMemoryContactForceModel.cpp
  Component/src/ObstacleImposedForce.cpp
  Component/src/Particle.cpp
  Component/src/ParticleKinematicsSphere.cpp
  Component/src/STLObstacle.cpp
  Component/src/Window.cpp
  Component/src/AllInsertionWindows.cpp
  Component/src/CompositeParticle.cpp
  Component/src/HertzContactForceModel.cpp
  Component/src/KinematicsBuilderFactory.cpp
  Component/src/ObstacleImposedVelocity.cpp
  Component/src/ParticleKinematics2D.cpp
  Component/src/QuadrilobeCylinder.cpp
  Component/src/Torsor.cpp
  Component/src/CollisionHistory.cpp
  Component/src/ContactBuilderFactory.cpp
  Component/src/HertzMemoryContactForceModel.cpp
  Component/src/ObstacleBuilderFactory.cpp
  Component/src/ObstacleKinematicsForce.cpp
  Component/src/ParticleKinematics3D.cpp
  Component/src/RoughWall.cpp             
  Component/src/TrilobeCylinder.cpp
  Component/src/Component.cpp
  Component/src/ContactForceModel.cpp
  Component/src/HookeContactForceModel.cpp
  Component/src/Obstacle.cpp
  Component/src/ObstacleKinematicsVelocity.cpp
  Component/src/ParticleKinematics.cpp
  Component/src/SimpleObstacle.cpp
  Component/src/TruncatedConicalShell.cpp
)

list(APPEND GRAINS_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Geometry/include>")


list(APPEND GRAINS_SOURCES
  Fluid/src/AppPRSHydroFT.cpp
  Fluid/src/GrainsCoupledWithFluid.cpp
)

list(APPEND GRAINS_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Fluid/include>")


list(APPEND GRAINS_SOURCES
  Geometry/src/BBox.cpp
  Geometry/src/Box.cpp
  Geometry/src/BVolume.cpp
  Geometry/src/Cone.cpp
  Geometry/src/ConvexBuilderFactory.cpp
  Geometry/src/Convex.cpp
  Geometry/src/Cylinder.cpp
  Geometry/src/Disc.cpp
  Geometry/src/GJK_SV2.cpp
  Geometry/src/GJK_SV.cpp
  Geometry/src/Group3.cpp
  Geometry/src/IndexArray.cpp
  Geometry/src/Matrix.cpp
  Geometry/src/OBB.cpp
  Geometry/src/OBC.cpp
  Geometry/src/Point3.cpp
  Geometry/src/PointC.cpp
  Geometry/src/PointContact.cpp
  Geometry/src/Polygon.cpp
  Geometry/src/Polyhedron.cpp
  Geometry/src/Polytope.cpp
  Geometry/src/Quaternion.cpp
  Geometry/src/Rectangle.cpp
  Geometry/src/RigidBody.cpp
  Geometry/src/RigidBodyWithCrust.cpp
  Geometry/src/Segment.cpp
  Geometry/src/Shape.cpp
  Geometry/src/Sphere.cpp
  Geometry/src/SpheroCylinder.cpp
  Geometry/src/SpheroCylindricalPrism.cpp
  Geometry/src/STLTriangle.cpp
  Geometry/src/STLVertex.cpp
  Geometry/src/Superquadric.cpp
  Geometry/src/Transform.cpp
  Geometry/src/TrapezoidalPrism.cpp
  Geometry/src/TruncatedCone.cpp
  Geometry/src/Vector3.cpp
  Geometry/src/VertexBase.cpp
)

list(APPEND GRAINS_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/PostProcessing/include>")

list(APPEND GRAINS_SOURCES
  PostProcessing/src/ParaviewPostProcessingWriter.cpp
  PostProcessing/src/ParaviewPostProcessingWriter_MPIIO.cpp
  PostProcessing/src/ParaviewPostProcessingWriter_Serial.cpp
  PostProcessing/src/PostProcessingWriterBuilderFactory.cpp
  PostProcessing/src/PostProcessingWriter.cpp
  PostProcessing/src/RawDataPostProcessingWriter.cpp
)

list(APPEND GRAINS_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/TimeIntegration/include>")

list(APPEND GRAINS_SOURCES 
  TimeIntegration/src/FirstOrderExplicit.cpp
  TimeIntegration/src/SecondOrderAdamsBashforth.cpp
  TimeIntegration/src/SecondOrderExplicit.cpp
  TimeIntegration/src/SecondOrderLeapFrog.cpp
  TimeIntegration/src/TimeIntegratorBuilderFactory.cpp
  TimeIntegration/src/TimeIntegrator.cpp
)

list(APPEND GRAINS_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Utility/include>")
list(APPEND GRAINS_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/Utility/include>")

list(APPEND GRAINS_SOURCES
  Utility/src/Data.cpp
)

list(APPEND GRAINS_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/XML3/include>")

list(APPEND GRAINS_SOURCES
  XML3/src/ReaderXML.cpp
  XML3/src/WriterXML.cpp
)

list(APPEND GRAINS_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/InterfaceCBasilisk/include>")

list(APPEND GRAINS_SOURCES
  InterfaceCBasilisk/src/InterfaceGrains3DBasilisk.cpp
)


#
# Build Grains shared object
#

# message(STATUS ">>> GRAINS_SOURCES = ${GRAINS_SOURCES}")

set(CMAKE_CXX_VISIBILITY_PRESET default)
set(CMAKE_VISIBILITY_INLINES_HIDDEN  OFF)

add_library(Grains SHARED ${GRAINS_SOURCES})
add_library(PacIFiC::Grains ALIAS Grains) 

set_target_properties(Grains
  PROPERTIES
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

set_target_properties(Grains PROPERTIES
  CXX_VISIBILITY_PRESET default
)
 
target_include_directories(Grains
  PUBLIC
    ${GRAINS_INCLUDE_DIRS}
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(Grains 
  PUBLIC
    MPI::MPI_C
    ZLIB::ZLIB
    XercesC::XercesC
)
 