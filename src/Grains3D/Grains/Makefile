CXX = $(PACIFIC_MPICXX)
CC = $(PACIFIC_MPICC)
CPPFLAGS += $(PACIFIC_XERCESC_CFLAGS) $(PACIFIC_ZLIB_CFLAGS)
LDLIBS   += $(PACIFIC_XERCESC_LFLAGS) $(PACIFIC_ZLIB_LFLAGS)

OBJDIR := $(LIBGRAINS_OBJDIR)
LIBDIR := $(LIBGRAINS_LIBDIR)
SHAREDIR := $(LIBGRAINS_SHAREDIR)

ifeq ($(strip $(SHAREDIR)),)
SHAREDIR := $(LIBDIR)/../share
endif

INCDIRS := \
  $(LIBGRAINS_SRCDIR)/Base/include \
  $(LIBGRAINS_SRCDIR)/Component/include \
  $(LIBGRAINS_SRCDIR)/Fluid/include \
  $(LIBGRAINS_SRCDIR)/Geometry/include \
  $(LIBGRAINS_SRCDIR)/PostProcessing/include \
  $(LIBGRAINS_SRCDIR)/TimeIntegration/include \
  $(LIBGRAINS_SRCDIR)/Utility/include \
  $(LIBGRAINS_SRCDIR)/XML3/include \
  $(LIBGRAINS_SRCDIR)/InterfaceCBasilisk/include

SRCS := \
  Base/src/AppBrownian.cpp  \
  Base/src/AppCollision.cpp  \
  Base/src/App.cpp  \
  Base/src/Cell.cpp  \
  Base/src/computingtime.cpp  \
  Base/src/Error.cpp  \
  Base/src/GrainsBuilderFactory.cpp  \
  Base/src/Grains.cpp  \
  Base/src/GrainsCRBFeatures.cpp  \
  Base/src/GrainsExec.cpp  \
  Base/src/GrainsMPI.cpp  \
  Base/src/GrainsMPIWrapper.cpp  \
  Base/src/GrainsParameters.cpp  \
  Base/src/GrainsPostProcessing.cpp  \
  Base/src/GrainsTestDev.cpp  \
  Base/src/LinkedCell.cpp  \
  Base/src/MPINeighbors.cpp  \
  Base/src/solvercomputingtime.cpp  \
  Component/src/AllComponents.cpp  \
  Component/src/CompositeObstacle.cpp  \
  Component/src/CylindricalShell.cpp  \
  Component/src/HookeMemoryContactForceModel.cpp  \
  Component/src/ObstacleImposedForce.cpp  \
  Component/src/Particle.cpp  \
  Component/src/ParticleKinematicsSphere.cpp  \
  Component/src/STLObstacle.cpp  \
  Component/src/Window.cpp  \
  Component/src/AllInsertionWindows.cpp  \
  Component/src/CompositeParticle.cpp  \
  Component/src/HertzContactForceModel.cpp  \
  Component/src/KinematicsBuilderFactory.cpp  \
  Component/src/ObstacleImposedVelocity.cpp  \
  Component/src/ParticleKinematics2D.cpp  \
  Component/src/QuadrilobeCylinder.cpp  \
  Component/src/Torsor.cpp  \
  Component/src/CollisionHistory.cpp  \
  Component/src/ContactBuilderFactory.cpp  \
  Component/src/HertzMemoryContactForceModel.cpp  \
  Component/src/ObstacleBuilderFactory.cpp  \
  Component/src/ObstacleKinematicsForce.cpp  \
  Component/src/ParticleKinematics3D.cpp  \
  Component/src/RoughWall.cpp  \
  Component/src/TrilobeCylinder.cpp  \
  Component/src/Component.cpp  \
  Component/src/ContactForceModel.cpp  \
  Component/src/HookeContactForceModel.cpp  \
  Component/src/Obstacle.cpp  \
  Component/src/ObstacleKinematicsVelocity.cpp  \
  Component/src/ParticleKinematics.cpp  \
  Component/src/SimpleObstacle.cpp  \
  Component/src/TruncatedConicalShell.cpp  \
  Fluid/src/AppPRSHydroFT.cpp  \
  Fluid/src/GrainsCoupledWithFluid.cpp  \
  Geometry/src/BBox.cpp  \
  Geometry/src/Box.cpp  \
  Geometry/src/BVolume.cpp  \
  Geometry/src/Cone.cpp  \
  Geometry/src/ConvexBuilderFactory.cpp  \
  Geometry/src/Convex.cpp  \
  Geometry/src/Cylinder.cpp  \
  Geometry/src/Disc.cpp  \
  Geometry/src/GJK_SV2.cpp  \
  Geometry/src/GJK_SV.cpp  \
  Geometry/src/Group3.cpp  \
  Geometry/src/IndexArray.cpp  \
  Geometry/src/Matrix.cpp  \
  Geometry/src/OBB.cpp  \
  Geometry/src/OBC.cpp  \
  Geometry/src/Point3.cpp  \
  Geometry/src/PointC.cpp  \
  Geometry/src/PointContact.cpp  \
  Geometry/src/Polygon.cpp  \
  Geometry/src/Polyhedron.cpp  \
  Geometry/src/Polytope.cpp  \
  Geometry/src/Quaternion.cpp  \
  Geometry/src/Rectangle.cpp  \
  Geometry/src/RigidBody.cpp  \
  Geometry/src/RigidBodyWithCrust.cpp  \
  Geometry/src/Segment.cpp  \
  Geometry/src/Shape.cpp  \
  Geometry/src/Sphere.cpp  \
  Geometry/src/SpheroCylinder.cpp  \
  Geometry/src/SpheroCylindricalPrism.cpp  \
  Geometry/src/STLTriangle.cpp  \
  Geometry/src/STLVertex.cpp  \
  Geometry/src/Superquadric.cpp  \
  Geometry/src/Transform.cpp  \
  Geometry/src/TrapezoidalPrism.cpp  \
  Geometry/src/TruncatedCone.cpp  \
  Geometry/src/Vector3.cpp  \
  Geometry/src/VertexBase.cpp  \
  PostProcessing/src/ParaviewPostProcessingWriter.cpp  \
  PostProcessing/src/ParaviewPostProcessingWriter_MPIIO.cpp  \
  PostProcessing/src/ParaviewPostProcessingWriter_Serial.cpp  \
  PostProcessing/src/PostProcessingWriterBuilderFactory.cpp  \
  PostProcessing/src/PostProcessingWriter.cpp  \
  PostProcessing/src/RawDataPostProcessingWriter.cpp  \
  TimeIntegration/src/FirstOrderExplicit.cpp  \
  TimeIntegration/src/SecondOrderAdamsBashforth.cpp  \
  TimeIntegration/src/SecondOrderExplicit.cpp  \
  TimeIntegration/src/SecondOrderLeapFrog.cpp  \
  TimeIntegration/src/TimeIntegratorBuilderFactory.cpp  \
  TimeIntegration/src/TimeIntegrator.cpp  \
  Utility/src/Data.cpp  \
  Utility/src/Exec.cpp  \
  XML3/src/ReaderXML.cpp  \
  XML3/src/WriterXML.cpp  \
  InterfaceCBasilisk/src/InterfaceGrains3DBasilisk.cpp

CPPFLAGS += $(addprefix -I,$(INCDIRS))

OBJS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

CXXFLAGS += -MMD -MP

# Compile rule
$(OBJDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

LIBSO := $(LIBDIR)/libGrains.so

$(LIBDIR):
	@mkdir -p $@

CXXFLAGS += -fPIC
$(LIBSO): $(OBJS) | $(LIBDIR)
	$(CXX) -shared -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)

.PHONY: all clean

all: $(LIBSO) dtd tools

clean:
	rm -rf $(OBJDIR)/* \
	rm -f $(LIBSO) \
  rm -rf $(SHAREDIR)/*

dtd:
	@mkdir -p $(SHAREDIR)/Dtd
	@cp -f \
		Dtd/Contact.dtd \
		Dtd/Convex.dtd \
		Dtd/Grains2D_InFluid.dtd \
		Dtd/Grains2D.dtd \
		Dtd/Grains3D_InFluid.dtd \
		Dtd/Grains3D.dtd \
		$(SHAREDIR)/Dtd/

tools:
	@mkdir -p $(SHAREDIR)/Tools/ExecScripts
	@cp -rf \
	  Tools/ExecScripts/*.exec \
		$(SHAREDIR)/Tools/ExecScripts/