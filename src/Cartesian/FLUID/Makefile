# Dependencies
CPPFLAGS += $(PACIFIC_XERCESC_CFLAGS) $(PACIFIC_ZLIB_CFLAGS)
LDLIBS   += $(PACIFIC_XERCESC_LFLAGS) $(PACIFIC_ZLIB_LFLAGS)

MAC_INCDIRS := \
  $(MAC_SRCDIR_ABS)/ArrayVector/include \
  $(MAC_SRCDIR_ABS)/Base/include \
  $(MAC_SRCDIR_ABS)/FiniteVolume/include \
  $(MAC_SRCDIR_ABS)/LinearAlgebra/include

CPPFLAGS += $(addprefix -I,$(MAC_INCDIRS))
LDLIBS += "-L$(MAC_LIBDIR) -lMAC"

GRAINS_INCDIRS := \
  $(GRAINS_SRCDIR_ABS)/Grains/Base/include \
  $(GRAINS_SRCDIR_ABS)/Grains/Component/include \
  $(GRAINS_SRCDIR_ABS)/Grains/Fluid/include \
  $(GRAINS_SRCDIR_ABS)/Grains/Geometry/include \
  $(GRAINS_SRCDIR_ABS)/Grains/PostProcessing/include \
  $(GRAINS_SRCDIR_ABS)/Grains/TimeIntegration/include \
  $(GRAINS_SRCDIR_ABS)/Grains/Utility/include \
  $(GRAINS_SRCDIR_ABS)/Grains/XML3/include \
  $(GRAINS_SRCDIR_ABS)/Grains/InterfaceCBasilisk/include

CPPFLAGS += $(addprefix -I,$(GRAINS_INCDIRS))
LDLIBS += "-L$(GRAINS_LIBDIR) -lGrains"

# Complilers
CXX = $(PACIFIC_MPICXX)
CC = $(PACIFIC_MPICC)

# Build directories
OBJDIR := $(FLUID_OBJDIR)
LIBDIR := $(FLUID_LIBDIR)

INCDIRS := \
  $(FLUID_SRCDIR_ABS)/FluidSolidInterface/include \
  $(FLUID_SRCDIR_ABS)/Misc/include \
  $(FLUID_SRCDIR_ABS)/PostProcessing/include \
  $(FLUID_SRCDIR_ABS)/Solvers/DLMFD/include \
  $(FLUID_SRCDIR_ABS)/Solvers/Viscoplastic_HeatTransfer/include 

SRCS := \
  FluidSolidInterface/src/FS_2Dcylinder.cc \
  FluidSolidInterface/src/FS_3Dbox.cc \
  FluidSolidInterface/src/FS_3Dcylinder.cc \
  FluidSolidInterface/src/FS_AllRigidBodies.cc \
  FluidSolidInterface/src/FS_GeneralPolyhedron.cc \
  FluidSolidInterface/src/FS_Grains3DPlugIn.cc \
  FluidSolidInterface/src/FS_RigidBody_BuilderFactory.cc \
  FluidSolidInterface/src/FS_RigidBody.cc \
  FluidSolidInterface/src/FS_SolidPlugIn_BuilderFactory.cc \
  FluidSolidInterface/src/FS_SolidPlugIn.cc \
  FluidSolidInterface/src/FS_Sphere.cc \
  Misc/src/PAC_computingtime.cc \
  Misc/src/PAC_solvercomputingtime.cc \
  Misc/src/PAC_Misc.cc \
  PostProcessing/src/PostProcessing.cc \
  Solvers/DLMFD/src/DLMFD_3Dbox.cc \
  Solvers/DLMFD/src/DLMFD_BoundaryMultiplierPoint.cc \
  Solvers/DLMFD/src/DLMFD_InteriorMultiplierPoint.cc \
  Solvers/DLMFD/src/DLMFD_RigidBody_BuilderFactory.cc \
  Solvers/DLMFD/src/DLMFD_3Dcylinder.cc \
  Solvers/DLMFD/src/DLMFD_FictitiousDomain.cc \
  Solvers/DLMFD/src/DLMFD_ParticlePoint.cc \
  Solvers/DLMFD/src/DLMFD_RigidBody.cc \
  Solvers/DLMFD/src/DLMFD_AllGeomBoundaries.cc \
  Solvers/DLMFD/src/DLMFD_GeneralPolyhedron.cc \
  Solvers/DLMFD/src/DLMFD_ProjectionNavierStokes.cc \
  Solvers/DLMFD/src/DLMFD_Sphere.cc \
  Solvers/DLMFD/src/DLMFD_AllRigidBodies.cc \
  Solvers/DLMFD/src/DLMFD_GeomBoundary.cc \
  Solvers/DLMFD/src/DLMFD_ProjectionNavierStokesSystem.cc \
  Solvers/DLMFD/src/DLMFD_XYZBoxPlan.cc \
  Solvers/Viscoplastic_HeatTransfer/src/VPH_HeatTransfer.cc \
  Solvers/Viscoplastic_HeatTransfer/src/VPH_HeatTransferSystem.cc \
  Solvers/Viscoplastic_HeatTransfer/src/VPH_Viscoplastic.cc \
  Solvers/Viscoplastic_HeatTransfer/src/VPH_Viscoplastic_Misc.cc \
  Solvers/Viscoplastic_HeatTransfer/src/VPH_ViscoplasticSystem.cc \

CPPFLAGS += $(addprefix -I,$(INCDIRS))

OBJS := $(patsubst %.cc,$(OBJDIR)/%.o,$(SRCS))
DEPS := $(OBJS:.o=.d)

CXXFLAGS += -MMD -MP

# Compile rule
$(OBJDIR)/%.o: %.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

LIBSO := $(LIBDIR)/libFLUID.so

$(LIBDIR):
	@mkdir -p $@

CXXFLAGS += -fPIC
$(LIBSO): $(OBJS) | $(LIBDIR)
	$(CXX) -shared -o $@ $(OBJS) $(LDFLAGS) $(LDLIBS)

.PHONY: all clean

all: $(LIBSO)

clean:
	rm -rf $(OBJDIR)
	rm -f $(LIBSO)
