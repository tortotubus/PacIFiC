# Compilation & linking flags
NVCC = $(PACIFIC_NVCC)
include ../mk/nvcc_link_flags.mk

NVCC_ARCHFLAGS += -arch=sm_${PACIFIC_NVCC_ARCH}
NVCCFLAGS += -O3 -g -Xcompiler -fPIC -x cu -rdc=true --extended-lambda --expt-relaxed-constexpr

CPPFLAGS += -I$(LIBGRAINS_GPU_INCDIR)
CPPFLAGS += -I$(LIBGRAINS_GPU_SRCDIR)/Base/include
CPPFLAGS += -I$(LIBGRAINS_GPU_SRCDIR)/CollisionDetection/include
CPPFLAGS += -I$(LIBGRAINS_GPU_SRCDIR)/Component/include
CPPFLAGS += -I$(LIBGRAINS_GPU_SRCDIR)/ContactForce/include
CPPFLAGS += -I$(LIBGRAINS_GPU_SRCDIR)/Geometry/include
CPPFLAGS += -I$(LIBGRAINS_GPU_SRCDIR)/PostProcessing/include
CPPFLAGS += -I$(LIBGRAINS_GPU_SRCDIR)/TimeIntegration/include
CPPFLAGS += -I$(LIBGRAINS_GPU_SRCDIR)/Utility/include
CPPFLAGS += -I$(LIBGRAINS_GPU_SRCDIR)/XML3/include
CPPFLAGS += $(PACIFIC_MPI_SKIP_CPPFLAGS) $(PACIFIC_XERCESC_CFLAGS) $(PACIFIC_ZLIB_CFLAGS) $(PACIFIC_CURAND_CFLAGS)
CPPFLAGS += $(PACIFIC_GTEST_CFLAGS) $(PACIFIC_GTEST_MAIN_CFLAGS)
CPPFLAGS += -DTEST_BUILD=1

LDLIBS += $(PACIFIC_GTEST_MAIN_LFLAGS) $(PACIFIC_GTEST_LFLAGS)
LDLIBS += $(PACIFIC_XERCESC_LFLAGS) $(PACIFIC_ZLIB_LFLAGS) $(PACIFIC_CUDART_LFLAGS) $(PACIFIC_CURAND_LFLAGS)
NVCC_LDLIBS := $(call to_nvcc_ldflags,$(LDLIBS))

# Target directories
TEST_GPU_BINDIR ?= $(GRAINS_GPU_BINDIR)
TEST_GPU_OBJDIR ?= $(GRAINS_GPU_OBJDIR)/tests

TARGET := $(TEST_GPU_BINDIR)/grains3dgpu_tests

# Sources
TEST_SRCS := \
	main.cpp \
	geometry/test_vector3.cpp \
	geometry/test_matrix3.cpp \
	geometry/test_quaternion.cpp \
	geometry/RotationMathTest.cpp \
	collision/test_gjk.cpp \
	collision/test_obb.cpp \
	collision/test_obc.cpp \
	collision/test_collision_detection.cpp \
	collision/test_neighborlist.cpp \
	contactForce/test_contact_table_host.cpp \
	contactForce/test_contact_table_device.cpp

TEST_SRCS := $(foreach f,$(TEST_SRCS),$(if $(wildcard $(f)),$(f),))
TEST_OBJS := $(patsubst %.cpp,$(TEST_GPU_OBJDIR)/%.o,$(TEST_SRCS))
TEST_DEPS := $(TEST_OBJS:.o=.d)

# Link directly with Grains object files to avoid CUDA TU/nvlink shared-boundary issues.
GRAINS_OBJS = $(filter-out $(LIBGRAINS_GPU_OBJDIR)/main.o,$(wildcard $(LIBGRAINS_GPU_OBJDIR)/*.o))

# -----------------------------------------------------------------------------
# Production rules
.PHONY: all run clean

all: $(TARGET)

$(TEST_GPU_BINDIR):
	@mkdir -p $@

$(TEST_GPU_OBJDIR):
	@mkdir -p $@

$(TEST_GPU_OBJDIR)/%.o: %.cpp | $(TEST_GPU_OBJDIR)
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCC_ARCHFLAGS) $(CPPFLAGS) $(NVCCFLAGS) -MMD -MP -MF $(@:.o=.d) -MT $@ -c $< -o $@

$(TARGET): $(TEST_OBJS) | $(TEST_GPU_BINDIR)
	$(NVCC) $(NVCC_ARCHFLAGS) -o $@ $(TEST_OBJS) $(GRAINS_OBJS) -L$(LIBGRAINS_GPU_LIBDIR) $(NVCC_LDLIBS)
	@echo -e '-----------------------\n''Tests grains3dgpu_tests built !\n''-----------------------'

run: $(TARGET)
	$(TARGET)

clean:
	rm -rf $(TEST_GPU_OBJDIR) $(TARGET)

-include $(TEST_DEPS)
