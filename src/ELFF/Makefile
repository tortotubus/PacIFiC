# Compilation and linking flags
CXX = $(if $(PACIFIC_MPICXX),$(PACIFIC_MPICXX),mpicxx)
CC = $(if $(PACIFIC_MPICC),$(PACIFIC_MPICC),mpicc)
RM ?= rm

ELFF_USE_MPI ?= 1
ELFF_USE_HDF5 ?= 1
ELFF_USE_EIGEN3 ?= 1
ELFF_PRECISION ?= double
ELFF_USE_EXCEPTIONS ?= 0
ELFF_USE_OPENMP ?= 0

ELFF_HDF5_CFLAGS ?= $(PACIFIC_HDF5_CFLAGS)
ELFF_HDF5_LFLAGS ?= $(PACIFIC_HDF5_LFLAGS)
ELFF_EIGEN_CFLAGS ?= $(PACIFIC_EIGEN_CFLAGS)

ifeq ($(strip $(ELFF_HDF5_CFLAGS)),)
ELFF_HDF5_CFLAGS := $(shell pkg-config --cflags hdf5 2>/dev/null)
endif

ifeq ($(strip $(ELFF_HDF5_LFLAGS)),)
ELFF_HDF5_LFLAGS := $(shell pkg-config --libs hdf5 2>/dev/null)
endif

ifeq ($(strip $(ELFF_HDF5_LFLAGS)),)
ELFF_HDF5_LFLAGS := -lhdf5_hl -lhdf5
endif

ifeq ($(strip $(ELFF_EIGEN_CFLAGS)),)
ELFF_EIGEN_CFLAGS := $(shell pkg-config --cflags eigen3 2>/dev/null)
endif

CPPFLAGS += $(PACIFIC_MPI_SKIP_CPPFLAGS)
CXXFLAGS += -O3 -g -MMD -MP -fPIC -std=c++17
LDFLAGS += -shared -Wl,-soname,libELFF.so

ifeq ($(filter 1 true TRUE yes YES on ON,$(ELFF_USE_OPENMP)),)
else
CPPFLAGS += -DELFF_WITH_OPENMP=1
CXXFLAGS += -fopenmp
LDLIBS += -fopenmp
endif

ifeq ($(filter 1 true TRUE yes YES on ON,$(ELFF_USE_HDF5)),)
else
CPPFLAGS += $(ELFF_HDF5_CFLAGS)
LDLIBS += $(ELFF_HDF5_LFLAGS)
endif

ifeq ($(filter 1 true TRUE yes YES on ON,$(ELFF_USE_EIGEN3)),)
else
CPPFLAGS += $(ELFF_EIGEN_CFLAGS)
endif

# Target directories
OBJDIR := $(ELFF_OBJDIR)
LIBDIR := $(ELFF_LIBDIR)
INCDIR := $(ELFF_INCDIR)
SHAREDIR := $(ELFF_SHAREDIR)

ifeq ($(strip $(OBJDIR)),)
OBJDIR := $(CURDIR)/build/obj
endif

ifeq ($(strip $(LIBDIR)),)
LIBDIR := $(CURDIR)/build/lib
endif

ifeq ($(strip $(INCDIR)),)
INCDIR := $(CURDIR)/build/include
endif

ifeq ($(strip $(SHAREDIR)),)
SHAREDIR := $(CURDIR)/build/share/ELFF
endif

LIBSO := libELFF.so
LIBSO_WP := $(LIBDIR)/$(LIBSO)
CONFIG_HEADER := $(INCDIR)/elff/config/_config.hpp

# Source directories. This mirrors the ELFF CMake library source set.
SRCDIRS := \
	c/general \
	c/models/beam \
	c/models/capsule \
	c/models/ibm \
	c/models/rigidbody \
	config \
	general \
	io/C \
	io/CXX \
	models/beam \
	models/ibm \
	models/rigidbody \
	mpi \
	fem

# Public header roots. Keep this broad so generated examples can include the
# same installed-style headers that CMake installs.
HEADER_ROOTS := \
	c \
	config \
	fem \
	general \
	io \
	models \
	mpi

SRCFILES := $(foreach dir,$(SRCDIRS),$(wildcard $(dir)/*.cpp))
HEADERFILES := \
	$(foreach dir,$(HEADER_ROOTS),$(shell find $(dir) -type f \( -name '*.h' -o -name '*.hpp' \) 2>/dev/null))
OBJS := $(SRCFILES:%.cpp=$(OBJDIR)/%.o)
DEPS := $(OBJS:.o=.d)
INSTALLED_HEADERS := $(HEADERFILES:%=$(INCDIR)/elff/%)

CPPFLAGS += -I$(INCDIR)

.PHONY: all headers config MESS clean FORCE

all: MESS $(LIBSO_WP) CPINC

MESS:
	@echo -e '************************************************************\n''ELFF compiled with' $(CXX) $(CPPFLAGS) $(CXXFLAGS)'\n************************************************************'

headers: CPINC
config: $(CONFIG_HEADER)

CPINC: $(CONFIG_HEADER) $(INSTALLED_HEADERS)

$(CONFIG_HEADER): FORCE | $(INCDIR)
	@mkdir -p $(dir $@)
	@tmp="$@.tmp"; \
	{ \
	  echo '#pragma once'; \
	  echo '#define ELFF_VERSION 1'; \
	  echo '#define ELFF_VERSION_STRING "0.0.1"'; \
	  echo '#define ELFF_VERSION_TYPE_RELEASE 0'; \
	  echo '#define ELFF_VERSION_TYPE_DEVELOPMENT 1'; \
	  echo '#define ELFF_VERSION_MAJOR ((ELFF_VERSION)/10000)'; \
	  echo '#define ELFF_VERSION_MINOR (((ELFF_VERSION)/100)%100)'; \
	  echo '#define ELFF_VERSION_PATCH ((ELFF_VERSION)%100)'; \
	  echo '#define ELFF_SOURCE_DIR "$(CURDIR)"'; \
	  echo '#define ELFF_INSTALL_DIR "$(INCDIR)/.."'; \
	  echo '#define ELFF_SHARED_BUILD'; \
	  if [ "$(ELFF_USE_MPI)" = "1" ] || [ "$(ELFF_USE_MPI)" = "true" ] || [ "$(ELFF_USE_MPI)" = "TRUE" ] || [ "$(ELFF_USE_MPI)" = "yes" ] || [ "$(ELFF_USE_MPI)" = "YES" ] || [ "$(ELFF_USE_MPI)" = "on" ] || [ "$(ELFF_USE_MPI)" = "ON" ]; then echo '#define ELFF_USE_MPI'; fi; \
	  if [ "$(ELFF_USE_OPENMP)" = "1" ] || [ "$(ELFF_USE_OPENMP)" = "true" ] || [ "$(ELFF_USE_OPENMP)" = "TRUE" ] || [ "$(ELFF_USE_OPENMP)" = "yes" ] || [ "$(ELFF_USE_OPENMP)" = "YES" ] || [ "$(ELFF_USE_OPENMP)" = "on" ] || [ "$(ELFF_USE_OPENMP)" = "ON" ]; then echo '#define ELFF_USE_OPENMP'; fi; \
	  if [ "$(ELFF_USE_HDF5)" = "1" ] || [ "$(ELFF_USE_HDF5)" = "true" ] || [ "$(ELFF_USE_HDF5)" = "TRUE" ] || [ "$(ELFF_USE_HDF5)" = "yes" ] || [ "$(ELFF_USE_HDF5)" = "YES" ] || [ "$(ELFF_USE_HDF5)" = "on" ] || [ "$(ELFF_USE_HDF5)" = "ON" ]; then echo '#define ELFF_USE_HDF5'; fi; \
	  if [ "$(ELFF_USE_EIGEN3)" = "1" ] || [ "$(ELFF_USE_EIGEN3)" = "true" ] || [ "$(ELFF_USE_EIGEN3)" = "TRUE" ] || [ "$(ELFF_USE_EIGEN3)" = "yes" ] || [ "$(ELFF_USE_EIGEN3)" = "YES" ] || [ "$(ELFF_USE_EIGEN3)" = "on" ] || [ "$(ELFF_USE_EIGEN3)" = "ON" ]; then echo '#define ELFF_USE_EIGEN3'; fi; \
	  if [ "$(ELFF_PRECISION)" = "single" ] || [ "$(ELFF_PRECISION)" = "Single" ] || [ "$(ELFF_PRECISION)" = "SINGLE" ]; then echo '#define ELFF_USE_SINGLE'; else echo '#define ELFF_USE_DOUBLE'; fi; \
	  if [ "$(ELFF_USE_EXCEPTIONS)" = "1" ] || [ "$(ELFF_USE_EXCEPTIONS)" = "true" ] || [ "$(ELFF_USE_EXCEPTIONS)" = "TRUE" ] || [ "$(ELFF_USE_EXCEPTIONS)" = "yes" ] || [ "$(ELFF_USE_EXCEPTIONS)" = "YES" ] || [ "$(ELFF_USE_EXCEPTIONS)" = "on" ] || [ "$(ELFF_USE_EXCEPTIONS)" = "ON" ]; then echo '#define ELFF_USE_EXCEPTIONS'; fi; \
	} > "$$tmp"; \
	if ! cmp -s "$$tmp" "$@"; then mv "$$tmp" "$@"; else rm -f "$$tmp"; fi

$(INCDIR)/elff/%: % | $(INCDIR)
	@mkdir -p $(dir $@)
	@cp -f $< $@

$(LIBSO_WP): $(OBJS) | $(LIBDIR)
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LDLIBS)
	@echo -e '-------------------------------\n''Library libELFF.so generated !\n''-------------------------------'

$(OBJDIR)/%.o: %.cpp | CPINC
	@mkdir -p $(dir $@)
	@echo 'Compiling '$<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(LIBDIR):
	@mkdir -p $@

$(INCDIR):
	@mkdir -p $@

clean:
	$(RM) -rf "$(OBJDIR)" "$(LIBSO_WP)" "$(INCDIR)/elff" "$(SHAREDIR)"

-include $(DEPS)
