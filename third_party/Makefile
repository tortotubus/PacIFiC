# Expects PACIFIC_BUILDDIR_ABS and PACIFIC_ROOT_ABS

# Derived config

MKDIR ?= mkdir
CP ?= cp
LN ?= ln

PACIFIC_THIRDPARTY_DIR := $(PACIFIC_ROOT_ABS)/third_party

.PHONY: all clean

all: third_party-basilisk third_party-hdf5 third_party-xercesc third_party-zlib

clean: 

third_party-builddir:
	@$(MKDIR) -p $(PACIFIC_THIRDPARTY_BUILDDIR_ABS)

third_party-installdir:
	@$(MKDIR) -p $(PACIFIC_THIRDPARTY_INSTALLDIR_ABS)

# =============================================================================
# BASILISK
# =============================================================================

PACIFIC_SUBMODULE_BASILISK := $(PACIFIC_THIRDPARTY_DIR)/basilisk/basilisk_submodule
PACIFIC_SUBMODULE_BASILISK_BINDIR := $(PACIFIC_THIRDPARTY_BUILDDIR_ABS)/basilisk
PACIFIC_SUBMODULE_BASILISK_INSTALLDIR := $(PACIFIC_THIRDPARTY_INSTALLDIR_ABS)
PACIFIC_SUBMODULE_BASILISK_CMAKE_DEFS ?= \
	-DCMAKE_INSTALL_PREFIX=$(PACIFIC_SUBMODULE_BASILISK_INSTALLDIR)

third_party-basilisk: third_party-basilisk-install

third_party-basilisk-configure: third_party-builddir
	@$(CMAKE) -S "$(PACIFIC_SUBMODULE_BASILISK)" -B "$(PACIFIC_SUBMODULE_BASILISK_BINDIR)" -G "Unix Makefiles" $(PACIFIC_SUBMODULE_BASILISK_CMAKE_DEFS)

third_party-basilisk-build: third_party-basilisk-configure
	@$(CMAKE) --build "$(PACIFIC_SUBMODULE_BASILISK_BINDIR)" --parallel

third_party-basilisk-install: third_party-basilisk-build third_party-installdir
	@$(CMAKE) --install "$(PACIFIC_SUBMODULE_BASILISK_BINDIR)"

# =============================================================================
# HDF5
# =============================================================================

PACIFIC_SUBMODULE_HDF5 := $(PACIFIC_THIRDPARTY_DIR)/hdf5/hdf5_submodule
PACIFIC_SUBMODULE_HDF5_BINDIR := $(PACIFIC_THIRDPARTY_BUILDDIR_ABS)/hdf5
PACIFIC_SUBMODULE_HDF5_INSTALLDIR := $(PACIFIC_THIRDPARTY_INSTALLDIR_ABS)
PACIFIC_SUBMODULE_HDF5_CMAKE_DEFS ?= \
	-DHDF5_ENABLE_PARALLEL=ON \
	-DBUILD_SHARED_LIBS=ON \
	-DHDF5_ALLOW_EXTERNAL_SUPPORT:STRING=GIT \
	-DHDF5_ENABLE_ZLIB_SUPPORT:BOOL=ON \
	-DHDF5_ENABLE_SZIP_SUPPORT:BOOL=ON \
	-DCMAKE_INSTALL_PREFIX=$(PACIFIC_SUBMODULE_HDF5_INSTALLDIR)

third_party-hdf5: third_party-hdf5-install

third_party-hdf5-configure: third_party-builddir
	@$(CMAKE) -S "$(PACIFIC_SUBMODULE_HDF5)" -B "$(PACIFIC_SUBMODULE_HDF5_BINDIR)" -G "Unix Makefiles" $(PACIFIC_SUBMODULE_HDF5_CMAKE_DEFS)

third_party-hdf5-build: third_party-hdf5-configure
	@$(CMAKE) --build "$(PACIFIC_SUBMODULE_HDF5_BINDIR)" --parallel

third_party-hdf5-install: third_party-hdf5-build third_party-installdir
	@$(CMAKE) --install "$(PACIFIC_SUBMODULE_HDF5_BINDIR)"

# =============================================================================
# ZLIB
# =============================================================================

PACIFIC_SUBMODULE_ZLIB := $(PACIFIC_THIRDPARTY_DIR)/zlib/zlib_submodule
PACIFIC_SUBMODULE_ZLIB_BINDIR := $(PACIFIC_THIRDPARTY_BUILDDIR_ABS)/zlib
PACIFIC_SUBMODULE_ZLIB_INSTALLDIR := $(PACIFIC_THIRDPARTY_INSTALLDIR_ABS)
PACIFIC_SUBMODULE_ZLIB_CMAKE_DEFS ?= \
	-DCMAKE_INSTALL_PREFIX=${PACIFIC_SUBMODULE_ZLIB_INSTALLDIR}

third_party-zlib: third_party-zlib-install

third_party-zlib-configure: third_party-builddir
	@$(CMAKE) -S "$(PACIFIC_SUBMODULE_ZLIB)" -B "$(PACIFIC_SUBMODULE_ZLIB_BINDIR)" -G "Unix Makefiles" $(PACIFIC_SUBMODULE_ZLIB_CMAKE_DEFS)

third_party-zlib-build: third_party-zlib-configure
	@$(CMAKE) --build "$(PACIFIC_SUBMODULE_ZLIB_BINDIR)" --parallel

third_party-zlib-install: third_party-zlib-build third_party-installdir
	@$(CMAKE) --install "$(PACIFIC_SUBMODULE_ZLIB_BINDIR)"


# =============================================================================
# XERCESC
# =============================================================================

PACIFIC_SUBMODULE_XERCESC := $(PACIFIC_THIRDPARTY_DIR)/xercesc/xercesc_submodule
PACIFIC_SUBMODULE_XERCESC_BINDIR := $(PACIFIC_THIRDPARTY_BUILDDIR_ABS)/xercesc
PACIFIC_SUBMODULE_XERCESC_INSTALLDIR := $(PACIFIC_THIRDPARTY_INSTALLDIR_ABS)
PACIFIC_SUBMODULE_XERCESC_CMAKE_DEFS ?= \
	-DCMAKE_INSTALL_PREFIX=${PACIFIC_SUBMODULE_XERCESC_INSTALLDIR}

third_party-xercesc: third_party-xercesc-install

third_party-xercesc-configure: third_party-builddir
	@$(CMAKE) -S "$(PACIFIC_SUBMODULE_XERCESC)" -B "$(PACIFIC_SUBMODULE_XERCESC_BINDIR)" -G "Unix Makefiles" $(PACIFIC_SUBMODULE_XERCESC_CMAKE_DEFS)

third_party-xercesc-build: third_party-xercesc-configure
	@$(CMAKE) --build "$(PACIFIC_SUBMODULE_XERCESC_BINDIR)" --parallel

third_party-xercesc-install: third_party-xercesc-build third_party-installdir
	@$(CMAKE) --install "$(PACIFIC_SUBMODULE_XERCESC_BINDIR)"
