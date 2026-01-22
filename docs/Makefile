# -------------------------------------------
# Makefile for building Doxygen + MkDocs docs
# -------------------------------------------

# Configuration
DOXYGEN       ?= doxygen
DOXYFILE      ?= DoxyGrains
MKDOCS        ?= mkdocs
SITE_DIR      ?= site
DOXY_OUTPUT   ?= docs/doxygen

# Optional intermediate directory (Doxygen HTML copied into MkDocs)
MKDOCS_DOXY_DIR := docs/api

# Default target
.PHONY: all
all: build

# -------------------------------------------
# Doxygen documentation
# -------------------------------------------

.PHONY: doxygen
doxygen:
	@echo ">>> Generating Doxygen documentation..."
	$(DOXYGEN) $(DOXYFILE)

build: clean doxygen 
	mkdocs build

develop: doxygen 
	mkdocs serve

clean: 
	rm -rf site pages/doxygen/grains_html