# .copr/Makefile

NAME      := pacific
TOPDIR    := $(CURDIR)/rpmbuild
SPECDIR   := $(TOPDIR)/SPECS
SOURCES   := $(TOPDIR)/SOURCES
SRPMS     := $(TOPDIR)/SRPMS

# Path to spec in your repo
SPEC_IN   := rpm/$(NAME).spec

# COPR runs: make -f .copr/Makefile srpm
.PHONY: srpm
srpm: prep spec sources
	@mkdir -p $(SRPMS)
	rpmbuild -bs \
	  --define "_topdir $(TOPDIR)" \
	  --define "_srcrpmdir $(SRPMS)" \
	  --define "_sourcedir $(SOURCES)" \
	  $(SPECDIR)/$(NAME).spec
	@echo "SRPMs in: $(SRPMS)"
	@ls -1 $(SRPMS)/*.src.rpm

.PHONY: prep
prep:
	@mkdir -p $(SPECDIR) $(SOURCES) $(SRPMS)

# Copy spec into rpmbuild tree (no templating)
.PHONY: spec
spec:
	@cp -f $(SPEC_IN) $(SPECDIR)/$(NAME).spec

# Create tarball matching Source0 = %{name}-%{version}.tar.gz
.PHONY: sources
sources:
	@set -euo pipefail; \
	ver="$$(rpmspec -q --qf '%{VERSION}\n' $(SPECDIR)/$(NAME).spec | head -n 1)"; \
	tarname="$(NAME)-$${ver}.tar.gz"; \
	topdir="$(NAME)-$${ver}"; \
	echo "ver=$${ver}"; \
	echo "tarname=$${tarname}"; \
	echo "topdir=$${topdir}"; \
	rm -rf "$(SOURCES)/$${topdir}"; \
	mkdir -p "$(SOURCES)/$${topdir}"; \
	rsync -a --delete \
	  --exclude '.git' \
	  --exclude 'rpmbuild' \
	  --exclude '*.swp' \
	  --exclude '*.swo' \
	  ./ "$(SOURCES)/$${topdir}/"; \
	( cd "$(SOURCES)" && tar -czf "$${tarname}" "$${topdir}" ); \
	test -f "$(SOURCES)/$${tarname}"; \
	echo "Wrote $(SOURCES)/$${tarname}"; \
	rm -rf "$(SOURCES)/$${topdir}"
