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
# inside .copr/Makefile
sources:
	ver=0.0.1; \
	tarname=pacific-$$ver.tar.gz; \
	topdir=pacific-$$ver; \
	rm -rf "$$topdir"; \
	mkdir -p "$$topdir"; \
	tar --exclude-vcs --exclude='./.git' -cf - . | tar -C "$$topdir" -xf - ; \
	tar -czf "$$tarname" "$$topdir"; \
	mv "$$tarname" "$(outdir)/"

