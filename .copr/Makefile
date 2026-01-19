.PHONY: srpm

install_buildreqs:
	dnf -y install rpmdevtools rpm-build git || true

test_srpm:
	mkdir -p test_output
	$(MAKE) -f .copr/Makefile srpm outdir="$(PWD)/test_output"

ensure_rpmbuild_dirs:
	@mkdir -p "$(HOME)/rpmbuild/SOURCES"
	@mkdir -p "$(HOME)/rpmbuild/SPECS"
	@mkdir -p "$(HOME)/rpmbuild/BUILD"
	@mkdir -p "$(HOME)/rpmbuild/RPMS"
	@mkdir -p "$(HOME)/rpmbuild/SRPMS"

srpm: install_buildreqs ensure_rpmbuild_dirs
	@if [ -z "$(outdir)" ]; then echo "ERROR: outdir is empty"; exit 2; fi
	./.copr/build_srpm.sh
	cp -v "$(HOME)/rpmbuild/SRPMS/"*.src.rpm "$(outdir)/"
	ls -lh "$(outdir)/"*.src.rpm
