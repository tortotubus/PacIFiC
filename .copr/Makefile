.PHONY: srpm

srpm:
	dnf -y install rpmdevtools rpm-build curl || true
	@if [ -z "$(outdir)" ]; then echo "ERROR: outdir is empty"; exit 2; fi
	mkdir -p "$(HOME)/rpmbuild/SRPMS"
	./.copr/build_srpm.sh
	cp -v "$(HOME)/rpmbuild/SRPMS/"*.src.rpm "$(outdir)/"
	ls -lh "$(outdir)/"*.src.rpm
