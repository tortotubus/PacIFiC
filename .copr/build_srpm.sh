#!/usr/bin/env bash
set -euo pipefail

mkdir -p "$HOME/rpmbuild/"{SPECS,SOURCES,SRPMS}

cp rpm/pacific.spec "$HOME/rpmbuild/SPECS/pacific.spec"

# fetch Source0 into SOURCES
spectool -g -R "$HOME/rpmbuild/SPECS/pacific.spec"

# build SRPM
rpmbuild -bs \
  --define "_topdir $HOME/rpmbuild" \
  --define "_sourcedir $HOME/rpmbuild/SOURCES" \
  --define "_srcrpmdir $HOME/rpmbuild/SRPMS" \
  "$HOME/rpmbuild/SPECS/pacific.spec"
