#!/usr/bin/env bash
set -euo pipefail

GIT_SHA="$(git rev-parse HEAD)"
GIT_DATE="$(git show -s --date=format:%Y%m%d --format=%cd HEAD)"

if git describe --tags --exact-match >/dev/null 2>&1; then
  SNAP=0
else
  SNAP=1
fi

mkdir -p "$HOME/rpmbuild/"{SPECS,SOURCES,SRPMS}

cp rpm/pacific.spec "$HOME/rpmbuild/SPECS/pacific.spec"

# fetch Source0 into SOURCES
spectool -g -R "$HOME/rpmbuild/SPECS/pacific.spec"

# build SRPM
rpmbuild -bs \
  --define "_topdir $HOME/rpmbuild" \
  --define "_sourcedir $HOME/rpmbuild/SOURCES" \
  --define "_srcrpmdir $HOME/rpmbuild/SRPMS" \
  --define "snap $SNAP" \
  --define "git_date $GIT_DATE" \
  --define "git_sha $GIT_SHA" \
  "$HOME/rpmbuild/SPECS/pacific.spec"
