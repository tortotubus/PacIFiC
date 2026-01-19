#!/usr/bin/env bash
set -euo pipefail

GIT_SHA="$(git rev-parse HEAD)"
GIT_DATE="$(git show -s --date=format:%Y%m%d --format=%cd HEAD)"

# Find latest semver tag reachable from HEAD (accepts v1.2.3 or 1.2.3)
SEMVER_TAG="$(
  git tag --list --sort=-v:refname \
    | grep -E '^(v)?[0-9]+\.[0-9]+\.[0-9]+$' \
    | head -n1 || true
)"

BASEVER="${SEMVER_TAG#v}"
BASEVER="${BASEVER:-0.0.1}"

# SNAP=0 only if HEAD is exactly at a semver tag
if git describe --tags --exact-match 2>/dev/null | grep -Eq '^(v)?[0-9]+\.[0-9]+\.[0-9]+$'; then
  SNAP=0
  VERSION="$(git describe --tags --exact-match | sed 's/^v//')"
else
  SNAP=1
  VERSION="$BASEVER"
fi

mkdir -p "$HOME/rpmbuild/"{SPECS,SOURCES,SRPMS}
cp rpm/pacific.spec "$HOME/rpmbuild/SPECS/pacific.spec"

git archive --format=tar.gz HEAD > "$HOME/rpmbuild/SOURCES/pacific-${VERSION}.tar.gz"

spectool -g -R "$HOME/rpmbuild/SPECS/pacific.spec"

rpmbuild -bs \
  --define "_topdir $HOME/rpmbuild" \
  --define "_sourcedir $HOME/rpmbuild/SOURCES" \
  --define "_srcrpmdir $HOME/rpmbuild/SRPMS" \
  --define "version $VERSION" \
  --define "snap $SNAP" \
  --define "git_date $GIT_DATE" \
  --define "git_sha $GIT_SHA" \
  "$HOME/rpmbuild/SPECS/pacific.spec"
