#!/usr/bin/env bash
set -euo pipefail

TOPDIR="${HOME}/rpmbuild"
SPEC_SRC="rpm/pacific.spec"
SPEC_DST="${TOPDIR}/SPECS/pacific.spec"
SOURCES="${TOPDIR}/SOURCES"
SRPMS="${TOPDIR}/SRPMS"

GIT_SHA="$(git rev-parse --short=12 HEAD)"
GIT_DATE="$(git show -s --date=format:%Y%m%d --format=%cd HEAD)"

# Latest semver tag (v1.2.3 or 1.2.3); fallback 0.0.1
SEMVER_TAG="$(
  git tag --list --sort=-v:refname \
    | grep -E '^(v)?[0-9]+\.[0-9]+\.[0-9]+$' \
    | head -n1 || true
)"
BASEVER="${SEMVER_TAG#v}"
BASEVER="${BASEVER:-0.0.1}"

# Determine whether HEAD is exactly at a semver tag
EXACT_TAG="$(git describe --tags --exact-match 2>/dev/null || true)"
if printf '%s' "$EXACT_TAG" | grep -Eq '^(v)?[0-9]+\.[0-9]+\.[0-9]+$'; then
  VERSION="${EXACT_TAG#v}"
  RELEASE="1%{?dist}"
else
  VERSION="$BASEVER"
  RELEASE="0.${GIT_DATE}git${GIT_SHA}%{?dist}"
fi

mkdir -p "${TOPDIR}"/{SPECS,SOURCES,SRPMS}
cp -f "$SPEC_SRC" "$SPEC_DST"

# Remove any previous injected block so reruns are idempotent
sed -i '/^# BEGIN AUTOGEN$/,/^# END AUTOGEN$/d' "$SPEC_DST"

# Prepend only version+release macros
tmp="$(mktemp)"
cat > "$tmp" <<EOF
# BEGIN AUTOGEN
%global version ${VERSION}
%global release ${RELEASE}
# END AUTOGEN

EOF
cat "$SPEC_DST" >> "$tmp"
mv -f "$tmp" "$SPEC_DST"

# Create Source0 tarball with the expected top-level directory
git archive --format=tar.gz --prefix="pacific-${VERSION}/" HEAD \
  > "${SOURCES}/pacific-${VERSION}.tar.gz"

# Build SRPM (no --define version/release needed)
rpmbuild -bs \
  --define "_topdir ${TOPDIR}" \
  --define "_sourcedir ${SOURCES}" \
  --define "_srcrpmdir ${SRPMS}" \
  "$SPEC_DST"
