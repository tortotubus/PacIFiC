#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: $0 <source.c>" >&2
  exit 1
fi

if [ -z "${PACIFIC_SRCDIR_ABS:-}" ]; then
  echo "Please source the pacific environment file first." >&2
  exit 1
fi

SRC="$1"
OUT="$(basename "${SRC%.*}")"

PACIFIC_OCTREE_ABS="${PACIFIC_SRCDIR_ABS}/Octree"

PACIFIC_THIRDPARTY_LIBDIR=${PACIFIC_THIRDPARTY_INSTALLDIR_ABS}/lib64

CC99="${PACIFIC_MPICC:-mpicc} -std=c99 -D_GNU_SOURCE=1" \
  "${PACIFIC_BASILISK_QCC}" \
  -disable-dimensions \
  -D_MPI=1 \
  -D_DEFAULT_SOURCE \
  -DDISPLAY=1 \
  -DMTRACE=2 \
  -DTRACE=2 \
  "$SRC" \
  -I"${PACIFIC_OCTREE_ABS}" \
  -L"${PACIFIC_THIRDPARTY_LIBDIR}" \
  -lfb_tiny \
  -lglutils \
  -ltinyrenderer \
  -lwsserver \
  -lm \
  ${PACIFIC_HDF5_CFLAGS} \
  ${PACIFIC_HDF5_LFLAGS} \
  ${PACIFIC_MPICC_LFLAGS} \
  -Wno-unused-result \
  -o "$OUT"
