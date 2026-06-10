/**
# Triangulated mesh of a capsule

This header is the integration point for the lagrangian capsule solver. Core
data structures and ownership utilities live in the split headers below.
*/

#pragma once

#include <string.h>
 
 

#ifndef CONSERVE_VOLUME
#define CONSERVE_VOLUME 1
#endif

// ============================================================================
// Globals
// ============================================================================

scalar Index_lagnode[];
vector Index_lag_id[];

#include "capsule-manager-ft.h"
#include "capsule-mesh-ft.h"
#include "capsule-mesh-prototype-ft.h"

// ============================================================================
// Macros
// ============================================================================

/* Capsule coordinates are material coordinates. The IBM layer owns wrapping
   node positions into the Eulerian periodic box for gather/spread support. */
#define ACROSS_PERIODIC(a, b) false
#define PERIODIC_1DIST(a, b) ((a) - (b))
#define GENERAL_1DIST(a, b) ((a) - (b))
#define PERIODIC_1DAVG(a, b) ((a) + (b))
#define GENERAL_1DAVG(a, b) ((a) + (b))
#define GENERAL_SQNORM(a, b)                                                   \
  (sq(GENERAL_1DIST(a.x, b.x)) + sq(GENERAL_1DIST(a.y, b.y)) +                 \
   sq(GENERAL_1DIST(a.z, b.z)))

#define cnorm(a) (sqrt(sq(a.x) + sq(a.y) + sq(a.z)))
#define cdot(a, b) (a.x * b.x + a.y * b.y + a.z * b.z)

#include "geometry-ft.h"

#if CONSERVE_VOLUME
#include "volume-conservation-ft.h"
#endif

#include "capsule-ibm-adapter.h"
#include "plugins-ft.h"
#include "dump-ft.h"

// ============================================================================
// Function Definitions
// ============================================================================

/*!
 * @brief
 */
trace void synchronize(scalar *list) {
  for (scalar s in list)
    s.dirty = true;
  boundary(list);
}
