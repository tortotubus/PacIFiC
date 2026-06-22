#include "ibm/IBConfig.h"
#include "ibm/IBNode.h"
#include "ibm/IBMacros.h"

/* Yang et al. (2009) smoothed 3-point kernel (IBM_stencil=11 in lagrangian_caps).
   Support ±1.5Δ. Normalised so that Σ w(rᵢ) = 1 in each dimension. */
static inline double weight_stencil_yang11 (double r) {
  if (r <= 0.5)
    return 0.75 - r*r;
  if (r <= 1.5)
    return 1.125 - 1.5*r + 0.5*r*r;
  return 0.;
}

static inline double ibm_kernel_distance_1d (double a, double b,
                                             bool periodic) {
  double d = a - b;
  if (periodic) {
    if (d > 0.5*L0)
      d -= L0;
    else if (d < -0.5*L0)
      d += L0;
  }
  return fabs (d);
}

macro peskin_cosine_kernel_gather_dimensionless (IBNode* node = node) {
  // bool ib_set_dirty = true;
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level (PESKIN_SUPPORT_RADIUS, node->depth, pos) {
    if (!is_local (cell))
      continue;
#else
  foreach_neighbor_coord (PESKIN_SUPPORT_RADIUS, pos) {
#endif
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x =
        ibm_kernel_distance_1d (d.x, cell_centre.x, Period.x) / Delta;
      // weight *= .25 * (1 + cos (.5 * pi * kernel_dist.x));  // Peskin cosine 4-pt
      weight *= weight_stencil_yang11 (kernel_dist.x);  // Yang smoothed 3-pt (IBM_stencil=11)
    }
    // clang-format off
    {...}
    // clang-format on
  }
}

macro peskin_cosine_kernel_spread_dimensionless (IBNode* node = node) {
  // bool ib_set_dirty = false;
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level (PESKIN_SUPPORT_RADIUS, node->depth, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x =
        ibm_kernel_distance_1d (d.x, cell_centre.x, Period.x) / Delta;
      // weight *= .25 * (1 + cos (.5 * pi * kernel_dist.x));  // Peskin cosine 4-pt
      weight *= weight_stencil_yang11 (kernel_dist.x);  // Yang smoothed 3-pt (IBM_stencil=11)
    }
    if (is_local (cell)) {
      // clang-format off
          {...}
      // clang-format on
    }
  }
#else
  foreach_neighbor_coord_nonlocal (PESKIN_SUPPORT_RADIUS, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x =
        ibm_kernel_distance_1d (d.x, cell_centre.x, Period.x) / Delta;
      // weight *= .25 * (1 + cos (.5 * pi * kernel_dist.x));  // Peskin cosine 4-pt
      weight *= weight_stencil_yang11 (kernel_dist.x);  // Yang smoothed 3-pt (IBM_stencil=11)
    }

    {
      // clang-format off
          {...}
      // clang-format on
    }
  }
#endif
}

macro peskin_cosine_kernel_gather (IBNode* node = node) {
  // bool ib_set_dirty = true;
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level (PESKIN_SUPPORT_RADIUS, node->depth, pos) {
    if (!is_local (cell))
      continue;
#else
  foreach_neighbor_coord (PESKIN_SUPPORT_RADIUS, pos) {
#endif
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = ibm_kernel_distance_1d (d.x, cell_centre.x, Period.x);
      if (kernel_dist.x <= Delta * PESKIN_SUPPORT_RADIUS) {
        weight *= (1. + cos (pi * kernel_dist.x / (Delta * PESKIN_SUPPORT_RADIUS))) / (2. * PESKIN_SUPPORT_RADIUS);
      } else {
        weight = 0;
      }
    }
    // clang-format off
    {...}
    // clang-format on
  }
}

macro peskin_cosine_kernel_spread (IBNode* node = node) {
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level (PESKIN_SUPPORT_RADIUS, node->depth, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = ibm_kernel_distance_1d (d.x, cell_centre.x, Period.x);
      if (kernel_dist.x <= Delta *  PESKIN_SUPPORT_RADIUS)
        weight *= (1. + cos (pi * kernel_dist.x / (Delta * PESKIN_SUPPORT_RADIUS))) / (2. * PESKIN_SUPPORT_RADIUS);
      else 
        weight = 0.;
    }

    // if (weight != 0.) {
      if (is_local (cell)) {
        // clang-format off
          {...}
        // clang-format on
      }
    // }
  }
#else
  foreach_neighbor_coord_nonlocal (PESKIN_SUPPORT_RADIUS, pos) {
    double weight = 1.0;
    coord kernel_dist = {0};
    coord cell_centre = {.x = x, .y = y, .z = z};
    foreach_dimension () {
      kernel_dist.x = ibm_kernel_distance_1d (d.x, cell_centre.x, Period.x);
      if (kernel_dist.x <= Delta *  PESKIN_SUPPORT_RADIUS)
        weight *= (1. + cos (pi * kernel_dist.x / (Delta * PESKIN_SUPPORT_RADIUS))) / (2. * PESKIN_SUPPORT_RADIUS);
      else 
        weight = 0.;
    }

    // if (weight != 0.) 
    {
      // clang-format off
          {...}
      // clang-format on
    }
  }
#endif
}
