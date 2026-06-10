#pragma once

#include "ibm/IBKernels.h"

// ============================================================================
// Globals
// ============================================================================

vector capsule_viscosity_grid_gradient[];
scalar capsule_viscosity_divergence[];

// ============================================================================
// Macros
// ============================================================================

#ifndef CAPSULE_VISCOSITY_STENCIL
#define CAPSULE_VISCOSITY_STENCIL 14
#endif

// ============================================================================
// Function Declarations
// ============================================================================

static double capsule_viscosity_weight_stencil(double dist);
static void capsule_viscosity_spread_node_area_normal(IBNode *node);
static void capsule_viscosity_set_wall_boundary_conditions(scalar indicator);
static void capsule_viscosity_add_ibmesh_grid_gradient(IBMesh *ibmesh);
static void capsule_viscosity_refresh_ibmesh_fields(IBMesh *ibmesh);
static void capsule_viscosity_construct_indicator(scalar indicator);
static void capsule_viscosity_set_from_indicator(face vector viscosity,
                                                 scalar indicator,
                                                 double outside_viscosity,
                                                 double inside_viscosity);

// ============================================================================
// Function Definitions
// ============================================================================

/*!
 * @brief
 */
static double capsule_viscosity_weight_stencil(double dist) {
  double adist = fabs(dist);

  if (CAPSULE_VISCOSITY_STENCIL == 11) {
    if (adist <= 0.5)
      return 3. / 4. - sq(dist);
    if (adist <= 1.5)
      return 9. / 8. - 3. * adist / 2. + sq(dist) / 2.;
    return 0.;
  }

  if (adist <= 0.5)
    return 3. / 8. + pi / 32. - sq(dist) / 4.;
  if (adist <= 1.5)
    return 1. / 4. +
           (1. - adist) * sqrt(-2. + 8. * adist - 4. * sq(dist)) / 8. -
           asin(sqrt(2.) * (adist - 1.)) / 8.;
  if (adist <= 2.5)
    return 17. / 16. - pi / 64. - 3. * adist / 4. + sq(dist) / 8. +
           (adist - 2.) * sqrt(-14. + 16. * adist - 4. * sq(dist)) / 16. +
           asin(sqrt(2.) * (adist - 2.)) / 16.;
  return 0.;
}

/*!
 * @brief
 */
static void capsule_viscosity_spread_node_area_normal(IBNode *node) {
  IBNODE_VARIABLES();
#if TREE
  foreach_neighbor_coord_level(3, node->depth, pos) {
    if (!is_local(cell))
      continue;
#else
  foreach_neighbor_coord_nonlocal(3, pos) {
#endif
    coord dist = {0};
    double weight = 1.;
    foreach_dimension() {
      dist.x = GENERAL_1DIST(x, pos.x);
      weight *= sq(dist.x) <= sq(2. * Delta)
                    ? capsule_viscosity_weight_stencil(dist.x / Delta)
                    : 0.;
    }
#if dimension == 2
    double kernel_volume = sq(Delta);
#else
    double kernel_volume = cube(Delta);
#endif
    foreach_dimension() capsule_viscosity_grid_gradient.x[] -=
        weight * ibval(capsule_viscosity_area_normal.x) / kernel_volume;
  }
}

/*!
 * @brief
 */
static void capsule_viscosity_set_wall_boundary_conditions(scalar indicator) {
  if (u.x.boundary[left] != periodic_bc) {
    indicator[left] = dirichlet(0.);
    indicator[right] = dirichlet(0.);
    capsule_viscosity_divergence[left] = dirichlet(0.);
    capsule_viscosity_divergence[right] = dirichlet(0.);
  }
#if dimension >= 2
  if (u.y.boundary[top] != periodic_bc) {
    indicator[top] = dirichlet(0.);
    indicator[bottom] = dirichlet(0.);
    capsule_viscosity_divergence[top] = dirichlet(0.);
    capsule_viscosity_divergence[bottom] = dirichlet(0.);
  }
#endif
#if dimension >= 3
  if (u.z.boundary[front] != periodic_bc) {
    indicator[front] = dirichlet(0.);
    indicator[back] = dirichlet(0.);
    capsule_viscosity_divergence[front] = dirichlet(0.);
    capsule_viscosity_divergence[back] = dirichlet(0.);
  }
#endif
}

/*!
 * @brief
 */
static void capsule_viscosity_add_ibmesh_grid_gradient(IBMesh *ibmesh) {
  if (!ibmesh || ibmesh->model.type != IB_MODEL_VELOCITY_COUPLED)
    return;

#if dimension == 3
  for (size_t i = 0; i < ibmesh->nodes.size; i++) {
    IBNode *node = ibmesh->nodes.ptrs[i];
    capsule_viscosity_spread_node_area_normal(node);
  }
#endif
}

/*!
 * @brief
 */
static void capsule_viscosity_refresh_ibmesh_fields(IBMesh *ibmesh) {
  if (!ibmesh || ibmesh->model.type != IB_MODEL_VELOCITY_COUPLED)
    return;

  if (ibmesh->sparse_managed)
    return;

  CapsuleIBMProxy *proxy = (CapsuleIBMProxy *)ibmesh->model.ctx;
  capsule_ibm_sync_viscosity_area_normals(proxy, ibmesh);
}

/*!
 * @brief
 */
static void capsule_viscosity_construct_indicator(scalar indicator) {
  foreach () {
    if (cm[] > 1.e-20) {
      capsule_viscosity_divergence[] = 0.;
      foreach_dimension() capsule_viscosity_grid_gradient.x[] = 0.;
    }
  }

  foreach_ibmesh() {
    capsule_viscosity_refresh_ibmesh_fields(mesh);
    capsule_viscosity_add_ibmesh_grid_gradient(mesh);
  }

  foreach_dimension() capsule_viscosity_grid_gradient.x.dirty = true;
  boundary((scalar *){capsule_viscosity_grid_gradient});

  foreach () {
    if (cm[] > 1.e-20) {
      capsule_viscosity_divergence[] = 0.;
      foreach_dimension() capsule_viscosity_divergence[] +=
          (capsule_viscosity_grid_gradient.x[1] -
           capsule_viscosity_grid_gradient.x[-1]) /
          (2. * Delta);
    }
  }

  poisson(indicator, capsule_viscosity_divergence);

  foreach ()
    if (cm[] > 1.e-20)
      indicator[] = clamp(indicator[], 0., 1.);

  indicator.dirty = true;
  boundary({indicator});
}

/*!
 * @brief
 */
static void capsule_viscosity_set_from_indicator(face vector viscosity,
                                                 scalar indicator,
                                                 double outside_viscosity,
                                                 double inside_viscosity) {
  foreach_face() if (fm.x[] > 1.e-20) viscosity.x[] =
      (outside_viscosity + (inside_viscosity - outside_viscosity) * 0.5 *
                               (indicator[] + indicator[-1])) *
      fm.x[];
}
