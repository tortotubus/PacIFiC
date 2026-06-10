#pragma once

#include "ibm/IBKernels.h"

// ============================================================================
// Globals
// ============================================================================

vector capsule_viscosity_grid_gradient[];
scalar capsule_viscosity_divergence[];

// ============================================================================
// Function Declarations
// ============================================================================

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
static void capsule_viscosity_spread_node_area_normal(IBNode *node) {
  peskin_cosine_kernel_spread_dimensionless(node) {
    // area_normal has units L^2 in 3D; weight/dv() gives the grid gradient
    // units 1/L.
    foreach_dimension() capsule_viscosity_grid_gradient.x[] -=
        weight * ibval(capsule_viscosity_area_normal.x) / dv();
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
