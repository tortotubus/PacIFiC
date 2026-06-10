/**
# Bending force for front-tracking membranes


*/
#define _BENDING_FT 1

#ifndef E_B
#define E_B 1.
#endif

#ifndef REF_CURV
#define REF_CURV 1
#endif

#ifndef GLOBAL_REF_CURV
#if REF_CURV
#define GLOBAL_REF_CURV 1
#else
#define GLOBAL_REF_CURV 0
#endif
#endif

#if GLOBAL_REF_CURV
#ifndef C0
/** $c_0/a = -2.09$ is the typical reference curvature of a red blood cell. */
// #define C0 (-2.09/RADIUS)
#endif
#endif

#ifndef LINEAR_BENDING
#define LINEAR_BENDING 0
#endif

#include "curvature-ft.h" 

static void comp_bending_force(CapsuleMesh* mesh) {
  if (!mesh || !mesh->isactive || !capsule_mesh_is_local_owner(mesh) ||
      !mesh->nodes)
    return;

  comp_curvature(mesh);
  for (int j = 0; j < mesh->nln; j++) {
#if (!LINEAR_BENDING)
    double curv = mesh->nodes[j].curv;
    double rcurv = mesh->nodes[j].ref_curv;
    double gcurv = mesh->nodes[j].gcurv;
#endif
    double lbcurv = laplace_beltrami(mesh, j, true);
    double Eb = ND_EB * mesh->cap_es * sq(mesh->cap_radius);
    double bending_surface_force =
        2 * Eb *
        (lbcurv
#if (!LINEAR_BENDING)
         + 2 * (curv - rcurv) * (sq(curv) - gcurv + rcurv * curv)
#endif
        );
    double area = compute_node_area(mesh, j);
    foreach_dimension() mesh->nodes[j].lagForce.x +=
        mesh->nodes[j].normal.x * bending_surface_force * area;
  }
}

#if CAPSULE_LEGACY_EVENTS
event acceleration(i++) {
  /*Compute borders of the current proc*/
  compute_proc_borders(&proc_max, &proc_min);

  for (int i = 0; i < allCaps.nbcaps; i++) {
    if (CAPS(i).isactive) {
      if (is_capsule_in_boundingbox(proc_max, proc_min, &CAPS(i))) {
        comp_bending_force(&CAPS(i));
      }
    }
  }
}
#endif
