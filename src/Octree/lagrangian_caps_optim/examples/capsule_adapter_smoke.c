#define RADIUS 0.1
#define MU 1.
#define E_S 1.

static int n_caps = 1;

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps_optim/capsule-ft.h"
#include "lagrangian_caps_optim/common-shapes-ft.h"

event init (i = 0) {
  CapsuleMeshPrototype* prototype = build_spherical_capsule_prototype(
    (_initialize_circular_capsule){
      .cap_es = E_S,
      .cap_radius = RADIUS,
      .level = 1
    });
  activate_spherical_capsule((_initialize_circular_capsule){
    .mesh = &CAPS(0),
    .prototype = prototype,
    .cap_es = E_S,
    .cap_radius = RADIUS,
    .cap_id = 0,
    .level = 1,
    .shift = {0.5, 0.5, 0.5}
  });
  capsule_mesh_prototype_release(prototype);
}

int main() {
  L0 = 1.;
  origin(0., 0., 0.);
  init_grid(8);
  capsule_manager_set_count(&allCaps, n_caps);
  run();
  return 0;
}
