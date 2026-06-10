#define RADIUS 0.05
#define MU 1.
#define E_S 1.

static int n_caps = 5000;
static int lag_level = 4;

#include <sys/resource.h>

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps_optim/capsule-ft.h"
#include "lagrangian_caps_optim/skalak-ft.h"
#include "lagrangian_caps_optim/common-shapes-ft.h"

static long rss_bytes(void) {
  struct rusage usage;
  getrusage(RUSAGE_SELF, &usage);
  return usage.ru_maxrss*1024L;
}

static long current_rss_bytes(void) {
  FILE* fp = fopen("/proc/self/statm", "r");
  if (!fp)
    return -1;
  long pages = 0;
  int matched = fscanf(fp, "%*s %ld", &pages);
  fclose(fp);
  if (matched != 1)
    return -1;
  return pages*4096L;
}

static size_t capsule_prototype_bytes(CapsuleMeshPrototype* prototype) {
  size_t bytes = sizeof(*prototype);
  bytes += (size_t) prototype->nln*sizeof(CapNode);
  bytes += (size_t) prototype->nle*(sizeof(EdgeTopology) + sizeof(EdgeState));
#if dimension > 2
  bytes += (size_t) prototype->nlt*(sizeof(TriangleTopology) + sizeof(TriangleState));
#endif
  return bytes;
}

static size_t capsule_state_bytes(CapsuleMesh* mesh) {
  size_t bytes = sizeof(*mesh);
  bytes += (size_t) mesh->nln*sizeof(CapNode);
  bytes += (size_t) mesh->nle*sizeof(EdgeState);
#if dimension > 2
  bytes += (size_t) mesh->nlt*sizeof(TriangleState);
#endif
  return bytes;
}

static size_t capsule_unshared_bytes(CapsuleMesh* mesh) {
  size_t bytes = sizeof(*mesh);
  bytes += (size_t) mesh->nln*sizeof(CapNode);
  bytes += (size_t) mesh->nle*sizeof(Edge);
#if dimension > 2
  bytes += (size_t) mesh->nlt*sizeof(Triangle);
#endif
  return bytes;
}

event init(i = 0) {
  long rss_before = rss_bytes();
  long current_before = current_rss_bytes();
  CapsuleMeshPrototype* prototype = build_spherical_capsule_prototype(
    (_initialize_circular_capsule){
      .cap_es = E_S,
      .cap_radius = RADIUS,
      .level = lag_level
    });
  long rss_after_prototype = rss_bytes();
  long current_after_prototype = current_rss_bytes();

  for (int i = 0; i < n_caps; i++) {
    coord shift = {
      0.15 + 0.08*i,
      0.5,
      0.5
    };
    activate_spherical_capsule((_initialize_circular_capsule){
      .mesh = &CAPS(i),
      .prototype = prototype,
      .cap_es = E_S,
      .cap_radius = RADIUS,
      .cap_id = i,
      .cap_type = 0,
      .level = lag_level,
      .shift = shift
    });
  }
  long rss_after_capsules = rss_bytes();
  long current_after_capsules = current_rss_bytes();

  assert(CAPS(0).prototype != NULL);
  for (int i = 1; i < n_caps; i++)
    assert(CAPS(i).prototype == CAPS(0).prototype);

  assert(CAPS(0).prototype->refcount == n_caps + 1);
  assert(CAPS(0).nodes != CAPS(1).nodes);
  assert(fabs(CAPS(0).nodes[0].pos.x - CAPS(1).nodes[0].pos.x) > 1.e-12);

  size_t proto_bytes = capsule_prototype_bytes(prototype);
  size_t state_bytes = capsule_state_bytes(&CAPS(0));
  size_t old_bytes = capsule_unshared_bytes(&CAPS(0));
  size_t shared_total = proto_bytes + (size_t) n_caps*state_bytes;
  size_t unshared_total = (size_t) n_caps*old_bytes;

  fprintf(stderr,
    "prototype smoke: n_caps=%d level=%d nodes=%d edges=%d triangles=%d refcount=%d\n",
    n_caps, lag_level, CAPS(0).nln, CAPS(0).nle, CAPS(0).nlt,
    CAPS(0).prototype->refcount);
  fprintf(stderr,
    "estimated bytes: prototype=%zu per_cap_state=%zu old_per_cap=%zu shared_total=%zu old_total=%zu saved=%zu\n",
    proto_bytes, state_bytes, old_bytes, shared_total, unshared_total,
    unshared_total > shared_total ? unshared_total - shared_total : 0);
  fprintf(stderr,
    "rss max bytes: before=%ld after_prototype=%ld after_capsules=%ld delta_capsules=%ld\n",
    rss_before, rss_after_prototype, rss_after_capsules,
    rss_after_capsules - rss_after_prototype);
  fprintf(stderr,
    "rss current bytes: before=%ld after_prototype=%ld after_capsules=%ld delta_capsules=%ld\n",
    current_before, current_after_prototype, current_after_capsules,
    current_after_capsules - current_after_prototype);

  capsule_mesh_prototype_release(prototype);
}

event end(i = 0) {
  fprintf(stderr, "rss before free_all_caps: max=%ld current=%ld\n",
    rss_bytes(), current_rss_bytes());
  free_all_caps(&allCaps);
  fprintf(stderr, "rss after free_all_caps: max=%ld current=%ld\n",
    rss_bytes(), current_rss_bytes());
}

int main() {
  L0 = 1.;
  origin(0., 0., 0.);
  init_grid(8);
  capsule_manager_set_count(&allCaps, n_caps);
  run();
  return 0;
}
