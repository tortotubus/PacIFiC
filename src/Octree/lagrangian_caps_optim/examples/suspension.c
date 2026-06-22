/**
# Rheology of a suspension of deformable capsules in shear flow

Optimized-capsule analogue of the Basilisk sandbox suspension case.  The case
uses prototype-backed spherical capsules, the IBM mesh manager, restart sidecar
checkpointing, HDF5 VTK output, and scalar rheology diagnostics.

Reference case:
https://basilisk.fr/sandbox/huet/cases/lagrangian_caps/suspension.c?raw
*/

#include <stdbool.h>
#include <stdlib.h>
#include <time.h>

#define PI 3.14159265358979323846
#define REF_CURV 1
#define GLOBAL_REF_CURV 1
#define ADVECT_LAG_RK2 0
#define LUBR_FORCE 0
#define LUBR_VEL 0
#define JACOBI 1

#ifndef RESTART_CASE
#define RESTART_CASE 0
#endif

#define DEFAULT_BASE_PATH "suspension_data"

static double domain_length = 1.;
static double volume_fraction = 0.05;
static double radius = 0.;
static double shear_rate = 1.;
static double reynolds_number = 10.;
static double fluid_viscosity = 1.;
static double capsule_inside_viscosity = 1.;
static double capillary_number = 0.1;
static double bending_ratio = 0.;
static double reference_curvature = 0.;
static double t_end = 100.;
static double dt_max = 1.e-3;
static double solver_tolerance = 1.e-6;
static double velocity_tolerance = 0.01;
static int min_level = 5;
static int max_level = 8;
static int lag_level = 4;
static int n_caps = 50;
static int random_seed = 1;
static int output_freq_stats = 5;
static int output_freq_vtk = 5;
static int output_freq_dump = 100;
static bool stokes_flow = false;
static bool capsule_viscosity_enabled = false;
char *base_path = DEFAULT_BASE_PATH;

#define RADIUS radius
#define SHEAR_RATE shear_rate
#define MU fluid_viscosity
#define RHO (reynolds_number * fluid_viscosity / (shear_rate * sq(radius)))
#define E_S (fluid_viscosity * radius * shear_rate / capillary_number)
#define ND_EB bending_ratio
#define C0 reference_curvature

#include "grid/octree.h"
// #include "grid/multigrid3D.h"

#include "ibm/IBMeshManager.h"
#include "ibm/navier-stokes/centered-ibm.h"
#include "ibm/IBOutput.h"

#include "lagrangian_caps_optim/capsule-ft.h"
#include "lagrangian_caps_optim/capsule-ibm-adapter.h"
#include "lagrangian_caps_optim/bending-ft.h"
#include "lagrangian_caps_optim/elasticity-ft.h"
#include "lagrangian_caps_optim/common-shapes-ft.h"
#include "lagrangian_caps_optim/capsule-viscosity-ratio.h"
#include "lagrangian_caps_optim/rheology-diagnostics-ft.h"

#include "navier-stokes/perfs.h"

#include "io/output-vtk-ibm.h"
#include "io/output-vtk.h"
#include "lagrangian_caps_optim/output-vtk-capsules.h"

#include "io/params/params-cli.h"

scalar myrho[];
scalar capsule_indicator[];
face vector mymu[];
face vector myalpha[];

static double suspension_radius_from_volume_fraction(void) {
  return cbrt(3. * volume_fraction * cube(domain_length) /
              (4. * pi * (double)n_caps));
}

static double suspension_rand_unit(void) {
  return (double)rand() / ((double)RAND_MAX + 1.);
}

static double suspension_rand_periodic_coord(void) {
  return (suspension_rand_unit() - 0.5) * domain_length;
}

static double suspension_rand_wall_coord(double margin) {
  return (suspension_rand_unit() - 0.5) * (domain_length - 2. * margin);
}

static double suspension_periodic_delta(double a, double b) {
  double d = a - b;
  if (d > 0.5 * domain_length)
    d -= domain_length;
  if (d < -0.5 * domain_length)
    d += domain_length;
  return d;
}

static bool suspension_centroid_is_clear(coord *centroids, int count,
                                         coord candidate, double min_dist) {
  for (int i = 0; i < count; i++) {
    double dx = suspension_periodic_delta(candidate.x, centroids[i].x);
    double dy = candidate.y - centroids[i].y;
    double dz = suspension_periodic_delta(candidate.z, centroids[i].z);
    if (sq(dx) + sq(dy) + sq(dz) < sq(min_dist))
      return false;
  }
  return true;
}

static void suspension_generate_centroids(coord *centroids) {
  srand((unsigned)random_seed);
  double min_delta = domain_length / (double)(1 << max_level);
  double min_initial_gap = 3. * min_delta;
  double margin = radius + min_initial_gap;
  double min_dist = 2. * radius + min_initial_gap;
  int max_attempts = 1000000;

  if (domain_length <= 2. * margin) {
    fprintf(stderr,
            "suspension: capsule radius/gap does not fit the wall-normal "
            "domain extent.\n");
    assert(false);
  }

  for (int c = 0; c < n_caps; c++) {
    bool placed = false;
    for (int attempt = 0; attempt < max_attempts; attempt++) {
      coord candidate = {suspension_rand_periodic_coord(),
                         suspension_rand_wall_coord(margin),
                         suspension_rand_periodic_coord()};
      if (suspension_centroid_is_clear(centroids, c, candidate, min_dist)) {
        centroids[c] = candidate;
        placed = true;
        if (pid() == 0)
          fprintf(stderr, "suspension: inserted capsule %d after %d attempts\n",
                  c, attempt + 1);
        break;
      }
    }

    if (!placed) {
      fprintf(stderr,
              "suspension: failed to insert capsule %d after %d attempts; "
              "reduce volume_fraction/count or increase max_level/domain.\n",
              c, max_attempts);
      assert(false);
    }
  }
}

static void suspension_capsule_forces(CapsuleMesh *mesh) {
  comp_elastic_stress(mesh);
  if (bending_ratio > 0.)
    comp_bending_force(mesh);
}

static void suspension_activate_capsules(void) {
  coord *centroids = (coord *)calloc((size_t)n_caps, sizeof(coord));
  assert(centroids);
  suspension_generate_centroids(centroids);

  CapsuleMeshPrototype *prototype =
      build_spherical_capsule_prototype((_initialize_circular_capsule){
          .cap_es = E_S, .cap_radius = RADIUS, .level = lag_level});

  for (int c = 0; c < n_caps; c++) {
    int owner = c % npe();
    activate_spherical_capsule(
        (_initialize_circular_capsule){.mesh = &CAPS(c),
                                       .prototype = prototype,
                                       .cap_es = E_S,
                                       .cap_radius = RADIUS,
                                       .level = lag_level,
                                       .cap_id = c,
                                       .cap_type = 0,
                                       .pid = owner,
                                       .shift = centroids[c]});
  }

  capsule_mesh_prototype_release(prototype);
  free(centroids);
}

static void suspension_prepare_ibm_models(void) {
  capsule_ibm_register_active_capsules();
  foreach_ibmesh()
    mesh->depth = max_level;
  foreach_ibnode_per_ibmesh() {
    node->depth = max_level;
  }
}

static void suspension_apply_fresh_initial_conditions(void) {
  foreach () {
    u.x[] = SHEAR_RATE * y;
    u.y[] = 0.;
    u.z[] = 0.;
    p[] = 0.;
    pf[] = 0.;
    capsule_indicator[] = 0.;
    capsule_viscosity_divergence[] = 0.;
    foreach_dimension() {
      ibmf.x[] = 0.;
      capsule_viscosity_grid_gradient.x[] = 0.;
    }
    Index_lagnode[] = -1;
    foreach_dimension() Index_lag_id.x[] = -1;
  }

  foreach_face() {
    mymu.x[] = fluid_viscosity;
    myalpha.x[] = 1. / RHO;
  }
  boundary({u, p, pf, capsule_indicator, capsule_viscosity_divergence,
            capsule_viscosity_grid_gradient, ibmf, Index_lagnode,
            Index_lag_id});

  ibmeshmanager_sync_velocity_coupled_model_outputs();

#if TREE
  adapt_wavelet_ibm(NULL,NULL, max_level, min_level, NULL, true);
#endif
}

int main(int argc, char *argv[]) {
  input_file_register_option_named("fluid", "length", domain_length,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("fluid", "shear_rate", shear_rate,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("fluid", "reynolds_number", reynolds_number,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("fluid", "fluid_viscosity", fluid_viscosity,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("fluid", "stokes_flow", stokes_flow,
                                   PARAM_VALUE_BOOL);
  input_file_register_option_named("fluid", "min_level", min_level,
                                   PARAM_VALUE_INT);
  input_file_register_option_named("fluid", "max_level", max_level,
                                   PARAM_VALUE_INT);
  input_file_register_option_named("fluid", "lag_level", lag_level,
                                   PARAM_VALUE_INT);

  input_file_register_option_named("capsule", "count", n_caps, PARAM_VALUE_INT);
  input_file_register_option_named("capsule", "volume_fraction",
                                   volume_fraction, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("capsule", "radius", radius,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("capsule", "inside_viscosity",
                                   capsule_inside_viscosity,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("capsule", "viscosity_enabled",
                                   capsule_viscosity_enabled, PARAM_VALUE_BOOL);
  input_file_register_option_named("capsule", "capillary_number",
                                   capillary_number, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("capsule", "bending_ratio", bending_ratio,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("capsule", "reference_curvature",
                                   reference_curvature, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("capsule", "random_seed", random_seed,
                                   PARAM_VALUE_INT);

  input_file_register_option_named("basilisk.solver", "dt_max", dt_max,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("basilisk.solver", "solver_tolerance",
                                   solver_tolerance, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("basilisk.solver", "velocity_tolerance",
                                   velocity_tolerance, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("basilisk.solver", "t_end", t_end,
                                   PARAM_VALUE_DOUBLE);

  input_file_register_option_named("output", "output_freq_vtk", output_freq_vtk,
                                   PARAM_VALUE_INT);
  input_file_register_option_named("output", "output_freq_dump",
                                   output_freq_dump, PARAM_VALUE_INT);
  input_file_register_option_named("output", "output_freq_stats",
                                   output_freq_stats, PARAM_VALUE_INT);
  input_file_register_option_named("output", "basepath", base_path,
                                   PARAM_VALUE_STRING);

  int input_file_parse_result = input_file_parse_cli(argc, argv);
  if (input_file_parse_result != 0)
    return input_file_parse_result > 0 ? 0 : input_file_parse_result;

  input_file_print_options();
  input_file_apply_options();

  if (n_caps <= 0) {
    fprintf(stderr, "suspension: capsule.count must be positive.\n");
    return 1;
  }
  if (radius <= 0.)
    radius = suspension_radius_from_volume_fraction();
  t_end = t_end > 0. ? t_end : 100. / shear_rate;

  L0 = domain_length;
  origin(-0.5 * domain_length, -0.5 * domain_length, -0.5 * domain_length);
  periodic(right);
  periodic(front);

#if TREE
  N = 1 << min_level;
#else 
  N = 1 << max_level;
#endif

  mu = mymu;
  rho = myrho;
  alpha = myalpha;
  TOLERANCE = solver_tolerance;
  stokes = stokes_flow;
  DT = dt_max;

  capsule_manager_set_count(&allCaps, n_caps);
  capsule_ibm_set_force_assembler(suspension_capsule_forces);
  capsule_viscosity_set_wall_boundary_conditions(capsule_indicator);

  checkpoint_configuration.checkpoint_on_sim_iter_iterations = output_freq_dump;
  checkpoint_configuration.checkpoint_on_sim_iter = true;
  checkpoint_configuration.checkpoint_on_wall_time = false;

  run();
}

u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
u.t[bottom] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.r[top] = dirichlet(shear_rate * y);
u.r[bottom] = dirichlet(shear_rate * y);
uf.n[top] = dirichlet(0);
uf.n[bottom] = dirichlet(0);

event init(i = 0) {
  suspension_activate_capsules();
  suspension_prepare_ibm_models();
  if (!restore_handler(base_path)) {
    suspension_apply_fresh_initial_conditions();
  } else {
    foreach () {
      p[] = 0.;
      pf[] = 0.;
    }
    boundary({p, pf});
  }
}

event properties(i++, last) {
  foreach ()
    myrho[] = RHO;

  if (capsule_viscosity_enabled) {
    capsule_viscosity_construct_indicator(capsule_indicator);
    capsule_viscosity_set_from_indicator(
        mymu, capsule_indicator, fluid_viscosity, capsule_inside_viscosity);
  } else {
    foreach_face() mymu.x[] = fluid_viscosity;
  }

  foreach_face() myalpha.x[] = 1. / RHO;
}

event vtk(i += output_freq_vtk) {
#if TREE
  output_hdf_htg({p, capsule_viscosity_divergence},
                 {u, ibmf, capsule_viscosity_grid_gradient}, base_path);
#else
  output_hdf_imagedata({p, capsule_viscosity_divergence},
                       {u, ibmf, capsule_viscosity_grid_gradient}, base_path);
#endif
  output_hdf_ibm(
      NULL, (IBvector[]){nforce, nvel, capsule_viscosity_area_normal, {{-1}}},
      base_path);
  output_hdf_capsules(base_path);
}

event progress_log(i++) {
  if (pid() == 0) {
    fprintf(ferr, "%d %g\n", i, t);
  }
}

event rheology_log(i += output_freq_stats) {
  if (i > 0)
    output_capsule_rheology(base_path, i, t, fluid_viscosity,
                            capsule_inside_viscosity, shear_rate,
                            capsule_viscosity_enabled);
}

#if TREE
event adapt(i++) {

  adapt_wavelet_ibm(
      {u},
      (double[]){velocity_tolerance, velocity_tolerance, velocity_tolerance},
      max_level, min_level,
      {u, p, pf, capsule_indicator, capsule_viscosity_divergence,
       capsule_viscosity_grid_gradient, ibmf, Index_lagnode, Index_lag_id});
}
#endif

event checkpoint_event(i++) {
  int status = checkpoint_handler(t, i, base_path);
  return status;
}

event stop(t = t_end) { return 1; }

event cleanup(t = end, last) {
  capsule_ibm_cleanup();
  ibmeshmanager_free();
  free_all_caps(&allCaps);
}
