/**
# Single optimized capsule with Langevin particles in shear flow

This example keeps the old particle/example/one_cap_old.c nanoparticle setup
but uses the optimized capsule/IBM implementation.
*/

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

#define REF_CURV 1
#define GLOBAL_REF_CURV 1
#define ADVECT_LAG_RK2 0 

#define TWO_WAY 0
#define LJ_FORCE 0
#define BROWNIAN 0
#define BROWNIAN_PE 0.34
#define BROWNIAN_SEED 1
#define STOKES_REFLECT_Y_WALLS 1
#define NANOPARTICLE_CAPSULE_INDICATOR 0

#define EPS_LJ 1e-9
#define SIGMA_LJ 8.908987e-3
#define RCUT_LJ (2.5 * SIGMA_LJ)

#define DEFAULT_BASE_PATH "one_cap_data"

static double domain_length = 1.;
static double capsule_radius = 0.1;
static double shear_rate = 1.;
static double reynolds_number = 10.;
static double fluid_viscosity = 1.;
static double capsule_inside_viscosity = 5.;
static double membrane_viscosity = 1.;
static double capillary_number = 0.2;
static double bending_ratio = .025;
static double reference_curvature = 0.;
static double area_dilatation_modulus = 50.;
static double dt_max = 1.e-3;
static double solver_tolerance = 1.e-6;
static double velocity_tolerance = 0.01;
static int min_level = 4;
static int max_level = 6;
static int lag_level = 4;
static int output_freq = 5;
static int log_freq = 100;
static int checkpoint_freq = 200;
static int particle_count = 10000;
static int particle_grid_level = 5;
static int particle_random_seed = 100;
static double particle_density = 2000.;
static double particle_radius = 0.5 * pow(2., 1. / 6.) * SIGMA_LJ;
static bool stokes_flow = false;
static bool capsule_viscosity_enabled = true;
char *base_path = DEFAULT_BASE_PATH;

#define RADIUS capsule_radius
#define SHEAR_RATE shear_rate
#define MU fluid_viscosity
#define RHO (reynolds_number * membrane_viscosity / (shear_rate * sq(RADIUS)))
#define E_S (membrane_viscosity * RADIUS * shear_rate / capillary_number)
#define ND_EB bending_ratio
#define E_B (ND_EB * E_S * sq(RADIUS))
#define C0 reference_curvature
#define AREA_DILATATION_MODULUS area_dilatation_modulus
#define TEND (30. / shear_rate)
#define PARTICLE_GRID_LEVEL particle_grid_level

#include "grid/octree.h"

#include "ibm/navier-stokes/centered-ibm.h"
#include "ibm/IBOutput.h"
#include "particle/navier-stokes/langevin.h"

#include "lagrangian_caps_optim/capsule-ft.h"
#include "lagrangian_caps_optim/capsule-ibm-adapter.h"
#include "lagrangian_caps_optim/bending-ft.h"
#include "lagrangian_caps_optim/skalak-ft.h"
#include "lagrangian_caps_optim/common-shapes-ft.h"
#include "lagrangian_caps_optim/capsule-viscosity-ratio.h"

#include "io/output-vtk-particles.h"
#include "io/output-vtk.h"
#include "lagrangian_caps_optim/output-vtk-capsules.h"
#include "io/params/params-cli.h"

scalar myrho[];
scalar capsule_indicator[];
face vector mymu[];
face vector myalpha[];

static FILE *particle_log_file = NULL;

static int overlaps_existing_particle(coord candidate, coord *placed,
                                      double min_sep, int placed_count) {
  for (int m = 0; m < placed_count; m++) {
    double r2 = 0.;
    domain_distance2(r2, candidate, placed[m]);

    if (r2 < sq(min_sep))
      return 1;
  }

  return 0;
}

static void one_cap_capsule_forces(CapsuleMesh *mesh) {
  comp_elastic_stress(mesh);
  comp_bending_force(mesh);
}

static void one_cap_activate_capsule(void) {
  CapsuleMeshPrototype *capsule_prototype =
      build_biconcave_capsule_prototype((_initialize_circular_capsule){
          .cap_es = E_S, .cap_radius = RADIUS, .level = lag_level});

  activate_biconcave_capsule(
      (_initialize_circular_capsule){.mesh = &CAPS(0),
                                     .prototype = capsule_prototype,
                                     .cap_es = E_S,
                                     .cap_radius = RADIUS,
                                     .level = lag_level,
                                     .cap_id = 0,
                                     .pid = 0});

  capsule_mesh_prototype_release(capsule_prototype);

  if (capsule_mesh_is_local_owner(&CAPS(0)))
    place_capsule_mesh_z_rotation(&CAPS(0), (coord){0., 0., 0.}, 0.5 * pi);
}

static void one_cap_prepare_ibm_models(void) {
  capsule_ibm_register_active_capsules();
  foreach_ibmesh() mesh->depth = max_level;
  foreach_ibnode_per_ibmesh() node->depth = max_level;
}

static void one_cap_initialize_particles(void) {
  particle_grid_delete_all_particles();

  const int np = particle_count;
  const double min_sep = 2. * particle_radius;
  const int max_attempts = 10000;
  coord *placed = (coord *)calloc((size_t)np, sizeof(*placed));
  assert(placed);

  for (int n = 0; n < np; n++) {
    coord candidate = {0., 0., 0.};
    int attempts = 0;

    do {
      const uint64_t candidate_index =
          particle_random_fold((uint64_t)n, (uint64_t)attempts);
      candidate = particle_random_domain_position(
          (uint64_t)particle_random_seed, PARTICLE_RANDOM_STREAM_PLACEMENT,
          candidate_index);
      attempts++;
    } while (overlaps_existing_particle(candidate, placed, min_sep, n) &&
             attempts < max_attempts);

    if (attempts >= max_attempts) {
      fprintf(stderr,
              "Could not place particle %d without overlap; reduce particle "
              "count/radius or use lattice initialization.\n",
              n);
      exit(1);
    }

    placed[n] = candidate;

    Particle *particle = particle_grid_add_particle(candidate);
    if (!particle)
      continue;

    pval(prho) = particle_density;
    pval(pradius) = particle_radius;
    foreach_dimension() pval(pvel.x) = 0.;
  }

  free(placed);

  particle_grid_update_cells();
}

static void one_cap_update_particle_grid_ownership(void) {
  particle_grid_update_cells();
#if _MPI
  particle_grid_update_pid();
#endif
}

static void one_cap_apply_fresh_initial_conditions(void) {
  foreach () {
    u.x[] = SHEAR_RATE * y;
    u.y[] = 0.;
    u.z[] = 0.;
    Index_lagnode[] = -1;
    foreach_dimension() Index_lag_id.x[] = -1;
  }

  ibmeshmanager_sync_velocity_coupled_model_outputs();
  one_cap_initialize_particles();

#if TREE
  adapt_wavelet_ibm(NULL, NULL, 0, 1, all, true);
  one_cap_update_particle_grid_ownership();
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

  input_file_register_option_named("capsule", "radius", capsule_radius,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("capsule", "membrane_viscosity",
                                   membrane_viscosity, PARAM_VALUE_DOUBLE);
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
  input_file_register_option_named("capsule", "area_dilatation_modulus",
                                   area_dilatation_modulus, PARAM_VALUE_DOUBLE);

  input_file_register_option_named("particles", "count", particle_count,
                                   PARAM_VALUE_INT);
  input_file_register_option_named("particles", "grid_level", particle_grid_level,
                                   PARAM_VALUE_INT);
  input_file_register_option_named("particles", "random_seed",
                                   particle_random_seed, PARAM_VALUE_INT);
  input_file_register_option_named("particles", "density", particle_density,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("particles", "radius", particle_radius,
                                   PARAM_VALUE_DOUBLE);

  input_file_register_option_named("basilisk.solver", "dt_max", dt_max,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("basilisk.solver", "solver_tolerance",
                                   solver_tolerance, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("basilisk.solver", "velocity_tolerance",
                                   velocity_tolerance, PARAM_VALUE_DOUBLE);

  input_file_register_option_named("fluid", "min_level", min_level,
                                   PARAM_VALUE_INT);
  input_file_register_option_named("fluid", "max_level", max_level,
                                   PARAM_VALUE_INT);
  input_file_register_option_named("fluid", "lag_level", lag_level,
                                   PARAM_VALUE_INT);

  input_file_register_option_named("output", "output_freq", output_freq,
                                   PARAM_VALUE_INT);
  input_file_register_option_named("output", "log_freq", log_freq,
                                   PARAM_VALUE_INT);
  input_file_register_option_named("output", "checkpoint_freq", checkpoint_freq,
                                   PARAM_VALUE_INT);
  input_file_register_option_named("output", "basepath", base_path,
                                   PARAM_VALUE_STRING);

  int input_file_parse_result = input_file_parse_cli(argc, argv);
  if (input_file_parse_result != 0)
    return input_file_parse_result > 0 ? 0 : input_file_parse_result;

  input_file_print_options();
  input_file_apply_options();

  L0 = domain_length;
  origin(-0.5 * L0, -0.5 * L0, -0.5 * L0);
  periodic(right);
  periodic(front);

  N = 1 << min_level;
  mu = mymu;
  rho = myrho;
  alpha = myalpha;
  TOLERANCE = solver_tolerance;
  stokes = stokes_flow;
  DT = dt_max;
  nG.y = 0.;

  capsule_manager_set_count(&allCaps, 1);
  capsule_ibm_set_force_assembler(one_cap_capsule_forces);
  capsule_viscosity_set_wall_boundary_conditions(capsule_indicator);

  checkpoint_configuration.checkpoint_on_sim_iter_iterations = checkpoint_freq;
  checkpoint_configuration.checkpoint_on_sim_iter = checkpoint_freq > 0;
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
  if (pid() == 0) {
    create_path(base_path);
    char fname[1024];
    snprintf(fname, sizeof(fname), "%s/particles.log", base_path);
    particle_log_file = fopen(fname, "a");
  }

  one_cap_activate_capsule();
  one_cap_prepare_ibm_models();

  if (!restore_handler(base_path))
    one_cap_apply_fresh_initial_conditions();
  else {
    foreach () {
      p[] = 0.;
      pf[] = 0.;
    }
    boundary({p, pf});
    one_cap_initialize_particles();
  }
}

event properties(i++, last) {
  foreach ()
    myrho[] = RHO;

  if (capsule_viscosity_enabled)
    capsule_viscosity_construct_indicator(capsule_indicator);

  foreach_face() myalpha.x[] = 1. / RHO;
}

event properties(i++, last) {
  if (capsule_viscosity_enabled)
    capsule_viscosity_set_from_indicator(
        mymu, capsule_indicator, fluid_viscosity, capsule_inside_viscosity);
  else
    foreach_face() mymu.x[] = fluid_viscosity;
}

event vtk(i += output_freq) {
#if TREE
  output_hdf_htg({p, capsule_indicator}, {u, ibmf}, base_path);
#else
  output_hdf_imagedata({p, capsule_indicator}, {u, ibmf}, base_path);
#endif
  output_hdf_particles((Pscalar[]){pradius, prho, {-1}},
                       (Pvector[]){ppos, pvel, {{-1}}}, base_path);
  output_hdf_capsules(base_path);
}

#if TREE
event adapt(i++) {
  adapt_wavelet_ibm(
      {u},
      (double[]){velocity_tolerance, velocity_tolerance, velocity_tolerance},
      minlevel = min_level, maxlevel = max_level);
  one_cap_update_particle_grid_ownership();
}
#endif

event console_log(i++) {
  if (pid() == 0)
    printf("%d %f\n", i, t);
}

event checkpoint_event(i++) {
  if (checkpoint_freq <= 0)
    return 0;
  return checkpoint_handler(t, i, base_path);
}

event end(t = TEND) {
  if (particle_log_file)
    fclose(particle_log_file);
  particle_grid_free();
  capsule_ibm_cleanup();
  ibmeshmanager_free();
  free_all_caps(&allCaps);
}
