/**
# Single capsule in shear flow

This is the particle-free version of the particle/example/one_cap.c setup. It
keeps the capsule mechanics and fluid/capsule coupling path, but does not use
the particle grid or nanoparticle output/logging code.
*/

#include <stdbool.h>

#define PI 3.14159265358979323846

#define REF_CURV 1
#define GLOBAL_REF_CURV 1 
#define ADVECT_LAG_RK2 0
  

#define DEFAULT_BASE_PATH "one_cap_data"

static double domain_length = 1.;
static double radius = 0.1;
static double shear_rate = 1.;
static double reynolds_number = 10.;
static double fluid_viscosity = 1.;
static double capsule_inside_viscosity = 5.;
static double membrane_viscosity = 1.;
static double capillary_number = 0.2;
static double bending_ratio = .025;
static double reference_curvature = 0.;
static double area_dilatation_modulus = 50.;
char *base_path = DEFAULT_BASE_PATH;
static int n_caps = 1;
static double dt_max = 1.e-3;
static double solver_tolerance = 1.e-6;
static double velocity_tolerance = 0.01;
static int min_level = 4;
static int max_level = 7;
static int lag_level = 4;
static int output_freq = 1;
static int output_freq_pv_dump = 5;
static bool stokes_flow = false;
static bool capsule_viscosity_enabled = true;

#define RADIUS radius
#define SHEAR_RATE shear_rate
#define MU fluid_viscosity
#define RHO (reynolds_number * membrane_viscosity / (shear_rate * sq(radius)))
#define E_S (membrane_viscosity * radius * shear_rate / capillary_number)
#define ND_EB bending_ratio
#define E_B (ND_EB * E_S * sq(radius))
#define C0 reference_curvature
#define AREA_DILATATION_MODULUS area_dilatation_modulus
#define TEND (30. / shear_rate)

#include "grid/octree.h"

#include "ibm/IBMeshManager.h"
#include "ibm/navier-stokes/centered-ibm.h"
#include "ibm/IBOutput.h"

#include "lagrangian_caps_optim/capsule-ft.h"
#include "lagrangian_caps_optim/capsule-ibm-adapter.h"
#include "lagrangian_caps_optim/bending-ft.h"
#include "lagrangian_caps_optim/skalak-ft.h"
#include "lagrangian_caps_optim/common-shapes-ft.h"
#include "lagrangian_caps_optim/capsule-viscosity-ratio.h"

#include "lambda2.h"
#include "navier-stokes/perfs.h"
#include "view.h"

#include "io/output-vtk-ibm.h"
#include "io/output-vtk.h"
#include "lagrangian_caps_optim/output-vtk-capsules.h"
#include "io/params/params-cli.h"

scalar myrho[];
scalar capsule_indicator[];
face vector mymu[];
face vector myalpha[];
vector capsule_lagforce_spread[];

static void one_cap_capsule_forces(CapsuleMesh *mesh) {
  comp_elastic_stress(mesh);
  comp_bending_force(mesh);
}

static void one_cap_spread_raw_lagforce(vector forcing, CapsuleMesh *lag) {
  CapsuleIBMProxy *proxy = capsule_ibm_register_mesh(lag);
  if (!proxy)
    return;

  IBMesh *ibmesh = &ibmm.meshes[proxy->mesh_id];
  capsule_ibm_model_sync(proxy, ibmesh);

  for (size_t i = 0; i < ibmesh->nodes.size; i++) {
    IBNode *node = ibmesh->nodes.ptrs[i];
    int node_id = node->node_lid >= 0 ? node->node_lid : (int)i;
    if (node_id < 0 || node_id >= lag->nln)
      continue;

    peskin_cosine_kernel_spread_dimensionless(node) {
      foreach_dimension()
        forcing.x[] += -weight * ibval(nforce.x) * ibval(nweight) / dv();
    }
  }
}

static void one_cap_update_lagforce_spread(void) {
  foreach() {
    foreach_dimension()
      capsule_lagforce_spread.x[] = 0.;
  }

  for (int i = 0; i < allCaps.nbcaps; i++)
    if (CAPS(i).isactive)
      one_cap_spread_raw_lagforce(capsule_lagforce_spread, &CAPS(i));

  foreach_dimension()
    capsule_lagforce_spread.x.dirty = true;
  boundary((scalar *){capsule_lagforce_spread});
}

static void one_cap_activate_capsules(void) {
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
}

static void one_cap_prepare_ibm_models(void) {
  capsule_ibm_register_active_capsules();
  foreach_ibmesh()
    mesh->depth = max_level;
  foreach_ibnode_per_ibmesh() {
    node->depth = max_level;
  }
}

static void one_cap_apply_fresh_initial_conditions(void) {
  foreach () {
    u.x[] = SHEAR_RATE * y;
    u.y[] = 0.;
    u.z[] = 0.;
    Index_lagnode[] = -1;
    foreach_dimension() Index_lag_id.x[] = -1;
  }

  if (capsule_mesh_is_local_owner(&CAPS(0)))
    place_capsule_mesh_z_rotation(&CAPS(0), (coord){0., 0., 0.}, 0.5 * pi);

  ibmeshmanager_sync_velocity_coupled_model_outputs();

#if TREE
  adapt_wavelet_ibm(NULL, NULL, 0, 1, all, true);
#endif
}

/* Here we do some basic initialization. */

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

  input_file_register_option_named("capsule", "membrane_viscosity",
                                   membrane_viscosity, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("capsule", "inside_viscosity",
                                   capsule_inside_viscosity,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("capsule", "viscosity_enabled",
                                   capsule_viscosity_enabled, PARAM_VALUE_BOOL);
  input_file_register_option_named("capsule", "count", n_caps, PARAM_VALUE_INT);
  input_file_register_option_named("capsule", "capillary_number",
                                   capillary_number, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("capsule", "bending_ratio", bending_ratio,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("capsule", "reference_curvature",
                                   reference_curvature, PARAM_VALUE_DOUBLE);
  input_file_register_option_named("capsule", "area_dilatation_modulus",
                                   area_dilatation_modulus, PARAM_VALUE_DOUBLE);

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
  input_file_register_option_named("output", "output_freq_pv_dump",
                                   output_freq_pv_dump, PARAM_VALUE_INT);
  input_file_register_option_named("output", "basepath", base_path,
                                   PARAM_VALUE_STRING);

  int input_file_parse_result = input_file_parse_cli(argc, argv);
  if (input_file_parse_result != 0)
    return input_file_parse_result > 0 ? 0 : input_file_parse_result;
  else {
    input_file_print_options();
    input_file_apply_options();
  }

  L0 = domain_length;
  origin(-0.5 * domain_length, -0.5 * domain_length, -0.5 * domain_length);
  periodic(right);
  periodic(front);

  N = 1 << min_level;
  mu = mymu;
  rho = myrho;
  alpha = myalpha;
  TOLERANCE = solver_tolerance;
  stokes = stokes_flow;
  DT = dt_max;

  capsule_manager_set_count(&allCaps, n_caps);
  capsule_ibm_set_force_assembler(one_cap_capsule_forces);
  capsule_viscosity_set_wall_boundary_conditions(capsule_indicator);

  checkpoint_configuration.checkpoint_on_sim_iter_iterations =
      output_freq_pv_dump;
  checkpoint_configuration.checkpoint_on_sim_iter = true;
  checkpoint_configuration.checkpoint_on_wall_time = false;

  run();
}

/*
 Here we set our boundary conditions on the fluid domain.
 */

u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
u.t[bottom] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.r[top] = dirichlet(shear_rate * y);
u.r[bottom] = dirichlet(shear_rate * y);
uf.n[top] = dirichlet(0);
uf.n[bottom] = dirichlet(0);

event init(i = 0) {
  one_cap_activate_capsules();
  one_cap_prepare_ibm_models();

  if (!restore_handler(base_path)) { // only if not a restore
    one_cap_apply_fresh_initial_conditions();
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

  foreach_face() { myalpha.x[] = 1. / RHO; }
}

event vtk(i += output_freq) {
  one_cap_update_lagforce_spread();

#if TREE
  output_hdf_htg({p, capsule_indicator, capsule_viscosity_divergence},
                 {u, ibmf, capsule_lagforce_spread,
                  capsule_viscosity_grid_gradient}, base_path);
#else
  output_hdf_imagedata({p, capsule_indicator, capsule_viscosity_divergence},
                       {u, ibmf, capsule_lagforce_spread,
                        capsule_viscosity_grid_gradient}, base_path);
#endif
  output_hdf_ibm(NULL, (IBvector[]){nforce, nvel, capsule_viscosity_area_normal, {{-1}}}, base_path);
  output_hdf_capsules(base_path);
}

#if TREE
event adapt(i++) {
  adapt_wavelet_ibm(
      {u},
      (double[]){velocity_tolerance, velocity_tolerance, velocity_tolerance},
      minlevel = min_level, maxlevel = max_level);
}
#endif

event console_log(i++) {
  if (pid() == 0) {
    printf("%d %f\n", i, t);
  }
}

event checkpoint_event(i++) { return checkpoint_handler(t, i, base_path); }

event end(i = 1000) {
  capsule_ibm_cleanup();
  ibmeshmanager_free();
  free_all_caps(&allCaps);
}
