#define TWO_WAY 0
#define LJ_FORCE 1
#define BROWNIAN 1

static int particle_grid_level = 5;
static double brownian_pe = 0.34;
static double epsilon_lj = 1e-9;
static double sigma_lj = 8.908987e-3;
static int brownian_seed = 1;

#define PARTICLE_GRID_LEVEL particle_grid_level
#define BROWNIAN_PE brownian_pe
#define EPS_LJ epsilon_lj
#define SIGMA_LJ sigma_lj
#define RCUT_LJ (2.5 * SIGMA_LJ)
#define BROWNIAN_SEED brownian_seed

#include "grid/octree.h"

#include "particle/ParticleRandom.h"
#include "particle/navier-stokes/langevin.h"

#include "io/output-vtk-particles.h"
#include "io/output-vtk.h"
#include "io/params/params-cli.h"

#define DEFAULT_BASE_PATH "2particles"

static int particle_count = 10000;
static int output_log_every = 5;
static double output_dt = 0.002;
static double end_time = 4.;
static double fluid_viscosity = 1e-3;
static double fluid_density = 1000.;
static double particle_density = 2000.;
static double domain_length = 1.;
static double dt_max = 1e-5;
static int min_level = 6;
char *base_path = DEFAULT_BASE_PATH;

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

static double two_particles_radius(void) {
  return 0.5 * pow(2., 1. / 6.) * SIGMA_LJ;
}

static void two_particles_place_particles(void) {
  const int np = particle_count;
  const double radius = two_particles_radius();
  const double min_sep = 2. * radius;
  const int max_attempts = 10000;
  coord *placed = (coord *)calloc((size_t)np, sizeof(*placed));
  assert(placed);

  for (int n = 0; n < np; n++) {
    coord candidate = {0., 0., 0.};
    int attempts = 0;

    do {
      const uint64_t candidate_index =
          particle_random_fold((uint64_t)n, (uint64_t)attempts);
      candidate =
          particle_random_domain_position((uint64_t)BROWNIAN_SEED,
                                          PARTICLE_RANDOM_STREAM_PLACEMENT,
                                          candidate_index);
      attempts++;
    } while (overlaps_existing_particle(candidate, placed, min_sep, n) &&
             attempts < max_attempts);

    if (attempts >= max_attempts) {
      fprintf(stderr,
              "Could not place particle %d without overlap; reduce count/radius "
              "or use lattice init.\n",
              n);
      exit(1);
    }

    placed[n] = candidate;

    Particle *particle = particle_grid_add_particle(candidate);
    if (particle) {
      pval(prho) = particle_density;
      pval(pradius) = radius;
      foreach_dimension() pval(pvel.x) = 0.;
    }
  }

  free(placed);
}

static void two_particles_output(int iter, double time) {
#if TREE
  output_hdf_htg(NULL, NULL, base_path, iter, time);
#else
  output_hdf_imagedata(NULL, NULL, base_path, iter, time);
#endif
  output_hdf_particles((Pscalar[]){pradius, prho, {-1}},
                       (Pvector[]){ppos, pvel, {{-1}}}, base_path, iter, time);
}

int main(int argc, char *argv[]) {
  input_file_register_option_named("particle", "count", particle_count,
                                   PARAM_VALUE_INT);
  input_file_register_option_named("particle", "density", particle_density,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("particle", "grid_level",
                                   particle_grid_level, PARAM_VALUE_INT);
  input_file_register_option_named("particle", "brownian_pe", brownian_pe,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("particle", "brownian_seed",
                                   brownian_seed, PARAM_VALUE_INT);
  input_file_register_option_named("particle", "lj_epsilon", epsilon_lj,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("particle", "lj_sigma", sigma_lj,
                                   PARAM_VALUE_DOUBLE);

  input_file_register_option_named("fluid", "density", fluid_density,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("fluid", "viscosity", fluid_viscosity,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("fluid", "length", domain_length,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("fluid", "min_level", min_level,
                                   PARAM_VALUE_INT);

  input_file_register_option_named("basilisk.solver", "dt_max", dt_max,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("output", "output_dt", output_dt,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("output", "log_every", output_log_every,
                                   PARAM_VALUE_INT);
  input_file_register_option_named("output", "end_time", end_time,
                                   PARAM_VALUE_DOUBLE);
  input_file_register_option_named("output", "basepath", base_path,
                                   PARAM_VALUE_STRING);

  int input_file_parse_result = input_file_parse_cli(argc, argv);
  if (input_file_parse_result != 0)
    return input_file_parse_result > 0 ? 0 : input_file_parse_result;
  else {
    input_file_print_options();
    input_file_apply_options();
  }

  const face vector muc[] = {fluid_viscosity, fluid_viscosity, fluid_viscosity};

  periodic(left);
  periodic(bottom);
  periodic(front);

  mu = muc;
  stokes = true;
  nG.y = 0.;

  L0 = domain_length;
  origin(-L0 / 2., -L0 / 2., -L0 / 2.);
  N = 1 << min_level;
  DT = dt_max;

  run();
}

event init(i = 0) {
  const scalar rhof[] = fluid_density;
  rho = rhof;

  two_particles_place_particles();
  particle_grid_update_cells();
}

#if TREE
event adapt(i++) {
  adapt_wavelet((scalar *){u}, (double[]){0.001, 0.001, 0.001}, 7, 4);
  particle_grid_update_cells();
#if _MPI
  particle_grid_update_pid();
#endif
}
#endif

event movie(t += output_dt) { two_particles_output(i, t); }

event log(i += output_log_every) {
  if (pid() == 0)
    fprintf(stderr, "i=%d t=%g np=%d radius=%g\n", i, t, particle_count,
            two_particles_radius());
}

event end(t = end_time);
