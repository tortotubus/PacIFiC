
#include "io/output-vtk-pd.h"
#include "io/output-vtk.h"
#include "particle/ParticleRandom.h"
#include "particle/navier-stokes/stokes.h"

u.x[bottom] = dirichlet(0.);
u.z[bottom] = dirichlet(0.);

#define MUZ 1e-3

#ifndef ASH_RANDOM_SEED
#define ASH_RANDOM_SEED PARTICLE_RANDOM_DEFAULT_SEED
#endif

double damp = 2.0;

static void ash_add_particle_at(coord pos, uint64_t particle_index) {
  Particle *particle = particle_grid_add_particle(pos);
  if (!particle)
    return;

  const double radius_noise =
      1. - 2. * particle_random_uniform01((uint64_t)ASH_RANDOM_SEED,
                                          PARTICLE_RANDOM_STREAM_INIT,
                                          particle_index);
  pval(nrho) = 2000.;
  pval(nradius) = 0.05e-3 + 0.01e-3 * radius_noise;
  foreach_dimension() pval(nvel.x) = 0.;
}

static void ash_place_particles(void) {
  int pni = 40;

  for (int ii = 0; ii < pni; ii++)
    for (int pj = 0; pj < pni; pj++)
      for (int pk = 0; pk < pni; pk++) {
        coord pos = {0., 0., 0.};
        pos.x = (double)ii / (100. * pni);
        pos.y = (double)pj / (100. * pni) + L0 / 4.;
        pos.z = (double)pk / (100. * pni);
        uint64_t particle_index =
            ((uint64_t)ii * (uint64_t)pni + (uint64_t)pj) * (uint64_t)pni +
            (uint64_t)pk;
        ash_add_particle_at(pos, particle_index);
      }

  particle_grid_update_cells();
}

int main() {
  const face vector muc[] = {MUZ, MUZ, MUZ};

  periodic(left);
  mu = muc;

  nG.y = -9.81;

  L0 = 0.1;
  origin(-L0 / 2., -L0 / 2., -L0 / 2.);
  N = 64;
  DT = 1e-3;

  run();
}

event init(i = 0) {
  const scalar rhof[] = 1000.;
  rho = rhof;
  ash_place_particles();
}

event ash_particle_bottom_wall(i++) {
  foreach_particle() {
    if (pval(npos.y) < Y0) {
      pval(npos.y) = Y0 + (Y0 - pval(npos.y));

      foreach_dimension() pval(nvel.x) /= damp;

      pval(nvel.y) *= -1.;
    }
  }

  // Both functions need to be called any time position (npos) is modified
  particle_grid_update_cells();
#if _MPI
  particle_grid_update_pid();
#endif
}

event adapt(i++) {
  adapt_wavelet((scalar *){u}, (double[]){0.001, 0.001, 0.001}, 7, 4);

// Since adapt_wavelet chagnes the ownership of fluid cells, the ownership of
// particles must also be updated, even though their position hasn't changed
#if ASH_ENABLE_PARTICLES && _MPI
  particle_grid_update_pid();
#endif
}

event movie(t += 0.02) {
  output_hdf_htg((scalar[]){p, {-1}}, (vector[]){u, {{-1}}}, "ash");
  output_hdf_pd((Pscalar[]){nradius, nrho, {-1}},
                (Pvector[]){npos, nvel, {{-1}}}, "ash");
}

event log(i++) { printf("%d %f\n", i, t); }

event end(t = 2);
