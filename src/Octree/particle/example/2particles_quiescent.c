#include "grid/octree.h"



#define TWO_WAY 0
#define LJ_FORCE 1
#define PARTICLE_GRID_LEVEL 5
#define BROWNIAN 1
#define BROWNIAN_PE 0.34

#define EPS_LJ 1e-9
#define SIGMA_LJ 8.908987e-3
#define RCUT_LJ (2.5*SIGMA_LJ)

#define BROWNIAN_SEED 1

#include "particle/navier-stokes/stokes.h"
#include "io/output-vtk-pd.h"
#include "io/output-vtk.h"

double muz = 1e-3;



static inline double rand01(void) {
  return (double) rand()/((double) RAND_MAX + 1.);
}

static int overlaps_existing_particle(coord candidate,
                                      double min_sep,
                                      int placed_count) {
  for (int m = 0; m < placed_count; m++) {
    Particle *particle = pg.pool.active.ptrs[m];
    double r2 = 0.;

    foreach_dimension() {
      double dr = candidate.x - pval(npos.x);
      if (Period.x)
        dr -= L0*nearbyint(dr/L0);
      r2 += sq(dr);
    }

    if (r2 < sq(min_sep))
      return 1;
  }

  return 0;
}

int main() {
  const face vector muc[] = {muz, muz, muz};

  periodic(left);
  periodic(bottom);
  periodic(front);
  
  mu = muc;

  nG.y = 0;

  L0 = 1.;
  origin(-L0 / 2., -L0 / 2., -L0 / 2.);
  N = 64;
  DT = 1e-5;

  run();
}

event init(i = 0) {
  const scalar rhof[] = 1000.;
  rho = rhof;

  int np = 10000;
  particle_grid_add_particles(np);

 double radius = 0.5*pow(2., 1./6.)*SIGMA_LJ; // about 0.0005

  const double min_sep = 2.*radius;
  const int max_attempts = 10000;

  for (int n = 0; n < np; n++) {
    Particle *particle = pg.pool.active.ptrs[n];

    pval(nrho) = 2000.;
    pval(nradius) = radius;


    coord candidate = {0., 0., 0.};
    int attempts = 0;

    do {
      candidate.x = X0 + rand01()*L0;
      candidate.y = Y0 + rand01()*L0;
      candidate.z = Z0 + rand01()*L0;
      attempts++;
    } while (overlaps_existing_particle(candidate, min_sep, n) &&
             attempts < max_attempts);

    if (attempts >= max_attempts) {
      fprintf(stderr,
              "Could not place particle %d without overlap; reduce np/radius or use lattice init.\n",
              n);
      exit(1);
    }

    pval(npos.x) = candidate.x;
    pval(npos.y) = candidate.y;
    pval(npos.z) = candidate.z;

    foreach_dimension()
      pval(nvel.x) = 0.;
  }

  particle_grid_update_cells();
}


#if TREE
event adapt(i++) {
  adapt_wavelet((scalar *){u}, (double[]){0.001, 0.001, 0.001}, 7, 4);
}
#endif

event movie(t+=0.002) {
#if TREE
  output_hdf_htg(NULL, NULL, "2particles");
#else 
  output_hdf_imagedata(NULL,NULL,"2particles");
#endif
output_hdf_pd((Pscalar[]){nradius, {-1}},(Pvector[]){npos, nvel, {{-1}}},"2particles");
  //output_hdf_pd(NULL, NULL, "2particles");
}

event log(i+=5) {
  int particle_index = 0;
  foreach_particle() {
    printf("%d %.12g %d %.12g %.12g %.12g %.12g %.12g %.12g\n",
           i, t, particle_index++,
           pval(npos.x), pval(npos.y), pval(npos.z),
           pval(nvel.x), pval(nvel.y), pval(nvel.z));
  }
}
event end(t = 4);
 
