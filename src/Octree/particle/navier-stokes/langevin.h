#include "particle/ParticleGrid.h"
#include "particle/ParticleRandom.h"
#if _MPI
#include "particle/ParticleGridMPI.h"
#endif

#ifndef TWO_WAY
#define TWO_WAY 1
#endif

#ifndef LJ_FORCE
#define LJ_FORCE 0
#endif

#ifndef EPS_LJ
#define EPS_LJ 0.
#endif

#ifndef SIGMA_LJ
#define SIGMA_LJ 1.
#endif

#ifndef RCUT_LJ
#define RCUT_LJ (2.5*SIGMA_LJ)
#endif

#ifndef RMIN_LJ
#define RMIN_LJ (0.8*SIGMA_LJ)
#endif

#ifndef PARTICLE_GRID_LEVEL
#define PARTICLE_GRID_LEVEL 0
#endif

#ifndef BROWNIAN
#define BROWNIAN 0
#endif

#ifndef BROWNIAN_PE
#define BROWNIAN_PE 1.
#endif

#ifndef BROWNIAN_SEED
#define BROWNIAN_SEED 1
#endif

#ifndef STOKES_REFLECT_Y_WALLS
#define STOKES_REFLECT_Y_WALLS 0
#endif

#ifndef NANOPARTICLE_CAPSULE_INDICATOR
#define NANOPARTICLE_CAPSULE_INDICATOR 0
#endif

#ifndef NANOPARTICLE_CAPSULE_INDICATOR_THRESHOLD
#define NANOPARTICLE_CAPSULE_INDICATOR_THRESHOLD 0.5
#endif

// extern face vector a;
// extern face vector mu;
// extern scalar rho;
// extern vector u;

// Particle quantities
Pscalar prho;      // particle density
Pscalar pradius;   // particle radius
Pscalar ptau;      // particle relaxation time
Pvector pgravity;  // buoyancy-corrected gravity acceleration
Pvector pforce;    // particle drag force on fluid

// Interpolated Eulerian quantities at particle positions
Pscalar peulrho;   // fluid density at particle
Pscalar peulmu;    // fluid viscosity at particle
Pvector peulvel;   // fluid velocity at particle

#if NANOPARTICLE_CAPSULE_INDICATOR
Pscalar pprevI;     // previous capsule indicator at particle
Pscalar pseenI;     // whether previous indicator has been initialized
Pscalar penteredI;  // whether this particle has ever entered
long ncap_entries = 0;
long ncap_unique_entries = 0;
long ncap_exits = 0;
long ncap_inside = 0;

double nanoparticle_capsule_indicator_value (coord pos);
#endif

double ndt = 0.;
double ndt_prev = 0.;

// Physical gravity vector
coord nG = {.x = 0., .y = -9.81, .z = 0.};

#if STOKES_REFLECT_Y_WALLS && dimension >= 2
static inline void stokes_reflect_y_walls (Particle* particle)
{
  if (Period.y)
    return;

  const double ymin = Y0;
  const double ymax = Y0 + L0;

  while (*_pval (ppos.y, particle) < ymin ||
         *_pval (ppos.y, particle) > ymax) {
    if (*_pval (ppos.y, particle) < ymin) {
      *_pval (ppos.y, particle) = 2.*ymin - *_pval (ppos.y, particle);
      *_pval (pvel.y, particle) *= -1.;
    }
    if (*_pval (ppos.y, particle) > ymax) {
      *_pval (ppos.y, particle) = 2.*ymax - *_pval (ppos.y, particle);
      *_pval (pvel.y, particle) *= -1.;
    }
  }
}
#endif

#if NANOPARTICLE_CAPSULE_INDICATOR
void nanoparticle_capsule_indicator_update (void)
{
  long entries = 0;
  long unique_entries = 0;
  long exits = 0;
  long inside = 0;

  foreach_particle () {
    coord pos = {0., 0., 0.};
    foreach_dimension()
      pos.x = pval (ppos.x);

    const double Ip = nanoparticle_capsule_indicator_value (pos);
    const int was_initialized = pval (pseenI) > 0.5;
    const int was_inside =
      was_initialized &&
      pval (pprevI) >= NANOPARTICLE_CAPSULE_INDICATOR_THRESHOLD;
    const int is_inside = Ip >= NANOPARTICLE_CAPSULE_INDICATOR_THRESHOLD;

    if (was_initialized) {
      if (!was_inside && is_inside) {
        entries++;
        if (pval (penteredI) < 0.5) {
          unique_entries++;
          pval (penteredI) = 1.;
        }
      }
      else if (was_inside && !is_inside)
        exits++;
    }

    if (is_inside)
      inside++;

    pval (pprevI) = Ip;
    pval (pseenI) = 1.;
  }

#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &entries, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &unique_entries, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &exits, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, &inside, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif

  ncap_entries += entries;
  ncap_unique_entries += unique_entries;
  ncap_exits += exits;
  ncap_inside = inside;
}
#endif

#if BROWNIAN
coord stokes_brownian_velocity_kick (Particle* particle, double dtau, int iter)
{
  const double rp = *_pval (pradius, particle);
  coord duB = {0., 0., 0.};

  if (dtau <= 0. || rp <= 0. || BROWNIAN_PE <= 0.)
    return duB;

  const double dp = 2.*rp;
  const double sigma_u = dp/sqrt (2.*BROWNIAN_PE*dtau);
  // Counter-based noise: independent of MPI rank and particle traversal order.
  const uint64_t index = particle_random_fold (particle->gid, (uint64_t) iter);
  coord G = particle_random_gaussian_vector ((uint64_t) BROWNIAN_SEED,
                                             PARTICLE_RANDOM_STREAM_BROWNIAN,
                                             index);
  foreach_dimension()
    duB.x = G.x*sigma_u;

  return duB;
}
#endif

#if LJ_FORCE
int stokes_lj_cell_radius (void)
{
  const double cell_width = L0/(1 << pg.level);
  int radius = cell_width > 0. ? (int) ceil (RCUT_LJ/cell_width) : 1;
  return radius > 1 ? radius : 1;
}

void stokes_lj_apply_pair_velocity_kick (Particle* particle,
                                         Particle* particle_other,
                                         double dtau)
{
  const double rp = *_pval (pradius, particle);
  const double rhop = *_pval (prho, particle);
  const double mp = rhop*(4./3.)*pi*cube (rp);

  const double rp_other = *_pval (pradius, particle_other);
  const double rhop_other = *_pval (prho, particle_other);
  const double mp_other = rhop_other*(4./3.)*pi*cube (rp_other);

  if (dtau <= 0. || mp <= 0. || mp_other <= 0.)
    return;

  double r2 = 0.;
  coord dr = {0., 0., 0.};
  coord pos = {0., 0., 0.};
  coord pos_other = {0., 0., 0.};
  foreach_dimension() {
    pos.x = *_pval (ppos.x, particle);
    pos_other.x = *_pval (ppos.x, particle_other);
  }
  domain_displacement (dr, pos, pos_other);
  foreach_dimension()
    r2 += sq (dr.x);

  if (r2 <= 1e-24) {
    dr.x = 1e-9;
    r2 = sq (dr.x);
  }

  if (r2 >= sq (RCUT_LJ) || r2 <= 1e-24)
    return;

  const double r2_lj = max (r2, sq (RMIN_LJ));
  const double inv_r2 = 1./r2_lj;
  const double sr2 = sq (SIGMA_LJ)*inv_r2;
  const double sr6 = cube (sr2);
  const double f_scalar = 24.*EPS_LJ*(2.*sq (sr6) - sr6)*inv_r2;

  foreach_dimension() {
    const double force = f_scalar*dr.x;
    *_pval (pvel.x, particle) += force*dtau/mp;
    *_pval (pvel.x, particle_other) -= force*dtau/mp_other;
  }
}
#endif

event defaults (i = 0)
{
  init_psolver ();

  new_pscalar (prho);
  new_pscalar (pradius);
  new_pscalar (ptau);
  new_pvector (pgravity);
#if TWO_WAY
  new_pvector (pforce);
#endif

  new_pscalar (peulrho);
  new_pscalar (peulmu);
  new_pvector (peulvel);

#if NANOPARTICLE_CAPSULE_INDICATOR
  new_pscalar (pprevI);
  new_pscalar (pseenI);
  new_pscalar (penteredI);
#endif

  particle_grid_init (PARTICLE_GRID_LEVEL);

#if TWO_WAY
  if (is_constant (a.x)) {
    a = new face vector;
    foreach_face()
      a.x[] = 0.;
    boundary ((scalar *) {a});
  }
#endif
}

event stokes_particles_1 (i += 2, last)
{
#if NANOPARTICLE_CAPSULE_INDICATOR
  nanoparticle_capsule_indicator_update ();
#endif

  // Interpolate fluid quantities to the particles
  foreach_particle () {
    if (!particle_is_local (particle))
      continue;

    double muc = 0.;
    double rhof = 0.;
#if dimension == 1
    foreach_dimension()
      pval (peulvel.x) = domain_interpolate_local (u.x, px);
    rhof = is_constant (rho) ? constant (rho) : domain_interpolate_local (rho, px);
    muc  = is_constant (mu.x) ? constant (mu.x) : domain_interpolate_local (mu.x, px);
#elif dimension == 2
    foreach_dimension()
      pval (peulvel.x) = domain_interpolate_local (u.x, px, py);
    rhof = is_constant (rho) ? constant (rho) : domain_interpolate_local (rho, px, py);
    muc  = is_constant (mu.x) ? constant (mu.x) : domain_interpolate_local (mu.x, px, py);
#else
    foreach_dimension()
      pval (peulvel.x) = domain_interpolate_local (u.x, px, py, pz);
    rhof = is_constant (rho) ? constant (rho) : domain_interpolate_local (rho, px, py, pz);
    muc  = is_constant (mu.x) ? constant (mu.x) : domain_interpolate_local (mu.x, px, py, pz);
#endif
    pval (peulrho) = rhof;
    pval (peulmu)  = muc;
  }

  // Match the inertial-particles/Stokes pattern:
  // first effective step uses dt, subsequent ones use 2*dt
  const double dt_now = dt > 0. ? dt : (ndt > 0. ? ndt : 1e-6);
  ndt_prev = (ndt == 0.) ? dt_now : ndt;
  ndt      = (ndt == 0.) ? dt_now : 2.*dt_now;

  // Compute tau, buoyancy-corrected gravity, and drag force
  foreach_particle () {
    if (!particle_is_local (particle))
      continue;

    const double rp = pval (pradius);
    const double rhop = pval (prho);
    const double muc = pval (peulmu);

    // Stokes relaxation time
    if (rp > 0. && muc > 0.)
      pval (ptau) = rhop*sq(2.*rp)/(18.*muc);
    else
      pval (ptau) = HUGE; // effectively uncoupled if invalid

    foreach_dimension() {
      pval (pgravity.x) =
        (rhop > 0.)
        ? nG.x*(rhop - pval(peulrho))/rhop
        : 0.;

#if TWO_WAY
      // equal-and-opposite drag force on fluid
      pval (pforce.x) =
        6.*pi*muc*rp*(pval(pvel.x) - pval(peulvel.x));
#endif
    }
  }

  // Implicit velocity relaxation
  foreach_particle () {
    if (!particle_is_local (particle))
      continue;

    const double it  = ndt > 0. ? 1./ndt : HUGE;
    const double ita = (pval(ptau) < HUGE) ? 1./pval(ptau) : 0.;

    foreach_dimension() {
      pval (pvel.x) =
        (it*pval(pvel.x) + ita*pval(peulvel.x) + pval(pgravity.x))/
        (it + ita);
    }
  }

#if LJ_FORCE
  int lj_cell_radius = stokes_lj_cell_radius ();
  foreach_particle_pair (lj_cell_radius) {
    stokes_lj_apply_pair_velocity_kick (particle, particle_other, ndt);
  }
#endif

#if BROWNIAN
  foreach_particle () {
    if (!particle_is_local (particle))
      continue;

    coord duB = stokes_brownian_velocity_kick (particle, ndt, i);
    foreach_dimension()
      pval (pvel.x) += duB.x;
  }
#endif
}

event stokes_particles_2 (i = 1; i += 2)
{
  const double dt_move = ndt + ndt_prev;
  if (dt_move <= 0.)
    return 0;

  // Update position
  foreach_particle () {
    if (!particle_is_local (particle))
      continue;

    foreach_dimension()
      pval (ppos.x) += pval (pvel.x)*dt_move;

    coord pos = {0., 0., 0.};
    foreach_dimension()
      pos.x = pval (ppos.x);
    domain_wrap_coord (pos);
    foreach_dimension()
      pval (ppos.x) = pos.x;
#if STOKES_REFLECT_Y_WALLS && dimension >= 2
    stokes_reflect_y_walls (particle);
#endif
  }

  particle_grid_update_cells ();
#if _MPI
  particle_grid_update_pid ();
#endif

#if NANOPARTICLE_CAPSULE_INDICATOR
  nanoparticle_capsule_indicator_update ();
#endif
}

#if TWO_WAY
event acceleration (i++)
{
  vector fp[];
  foreach()
    foreach_dimension()
      fp.x[] = 0.;

  // Deposit particle drag force to the cell containing each particle
  foreach_particle () {
    if (!particle_is_local (particle))
      continue;

    if (point.level >= 0) {
      foreach_dimension()
        fp.x[] += pval(pforce.x);
    }
  }

  boundary ((scalar *) {fp});

  // Convert centered force to face acceleration
  foreach_face() {
    double rhof = face_value (rho, 0);
    if (rhof > 0.)
      a.x[] += face_value (fp.x, 0)/(dv()*rhof);
  }

  boundary ((scalar *) {a});
}
#endif
