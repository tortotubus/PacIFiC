#include "navier-stokes/centered.h"
#include "particle/ParticleGrid.h"
#if _MPI
#include "particle/ParticleGridMPI.h"
#endif

#ifndef PARTICLE_GRID_LEVEL
#define PARTICLE_GRID_LEVEL 0
#endif

#ifndef TWO_WAY
#define TWO_WAY 1
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
#if TWO_WAY
Pvector pforce;    // particle drag force on fluid
#endif

// Interpolated Eulerian quantities at particle positions
Pscalar peulrho;   // fluid density at particle
Pscalar peulmu;    // fluid viscosity at particle
Pvector peulvel;   // fluid velocity at particle

double ndt = 0.;
double ndt_prev = 0.;

// Physical gravity vector
coord nG = {.x = 0., .y = -9.81, .z = 0.};

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
  }

  particle_grid_update_cells ();
#if _MPI
  particle_grid_update_pid ();
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
