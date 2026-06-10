/**
# Rheology of a suspension of deformable capsules in a shear flow

## Definition of the relevant parameters

We define the following physical quantities:

* the size of the box $L_0$,
* the radius $a$,
* the shear rate $\dot{\gamma}$,
* the Reynolds number $Re$,
* the viscosity $\mu$,
* the density $\rho = \frac{Re \, \mu}{\dot{\gamma} a^2}$,
* the Capillary number defined $-$ in the case of capsules $-$ as the ratio of
viscous forces over elastic forces: $Ca = \frac{\mu a \dot{\gamma}}{E_S}$
* the elastic modulus $E_S$
* the bending modulus $E_b$
*/

#define L0 1.
#define PI 3.14159265358979323846

#define NCAPS 1 
/*c1 volume fraction 1.25e-6, c2 volume fraction 1.5625e-7*/
#define RADIUS 0.1

#define PSI_DEG 0
#define XI_DIST 1
#define THETAZX ((180. - PSI_DEG)*PI/180.) // 0-pi
#define DIST_cap01 (RADIUS + RADIUS*XI_DIST/8.)

#define SHEAR_RATE 1.
#ifndef RE
  #define RE 10
#endif
#define MU 1.
#define MUP 1.
#define MUC 5.

//Viscosity ratio

#define RHO (RE*MUP/(SHEAR_RATE*sq(RADIUS)))
#ifndef CA
  #define CA 0.2
#endif
#define E_S (MUP*RADIUS*SHEAR_RATE/CA)
#define ND_EB .025
#define E_B (ND_EB*E_S*sq(RADIUS))
#define REF_CURV 1
#define GLOBAL_REF_CURV 1
// #define C0 (-2.09/RADIUS)
#define C0 (0.0)
#define AREA_DILATATION_MODULUS 50.

#define STR(s) #s
#define STRING(s) STR(s)

//paralle test
#define ADVECT_LAG_RK2 0
#define LUBR_FORCE 0
#define LUBR_VEL 1 

#ifndef RESTART_CASE
  #define RESTART_CASE 0
#endif
/**
We also define some solver-ralated quantities:

* the non-dimensional duration time of the simulation TEND
* the minimum and maximum refinement levels of the Eulerian grid
* the number of refinement levels LAGLEVEL of the Lagrangian mesh of the capsule
* the maximum time-step DT_MAX, which appears to be dependant on the Capillary number
* the tolerance of the Poisson solver (maximum admissible residual)
* the tolerance of the wavelet adaptivity algorithm for the velocity
* an frequency of post-processing output (here, every 10 iterations)
* the $stokes$ boolean, in order to ignore the convective term in the
Navier-Stokes solver
* the Jacobi preconditionner, which we switch on for this case.
*/

#ifndef TEND
  #define TEND (30./SHEAR_RATE)
#endif
#ifndef MINLEVEL
  #define MINLEVEL 4
#endif
#ifndef MAXLEVEL
  #define MAXLEVEL 6
#endif
#ifndef LAG_LEVEL
  #define LAG_LEVEL 4 
#endif
#ifndef DT_MAX
  #define DT_MAX 1.e-3
  #endif
#ifndef MY_TOLERANCE
  #define MY_TOLERANCE 1.e-6
#endif
#ifndef U_TOL
  #define U_TOL 0.01
#endif
#ifndef OUTPUT_FREQ
  #define OUTPUT_FREQ 100
#endif
#ifndef OUTPUT_FREQ_PV_DUMP
  #define OUTPUT_FREQ_PV_DUMP 200
#endif
#ifndef STOKES
  #define STOKES false
#endif
#define JACOBI 1
#define PARAVIEW_CAPSULE 1
#define PARAVIEW_FLOW_FIELD 0
#define OUTPUT_CAPS_NODE_TRI_INFO 0

//Have to fix the error for HDF5 writer 
#ifndef OUTPUT_NP_FLOW_FIELD
  #define OUTPUT_NP_FLOW_FIELD 0
#endif
#ifndef OUTPUT_NP_PARTICLES
  #define OUTPUT_NP_PARTICLES 0
#endif


#ifndef LOG_NP_PARTICLES
  #define LOG_NP_PARTICLES 1
#endif

#define STENCIL_TYPE 3


//BS I have to improve after testing
#define TWO_WAY 0
#define LJ_FORCE 0
#define BROWNIAN 0
#define BROWNIAN_PE 0.34
#define STOKES_REFLECT_Y_WALLS 1

#define EPS_LJ 1e-9
#define SIGMA_LJ 8.908987e-3
#define RCUT_LJ (2.5*SIGMA_LJ)

#define BROWNIAN_SEED 1


#define NANOPARTICLE_CAPSULE_INDICATOR 0
/////////////////////////////////////////



/**
## Simulation setup

We import the octree grid, the centered Navier-Stokes solver, the Lagrangian
mesh, the neo-Hookean elasticity, a header file containing functions to
mesh a sphere, and the Basilisk viewing functions supplemented by a custom
function $draw\_lag$ useful to visualize the front-tracking interface.
*/

#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "particle/navier-stokes/stokes.h"
#include "io/output-vtk.h"
#include "io/output-vtk-ibm.h"

//#include "lagrangian_caps/my_multigrid3D.h"
#include "lagrangian_caps/capsule-ft.h"
//#include "lagrangian_caps/neo-hookean-ft.h"
#include "lagrangian_caps/skalak-ft.h"
#include "lagrangian_caps/bending-ft.h"
#include "lagrangian_caps/viscosity-ft.h"
#include "lagrangian_caps/common-shapes-ft.h"
#include "lagrangian_caps/view-ft.h"
#if PARAVIEW_FLOW_FIELD
  #include "lagrangian_caps/vtu_output.h"
#endif
#include "navier-stokes/perfs.h"
# include "lambda2.h"
# include "view.h"

FILE* fvelocity = NULL; //velpro
FILE* foutput = NULL;
FILE* fperf = NULL;
FILE* foutput_np = NULL;
FILE* fcrossing_caps = NULL;

int init_cycle_number = -1;
int nb_dump = 0;
int nb_pic = 0;

const scalar myrho[] = RHO;
const face vector mymu[] = {MU, MU, MU};
const face vector myalpha[] = {1./RHO, 1./RHO, 1./RHO};

static inline void sync_caps_grid_depth(void) {
#if TREE
  grid->maxdepth = depth();
#endif
}

#if NANOPARTICLE_CAPSULE_INDICATOR //to track if NP crossed cap membrane via the indicator function I(x)
double nanoparticle_capsule_indicator_value (coord pos)
{
#if dimension == 1
  return interpolate (I, pos.x);
#elif dimension == 2
  return interpolate (I, pos.x, pos.y);
#else
  return interpolate (I, pos.x, pos.y, pos.z);
#endif
}
#endif

int main(int argc, char* argv[]) {
  origin(-0.5*L0, -0.5*L0, -0.5*L0);
  /** We set periodic boundary conditions on the non-horizontal walls.*/
  periodic(right);
  periodic(front);

  N = 1 << MINLEVEL;


  nG.y = 0;

  mu = mymu;
  rho = myrho;
  alpha = myalpha;
  TOLERANCE = MY_TOLERANCE;
  /** We don't need to compute the convective term in this case, so we set the
  boolean $stokes$ to false. However it is still important to choose $Re \ll 1$
  since we are solving the unsteady Stokes equation. */
  stokes = STOKES;
  DT = DT_MAX;
  srand (100);
  run();
}


/** We impose shear-flow boundary conditions... */
u.n[bottom] = dirichlet(0);
u.n[top] = dirichlet(0);
u.t[bottom] = dirichlet(0.);
u.t[top] = dirichlet(0.);
u.r[top] = dirichlet(SHEAR_RATE*y);
u.r[bottom] =  dirichlet(SHEAR_RATE*y);
uf.n[top] = dirichlet(0);
uf.n[bottom] = dirichlet(0);

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

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

event init (i = 0) {
  sync_caps_grid_depth();


  fvelocity = fopen("velocity.txt", "a"); //velpro
  foutput  = fopen("output.txt","a");
  fperf = fopen("perf.csv", "a");
  if (pid() == 0) {
    foutput_np = fopen("output_np.txt", "a");
    fcrossing_caps = fopen("crossing_caps.txt", "a");
  }
  fprintf(fperf, "total_nb_cells, nb_iter_viscous, resb_viscous, resa_viscous, nb_relax_viscous, nb_iter_pressure, resb_pressure, resa_pressure, nb_relax_pressure\n");
 
  // Check if we are in a restart mode
  init_cycle_number = RESTART_CASE ? reinitialize_restart() : -1;
  if(init_cycle_number >= 0)
  {
      nb_pic = (int)(init_cycle_number*OUTPUT_FREQ_PV_DUMP)/OUTPUT_FREQ + 1;
      nb_dump = init_cycle_number + 1;
      // read_t_restart(dump_dir, &t, &dt);
    
    /**1_ We activate all the capsules.*/
    foreach() 
    {
      Index_lagnode[] = -1;
      foreach_dimension() Index_lag_id.x[] = -1;
    }

  }
  else
  {
    /** If not, we initialize the flow field to that of an undisturbed,
    fully-developed shear.*/
    foreach() 
    {
      u.x[] = SHEAR_RATE*y;
      u.y[] = 0.;
      u.z[] = 0.;
      Index_lagnode[] = -1;
      foreach_dimension() Index_lag_id.x[] = -1;
    }


    /**1_ We activate all the capsules from large to small.*/
    // for(int k = 0; k < NCAPS; k++)
    //   activate_spherical_capsule(&CAPS(k), cap_es = E_S_L, cap_radius = RADIUS, level = LAG_LEVEL_L, cap_id = k);
    // for(int k = NCAPS; k < NCAPS + NCAPS_S; k++)
    //   activate_spherical_capsule(&CAPS(k), cap_es = E_S_S, cap_radius = RADIUS_S, level = LAG_LEVEL_S, cap_id = k);

    /**2_ Generate random array of capsules, NCAPS should be around or less than 50 **/
     //generate_random_capsules();
    /**2_ Or generate capsules from file**/
    //generate_bidisperse_capsules_from_input();

//Generate one spherical capsule at center 
    activate_biconcave_capsule((_initialize_circular_capsule){.mesh=&CAPS(0), .cap_es = E_S, .cap_radius = RADIUS, .level = LAG_LEVEL, .cap_id = 0});
//    activate_spherical_capsule(&CAPS(1), cap_es = E_S_S, cap_radius = RADIUS_S, level = LAG_LEVEL_S, cap_id = 1);

      CAPS(0).centroid.x = 0.0;
      CAPS(0).centroid.y = 0.0;
      CAPS(0).centroid.z = 0.0;
      for(int i=0; i<CAPS(0).nln; i++)
        foreach_dimension() CAPS(0).nodes[i].pos.x += CAPS(0).centroid.x;

      for (int i = 0; i < CAPS(0).nln; i++) {
          double x = CAPS(0).nodes[i].pos.x - CAPS(0).centroid.x;
          double y = CAPS(0).nodes[i].pos.y - CAPS(0).centroid.y;
          double z = CAPS(0).nodes[i].pos.z - CAPS(0).centroid.z;
      
          // Apply 90° rotation about z-axis
          double x_new = -y;
          double y_new = x;
          double z_new = z;
          // double x_new = x;
          // double y_new = -z;
          // double z_new = y;
      
          CAPS(0).nodes[i].pos.x = x_new + CAPS(0).centroid.x;
          CAPS(0).nodes[i].pos.y = y_new + CAPS(0).centroid.y;
          CAPS(0).nodes[i].pos.z = z_new + CAPS(0).centroid.z;
      }
      
/*
      CAPS(1).centroid.x = -1.*cos(THETAZX)*DIST_cap01;
      CAPS(1).centroid.y = 0.0;
      CAPS(1).centroid.z = sin(THETAZX)*DIST_cap01;
      for(int i=0; i<CAPS(1).nln; i++)
        foreach_dimension() CAPS(1).nodes[i].pos.x += CAPS(1).centroid.x;
*/
// generate the stencils around the capsules
  for(int k=0; k<NCAPS; k++) correct_lag_pos(&CAPS(k));
    sync_caps_grid_depth();
    generate_lag_stencils(true);  


    #ifdef CAPS_VISCOSITY  
    
    foreach()
    {
        coord loc;
        loc.x = x;
        loc.y = y;
        loc.z = z;

        for(int k = 0; k < NCAPS; k++)
        {
          if(I[] < 1.e-6)
          	I[] = GENERAL_SQNORM(loc, CAPS(k).centroid) - sq(CAPS(k).cap_radius) > 0 ? 0 : 1;
        } //ggd 
    }
    foreach() if (I[] > 1.e-6) prevI[] = I[];
 
    #endif

    const scalar rhof[] = 1000.;
  rho = rhof;

  int np = 100;
  particle_grid_add_particles(np);

 double radius = 0.5*pow(2., 1./6.)*SIGMA_LJ; // about 0.0005

  const double min_sep = 2.*radius;
  const int max_attempts = 1000;

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

    //We initialize the viscosity field
    //fraction(mu,sq(x) + sq(y) + sq(z) - sq(RADIUS));

  }

}

#if !(MULT_GRID == 1)   
  event adapt (i++) {
    sync_caps_grid_depth();
    adapt_wavelet({u}, (double []){U_TOL, U_TOL, U_TOL},
      minlevel = MINLEVEL, maxlevel = MAXLEVEL);
    sync_caps_grid_depth();
    generate_lag_stencils(true);  
  }
#endif


event movie(t+=0.002) {
#if OUTPUT_NP_FLOW_FIELD
  #if TREE
    output_hdf_htg(NULL, NULL, "one_cap");
  #else 
    output_hdf_imagedata(NULL,NULL,"one_cap");
  #endif
#endif
#if OUTPUT_NP_PARTICLES
  output_hdf_pd((Pscalar[]){nradius, {-1}},(Pvector[]){npos, nvel, {{-1}}},"one_cap");
#endif
  //output_hdf_pd(NULL, NULL, "one_cap");
}

#if LOG_NP_PARTICLES
event log(i+=100) {
  if (pid() == 0 && foutput_np) {
    for (int particle_index = 0; particle_index < pg.pool.active.size; particle_index++) {
      Particle *particle = pg.pool.active.ptrs[particle_index];
      fprintf(foutput_np,
              "%d %.12g %d %.12g %.12g %.12g %.12g %.12g %.12g\n",
              i, t, particle_index,
              *_pval(npos.x, particle), *_pval(npos.y, particle),
              *_pval(npos.z, particle), *_pval(nvel.x, particle),
              *_pval(nvel.y, particle), *_pval(nvel.z, particle));
    }
    fflush(foutput_np);
  }

  
   #if NANOPARTICLE_CAPSULE_INDICATOR   ///to track NPs crossing cap membrane

  if (pid() == 0 && fcrossing_caps) {
    fprintf(fcrossing_caps, "%d %.12g %ld %ld %ld %ld\n",
            i, t, ncap_unique_entries, ncap_entries,
            ncap_exits, ncap_inside);
    fflush(fcrossing_caps);
  }
  #endif

}
#endif
event end(t = 4);
 


/*We output the field data for analysis*/
event output_and_dump_caps (i+=OUTPUT_FREQ_PV_DUMP) {
  #if OUTPUT_CAPS_NODE_TRI_INFO
    output_caps_node_tri(); //dump node and triangle infos
  #endif
  #if PARAVIEW_CAPSULE
    pv_output_ascii(nb_dump); // dump capsules in vtk
  #endif
    #if PARAVIEW_FLOW_FIELD
      scalar * list = {p};
      vector * vlist = {u};
      save_data(list, vlist, t, nb_dump);
    #endif

    char fname[62];
    sprintf(fname, "%s/flow_%d.dump", dump_dir, nb_dump);
    //scalar * dump_list =  (scalar *){u, p, I};
    //dump(fname, dump_list);  
    dump(fname);  
    sprintf(fname, "%s/caps_%d.dump", dump_dir, nb_dump);
    dump_capsules(fname, NULL); // dump lagrangian capsules

    /*We dump and update the cycle number for the next restart*/
    dump_cycle_number(nb_dump); 
    save_t_dt_restart(dump_dir, t, dt);
    nb_dump++;
}

event logfile (i += 25) {

  double top_visc_stress = 0;
  int top_nb_cells = 0;
  foreach_boundary(top, reduction(+:top_visc_stress) reduction(+:top_nb_cells)) {
    top_nb_cells++;
    top_visc_stress += (u.x[0, 1] - u.x[])*Delta +
      (u.y[1] - u.y[-1])*.5*Delta;
  }
  top_visc_stress *= MU/(sq(L0));

  double bottom_visc_stress = 0;
  int bottom_nb_cells = 0;
  foreach_boundary(bottom, reduction(+:bottom_visc_stress) reduction(+:bottom_nb_cells)) {
    bottom_nb_cells++;
    bottom_visc_stress += (u.x[0, 1] - u.x[])*Delta +
      (u.y[1] - u.y[-1])*.5*Delta;
  }
  bottom_visc_stress *= MU/(sq(L0));

  double fluid_visc_stress = (top_visc_stress + bottom_visc_stress) / 2.;

////////////////////////////////////// Particle Stresslet

  double* send_stress_pack = (double*)calloc(NCAPS*4, sizeof(double));
  double* recv_stress_pack = (double*)calloc(NCAPS*4, sizeof(double));

  double pN1 = 0.;
  double pN2 = 0.;
  double pmu = 0.;
  double ppres = 0.;

  for(int k = 0; k < NCAPS; k++)
  {
    double sigmaxy = 0.;
    double sigmaxx = 0.;
    double sigmayy = 0.;
    double sigmazz = 0.;

    Point point = locate(CAPS(k).centroid.x, CAPS(k).centroid.y, CAPS(k).centroid.z);

    if(point.level > -1)
    {

	for(int i=0; i<CAPS(k).nln; i++) 
	{
	/** The post-processing is only carried out if we are in the shear plane */ 
	double rx, ry, rz;


        rx = CAPS(k).centroid.x + GENERAL_1DIST(CAPS(k).nodes[i].pos.x, CAPS(k).centroid.x);
        ry = CAPS(k).centroid.y + GENERAL_1DIST(CAPS(k).nodes[i].pos.y, CAPS(k).centroid.y);
        rz = CAPS(k).centroid.z + GENERAL_1DIST(CAPS(k).nodes[i].pos.z, CAPS(k).centroid.z);

	#ifndef CAPS_VISCOSITY
	sigmaxx += - CAPS(k).nodes[i].lagForce.x * rx;
	sigmayy += - CAPS(k).nodes[i].lagForce.y * ry;
	sigmazz += - CAPS(k).nodes[i].lagForce.z * rz;
	sigmaxy += - (CAPS(k).nodes[i].lagForce.x * ry + CAPS(k).nodes[i].lagForce.y * rx) / 2.;
	#else 
	double visc_ratio = 1./MUP * MUC;
          /** We now have to compute the area associated with each node */
        double nodal_area = compute_node_area(&(CAPS(k)), i);
        
	sigmaxx += - CAPS(k).nodes[i].lagForce.x * rx + 2.*MU*(visc_ratio - 1)*(CAPS(k).nodes[i].lagVel.x*CAPS(k).nodes[i].normal.x)*nodal_area;
	sigmayy += - CAPS(k).nodes[i].lagForce.y * ry + 2.*MU*(visc_ratio - 1)*(CAPS(k).nodes[i].lagVel.y*CAPS(k).nodes[i].normal.y)*nodal_area;
	sigmazz += - CAPS(k).nodes[i].lagForce.z * rz + 2.*MU*(visc_ratio - 1)*(CAPS(k).nodes[i].lagVel.z*CAPS(k).nodes[i].normal.z)*nodal_area;
	sigmaxy += - (CAPS(k).nodes[i].lagForce.x * ry + CAPS(k).nodes[i].lagForce.y * rx) / 2. 
	           + MU*(visc_ratio - 1)*(CAPS(k).nodes[i].lagVel.x*CAPS(k).nodes[i].normal.y + CAPS(k).nodes[i].lagVel.y*CAPS(k).nodes[i].normal.x)*nodal_area;
	#endif	

	}

    pN1 = (sigmaxx - sigmayy);
    pN2 = (sigmayy - sigmazz);
    pmu = sigmaxy;
    ppres = -(sigmaxx + sigmayy + sigmazz)/3.;

    send_stress_pack[k*4] = pN1;
    send_stress_pack[k*4 + 1] = pN2;
    send_stress_pack[k*4 + 2] = pmu;
    send_stress_pack[k*4 + 3] = ppres;
    }

  pN1 = 0.;
  pN2 = 0.;
  pmu = 0.;
  ppres = 0.;
 }

MPI_Reduce(send_stress_pack, recv_stress_pack, 4*NCAPS, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

if (pid() == 0) 
{
    pN1 = 0.;
    pN2 = 0.;
    pmu = 0.;
    ppres = 0.;

   for(int k = 0; k < NCAPS; k++)
   {
	pN1 += recv_stress_pack[k*4];
	pN2 += recv_stress_pack[k*4 + 1];
	pmu += recv_stress_pack[k*4 + 2];
	ppres += recv_stress_pack[k*4 + 3];
   }

}
free(send_stress_pack);
free(recv_stress_pack);

//////////////////////////////////////////////////////////////////////////////////////


  if (pid() == 0) 
  {
    double avg_ncaps_area = 0;
    double avg_ncaps_volume = 0;
    for(int k=0; k<NCAPS; k++) {
      avg_ncaps_volume += CAPS(k).volume/CAPS(k).initial_volume;
      for(int i=0; i<CAPS(k).nlt; i++) avg_ncaps_area += CAPS(k).triangles[i].area/(4*pi*sq(CAPS(k).cap_radius));
    }
    avg_ncaps_area /= NCAPS;
    avg_ncaps_volume /= NCAPS;

  /*Compute average Taylor deformation and angular velocity*/
    double avg_TDmaxmin = 0;
    double avg_TDang = 0;
    double TDmaxmin = 0;
    double TDang = 0;
    double avg_taylor_deform = 0;
    double avg_inclin_angle = 0;
    double taylor_deform = 0;
    double inclin_angle = 0;
    double avg_ang_vel = 0;
    coord avg_rs = {0., 0., 0.};
    coord rs = {0., 0., 0.};
    for(int k=0; k<NCAPS; k++) {
      compute_taylor_factor(&CAPS(k), &taylor_deform, &inclin_angle, &rs, &TDmaxmin, &TDang);
      foreach_dimension() avg_rs.x += rs.x/CAPS(k).cap_radius; 
      avg_TDmaxmin += TDmaxmin;
      avg_TDang += TDang;
      avg_taylor_deform += taylor_deform;
      avg_inclin_angle += inclin_angle;
      avg_ang_vel += CAPS(k).ang_vel.z;
    }
    
    foreach_dimension() avg_rs.x /= NCAPS; 
    avg_TDmaxmin /= NCAPS;
    avg_TDang /= NCAPS;
    avg_inclin_angle /= NCAPS;
    avg_taylor_deform /= NCAPS;
    avg_ang_vel /= NCAPS;

    fprintf(foutput, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", i, t, 
      fluid_visc_stress, ppres, pmu, pN1, pN2, 
      avg_taylor_deform, avg_inclin_angle, avg_TDmaxmin, avg_TDang,
      avg_rs.x, avg_rs.y, avg_rs.z, avg_ang_vel, avg_ncaps_area, avg_ncaps_volume);
    fflush(foutput);
}

  // for(int k=0; k<NCAPS; k++) {
  //    output_physics(k, i);
  // }

  if (pid() == 0) 
  {
    fprintf(fperf, "%ld %d %g %g %d %d %g %g %d\n", grid->tn, mgu.i, mgu.resb,
      mgu.resa, mgu.nrelax, mgp.i, mgp.resb, mgp.resa, mgp.nrelax);
    fflush(fperf);
  }

}

/*
event profiling (i += 20)
{
  static FILE * fp = fopen ("profiling", "a"); // In case of restart
  trace_print (fp, 1); // Display functions taking more than 1% of runtime.
}
*/

/** We also output a movie frame every OUTPUT_FREQ iteration */
// int nb_pic = 0;

void draw_nanoparticles(float s = 3,
                        float pc[3] = {0.4, 0.0, 0.1}) {
  bview * view = draw();
#if dimension == 2
  glTranslatef (0., 0., view->lc*view->fov/24.);
#endif
  glColor3f(pc[0], pc[1], pc[2]);
  glPointSize(s*view->samples);
  glBegin(GL_POINTS);
  foreach_particle() {
#if dimension == 2
    glvertex2d(view, pos.x, pos.y);
#else
    glvertex3d(view, pos.x, pos.y, pos.z);
#endif
  }
  glEnd();
  view->ni++;
}


event pictures (i+=OUTPUT_FREQ) {
  char fname[32];
  view(fov = 25, bg = {1,1,1}, camera = "front", tx = - CAPS(0).centroid.x/L0, ty = -(CAPS(0).centroid.y/L0));
  clear();
  // cells(n = {0,0,1}, alpha = -.49*L0);
  cells(n = {0,0,1}, alpha = CAPS(0).centroid.z/L0);
  draw_nanoparticles();
  //squares("I", n = {0,0,1}, min=-.5, max=1.5);
  //squares("mu.x", n = {0,0,1}, min=1., max=5.);
  // squares("u.x", n = {0,0,1}, alpha = -.49*L0, map = cool_warm);
  draw_lags((_draw_lag){.lw = .5, .edges = true, .facets = true});
  sprintf(fname, "ux_%d.png", nb_pic);
  save(fname);
  nb_pic++;
}

event end (t = TEND) {
// event end (i = 100) {
  fclose(foutput);
  if (foutput_np)
    fclose(foutput_np);
  if (fcrossing_caps)
    fclose(fcrossing_caps);
  return 0.;
}

   
