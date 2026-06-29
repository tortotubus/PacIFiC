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
#define L0_X (1*L0)
#define L0_Y (6*L0)
#define L0_Z (1*L0) 
#define VELPRO_Y_MIN (-0.5*L0_Y)
#define VELPRO_Y_MAX (0.5*L0_Y)
#define VELPRO_N_INTV ((int)((L0_Y/L0)*pow(2, MAXLEVEL)))
#define PI 3.14159265358979323846



#define NCAPS 276
#define NCAPS_L 276
#define NCAPS_S 0
#define VOLUME_FRACTION_L 0.3 //corresponding to radius 0.1
#define VOLUME_FRACTION_S 0.0 //corresponding to radius 0.05, R_a=a_L/a_s=2
#define VOLUME_FRACTION (VOLUME_FRACTION_L + VOLUME_FRACTION_S)
#define RADIUS_L (pow(NCAPS_L*3.1415926535*4/(3*VOLUME_FRACTION_L*6), -1./3))
#define RADIUS_S (0)
#define RADIUS ((pow(RADIUS_L*2., 3)*NCAPS_L + pow(RADIUS_S*2., 3)*NCAPS_S)/(pow(RADIUS_L*2., 2)*NCAPS_L + pow(RADIUS_S*2., 2)*NCAPS_S)/2.) //Sauter-Mean diameter

#define SHEAR_RATE 1.
#ifndef RE
  #define RE 10
#endif
#define MU 1.
#define RHO (RE*MU/(SHEAR_RATE*sq(RADIUS)))


#define PSI_DEG 0
#define XI_DIST 1
#define THETAZX ((180. - PSI_DEG)*PI/180.) // 0-pi
#define DIST_cap01 (RADIUS + RADIUS*XI_DIST/8.)

//Elastic force
#define CA_L 0.2
#define CA_S 0.2 //smaller ones more deformable 
#define E_S_L (MU*RADIUS_L*SHEAR_RATE/CA_L) //E_S_L/E_S_S=1 when we use RADIUS
#define E_S_S (MU*RADIUS_S*SHEAR_RATE/CA_S)

#define ND_EB .025
#define E_B (ND_EB*E_S*sq(RADIUS))
#define REF_CURV 1

#define AVG_VELPRO 0 //velpro

#define STR(s) #s
#define STRING(s) STR(s)

//paralle test
#define ADVECT_LAG_RK2 0
#define LUBR_FORCE 0
#define LUBR_VEL 1 

#define RESTART_CASE 1
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

//#ifndef TEND
//  #define TEND (100./SHEAR_RATE)
//#endif


#ifndef TEND
  #define TEND (100./SHEAR_RATE)
#endif
#ifndef MINLEVEL
  #define MINLEVEL 5
#endif
#ifndef MAXLEVEL
  #define MAXLEVEL 5
#endif
#ifndef LAG_LEVEL_L
  #define LAG_LEVEL_L 3
#endif
#ifndef LAG_LEVEL_S
  #define LAG_LEVEL_S 3
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
  #define OUTPUT_FREQ_PV_DUMP 100
#endif
#ifndef STOKES
  #define STOKES false
#endif
#define JACOBI 1
#define PARAVIEW_CAPSULE 0
#define PARAVIEW_FLOW_FIELD 0
#define OUTPUT_CAPS_NODE_TRI_INFO 0


#define STENCIL_TYPE 3


/**
## Simulation setup

We import the octree grid, the centered Navier-Stokes solver, the Lagrangian
mesh, the neo-Hookean elasticity, a header file containing functions to
mesh a sphere, and the Basilisk viewing functions supplemented by a custom
function $draw\_lag$ useful to visualize the front-tracking interface.
*/
//#include "grid/octree.h"
#include "lagrangian_caps_optim_AM/my_multigrid3D.h"
#include "navier-stokes/centered.h"
#include "lagrangian_caps_optim_AM/capsule-ft.h"
//#include "lagrangian_caps/neo-hookean-ft.h"
#include "lagrangian_caps_optim_AM/skalak-ft.h"
//#include "lagrangian_caps/bending-ft.h"
// #include "lagrangian_caps/viscosity-ft.h"
#include "lagrangian_caps_optim_AM/common-shapes-ft.h"
#include "lagrangian_caps_optim_AM/view-ft.h"
#if PARAVIEW_FLOW_FIELD
  #include "lagrangian_caps_optim_AM/vtu_output.h"
#endif
#include "navier-stokes/perfs.h"
# include "lambda2.h"
# include "view.h"

FILE* fvelocity = NULL; //velpro
FILE* foutput = NULL;
FILE* foutput_fluc = NULL; //velpro
FILE* fperf = NULL;
int init_cycle_number = -1;
int nb_dump = 0;
int nb_pic = 0;

const scalar myrho[] = RHO;
const face vector mymu[] = {MU, MU, MU};
const face vector myalpha[] = {1./RHO, 1./RHO, 1./RHO};


int main(int argc, char* argv[]) {
  /* Set the origin*/
  origin(-0.5*L0_X, -0.5*L0_Y, -0.5*L0_Z); 
  
  /*Set the mpi dimensions*/
  dimensions(nx=1, ny=6, nz=1);

  /** We set periodic boundary conditions on the non-horizontal walls.*/
  periodic(right);
  periodic(front);
  N = 1 << MINLEVEL;

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
u.t[bottom] = dirichlet(0);
u.t[top] = dirichlet(0);
u.r[top] = dirichlet(SHEAR_RATE*y);
u.r[bottom] =  dirichlet(SHEAR_RATE*y);
//uf.n[top] = dirichlet(0);
//uf.n[bottom] = dirichlet(0);

#include <stdio.h>
#include <stdbool.h>
#include <math.h>

int read_average_velocity_profile(const char *filename, double *time,
                                  double *velocities, int max_velocities,
                                  int *vel_count);

event init (i = 0) {


  fvelocity = fopen("velocity.txt", "a"); //velpro
  foutput  = fopen("output.txt","a");
#if (AVG_VELPRO == 1)
  foutput_fluc = fopen("output_fluc.txt", "a"); //velpro
#endif
  fperf = fopen("perf.csv", "a");
  fprintf(fperf, "total_nb_cells, nb_iter_viscous, resb_viscous, resa_viscous, nb_relax_viscous, nb_iter_pressure, resb_pressure, resa_pressure, nb_relax_pressure\n");
 
  // Check if we are in a restart mode
  init_cycle_number = reinitialize_restart();
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


    //Read the particle positions from file
    size_t count = NCAPS;
    const char* ini_pos = "./capsule_positions.txt";
    Cap_read* cap_read;
    cap_read = read_bidisperse_positions(ini_pos, count);


    for (size_t k = 0; k < NCAPS; k++) {
        // printf("Capsule %zu: type=%d, x=%.2f, y=%.2f, z=%.2f\n",
        //        i, cap_read[i].type, cap_read[i].x, cap_read[i].y, cap_read[i].z);
      if(cap_read[k].type == 0)
      {
        //activate the cap mesh
        activate_spherical_capsule((_initialize_circular_capsule){.mesh=&CAPS(k), .cap_es = E_S_L, .cap_radius = RADIUS_L, .level = LAG_LEVEL_L, .cap_id = k, .cap_type = cap_read[k].type});
        //activate_spherical_capsule(&CAPS(k), cap_es = E_S_L, cap_radius = RADIUS_L, level = LAG_LEVEL_L, cap_id = k, cap_type = cap_read[k].type);
        
        CAPS(k).centroid.x = cap_read[k].x;
        CAPS(k).centroid.y = cap_read[k].y;
        CAPS(k).centroid.z = cap_read[k].z;
        for(int i=0; i<CAPS(k).nln; i++)
          foreach_dimension() CAPS(k).nodes[i].pos.x += CAPS(k).centroid.x;
      }
      else if (cap_read[k].type == 1)
      {
        //activate the cap mesh
        //activate_spherical_capsule(&CAPS(k), cap_es = E_S_S, cap_radius = RADIUS_S, level = LAG_LEVEL_S, cap_id = k, cap_type = cap_read[k].type);
        activate_spherical_capsule((_initialize_circular_capsule){.mesh=&CAPS(k), .cap_es = E_S_S, .cap_radius = RADIUS_S, .level = LAG_LEVEL_S, .cap_id = k, .cap_type = cap_read[k].type});

        CAPS(k).centroid.x = cap_read[k].x;
        CAPS(k).centroid.y = cap_read[k].y;
        CAPS(k).centroid.z = cap_read[k].z;
        for(int i=0; i<CAPS(k).nln; i++)
          foreach_dimension() CAPS(k).nodes[i].pos.x += CAPS(k).centroid.x;
      }
      else
      {
        fprintf(stderr,"Unknow capsule type!! Simulation terminated!\n");
        exit(0); 
      }
    }


// generate the stencils around the capsules
  for(int k=0; k<NCAPS; k++) correct_lag_pos(&CAPS(k));
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
    /*In the octree case, make sure the mesh is max refined on the 
    capsule surface */
    #if !(MULT_GRID == 1)   
      astats ss;
      int ic = 0;
      do {
        ic++;
        tag_ibm_stencils();
        ss = adapt_wavelet ({stencils}, (double[]) {1.e-30}, maxlevel = MAXLEVEL, minlevel = MINLEVEL);
        generate_lag_stencils(no_warning = true);
        // rearrange_lag_stencils(no_warning = true);
        fprintf(stderr, "Refine initial mesh: step %d\n", ic);
      } while ((ss.nf || ss.nc) && ic < 100);
    #endif
    
  }

}

/**
The Eulerian mesh is adapted at every time step, according to two criteria:

* first, the 5x5x5 stencils neighboring each Lagrangian node need to be at a
constant level. For this purpose we tag them in the $stencil$ scalar, which is
fed to the $adapt\_wavelet$ algorithm. We also force the top and bottom walls to be at a fine level.
* second, we adapt according to the velocity field.
**/

#if !(MULT_GRID == 1)   
  event adapt (i++) {
    tag_ibm_stencils();
    adapt_wavelet({stencils, u}, (double []){1.e-2, U_TOL, U_TOL, U_TOL},
      minlevel = MINLEVEL, maxlevel = MAXLEVEL);
    generate_lag_stencils(true);  
  }
#endif


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

event logfile (i += 10) {

  int N_pops = 2;

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


#if (AVG_VELPRO == 0)
////////////////////////////////////// Flow velocity profile

   int N_intv = VELPRO_N_INTV;
   double delta_cells = (VELPRO_Y_MAX - VELPRO_Y_MIN) / N_intv;
   double* uy_avg = (double*)calloc(N_intv, sizeof(double));
   double* uy_vol = (double*)calloc(N_intv, sizeof(double));
   foreach(reduction(+:uy_avg[:N_intv]) reduction(+:uy_vol[:N_intv])){
        double Delta3 = (L0*_DELTA)*(L0*_DELTA)*(L0*_DELTA);
        for (int k = 0; k < N_intv; k++){
            if( (y > VELPRO_Y_MIN + k*delta_cells )&&(y <= VELPRO_Y_MIN + (k+1)*delta_cells)){
                uy_avg[k] += u.x[]*Delta3;
                uy_vol[k] += Delta3;
                }
        }
   }

    if (pid() == 0) {
    // Write the position of the averaging layer at t = 0
        if (i == 0) {
          fprintf(fvelocity, "t ");
          for (int k = 0; k < N_intv; k++)
                  fprintf(fvelocity, "%g ", VELPRO_Y_MIN + (k + 0.5)*delta_cells);
         fprintf(fvelocity, " \n");
        }
      fprintf(fvelocity, "%g ", t);
      for (int k = 0; k < N_intv; k++)
        fprintf(fvelocity, "%g ", uy_avg[k]/uy_vol[k]);
      fprintf(fvelocity, " \n");
      fflush(fvelocity);
    }
#else
////////////////////////////////////// Flow Fluctuation Stress //velpro

   int N_intv = VELPRO_N_INTV;
   double delta_cells = (VELPRO_Y_MAX - VELPRO_Y_MIN) / N_intv;
   double intstress_xy = 0.;
   double intstress_xx = 0.;
   double intstress_yy = 0.;
   double intstress_zz = 0.;

   double* velocities = (double*)calloc(N_intv, sizeof(double));
   double time_rd = 0;
   int velocity_count = 0;
   if (read_average_velocity_profile("avgvelprofile.txt", &time_rd, velocities, N_intv, &velocity_count) == 0)
   {}
   else
   {
        return 1;
   }
   if (velocity_count != N_intv) {
        fprintf(stderr, "avgvelprofile.txt has %d velocity bins, expected %d for L0_Y = %g.\n",
                velocity_count, N_intv, L0_Y);
        return 1;
   }

  foreach (reduction(+:intstress_xy) reduction(+:intstress_xx) reduction(+:intstress_yy) reduction(+:intstress_zz))
  {

          double Delta3 = (L0*_DELTA)*(L0*_DELTA)*(L0*_DELTA);

        for (int k = 0; k < N_intv; k++){
            if( (y > VELPRO_Y_MIN + k*delta_cells )&&(y <= VELPRO_Y_MIN + (k+1)*delta_cells)){
                intstress_xy += (u.x[] - velocities[k])*(u.y[] - 0.)*Delta3*RHO;
                intstress_xx += (u.x[] - velocities[k])*(u.x[] - velocities[k])*Delta3*RHO;
                intstress_yy += (u.y[] - 0.)*(u.y[] - 0.)*Delta3*RHO;
          intstress_zz += (u.z[] - 0.)*(u.z[] - 0.)*Delta3*RHO;
            }
        }
  }
#endif




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
	double visc_ratio = 1.;
	double nodal_area = 1.;


        rx = CAPS(k).centroid.x + GENERAL_1DIST(CAPS(k).nodes[i].pos.x, CAPS(k).centroid.x, L0*L0_ratio.x);
        ry = CAPS(k).centroid.y + GENERAL_1DIST(CAPS(k).nodes[i].pos.y, CAPS(k).centroid.y, L0*L0_ratio.y);
        rz = CAPS(k).centroid.z + GENERAL_1DIST(CAPS(k).nodes[i].pos.z, CAPS(k).centroid.z,L0*L0_ratio.z);

	#ifndef CAPS_VISCOSITY
	sigmaxx += - CAPS(k).nodes[i].lagForce.x * rx;
	sigmayy += - CAPS(k).nodes[i].lagForce.y * ry;
	sigmazz += - CAPS(k).nodes[i].lagForce.z * rz;
	sigmaxy += - (CAPS(k).nodes[i].lagForce.x * ry + CAPS(k).nodes[i].lagForce.y * rx) / 2.;
	#else 
	visc_ratio = 1./MUP * MUC;
          /** We now have to compute the area associated with each node */
        nodal_area = compute_node_area(&(CAPS(k)), i);
        
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


int* pop_count = (int*)calloc(N_pops, sizeof(int));
int tot_count = 0;
double* pN1_p = (double*)calloc(N_pops, sizeof(double));
double* pN2_p = (double*)calloc(N_pops, sizeof(double));
double* pmu_p = (double*)calloc(N_pops, sizeof(double));
double* ppres_p = (double*)calloc(N_pops, sizeof(double));

if (pid() == 0) 
{
  tot_count = 0;
  pN1 = 0.;
  pN2 = 0.;
  pmu = 0.;
  ppres = 0.;

   for(int k = 0; k < NCAPS; k++)
   {
    for (int j = 0; j < N_pops; j++)
    { 
      if(CAPS(k).cap_type == j)
      {
        pN1_p[j] += recv_stress_pack[k*4];
        pN2_p[j] += recv_stress_pack[k*4 + 1];
        pmu_p[j] += recv_stress_pack[k*4 + 2];
        ppres_p[j] += recv_stress_pack[k*4 + 3];
        pop_count[j] ++;
        tot_count ++;
      }
    }
   }

   assert(tot_count == NCAPS && "Error pop count is different from NCAPS!\n");
 
}
free(send_stress_pack);
free(recv_stress_pack);

// We compute the averaged quantities by population 
//////////////////////////////////////////////////////////////////////////////////////
  for(int pop_type = 0; pop_type < N_pops; pop_type++)
  {
    if (pid() == 0) 
    {
      double avg_ncaps_area = 0;
      double avg_ncaps_volume = 0;
      for(int k = 0; k < NCAPS; k++) {
        if(CAPS(k).cap_type == pop_type)
        {
          avg_ncaps_volume += CAPS(k).volume/LAG_INITIAL_VOLUME(&CAPS(k));
          for(int i=0; i<CAPS(k).nlt; i++) avg_ncaps_area += CAPS(k).triangles[i].area/(4*pi*sq(CAPS(k).cap_radius));
        }
      }
      
      avg_ncaps_area /= pop_count[pop_type];
      avg_ncaps_volume /= pop_count[pop_type];

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
        if(CAPS(k).cap_type == pop_type)
        {
          compute_taylor_factor(&CAPS(k), &taylor_deform, &inclin_angle, &rs, &TDmaxmin, &TDang);
          foreach_dimension() avg_rs.x += rs.x/CAPS(k).cap_radius; 
          avg_TDmaxmin += TDmaxmin;
          avg_TDang += TDang;
          avg_taylor_deform += taylor_deform;
          avg_inclin_angle += inclin_angle;
          avg_ang_vel += CAPS(k).ang_vel.z;
        }
      }
      
      foreach_dimension() avg_rs.x /= pop_count[pop_type]; 
      avg_TDmaxmin /= pop_count[pop_type];
      avg_TDang /= pop_count[pop_type];
      avg_inclin_angle /= pop_count[pop_type];
      avg_taylor_deform /= pop_count[pop_type];
      avg_ang_vel /= pop_count[pop_type];

      char name[128];
      char default_name[30];
      sprintf(default_name, "%s", result_dir );
      strcat(default_name, "/" );
      strcat(default_name, "Individual_outputs\0" );
      // char* prefix = p.name ? p.name : default_name;
      char* prefix = default_name;
      char suffix[64];
      sprintf(suffix, "_cap%d.txt", pop_type);
      sprintf(name, "%s%s", prefix, suffix);
      FILE* foutput_physics = fopen(name, "a+");
      assert(foutput_physics);

      fprintf(foutput_physics, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", iter, t, 
        ppres_p[pop_type], pmu_p[pop_type], pN1_p[pop_type], pN2_p[pop_type], 
        avg_taylor_deform, avg_inclin_angle, avg_TDmaxmin, avg_TDang,
        avg_rs.x, avg_rs.y, avg_rs.z, avg_ang_vel, avg_ncaps_area, avg_ncaps_volume);
      fflush(foutput_physics);

      fclose(foutput_physics);
    }
  }


// We compute the ensemble average quantities
//////////////////////////////////////////////////////////////////////////////////////
  if (pid() == 0) 
    {
      double avg_ncaps_area = 0;
      double avg_ncaps_volume = 0;
      for(int k = 0; k < NCAPS; k++) {
          avg_ncaps_volume += CAPS(k).volume/LAG_INITIAL_VOLUME(&CAPS(k));
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

      char name[128];
      sprintf(name, "output.txt");
      FILE* fouttot_physics = fopen(name, "a+");
      assert(fouttot_physics);

      fprintf(fouttot_physics, "%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", iter, t, 
        fluid_visc_stress, ppres_p[0]+ppres_p[1], pmu_p[0]+pmu_p[1], pN1_p[0]+pN1_p[1], pN2_p[0]+pN2_p[1], 
        avg_taylor_deform, avg_inclin_angle, avg_TDmaxmin, avg_TDang,
        avg_rs.x, avg_rs.y, avg_rs.z, avg_ang_vel, avg_ncaps_area, avg_ncaps_volume);
      fflush(fouttot_physics);

      fclose(fouttot_physics);

#if (AVG_VELPRO==1)    //velpro
    fprintf(foutput_fluc, "%d %lf %lf %lf %lf \n", i, t,
    intstress_xy, (intstress_xx - intstress_yy), (intstress_yy - intstress_zz));
    fflush(foutput_fluc);
#endif

}
 

free(pop_count);
free(pN1_p);
free(pN2_p);
free(pmu_p);
free(ppres_p);
}

// event logfile2(i++)
// {
//   int N_pops = 2;
//   output_bidispse_physics(N_pops, i);
// }

/*
event profiling (i += 20)
{
  static FILE * fp = fopen ("profiling", "a"); // In case of restart
  trace_print (fp, 1); // Display functions taking more than 1% of runtime.
}
*/

// /** We also output a movie frame every OUTPUT_FREQ iteration */
// // int nb_pic = 0;
// event pictures (i+=OUTPUT_FREQ) {
//   char fname[32];
//   view(fov = 20, bg = {1,1,1}, camera = "front");
//   clear();
//   // cells(n = {0,0,1}, alpha = -.49*L0);
//   cells(n = {0,0,1});
//   //squares("I", n = {0,0,1}, min=-.5, max=1.5);
//   //squares("mu.x", n = {0,0,1}, min=1., max=5.);
//   // squares("u.x", n = {0,0,1}, alpha = -.49*L0, map = cool_warm);
//   draw_lags((_draw_lag){.lw = .5, .edges = true, .facets = true});
//   sprintf(fname, "ux_%d.png", nb_pic);
//   save(fname);
//   nb_pic++;
// }

event end (t = TEND) {
  fclose(foutput);
#if (AVG_VELPRO == 1)
  fclose(foutput_fluc);
#endif
  return 0.;
}

int read_average_velocity_profile(const char *filename, double *time, double *velocities, int max_velocities, int *vel_count) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        perror("Error opening file");
        return 1;
    }

    char buffer[20000];
    char last_line[20000];

    last_line[0] = '\0';
    while (fgets(buffer, sizeof(buffer), fp) != NULL) {
        strcpy(last_line, buffer);
    }

    fclose(fp);

    if (last_line[0] == '\0') {
        fprintf(stderr, "The file seems to be empty or no line read.\n");
        return 2;
    }

    // Parse the line
    char *token = strtok(last_line, " \t\n");
    if (token == NULL) {
        fprintf(stderr, "No data found in the last line.\n");
        return 3;
    }

    // First token is time
    *time = atof(token);

    // Reset velocity count
    *vel_count = 0;
    token = strtok(NULL, " \t\n");
    while (token != NULL && *vel_count < max_velocities) {
        velocities[*vel_count] = atof(token);
        (*vel_count)++;
        token = strtok(NULL, " \t\n");
    }
    return 0;
}
