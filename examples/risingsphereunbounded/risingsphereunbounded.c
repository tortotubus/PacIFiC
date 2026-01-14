/**
One lighter sphere in a large box
*/

# include "grid/octree.h"
# define DLM_ALPHA_COUPLING 1
# define DLMFD_OPT 1
# define DLMFD_PROB_AFTER_NAVIERSTOKES 0
# define INITIALGRIDADAPTIVE_NEWMETHOD 1
# define B_SPLIT_EXPLICIT_ACCELERATION 1
# define DLM_UZAWA_TOL 1.e-5


/* Parameters related to the spatial resolution */
# define LEVEL 2
# define ADAPTIVE 1 

/* if adaptivity is enabled define the maximum level of refinement MAXLEVEL */
# if ADAPTIVE
#   define MAXLEVEL (LEVEL + 8)
#   define FLAG_ADAPT_CRIT (1.E-9)
#   define UX_ADAPT_CRIT (1.E-3)
#   define UY_ADAPT_CRIT (1.E-3)
#   define UZ_ADAPT_CRIT (5.E-4)
# endif

/* Physical parameters */
# define FLUID_DENSITY 1000. // fluid density
# define FLUID_VISCOSITY 1.771779e-03 // fluid dynamic viscosity
# define GRAVITY_VECTOR ((coord){0.,0.,-9.81})
# define rho_solid 900. // rigid body density
# define Deq 0.002 // equivalent diameter
# define Density_ratio ( rho_solid / FLUID_DENSITY ) // density ratio
# define U_c ( sqrt( 4. * fabs( Density_ratio - 1. ) * 9.81 * Deq / 3. ) )
# define T_c ( Deq / U_c ) // characteristic time scale
# define Re ( FLUID_DENSITY * Deq * U_c / FLUID_VISCOSITY )

/* output and numerical parameters */
# define MAXDT (4.e-5) // time-step
# define TINTERVALOUTPUT (1.E-2) // time interval between outputs
# define SIMUTIMEINTERVAL (1*TINTERVALOUTPUT) // simulation time interval

/* Post-processing */
# define VORTICITY 1
# define PARAVIEW_VTU 1
# define PARAVIEW_HTG 1
# define PARAVIEW_BINFILE 1
# define PARAVIEW_VTU_MPIIO_WRITER 1
# define PARAVIEW_SCALAR_LIST p,DLM_Flag,DLM_FlagMesh
# define PARAVIEW_VECTOR_LIST u,DLM_Index
# define PARAVIEW_DLMFD_BNDPTS 1
# define PARAVIEW_DLMFD_INTPTS 1
# define BVIEW 0

/* 1st order Marshuk-Yanenko implementation of the DLM-FD method */
# include "DLMFD_Grains3D.h"


int main () {

  origin (0., 0., 0.);

  L0 = 0.2; // length of the fluid-domain

  /* set time step */
  DT = MAXDT;

  /* initialise a uniform grid */
  init_grid (1 << LEVEL);

  /* boundary conditions */
  u.t[left]   = dirichlet(0.); //v
  u.r[left]   = dirichlet(0.); //w
  u.n[left]   = dirichlet(0.); //u

  u.t[right]  = dirichlet(0.); //v
  u.r[right]  = dirichlet(0.); //w
  u.n[right]  = dirichlet(0.); //u

  u.n[bottom] = dirichlet(0.); //v
  u.t[bottom] = dirichlet(0.); //w
  u.r[bottom] = dirichlet(0.); //u

  u.n[top]    = dirichlet(0.); //v
  u.t[top]    = dirichlet(0.); //w
  u.r[top]    = dirichlet(0.); //u

  u.n[back]   = dirichlet(0.); //w 
  u.t[back]   = dirichlet(0.); //u 
  u.r[back]   = dirichlet(0.); //v 

  u.n[front]   = dirichlet(0.); //w 
  u.t[front]   = dirichlet(0.); //u 
  u.r[front]   = dirichlet(0.); //v

  // Convergence criteria for the N-S solver
  TOLERANCE = 1e-5;

  run();
}


event init (i = 0) {
  if ( pid() == 0 ) 
  {
    printf( "Time scale: %8.5f \n", T_c );  
    printf( "Velocity scale: %8.5f \n", U_c );
    printf( "Reynolds number: %8.5f \n", Re ); 
  }   
}
