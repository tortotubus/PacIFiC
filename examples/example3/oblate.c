/**
One flat oblate in a tall columnar box
*/

# include "grid/octree.h"
# define ADAPTIVE 1

//# include "DLMFD_multigrid3D.h"
//# define ADAPTIVE 0

# define DLM_ALPHA_COUPLING 1
# define DLMFD_OPT 1
# define DLMFD_PROB_AFTER_NAVIERSTOKES 0
# define INITIALGRIDADAPTIVE_NEWMETHOD 1
# define RIGIDBODY_VERBOSE 0
# define INTERBPCOEF 1.75
# define VELOCITYCHANGE_COMPONENT z
# define PRODUCTION_LOG 1
# define SINGLE_SETTLING_SHIFT_TO_ZERO_Z 1


/* Physical parameters */
# define FLUID_DENSITY 1.225 // fluid density
# define GRAVITY_VECTOR ((coord){0.,0.,-9.81})
# define SOLID_DENSITY 122.5 // rigid body density
# define SOLID_VOLUME 6.2832E-07 // rigid body volume
# define LC 0.02 // characteristic length
# define SOLID_DIMLESS_VOLUME ( SOLID_VOLUME / pow( LC, 3. ) ) 
# define DENSITY_RATIO ( SOLID_DENSITY / FLUID_DENSITY ) // density ratio
# define SOLID_DIMLESS_MASS ( DENSITY_RATIO * SOLID_DIMLESS_VOLUME )
# define GALILEO 50. // Galileo number
# define UCGA ( sqrt( 8. * 9.81 * LC * \
	fabs( SOLID_DIMLESS_MASS - SOLID_DIMLESS_VOLUME ) / pi ) ) \
	// ref velocity
# define FLUID_VISCOSITY ( FLUID_DENSITY * UCGA * LC / GALILEO ) // fluid \
	// dynamic viscosity
# define TC ( LC / UCGA ) // characteristic time scale


/* Parameters related to the spatial resolution 
If adaptivity is enabled define the maximum level of refinement MAXLEVEL */
# define LEVEL 1
# if ADAPTIVE
#   define MAXLEVEL (LEVEL + 12)
#   define LEVELDIFF_FLAG_U 1
#   define UX_ADAPT_CRIT (2.E-3 * UCGA)
#   define UY_ADAPT_CRIT (2.E-3 * UCGA)
#   define UZ_ADAPT_CRIT (2.E-3 * UCGA)
# else
#   define MAXLEVEL (LEVEL)
# endif


/* output and numerical parameters */
# define MAXDT ( TC / 2000. ) // time-step
# define TINTERVALOUTPUT (2.5E-1) // time interval between outputs
# define SIMUTIMEINTERVAL (1*TINTERVALOUTPUT) // simulation time interval


/* Post-processing */
# define VORTICITY 1
# define LAMBDA2 1
# define PARAVIEW_VTU 0
# define PARAVIEW_HTG 1
# define PARAVIEW_BINFILE 1
# define PARAVIEW_VTU_MPIIO_WRITER 1
# define PARAVIEW_SCALAR_LIST p
# define PARAVIEW_VECTOR_LIST u
# define PARAVIEW_DLMFD_BNDPTS 1
# define PARAVIEW_DLMFD_INTPTS 1
# define BVIEW 0


/* 1st order Marshuk-Yanenko implementation of the DLM-FD method */
# include "DLMFD_Grains3D.h"


int main () {

  /* Define fluid domain */
  origin( 0.2, 0.2, 0.0 );
  size( 0.8 ); 
  dimensions( nx = 1, ny = 1, nz = 25 );


  /* set time step */
  DT = MAXDT;

  /* initialise a uniform grid */
  init_grid( 1 << LEVEL );

  /* boundary conditions */
  p[front]     = dirichlet(0.);
  pf[front]    = dirichlet(0.);  
  u.n[front]   = neumann(0.);
  u.t[front]   = dirichlet(0.); 
  u.r[front]   = dirichlet(0.);

  u.n[back]   = dirichlet(0.); //w 
  u.t[back]   = dirichlet(0.); //u 
  u.r[back]   = dirichlet(0.); //v

  periodic(right);
  periodic(top);

  // Convergence criteria for the N-S solver
  TOLERANCE = 1e-5;

  run();
}


event init (i = 0) {
  if ( pid() == 0 ) 
  {
    printf( "Particle density: %8.5f \n", SOLID_DENSITY );
    printf( "Fluid density: %8.5f \n", FLUID_DENSITY );         
    printf( "Fluid dynamic viscosity: %8.5f \n", FLUID_VISCOSITY );    
    printf( "Length scale: %8.5f \n", LC );  
    printf( "Velocity scale: %8.5f \n", UCGA );    
    printf( "Time scale: %8.5f \n", TC );
    printf( "Density ratio: %8.5f \n", DENSITY_RATIO );
    printf( "Dimensionless volume: %8.5f \n", SOLID_DIMLESS_VOLUME );
    printf( "Dimensionless mass: %8.5f \n", SOLID_DIMLESS_MASS );           
    printf( "Galileo number: %8.5f \n", GALILEO );       
    printf( "Points per ellipsoid diameter: %u \n", 
    	(unsigned int) ( pow(2.,MAXLEVEL) * LC / L0 ) );     
  } 
}
