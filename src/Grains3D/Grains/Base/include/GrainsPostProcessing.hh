#ifndef _GRAINSPOSTPROCESSING_HH_
#define _GRAINSPOSTPROCESSING_HH_

#include "Grains.hh"
#include "Xml.hh"
#include "Basic.hh"
#include <list>
#include <string>
using namespace std;


/** @brief Global porosity */
struct GlobalPorosity
{
  Window domain; /**< porosity domain */
  size_t nintervals[3]; /**< Number of intervals in each direction */   
};


/** @brief FreeSurface */
struct FreeSurface
{
  Window domain; /**< porosity domain */
  double dl; /**< Ray tracing step */
  double min_height; /**< Minimal height if ray tracing fails */
  size_t n0; /**< number of grid points in direction 0 */
  size_t n1; /**< number of grid points in direction 1 */ 
  string filename; /**< VTU output file name */
  size_t nlevels; /**< Number of times (minus 1 ) the step is divided by 10 */   
};


/** @brief The class GrainsPostProcessing.

    Grains3D application to postprocess the granular medium at a given time
    through reloading the system as in a restarted simulation.

    @author A.WACHS - 2022 - Major cleaning & refactoring */
// ============================================================================
class GrainsPostProcessing : virtual public Grains
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor 
    @param fluid_density_ fluid density */
    GrainsPostProcessing();

    /** @brief Destructor */
    virtual ~GrainsPostProcessing();
    //@}


    /** @name Core methods */
    //@{
    /** @brief Tasks to perform before time-stepping 
    @param rootElement XML root */
    virtual void do_before_time_stepping( grains::xml::OptionalElement rootElement );

    /** @brief Tasks to perform after time-stepping */
    virtual void do_after_time_stepping();

    /** @brief Runs the simulation over the prescribed time interval
    @param time_interval runs the simulation on a time interval different than 
    	m_tend - m_tstart  */
    virtual void Simulation( double time_interval = 0. ) ; 	
    //@}


    /** @name Methods */
    //@{  
    //@}
    
    
    /** @name I/O methods */
    //@{  
    /** @brief Writes an initial message in the standard output only
    on the process ranked 0 */
    virtual void initialOutputMessage(); 	          
    //@}
    

  protected:
    /** @name Parameters */
    //@{ 
    struct GlobalPorosity* m_global_porosity;
    struct FreeSurface* m_free_surface;    	       
    //@}


    /**@name Methods */
    //@{
    /** @brief Construction of the simulation: linked cell, particles &
    obstacles, domain decomposition 
    @param rootElement XML root */
    virtual void Construction( grains::xml::OptionalElement rootElement );

    /** @brief Additional features of the simulation: time features, insertion,
    post-processing
    @param rootElement XML root */
    virtual void AdditionalFeatures( grains::xml::OptionalElement rootElement );

    /** @brief External force definition
    @param rootElement XML root */
    virtual void Forces( grains::xml::OptionalElement rootElement );
    //@}    
    
    
  private:
};

#endif
