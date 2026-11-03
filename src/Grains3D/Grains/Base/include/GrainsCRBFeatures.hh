#ifndef _GRAINSCRBFEATURES_HH_
#define _GRAINSCRBFEATURES_HH_

#include "Grains.hh"
#include "AllComponents.hh"
#include "App.hh"
#include "LinkedCell.hh"

#include <list>
#include <vector>
#include <string>
#include <iostream>
using namespace std;

#include "Xml.hh"
#include "Basic.hh"


/** @brief The class GrainsCRBFeatures.

    Application that computes the volume, center of mass coordinates and moment
    of inertia tensor of a composite rigid body.

    @author A.WACHS - Institut Francais du Petrole - 2010 - Creation
    @author A.WACHS - 2021 - Major cleaning & refactoring */
// ============================================================================
class GrainsCRBFeatures : public Grains
{
  public:
    /** @name Contructors & Destructor */
    //@{
    /** @brief Default constructor */
    GrainsCRBFeatures();

    /** @brief Destructor */
    virtual ~GrainsCRBFeatures();
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
    
    
    /** @name I/O methods */
    //@{  
    /** @brief Writes an initial message in the standard output only
    on the process ranked 0 */
    virtual void initialOutputMessage();        
    //@}       

  
  private:
    /** @name Parameters */
    //@{  
    string m_outputfilename; /**< output file name */
    int m_geomtypeID; /**< geometric type number of the composite particle */
    double m_eqsphrad; /**< equivalent sphere radius of the composite 
    	particle */
    double m_dimensionless_delta; /**< dimensionless grid cell size */
    size_t m_Nmax; /**< maximum number of cells of the grid in all 
    	directions */         
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
};

#endif
  
