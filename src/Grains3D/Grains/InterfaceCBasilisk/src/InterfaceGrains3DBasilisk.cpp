#include "SimulationXml.hh"
/**
# Coupling interface for Grains3D and Basilisk 
*/

#include "GrainsBuilderFactory.hh"
#include "InterfaceGrains3DBasilisk.h"
#include "Exec.hh"
#include <iostream>
#include <cstring>
#include <string>

#ifdef __cplusplus
extern "C" {
#endif
  
  static GrainsCoupledWithFluid* coupledGrains = NULL;

  
  void Init_Grains ( char const* inputfile, 
  	double fluid_rho, const bool b_restart,
	const bool b_fluidcorrectedacc ) 
  {
    string simulation_file( inputfile );
    string simulation_file_exe = GrainsBuilderFactory::init( simulation_file, 
    	0, 1 );
  
    grains::xml::OptionalElement rootNode = grains::xml::readSimulationFile(simulation_file_exe).root();
    
    coupledGrains = GrainsBuilderFactory::createCoupledWithFluid( rootNode,
    	fluid_rho );
    if ( b_restart ) coupledGrains->setReloadSame();
    coupledGrains->do_before_time_stepping( rootNode );
    // REPLACED_EXEC_SCRIPTS
    
    // string cmd = "/bin/rm " + simulation_file_exe;
    // GrainsExec::m_return_syscmd = system( cmd.c_str() );
    
    std::error_code ec;
    const bool removed = std::filesystem::remove(std::filesystem::path(simulation_file_exe), ec);
    
    if ( !b_fluidcorrectedacc )
      coupledGrains->setFluidCorrectedAcceleration( b_fluidcorrectedacc );
         
    cout << "Construction of Grains completed" << endl;
  }




  void Simu_Grains( const double dt_fluid ) 
  {
    coupledGrains->Simulation( dt_fluid );
  }

 
 
  
  char* GrainsToBasilisk( int* pstrsize )
  {
    // We use the interface function of PeliGRIFF
    istringstream iss;
    coupledGrains->GrainsToFluid2( iss );

    // We remove the formatting and separate each entry by " "
    string buff;
    *pstrsize = int(iss.str().size()) + 1;    
    char* pstr = new char [*pstrsize];
    int pos = 0;
    while ( iss >> buff )
    {
      std::strcpy( &pstr[pos], buff.c_str() );
      pos += int(buff.size());
      std::strcpy( &pstr[pos], " " ); 
      pos += 1;     
    }    
    
    // Return the pointer to the char
    return pstr;
  } 




  void SaveResults_Grains() 
  {
    static unsigned int ppcounter = 0;
    if ( !ppcounter ) 
      coupledGrains->InitialPostProcessing( 6 );
    else
      coupledGrains->doPostProcessing( 6 );
    ++ppcounter; 
  }  




  void checkParaviewPostProcessing_Grains( char* solid_resDir )
  {    
    coupledGrains->checkParaviewPostProcessing( "grains", solid_resDir, true );
  }



  
  void UpdateVelocityGrains( double** arrayv, const int m ) 
  {
    // Transfer into a vector< vector<double> >
    vector<double> buf( 6, 0.);
    vector< vector<double> > vecv( m, buf );
    for (size_t i=0;i<size_t(m);++i)
      for (size_t j=0;j<6;++j)
        vecv[i][j] = arrayv[i][j];

    // We use the interface function of PeliGRIFF
    coupledGrains->updateParticlesVelocity( vecv, false );
  }  



  
  void SetInitialCycleNumber( int cycle0 ) 
  {
    coupledGrains->setInitialCycleNumber( cycle0 );
  }



  
  void SetInitialTime( double tinit ) 
  {
    coupledGrains->setInitialTime( tinit );
  } 
  
  
  
     
  void NumberOfRigidBodiesInBasilisk( size_t* nparticles, size_t* nobstacles )
  {
    coupledGrains->numberOfRBToFluid( nparticles, nobstacles );
  }
  
  
  
  
  void NumberOfReferenceRigidBodiesInBasilisk( size_t* nrefRB )
  {
    coupledGrains->numberOfReferenceRBToFluid( nrefRB );
  }



  
  char* GrainsToBasiliskReference( int* pstrsize )
  {
    // We use the interface function of PeliGRIFF
    istringstream iss;
    coupledGrains->GrainsToFluidReference( iss );

    // We remove the formatting and separate each entry by " "
    string buff;
    *pstrsize = int(iss.str().size()) + 1;    
    char* pstr = new char [*pstrsize];
    int pos = 0;
    while ( iss >> buff )
    {
      std::strcpy( &pstr[pos], buff.c_str() );
      pos += int(buff.size());
      std::strcpy( &pstr[pos], " " ); 
      pos += 1;     
    }    
    
    // Return the pointer to the char
    return pstr;
  }  
      
#ifdef __cplusplus
}
#endif
