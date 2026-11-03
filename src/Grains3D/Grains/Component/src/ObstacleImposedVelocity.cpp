#include "ObstacleImposedVelocity.hh"
#include "Obstacle.hh"
#include "GrainsExec.hh"
#include "LinkedCell.hh"
#include <stdlib.h>


// ----------------------------------------------------------------------------
// Default constructor
ObstacleImposedVelocity::ObstacleImposedVelocity()
{
  m_ObstacleName = "Undefined";
  m_type = "Undefined";
  m_tstart = 0.; 
  m_tend = 0.;
  m_translationalVelocity = Vector3Null;
  m_previous_translationalVelocity = Vector3Null;
  m_angularVelocity = Vector3Null;
  m_rotationCenterIsCenterOfMass = true;
  m_rotationCenter = Point3Null;
  m_amplitude = 0.;
  m_period = 0.;
  m_Sin_phase_shift = 0.;  
  m_unit_vitRef = Vector3Null;
  m_MultiSin_period = Vector3Null;
  m_MultiSin_amplitude = Vector3Null;
  m_MultiSin_phase_shift = Vector3Null; 
  m_stress_max = 0.;
  m_stress = 0.;
  m_previous_stress = 0.; 
}




// ----------------------------------------------------------------------------
// Constructor with an XML node as input parameter
ObstacleImposedVelocity::ObstacleImposedVelocity( grains::xml::OptionalElement root,
	double dt, int rank, size_t& error )
{
  m_type = "Undefined";
  m_translationalVelocity = Vector3Null;
  m_previous_translationalVelocity = Vector3Null;  
  m_angularVelocity = Vector3Null;
  m_rotationCenterIsCenterOfMass = true;
  m_rotationCenter = Point3Null;
  m_amplitude = 0.;
  m_period = 0.;
  m_Sin_phase_shift = 0.;  
  m_unit_vitRef = Vector3Null;
  m_MultiSin_period = Vector3Null;
  m_MultiSin_amplitude = Vector3Null;
  m_MultiSin_phase_shift = Vector3Null;
  m_stress_max = 0.;
  m_stress = 0.;
  m_previous_stress = 0.; 
    
  error = 0;   
    
  m_ObstacleName = grains::xml::requireElement(root).attributeOr<std::string>("ObstacleName", std::string{});
  
  grains::xml::OptionalElement nTimeInterval = grains::xml::requireElement(root).getChild("TimeInterval");
  m_tstart = grains::xml::requireElement(nTimeInterval).requireAttribute<double>("Start");
  m_tend = grains::xml::requireElement(nTimeInterval).requireAttribute<double>("End");

  // Constant translation
  if ( grains::xml::requireElement(root).getChild("ConstantTranslation") )
  {
    m_type = "ConstantTranslation";
    grains::xml::OptionalElement nVTranslation = grains::xml::requireElement(root).getChild("ConstantTranslation");
    grains::xml::OptionalElement nTV = grains::xml::requireElement(nVTranslation).getChild("TranslationalVelocity");
    m_translationalVelocity[X] = grains::xml::requireElement(nTV).requireAttribute<double>("VX");
    m_translationalVelocity[Y] = grains::xml::requireElement(nTV).requireAttribute<double>("VY");
    m_translationalVelocity[Z] = grains::xml::requireElement(nTV).requireAttribute<double>("VZ");
    if ( rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Obstacle name = " << m_ObstacleName 
      	<< endl;      
      cout << GrainsExec::m_shift12 << "Time interval = [" 
      	<< m_tstart << "," << m_tend << "]" << endl;
      cout << GrainsExec::m_shift12 << "Type = " << m_type << endl;
      cout << GrainsExec::m_shift12 << "Translational velocity = " << 
      	m_translationalVelocity << endl;
    }   
  }
  // Sinusoidal translation  
  else if ( grains::xml::requireElement(root).getChild("SinTranslation") )
  {
    m_type = "SinTranslation";
    grains::xml::OptionalElement nVTranslation = grains::xml::requireElement(root).getChild("SinTranslation");
    grains::xml::OptionalElement nTV = grains::xml::requireElement(nVTranslation).getChild("Direction");
    m_unit_vitRef[X] = grains::xml::requireElement(nTV).requireAttribute<double>("VX");
    m_unit_vitRef[Y] = grains::xml::requireElement(nTV).requireAttribute<double>("VY");
    m_unit_vitRef[Z] = grains::xml::requireElement(nTV).requireAttribute<double>("VZ");
    m_unit_vitRef.normalize();
    grains::xml::OptionalElement nPar = grains::xml::requireElement(nVTranslation).getChild("Parameters");
    m_amplitude = grains::xml::requireElement(nPar).requireAttribute<double>("Amplitude");
    m_period = grains::xml::requireElement(nPar).requireAttribute<double>("Period");
    m_Sin_phase_shift = grains::xml::requireElement(nPar).requireAttribute<double>("PhaseShift")
    	* PI / 180.;
    if ( rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Obstacle name = " << m_ObstacleName 
      	<< endl;      
      cout << GrainsExec::m_shift12 << "Time interval = [" 
      	<< m_tstart << "," << m_tend << "]" << endl;
      cout << GrainsExec::m_shift12 << "Type = " << m_type << endl;
      cout << GrainsExec::m_shift12 << "Unit direction vector = " << 
      	m_unit_vitRef << endl;
      cout << GrainsExec::m_shift12 << "Amplitude = " << 
      	m_amplitude << endl;
      cout << GrainsExec::m_shift12 << "Period = " << 
      	m_period << endl;
      cout << GrainsExec::m_shift12 << "Phase shift in rad = " << 
      	m_Sin_phase_shift << endl;	
      cout << GrainsExec::m_shift12 << "Maximum motion = " << 
      	m_amplitude * m_period / ( 2. * PI ) << endl;
      cout << GrainsExec::m_shift12 << "Maximum acceleration = " << 
      	m_amplitude * 2. * PI / m_period << endl; 
    }             
  }
  // Sinusoidal translation  
  else if ( grains::xml::requireElement(root).getChild("CrenelTranslation") )
  {
    m_type = "CrenelTranslation";
    grains::xml::OptionalElement nVTranslation = grains::xml::requireElement(root).getChild("CrenelTranslation");
    grains::xml::OptionalElement nTV = grains::xml::requireElement(nVTranslation).getChild("Direction");
    m_unit_vitRef[X] = grains::xml::requireElement(nTV).requireAttribute<double>("VX");
    m_unit_vitRef[Y] = grains::xml::requireElement(nTV).requireAttribute<double>("VY");
    m_unit_vitRef[Z] = grains::xml::requireElement(nTV).requireAttribute<double>("VZ");
    m_unit_vitRef.normalize();
    grains::xml::OptionalElement nPar = grains::xml::requireElement(nVTranslation).getChild("Parameters");
    m_amplitude = grains::xml::requireElement(nPar).requireAttribute<double>("Amplitude");
    m_period = grains::xml::requireElement(nPar).requireAttribute<double>("Period");
    if ( rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Obstacle name = " << m_ObstacleName 
      	<< endl;      
      cout << GrainsExec::m_shift12 << "Time interval = [" 
      	<< m_tstart << "," << m_tend << "]" << endl;
      cout << GrainsExec::m_shift12 << "Type = " << m_type << endl;
      cout << GrainsExec::m_shift12 << "Unit direction vector = " << 
      	m_unit_vitRef << endl;
      cout << GrainsExec::m_shift12 << "Amplitude = " << 
      	m_amplitude << endl;
      cout << GrainsExec::m_shift12 << "Period = " << 
      	m_period << endl;	
      cout << GrainsExec::m_shift12 << "Maximum motion = " << 
      	m_amplitude * m_period << endl;
    }             
  } 
  // Shear stress dependent cyclic translation  
  else if ( grains::xml::requireElement(root).getChild("SSDCyclicTranslation") )
  {
    m_type = "SSDCyclicTranslation";
    grains::xml::OptionalElement nVTranslation = grains::xml::requireElement(root).getChild("SSDCyclicTranslation");
    grains::xml::OptionalElement nTV = grains::xml::requireElement(nVTranslation).getChild("Direction");
    string sfirst = grains::xml::requireElement(nTV).attributeOr<std::string>("K", std::string{});
    if ( sfirst == "X" ) m_stressIndices.first = 0;
    else if ( sfirst == "Y" ) m_stressIndices.first = 1;
    else if ( sfirst == "Z" ) m_stressIndices.first = 2;
    else
      cout << "WARNING: unknown direction in SSDCyclicTranslation associated"
      	<< " to " << m_ObstacleName << endl;
    m_unit_vitRef[m_stressIndices.first] = 1.;
    string ssecond = grains::xml::requireElement(nTV).attributeOr<std::string>("L", std::string{});
    if ( ssecond == "X" ) m_stressIndices.second = 0;
    else if ( ssecond == "Y" ) m_stressIndices.second = 1;
    else if ( ssecond == "Z" ) m_stressIndices.second = 2;
    else
      cout << "WARNING: unknown direction in SSDCyclicTranslation associated"
      	<< " to " << m_ObstacleName << endl;    
    grains::xml::OptionalElement nPar = grains::xml::requireElement(nVTranslation).getChild("Parameters");
    m_amplitude = grains::xml::requireElement(nPar).requireAttribute<double>("Amplitude");
    m_stress_max = grains::xml::requireElement(nPar).requireAttribute<double>("MaxStress");
    m_translationalVelocity[m_stressIndices.first] = m_amplitude;
    if ( rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Obstacle name = " << m_ObstacleName 
      	<< endl;      
      cout << GrainsExec::m_shift12 << "Time interval = [" 
      	<< m_tstart << "," << m_tend << "]" << endl;
      cout << GrainsExec::m_shift12 << "Type = " << m_type << endl;
      cout << GrainsExec::m_shift12 << "Unit direction vector = " << 
      	m_unit_vitRef << endl;
      cout << GrainsExec::m_shift12 << "Amplitude = " << 
      	m_amplitude << endl;
      cout << GrainsExec::m_shift12 << "MaxStress = " << 
      	m_stress_max << endl;	
      cout << GrainsExec::m_shift12 << "Stress component = " << 
      	sfirst << ssecond << endl;
    }             
  }    
  // Constant rotation
  else if ( grains::xml::requireElement(root).getChild("ConstantRotation") )
  {
    m_type = "ConstantRotation";
    grains::xml::OptionalElement nVRotation = grains::xml::requireElement(root).getChild("ConstantRotation");
    grains::xml::OptionalElement nRV = grains::xml::requireElement(nVRotation).getChild("AngularVelocity");
    m_angularVelocity[X] = grains::xml::requireElement(nRV).requireAttribute<double>("WX");
    m_angularVelocity[Y] = grains::xml::requireElement(nRV).requireAttribute<double>("WY");
    m_angularVelocity[Z] = grains::xml::requireElement(nRV).requireAttribute<double>("WZ");
    grains::xml::OptionalElement nRC = grains::xml::requireElement(nVRotation).getChild("RotationCenter");
    if ( nRC )
    {
      m_rotationCenterIsCenterOfMass = false;
      m_rotationCenter[X] = grains::xml::requireElement(nRC).requireAttribute<double>("CX");
      m_rotationCenter[Y] = grains::xml::requireElement(nRC).requireAttribute<double>("CY");
      m_rotationCenter[Z] = grains::xml::requireElement(nRC).requireAttribute<double>("CZ");
    }      
    if ( rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Obstacle name = " << m_ObstacleName 
      	<< endl;      
      cout << GrainsExec::m_shift12 << "Time interval = [" 
      	<< m_tstart << "," << m_tend << "]" << endl;
      cout << GrainsExec::m_shift12 << "Type = " << m_type << endl;
      cout << GrainsExec::m_shift12 << "Angular velocity = " << 
      	m_angularVelocity << endl;
      cout << GrainsExec::m_shift12 << "Center of rotation = ";
      if ( m_rotationCenterIsCenterOfMass ) cout << "center of mass";
      else cout << m_rotationCenter[X] << " " << m_rotationCenter[Y] <<
	" " <<  m_rotationCenter[Z]; 
      cout << endl;	
    }   
  }
  // Multi-dimensional sinusoidal translation
  else if ( grains::xml::requireElement(root).getChild("MultiSinTranslation") )
  {
    m_type = "MultiSinTranslation";
    grains::xml::OptionalElement nMS = grains::xml::requireElement(root).getChild("MultiSinTranslation");
    grains::xml::OptionalElement nTV = grains::xml::requireElement(nMS).getChild("Direction");
    m_unit_vitRef[X] = grains::xml::requireElement(nTV).requireAttribute<double>("VX");
    m_unit_vitRef[Y] = grains::xml::requireElement(nTV).requireAttribute<double>("VY");
    m_unit_vitRef[Z] = grains::xml::requireElement(nTV).requireAttribute<double>("VZ");
    m_unit_vitRef.normalize();    
    grains::xml::OptionalElement nPer = grains::xml::requireElement(nMS).getChild("Period");
    m_MultiSin_period[X] = grains::xml::requireElement(nPer).requireAttribute<double>("PX");
    m_MultiSin_period[Y] = grains::xml::requireElement(nPer).requireAttribute<double>("PY");
    m_MultiSin_period[Z] = grains::xml::requireElement(nPer).requireAttribute<double>("PZ");
    grains::xml::OptionalElement nPhaseShift = grains::xml::requireElement(nMS).getChild("PhaseShift");
    m_MultiSin_phase_shift[X] = grains::xml::requireElement(nPhaseShift).requireAttribute<double>("PhiX") * PI / 180.;
    m_MultiSin_phase_shift[Y] = grains::xml::requireElement(nPhaseShift).requireAttribute<double>("PhiY") * PI / 180.;
    m_MultiSin_phase_shift[Z] = grains::xml::requireElement(nPhaseShift).requireAttribute<double>("PhiZ") * PI / 180.;
    grains::xml::OptionalElement nAmp = grains::xml::requireElement(nMS).getChild("Amplitude");
    m_MultiSin_amplitude[X] = grains::xml::requireElement(nAmp).requireAttribute<double>("AX")
    	* m_unit_vitRef[X];
    m_MultiSin_amplitude[Y] = grains::xml::requireElement(nAmp).requireAttribute<double>("AY")
    	* m_unit_vitRef[Y];
    m_MultiSin_amplitude[Z] = grains::xml::requireElement(nAmp).requireAttribute<double>("AZ")
    	* m_unit_vitRef[Z];
    if ( rank == 0 )
    {
      cout << GrainsExec::m_shift12 << "Obstacle name = " << m_ObstacleName 
      	<< endl;      
      cout << GrainsExec::m_shift12 << "Time interval = [" 
      	<< m_tstart << "," << m_tend << "]" << endl;
      cout << GrainsExec::m_shift12 << "Type = " << m_type << endl;
      cout << GrainsExec::m_shift12 << "Amplitude = " << 
      	 m_MultiSin_amplitude[X] << " " << m_MultiSin_amplitude[Y] << " " << 
	 m_MultiSin_amplitude[Z] << endl;
      cout << GrainsExec::m_shift12 << "Period = " << 
      	m_MultiSin_period[X] << " " << m_MultiSin_period[Y] << " " <<
	m_MultiSin_period[Z] << endl;
      cout << GrainsExec::m_shift12 << "Phase shift in rad = " << 
      	m_MultiSin_phase_shift[X] << " " << m_MultiSin_phase_shift[Y] << " " 
	<< m_MultiSin_phase_shift[Z] << endl;	 
    }   
  }    
  else 
  {
    error = 1;
    if ( rank == 0 ) cout << GrainsExec::m_shift6 << 
	"Unknown or missing obstacle imposed velocity node !!" << endl;
  }
}




// ----------------------------------------------------------------------------
// Destructor
ObstacleImposedVelocity::~ObstacleImposedVelocity()
{}




// ----------------------------------------------------------------------------
// Returns obstacle name
string ObstacleImposedVelocity::getObstacleName() const
{
  return ( m_ObstacleName );
}




// ----------------------------------------------------------------------------
// Returns whether the imposed motion is activ over the time interval [td,te] 
// and if it is the sub-interval length within [td,te] when it is actually 
// active
bool ObstacleImposedVelocity::isActif( double const& td, double const& te, 
	double const& dt, double& subinterval ) const 
{
  subinterval = te - td;

  if ( td < m_tstart ) subinterval -= ( m_tstart - td );
  if ( m_tend < te ) subinterval -= ( te - m_tend );

  return ( subinterval > dt * 1.e-5 );
}




// ----------------------------------------------------------------------------
// Returns whether the imposed motion is completed at time t
bool ObstacleImposedVelocity::isCompleted( double t, double dt ) const 
{
  return ( t > m_tend + dt * 1.e-5 );
}




// ----------------------------------------------------------------------------
// Returns the translational velocity at time t 
Vector3 const* ObstacleImposedVelocity::translationalVelocity( double time, 
	double dt, Point3 const& cg )
{
  if ( isCompleted( time, dt ) ) m_translationalVelocity = 0.;
  else
  {
    if ( m_type == "SinTranslation" )
      m_translationalVelocity = m_amplitude * 
    	sin( 2. * PI * ( time - m_tstart ) / m_period + m_Sin_phase_shift )
	* m_unit_vitRef ;
    else if ( m_type == "CrenelTranslation" )
    {
      m_translationalVelocity = m_amplitude * 
      	( int( ( time - m_tstart ) / m_period ) % 2 == 0 ? 1 : -1 )
	* m_unit_vitRef;  
    }
    else if ( m_type == "SSDCyclicTranslation" )
    {
      m_previous_translationalVelocity = m_translationalVelocity;
      if ( m_stress > m_stress_max && m_stress > m_previous_stress &&
      	m_translationalVelocity[m_stressIndices.first] > 0. )
	m_translationalVelocity[m_stressIndices.first] = - m_amplitude;
      else if ( m_stress < - m_stress_max && m_stress < m_previous_stress &&
      	m_translationalVelocity[m_stressIndices.first] < 0. )
	m_translationalVelocity[m_stressIndices.first] = m_amplitude;
    }      
    else if ( m_type == "MultiSinTranslation" )
    {
      for (size_t i=0;i<3;++i)
        m_translationalVelocity[i] = m_MultiSin_amplitude[i] * 
		sin( 2. * PI * ( time - m_tstart ) / m_MultiSin_period[i] 
		+ m_MultiSin_phase_shift[i] ) * m_unit_vitRef[i] ;	    	
    }
    else if ( m_type == "ConstantRotation" && !m_rotationCenterIsCenterOfMass )
      m_translationalVelocity = m_angularVelocity ^ ( cg - m_rotationCenter );
  }

  return ( &m_translationalVelocity );
}




// ----------------------------------------------------------------------------
// Returns the angular velocity at time t 
Vector3 const* ObstacleImposedVelocity::angularVelocity( double time, 
	double dt )
{    
  if ( isCompleted( time, dt ) ) m_angularVelocity = 0.;  
  return ( &m_angularVelocity );
}




// ----------------------------------------------------------------------------
// Returns the translational motion over dt at time t
Vector3 ObstacleImposedVelocity::translationalMotion( double time, 
	double dt, double const& subinterval, Point3 const& cg ) 
{
  // We integrate analytically the velocity over the subinterval of 
  // [time-dt,time] where the imposed motion is active
  Vector3 translation;
  double ti = time - dt < m_tstart ? m_tstart : time - dt;
  double te = time > m_tend ? m_tend : time;   
  
  if ( m_type == "ConstantTranslation" ) 
    translation = m_translationalVelocity * subinterval;
  else if ( m_type == "CrenelTranslation" )
    translation = m_amplitude * 
      	( int( ( 0.5 * ( ti + te ) - m_tstart ) / m_period ) % 2 == 0 ? 1 : -1 )
	* subinterval * m_unit_vitRef;
  else if ( m_type == "SSDCyclicTranslation" )
  {
    translation = 0.5 * ( m_translationalVelocity
    	+ m_previous_translationalVelocity ) * subinterval ;    
  }
  else if ( m_type == "SinTranslation" )
    translation = ( m_amplitude * m_period / ( 2. * PI ) )
    	* ( cos( 2. * PI * ( ti - m_tstart ) / m_period 
		+ m_Sin_phase_shift ) 
	- cos( 2. * PI * ( te - m_tstart ) / m_period 
		+ m_Sin_phase_shift ) ) * m_unit_vitRef ;
  else if ( m_type == "MultiSinTranslation" )
  {
    for (size_t i=0;i<3;++i)
      translation[i] = 
      	( m_MultiSin_amplitude[i] * m_MultiSin_period[i] / ( 2. * PI ) )
    	* ( cos( 2. * PI * ( ti - m_tstart ) / m_MultiSin_period[i] 
		+ m_MultiSin_phase_shift[i] ) 
	- cos( 2. * PI * ( te - m_tstart ) / m_MultiSin_period[i] 
		+ m_MultiSin_phase_shift[i] ) ) * m_unit_vitRef[i] ;	    	
  }
  else if ( m_type == "ConstantRotation" && !m_rotationCenterIsCenterOfMass )
  {
    Vector3 rota = angularMotion( time, dt, subinterval );
    double d = Norm(rota);
    Quaternion q;

    if ( d != 0. ) 
    {
      Vector3 vect = ( sin( d / 2. ) / d ) * rota;
      q = Quaternion( vect, cos( d / 2. ) );
    }
    else 
      q = Quaternion( 0., 0., 0., 1. );
    translation = q.rotateVector( cg - m_rotationCenter ) 
    	- ( cg - m_rotationCenter );    
  }
    
  return ( translation );
}  




// ----------------------------------------------------------------------------
// Returns the angular motion over dt at time t 
Vector3 ObstacleImposedVelocity::angularMotion( double time, double dt, 
	double const& subinterval )
{
  // We only consider contant angular velocity so far
  return ( m_angularVelocity * subinterval );
}  




// ----------------------------------------------------------------------------
// Operator == based on the object address
bool ObstacleImposedVelocity::operator == ( 
	ObstacleImposedVelocity const& other ) const
{
  return ( this == &other );
}




// ----------------------------------------------------------------------------
// Operator < based on the start time of the imposed motion.
// Returns true if c0.tdebut < c1.tdebut
bool operator < ( ObstacleImposedVelocity const& c0,
	ObstacleImposedVelocity const& c1 )
{
  return ( c0.m_tstart < c1.m_tstart );
}




// ----------------------------------------------------------------------------
// Debug
void ObstacleImposedVelocity::debug( char* c )
{
  cout << m_tstart << '-' << m_tend << '\t';
}




// ----------------------------------------------------------------------------
// Returns the imposed motion type
string ObstacleImposedVelocity::getType() const
{
  return ( m_type );
}



 
// ----------------------------------------------------------------------------
// Updates imposed velocity based on a stress criterion (for cyclic shearing) 
void ObstacleImposedVelocity::updateImposedVelocity( LinkedCell const* LC )
{
  if ( m_type == "SSDCyclicTranslation" )
  {
    double tau = LC->getStressTensorComponent( m_stressIndices.first,
    	m_stressIndices.second );

    // Check whether the stress component was updated over the last time step
    if ( fabs( tau - m_stress ) > EPSILON )
    {
      // Update the stress
      m_previous_stress = m_stress;
      m_stress = tau;
    }
  }
}
