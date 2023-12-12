#include <FS_Disc.hh>
#include <MAC.hh>
#include <math.h>
using std::endl;


//---------------------------------------------------------------------------
FS_Disc:: FS_Disc()
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_Disc:: FS_Disc" ) ;

  m_space_dimension = 2;
  m_shape_type = GEOM_DISC;
  m_agp_Disc.radius = 0.;
}




//---------------------------------------------------------------------------
FS_Disc:: FS_Disc( istream& in, size_t& id_ )
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_Disc:: FS_Disc" ) ;

  // Default parameter
  m_space_dimension = 2;
  m_Id = id_;
  m_shape_type = GEOM_DISC;

  // Resize parameters
  m_gravity_center.resize(3);
  m_translational_velocity.resize(3);
  m_angular_velocity.resize(3);
  m_hydro_force.resize(3);
  m_hydro_torque.resize(3);
  m_heat_flux.resize(3);

  // Set the rigid body features from the input stream
  set( in );

}




//---------------------------------------------------------------------------
FS_Disc:: ~FS_Disc()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Disc:: ~FS_Disc" ) ;

}




//---------------------------------------------------------------------------
void FS_Disc:: update( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Disc:: update" ) ;

  // Set the rigid body features from the input stream
  set( in );

}




//---------------------------------------------------------------------------
void FS_Disc:: set( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Disc:: set" ) ;

  size_t ncorners, i, nfaces, nper ;

  // Read the input stream
  in >> m_type >> m_translational_velocity >>
  	m_angular_velocity >> m_density >> m_mass;
  in >> m_inertia[0][0];
  in >> m_inertia[0][1];
  in >> m_inertia[0][2];
  in >> m_inertia[1][1];
  in >> m_inertia[1][2];
  in >> m_inertia[2][2];
  m_inertia[1][0] = m_inertia[0][1];
  m_inertia[2][0] = m_inertia[0][2];
  m_inertia[2][1] = m_inertia[1][2];

  in >> m_rotation_matrix[0][0];
  in >> m_rotation_matrix[0][1];
  in >> m_rotation_matrix[0][2];
  in >> m_rotation_matrix[1][0];
  in >> m_rotation_matrix[1][1];
  in >> m_rotation_matrix[1][2];
  in >> m_rotation_matrix[2][0];
  in >> m_rotation_matrix[2][1];
  in >> m_rotation_matrix[2][2];

  in >> m_gravity_center;

  if ( m_periodic_directions )
  {
    delete m_periodic_directions;
    m_periodic_directions = NULL;
  }

  if ( m_type == "PP" )
  {
    geomVector pv( 3 );
    in >> nper;
    m_periodic_directions = new vector<geomVector>( nper, pv ) ;
    for (i=0;i<nper;++i)
      in >> (*m_periodic_directions)[i](0) >> (*m_periodic_directions)[i](1)
      	>> (*m_periodic_directions)[i](2);
  }
  in >> m_circumscribed_radius >> ncorners;
  in >> m_gravity_center;

  // Set volume
  m_volume = m_mass / m_density ;

  // Set radius
  m_agp_Disc.radius = m_circumscribed_radius;

  // Force moi of a 2D disk
  // m_inertia[0][0] = (1./2.) * m_mass
  //                           * m_circumscribed_radius * m_circumscribed_radius;
  // m_inertia[1][1] = m_inertia[0][0];
  // m_inertia[2][2] = m_inertia[0][0];

}




//---------------------------------------------------------------------------
void FS_Disc:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Disc:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;
  out << space << "Shape type = " <<
  	FS_RigidBody::GEOMETRICSHAPE_name[m_shape_type] << endl;
  out << space << "Specific attributes" << endl;
  out << space << three << "Radius = " << m_agp_Disc.radius << endl;
  out << space << "General attributes" << endl;

  display_general( out, indent_width + 3 );

}




//---------------------------------------------------------------------------
bool FS_Disc:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Disc:: isIn(pt)" ) ;

  bool status = ( m_gravity_center.calcDist( pt ) <= m_agp_Disc.radius);

  if (m_periodic_directions) {
     for (size_t i = 0; i < m_periodic_directions->size(); ++i) {
        if (status) break;
        status = (m_gravity_center + (*m_periodic_directions)[i]).calcDist( pt )
                  <= m_agp_Disc.radius;
     }
  }

  return (status);

  // return ( m_gravity_center.calcDist( pt ) <= m_agp_Disc.radius );

}




//---------------------------------------------------------------------------
bool FS_Disc:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Disc:: isIn(x,y,z)" ) ;

  geomVector pt(x, y, z);

  bool status = ( m_gravity_center.calcDist( pt ) <= m_agp_Disc.radius);

  if (m_periodic_directions) {
     for (size_t i = 0; i < m_periodic_directions->size(); ++i) {
        if (status) break;
        status = (m_gravity_center + (*m_periodic_directions)[i]).calcDist( pt )
                  <= m_agp_Disc.radius;
     }
  }

  return (status);

  // return ( m_gravity_center.calcDist( x, y, z ) <= m_agp_Disc.radius );

}




//---------------------------------------------------------------------------
double FS_Disc:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Disc:: level_set_value(pt)" ) ;

  double value = ( m_gravity_center.calcDist( pt ) - m_agp_Disc.radius);

  if (m_periodic_directions) {
     for (size_t i = 0; i < m_periodic_directions->size(); ++i) {
        double temp = (m_gravity_center +
                     (*m_periodic_directions)[i]).calcDist( pt )
                     - m_agp_Disc.radius;
        value = MAC::min(temp,value);
     }
  }

  return (value);

  // return ( m_gravity_center.calcDist( pt ) - m_agp_Disc.radius );

}




//---------------------------------------------------------------------------
double FS_Disc:: level_set_value( double const& x
                                , double const& y
                                , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_Disc:: level_set_value(x,y,z)" ) ;

  geomVector pt(x, y, z);

  double value = ( m_gravity_center.calcDist( pt ) - m_agp_Disc.radius);

  if (m_periodic_directions) {
     for (size_t i = 0; i < m_periodic_directions->size(); ++i) {
        double temp = (m_gravity_center +
                     (*m_periodic_directions)[i]).calcDist( pt )
                     - m_agp_Disc.radius;
        value = MAC::min(temp,value);
     }
  }

  return (value);

  // return ( m_gravity_center.calcDist( x, y, z ) - m_agp_Disc.radius );

}

//---------------------------------------------------------------------------
double FS_Disc::analytical_distanceTo(geomVector const &source,
                                      geomVector const &rayDir) const
//---------------------------------------------------------------------------
{
  MAC_LABEL("FS_Disc:: analytical_distanceTo");

  double value = analytical_distanceTo_nonPeriodic(m_gravity_center,
                                                  source,
                                                  rayDir);

  if (m_periodic_directions) {
    for (size_t i = 0; i < m_periodic_directions->size(); ++i) {
      double temp = analytical_distanceTo_nonPeriodic
                      (m_gravity_center + (*m_periodic_directions)[i]
                      , source, rayDir);
      value = MAC::min(temp, value);
    }
  }

  return (value);
}


//---------------------------------------------------------------------------
double FS_Disc::analytical_distanceTo_nonPeriodic(geomVector const &gc,
                                                  geomVector const &source,
                                                  geomVector const &rayDir) const
//---------------------------------------------------------------------------
{
  MAC_LABEL("FS_Disc:: analytical_distanceTo");

  double value = 0.;

  // Intersection distance in x direction
  if (rayDir(0) != 0)
  {
    double dx = 0.;
    double dy = source(1) - gc(1);
    double x1 = 0., x2 = 0.;
    if (MAC::abs(dy) < m_agp_Disc.radius)
    {
      dx = MAC::sqrt(MAC::pow(m_agp_Disc.radius, 2) - MAC::pow(dy, 2));
      x1 = gc(0) + dx;
      x2 = gc(0) - dx;
    }

    value = MAC::min(MAC::abs(source(0) - x1), MAC::abs(source(0) - x2));
  }
  else if (rayDir(1) != 0)
  {
    double dx = source(0) - gc(0);
    double dy = 0.;
    double x1 = 0., x2 = 0.;
    if (MAC::abs(dx) < m_agp_Disc.radius)
    {
      dy = MAC::sqrt(MAC::pow(m_agp_Disc.radius, 2) - MAC::pow(dx, 2));
      x1 = gc(1) + dy;
      x2 = gc(1) - dy;
    }

    value = MAC::min(MAC::abs(source(1) - x1), MAC::abs(source(1) - x2));
  }

  return (value);
}

//---------------------------------------------------------------------------
struct FS_Disc_Additional_Param const *FS_Disc::
    get_ptr_FS_Disc_Additional_Param() const
//---------------------------------------------------------------------------
{
  MAC_LABEL("FS_Disc:: get_ptr_FS_Disc_Additional_Param");

  return (&m_agp_Disc);
}

//---------------------------------------------------------------------------
void FS_Disc::update_additional_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL("FS_Disc:: update_additional_parameters( )");
}
