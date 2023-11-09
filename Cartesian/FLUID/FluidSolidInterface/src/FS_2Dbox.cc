#include <FS_2Dbox.hh>
#include <math.h>
using std::endl;
#define THRESHOLD 1.e-14


//---------------------------------------------------------------------------
FS_2Dbox:: FS_2Dbox()
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_2Dbox:: FS_2Dbox" ) ;

  m_space_dimension = 2;
  m_shape_type = GEOM_2DBOX;

}




//---------------------------------------------------------------------------
FS_2Dbox:: FS_2Dbox( istream& in, size_t& id_ )
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_2Dbox:: FS_2Dbox" ) ;

  // Default parameter
  m_space_dimension = 2;
  m_Id = id_;
  m_shape_type = GEOM_2DBOX;

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
FS_2Dbox:: ~FS_2Dbox()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dbox:: ~FS_2Dbox" ) ;

}




//---------------------------------------------------------------------------
void FS_2Dbox:: update( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dbox:: update" ) ;

  // Set the rigid body features from the input stream
  set( in );

}




//---------------------------------------------------------------------------
void FS_2Dbox:: set( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dbox:: set" ) ;

  size_t ncorners, i, nper ;
  geomVector node(3);

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

  // Read the corners and faces numbering of the 2D box

  // Build the polyhedron corners
  if ( m_agp_2dbox.corners.empty() ) {
    m_agp_2dbox.corners.reserve( ncorners );
    for (i = 0; i < ncorners; i++) {
      in >> node;
      m_agp_2dbox.corners.push_back(node);
    }
  } else {
    for (i = 0; i < ncorners; i++) in >> m_agp_2dbox.corners[i];
  }

  // Build the reference polyhedron corners
  if ( m_agp_2dbox.ref_corners.empty() ) {
    m_agp_2dbox.ref_corners.reserve( ncorners );
    for (i = 0; i < ncorners; i++) {
      m_agp_2dbox.ref_corners.push_back(node);
    }
  }

  compute_reverseTransformationOfCorners( );

  // Set volume
  m_volume = m_mass / m_density ;

  // display(std::cout, '\t');
}




//---------------------------------------------------------------------------
void FS_2Dbox:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dbox:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;
  out << space << "Shape type = " <<
  	FS_RigidBody::GEOMETRICSHAPE_name[m_shape_type] << endl;
  out << space << "Specific attributes" << endl;
  for (int i = 0; i < (int) m_agp_2dbox.corners.size(); i++)
     out << space << three << i+1 << "st Corner = " << m_agp_2dbox.corners[i](0)
                                           << three << m_agp_2dbox.corners[i](1)
                                           << three << m_agp_2dbox.corners[i](2)
                                           << endl;

  out << space << "General attributes" << endl;
  display_general( out, indent_width + 3 );

}




//---------------------------------------------------------------------------
bool FS_2Dbox:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dbox:: isIn(pt)" ) ;

  // Inspired from the following website:
  // www.jeffreythompson.org/collision-detection/line-line.php

  bool b_isIn = true;

  geomVector p1(pt);
  geomVector p2(m_gravity_center);

  for (int i = 0; i < (int)m_agp_2dbox.corners.size() && b_isIn; i++) {

    geomVector p3(3), p4(3);
    if (i < (int)m_agp_2dbox.corners.size() - 1) {
      p3 = m_agp_2dbox.corners[i];
      p4 = m_agp_2dbox.corners[i+1];
    } else {
      p3 = m_agp_2dbox.corners[i];
      p4 = m_agp_2dbox.corners[0];
    }

    double t = ((p1(0) - p3(0)) * (p3(1) - p4(1)) - (p1(1) - p3(1)) * (p3(0) - p4(0))) 
             / ((p1(0) - p2(0)) * (p3(1) - p4(1)) - (p1(1) - p2(1)) * (p3(0) - p4(0)));

    double u = ((p1(0) - p3(0)) * (p1(1) - p2(1)) - (p1(1) - p3(1)) * (p1(0) - p2(0))) 
             / ((p1(0) - p2(0)) * (p3(1) - p4(1)) - (p1(1) - p2(1)) * (p3(0) - p4(0)));

    if ((t >= 0.) && (t <= 1.) && (u >= 0.) && (u <= 1.)) b_isIn = false;
  }

  return ( b_isIn );

}




//---------------------------------------------------------------------------
bool FS_2Dbox:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dbox:: isIn(x,y,z)" ) ;

  geomVector pt( x, y, z );

  return ( isIn(pt) );

}




//---------------------------------------------------------------------------
double FS_2Dbox:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dbox:: level_set_value(pt)" ) ;

  double value = 1.;

  if (isIn(pt)) value = -1.;



  return(value);

}




//---------------------------------------------------------------------------
double FS_2Dbox:: level_set_value( double const& x
                                 , double const& y
                                 , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dbox:: level_set_value(x,y,z)" ) ;

  geomVector pt( x, y, z );

  return (level_set_value(pt));
}




//---------------------------------------------------------------------------
double FS_2Dbox::analytical_distanceTo(geomVector const &source,
                                       geomVector const &rayDir) const
//---------------------------------------------------------------------------
{
  MAC_LABEL("FS_2Dbox:: analytical_distanceTo");

  geomVector p1(source);
  geomVector p2(source + m_circumscribed_radius*rayDir);
  vector<double> distance;

  for (int i = 0; i < (int)m_agp_2dbox.corners.size(); i++)
  {

    geomVector p3(3), p4(3);
    if (i < (int)m_agp_2dbox.corners.size() - 1)
    {
      p3 = m_agp_2dbox.corners[i];
      p4 = m_agp_2dbox.corners[i + 1];
    }
    else
    {
      p3 = m_agp_2dbox.corners[i];
      p4 = m_agp_2dbox.corners[0];
    }

    double t = ((p1(0) - p3(0)) * (p3(1) - p4(1)) - (p1(1) - p3(1)) * (p3(0) - p4(0))) 
             / ((p1(0) - p2(0)) * (p3(1) - p4(1)) - (p1(1) - p2(1)) * (p3(0) - p4(0)));

    double u = ((p1(0) - p3(0)) * (p1(1) - p2(1)) - (p1(1) - p3(1)) * (p1(0) - p2(0))) 
             / ((p1(0) - p2(0)) * (p3(1) - p4(1)) - (p1(1) - p2(1)) * (p3(0) - p4(0)));

    if ((t >= 0.) && (t <= 1.) && (u >= 0.) && (u <= 1.)) {
      distance.push_back(t * (p1 - p2).calcNorm());
    }
  }

  if (distance.size() == 0) {
    return (0.);
  } else if (distance.size() == 1) {
    return (distance[0]);
  } else if (distance.size() == 2) {
    return (MAC::min(distance[0], distance[1]));
  }

    return (0.);
}

//---------------------------------------------------------------------------
struct FS_2Dbox_Additional_Param const* FS_2Dbox::
	get_ptr_FS_2Dbox_Additional_Param() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dbox:: get_ptr_FS_2Dbox_Additional_Param()" ) ;

  return( &m_agp_2dbox );

}




//---------------------------------------------------------------------------
void FS_2Dbox::compute_reverseTransformationOfCorners( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dbox:: compute_reverseTransformationOfCorners()" ) ;

  for (int i = 0; i < (int) m_agp_2dbox.ref_corners.size(); i++) {
     geomVector pt(m_agp_2dbox.corners[i]);

     m_agp_2dbox.ref_corners[i](0) = (pt(0) - m_gravity_center(0))*m_rotation_matrix[0][0]
                                   + (pt(1) - m_gravity_center(1))*m_rotation_matrix[1][0]
                                   + (pt(2) - m_gravity_center(2))*m_rotation_matrix[2][0];
     m_agp_2dbox.ref_corners[i](1) = (pt(0) - m_gravity_center(0))*m_rotation_matrix[0][1]
                                   + (pt(1) - m_gravity_center(1))*m_rotation_matrix[1][1]
                                   + (pt(2) - m_gravity_center(2))*m_rotation_matrix[2][1];
     m_agp_2dbox.ref_corners[i](2) = (pt(0) - m_gravity_center(0))*m_rotation_matrix[0][2]
                                   + (pt(1) - m_gravity_center(1))*m_rotation_matrix[1][2]
                                   + (pt(2) - m_gravity_center(2))*m_rotation_matrix[2][2];
  }

}




//---------------------------------------------------------------------------
void FS_2Dbox::compute_TransformationOfCorners( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dbox:: compute_TransformationOfCorners()" ) ;

  for (int i = 0; i < (int) m_agp_2dbox.corners.size(); i++) {
     geomVector pt(m_agp_2dbox.ref_corners[i]);

     m_agp_2dbox.corners[i](0) = pt(0)*m_rotation_matrix[0][0]
                               + pt(1)*m_rotation_matrix[0][1]
                               + pt(2)*m_rotation_matrix[0][2]
                               + m_gravity_center(0);
     m_agp_2dbox.corners[i](1) = pt(0)*m_rotation_matrix[1][0]
                               + pt(1)*m_rotation_matrix[1][1]
                               + pt(2)*m_rotation_matrix[1][2]
                               + m_gravity_center(1);
     m_agp_2dbox.corners[i](2) = pt(0)*m_rotation_matrix[2][0]
                               + pt(1)*m_rotation_matrix[2][1]
                               + pt(2)*m_rotation_matrix[2][2]
                               + m_gravity_center(2);
  }

}




//---------------------------------------------------------------------------
void FS_2Dbox::update_additional_parameters( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_2Dbox:: update_additional_parameters( )" ) ;

  compute_TransformationOfCorners();
  // display (std::cout, '\t');

}
