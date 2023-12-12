#include <DS_Disc.hh>
#include <FS_RigidBody.hh>
#include <FS_Disc.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
using std::endl;


//---------------------------------------------------------------------------
DS_Disc:: DS_Disc()
//---------------------------------------------------------------------------
  : DS_RigidBody()
{
  MAC_LABEL( "DS_Disc:: DS_Disc" ) ;

}




//---------------------------------------------------------------------------
DS_Disc:: DS_Disc( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
  : DS_RigidBody( pgrb )
{
  MAC_LABEL( "DS_RigidBody:: DS_RigidBody" ) ;

}




//---------------------------------------------------------------------------
DS_Disc:: ~DS_Disc()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: ~DS_Disc" ) ;

}




//---------------------------------------------------------------------------
void DS_Disc:: update()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: update" ) ;

}




//---------------------------------------------------------------------------
void DS_Disc:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;

  out << space << "Geometric rigid body features" << endl;
  m_geometric_rigid_body->display( out, indent_width + 3 );
  out << space << "Direction splitting specific features" << endl;
  out << space << three << "None so far" << endl;

}




//---------------------------------------------------------------------------
void DS_Disc:: compute_rigid_body_halozone( double const& dx )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: compute_rigid_body_halozone" ) ;

  geomVector const* pgc = dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre();

  double r_equi = get_circumscribed_radius() + dx;

  geomVector delta(r_equi, r_equi, r_equi);

  m_halo_zone[0]->operator=(*pgc);
  m_halo_zone[1]->operator=(*pgc);

  m_halo_zone[0]->operator-=(delta);
  m_halo_zone[1]->operator+=(delta);

}




//---------------------------------------------------------------------------
bool DS_Disc:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: isIn(pt)" ) ;

  return ( m_geometric_rigid_body->isIn( pt ) );

}




//---------------------------------------------------------------------------
bool DS_Disc:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: isIn(x,y,z)" ) ;

  return ( m_geometric_rigid_body->isIn( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_Disc:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: level_set_value(pt)" ) ;

  return ( m_geometric_rigid_body->level_set_value( pt ) );

}




//---------------------------------------------------------------------------
double DS_Disc:: level_set_value( double const& x
                                     , double const& y
                                     , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: level_set_value(x,y,z)" ) ;

  return ( m_geometric_rigid_body->level_set_value( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_Disc:: get_distanceTo( geomVector const& source,
                                      geomVector const& rayDir,
                                      double const& delta ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: get_distanceTo" ) ;

  return (m_geometric_rigid_body->distanceTo(source, rayDir, delta));
  // return (m_geometric_rigid_body->analytical_distanceTo(source, rayDir));
}




//---------------------------------------------------------------------------
geomVector DS_Disc:: get_rigid_body_velocity( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: rigid_body_velocity(pt)" ) ;

  return (m_geometric_rigid_body->rigid_body_velocity(pt));

}




//---------------------------------------------------------------------------
geomVector DS_Disc:: get_rigid_body_angular_velocity( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: rigid_body_angular_velocity()" ) ;

  return (m_geometric_rigid_body->rigid_body_angular_velocity());

}




//---------------------------------------------------------------------------
std::tuple<double,double,double> DS_Disc:: get_mass_and_density_and_moi() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: get_mass_and_density()" ) ;

  return ( m_geometric_rigid_body->get_mass_and_density_and_moi() );

}




//---------------------------------------------------------------------------
double DS_Disc:: get_circumscribed_radius( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: get_circumscribed_radius()" ) ;

  return (m_geometric_rigid_body->get_circumscribed_radius());

}




//---------------------------------------------------------------------------
geomVector const* DS_Disc:: get_ptr_to_gravity_centre( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: get_ptr_to_gravity_centre( )" ) ;

  return (dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre());

}




//---------------------------------------------------------------------------
void DS_Disc:: update_RB_position_and_velocity(geomVector const& pos,
                                                    geomVector const& vel,
                                                    geomVector const& ang_vel,
                                   vector<geomVector> const& periodic_directions,
                                   double const& time_step)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: update_RB_position_and_velocity" ) ;

  return (m_geometric_rigid_body->update_RB_position_and_velocity(pos,vel
                                    ,ang_vel,periodic_directions, time_step));

}




//---------------------------------------------------------------------------
void DS_Disc:: update_additional_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: update_additional_parameters" ) ;

  m_geometric_rigid_body->update_additional_parameters();

}




//---------------------------------------------------------------------------
void DS_Disc:: compute_surface_points( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: compute_surface_points" ) ;

  // Pointers to location and additional parameters
  struct FS_Disc_Additional_Param const* pagp =
   dynamic_cast<FS_Disc*>(m_geometric_rigid_body)
      ->get_ptr_FS_Disc_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
                            ->get_ptr_to_gravity_centre();

  size_t Npoints = m_surface_area.size();
  double d_theta = 2.*MAC::pi()/((double)Npoints);
  double theta = 0.5*d_theta;

  for (size_t i = 0; i < Npoints; i++) {
     theta = theta + d_theta;

     geomVector point( pagp->radius*MAC::cos(theta)
                     , pagp->radius*MAC::sin(theta)
                     , 0. );

     m_surface_points[i]->operator=(point);
     m_surface_area[i]->operator()(0) = pagp->radius*d_theta;

     // Create surface normal vectors
     m_surface_normal[i]->operator=(point);
  }

  // Translate and rotate
  for (size_t i = 0; i < m_surface_area.size(); i++) {
//     m_geometric_rigid_body->rotate(m_surface_points[i]);
//     m_geometric_rigid_body->rotate(m_surface_normal[i]);
     m_surface_points[i]->translate(*pgc);
  }


}




//---------------------------------------------------------------------------
void DS_Disc:: compute_number_of_surface_variables(
                                          double const& surface_cell_scale
                                        , double const& dx)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_Disc:: compute_number_of_surface_variables" ) ;

  struct FS_Disc_Additional_Param const* pagp =
   dynamic_cast<FS_Disc*>(m_geometric_rigid_body)
      ->get_ptr_FS_Disc_Additional_Param();

  size_t temp = (size_t) ((1./surface_cell_scale)
               *(2.*MAC::pi()*pagp->radius)
               /(dx));

  // Getting the nearest even number
  Ntot = (size_t) (round((double)temp * 0.5) * 2.);

}
