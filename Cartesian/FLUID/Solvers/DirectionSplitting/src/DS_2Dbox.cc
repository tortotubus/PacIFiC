#include <DS_2Dbox.hh>
#include <FS_RigidBody.hh>
#include <FS_2Dbox.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
using std::endl;


//---------------------------------------------------------------------------
DS_2Dbox:: DS_2Dbox()
//---------------------------------------------------------------------------
  : DS_RigidBody()
{
  MAC_LABEL( "DS_2Dbox:: DS_2Dbox" ) ;

}




//---------------------------------------------------------------------------
DS_2Dbox:: DS_2Dbox( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
  : DS_RigidBody( pgrb )
{
  MAC_LABEL( "DS_RigidBody:: DS_RigidBody" ) ;

}




//---------------------------------------------------------------------------
DS_2Dbox:: ~DS_2Dbox()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: ~DS_2Dbox" ) ;

}




//---------------------------------------------------------------------------
void DS_2Dbox:: update()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: update" ) ;

}




//---------------------------------------------------------------------------
void DS_2Dbox:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;

  out << space << "Geometric rigid body features" << endl;
  m_geometric_rigid_body->display( out, indent_width + 3 );
  out << space << "Direction splitting specific features" << endl;
  out << space << three << "None so far" << endl;

}




//---------------------------------------------------------------------------
void DS_2Dbox:: compute_rigid_body_halozone( double const& dx )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: compute_rigid_body_halozone" ) ;

  geomVector const* pgc = dynamic_cast<FS_2Dbox*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre();

  geomVector delta(3);

  double r_equi = get_circumscribed_radius() + dx;

  delta(0) = r_equi;
  delta(1) = r_equi;
  delta(2) = r_equi;

  m_halo_zone[0]->operator=(*pgc);
  m_halo_zone[1]->operator=(*pgc);

  m_halo_zone[0]->operator-=(delta);
  m_halo_zone[1]->operator+=(delta);

}




//---------------------------------------------------------------------------
bool DS_2Dbox:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: isIn(pt)" ) ;

  return ( m_geometric_rigid_body->isIn( pt ) );

}




//---------------------------------------------------------------------------
bool DS_2Dbox:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: isIn(x,y,z)" ) ;

  return ( m_geometric_rigid_body->isIn( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_2Dbox:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: level_set_value(pt)" ) ;

  return ( m_geometric_rigid_body->level_set_value( pt ) );

}




//---------------------------------------------------------------------------
double DS_2Dbox:: level_set_value( double const& x
                                     , double const& y
                                     , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: level_set_value(x,y,z)" ) ;

  return ( m_geometric_rigid_body->level_set_value( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_2Dbox:: get_distanceTo( geomVector const& source,
                                      geomVector const& rayDir,
                                      double const& delta ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: get_distanceTo" ) ;

  // return (m_geometric_rigid_body->distanceTo(source, rayDir, delta));
  return (m_geometric_rigid_body->analytical_distanceTo(source, rayDir));
}




//---------------------------------------------------------------------------
geomVector DS_2Dbox:: get_rigid_body_velocity( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: rigid_body_velocity(pt)" ) ;

  return (m_geometric_rigid_body->rigid_body_velocity(pt));

}




//---------------------------------------------------------------------------
geomVector DS_2Dbox:: get_rigid_body_angular_velocity( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: rigid_body_angular_velocity()" ) ;

  return (m_geometric_rigid_body->rigid_body_angular_velocity());

}




//---------------------------------------------------------------------------
std::tuple<double,double,double> DS_2Dbox:: get_mass_and_density_and_moi() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: get_mass_and_density()" ) ;

  return ( m_geometric_rigid_body->get_mass_and_density_and_moi() );

}




//---------------------------------------------------------------------------
double DS_2Dbox:: get_circumscribed_radius( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: get_circumscribed_radius()" ) ;

  return (m_geometric_rigid_body->get_circumscribed_radius());

}




//---------------------------------------------------------------------------
geomVector const* DS_2Dbox:: get_ptr_to_gravity_centre( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: get_ptr_to_gravity_centre( )" ) ;

  return (dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre());

}




//---------------------------------------------------------------------------
void DS_2Dbox:: update_RB_position_and_velocity(geomVector const& pos,
                                                    geomVector const& vel,
                                                    geomVector const& ang_vel,
                                vector<geomVector> const& periodic_directions,
                                   double const& time_step)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: update_RB_position_and_velocity" ) ;

  return (m_geometric_rigid_body->update_RB_position_and_velocity(pos,vel
                                                                  ,ang_vel
                                                         ,periodic_directions
                                                         , time_step));

}




//---------------------------------------------------------------------------
void DS_2Dbox:: update_additional_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: update_additional_parameters" ) ;

  m_geometric_rigid_body->update_additional_parameters();

}




//---------------------------------------------------------------------------
void DS_2Dbox:: compute_surface_points(  )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: compute_surface_points" ) ;

  doubleVector di(2,0.);

  di(0) = (box_max[0] - box_min[0])/(double)Npoints[0];
  di(1) = (box_max[1] - box_min[1])/(double)Npoints[1];

  size_t cntr = 0;

  // x face
  for (size_t i = 0; i < Npoints[1]; i++) {
      geomVector point_min (box_min[0]
                          , box_min[1] + (0.5+(double)i)*di(1)
                          , 0. );
      geomVector point_max (box_max[0]
                          , box_min[1] + (0.5+(double)i)*di(1)
                          , 0. );
      m_surface_points[cntr]->operator=(point_min);
      m_surface_points[Npoints[1] + cntr]->operator=(point_max);

      m_surface_area[cntr]->operator()(0) = di(1);
      m_surface_area[Npoints[1]+cntr]->operator()(0) = di(1);

      geomVector normal(-1., 0., 0.);
      m_surface_normal[cntr]->operator=(normal);
      m_surface_normal[Npoints[1] + cntr]->operator=(-1.*normal);

      cntr++;
  }

  cntr = cntr + Npoints[1];

  // y face
  for (size_t i = 0; i < Npoints[0]; i++) {
      geomVector point_min (box_min[0] + (0.5+(double)i)*di(0)
                          , box_min[1]
                          , 0. );
      geomVector point_max (box_min[0] + (0.5+(double)i)*di(0)
                          , box_max[1]
                          , 0. );
      m_surface_points[cntr]->operator=(point_min);
      m_surface_points[Npoints[0] + cntr]->operator=(point_max);

      m_surface_area[cntr]->operator()(0) = di(0);
      m_surface_area[Npoints[0] + cntr]->operator()(0) = di(0);

      geomVector normal(0., -1., 0.);
      m_surface_normal[cntr]->operator=(normal);
      m_surface_normal[Npoints[0] + cntr]->operator=(-1.*normal);

      cntr++;
  }

  // Translate and rotate
  for (size_t i = 0; i < m_surface_area.size(); i++) {
     m_geometric_rigid_body->rotate(m_surface_points[i]);
     m_geometric_rigid_body->rotate(m_surface_normal[i]);
     m_geometric_rigid_body->translate(m_surface_points[i]);
  }

}




//---------------------------------------------------------------------------
void DS_2Dbox:: compute_number_of_surface_variables(
                                          double const& surface_cell_scale
                                        , double const& dx )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2Dbox:: compute_number_of_surface_variables" ) ;

  struct FS_2Dbox_Additional_Param const* pagp =
                        dynamic_cast<FS_2Dbox*>(m_geometric_rigid_body)
                           ->get_ptr_FS_2Dbox_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_2Dbox*>(m_geometric_rigid_body)
                            ->get_ptr_to_gravity_centre();

  box_min.assign(2,1.e14);
  box_max.assign(2,-1.e14);

  for (int i = 0; i < (int) pagp->ref_corners.size(); i++) {
     for (int dir = 0; dir < 2; dir++) {
        if (pagp->ref_corners[i](dir) < box_min[dir])
           box_min[dir] = pagp->ref_corners[i](dir);

        if (pagp->ref_corners[i](dir) > box_max[dir])
           box_max[dir] = pagp->ref_corners[i](dir);
     }
  }

  doubleVector delta(2,0.);

  delta(0) = (box_max[0]-box_min[0]);
  delta(1) = (box_max[1]-box_min[1]);

  double scale = 1. / sqrt(surface_cell_scale) / dx;

  if (Npoints.empty()) Npoints.reserve(2);

  Npoints[0] = (size_t) round(scale * delta(0) );
  Npoints[1] = (size_t) round(scale * delta(1) );

  Ntot = 2 * (Npoints[0] + Npoints[1]);

}
