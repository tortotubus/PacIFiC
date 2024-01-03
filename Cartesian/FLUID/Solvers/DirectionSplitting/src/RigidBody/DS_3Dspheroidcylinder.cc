#include <DS_3Dspheroidcylinder.hh>
#include <FS_RigidBody.hh>
#include <FS_3Dspheroidcylinder.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
using std::endl;


//---------------------------------------------------------------------------
DS_3Dspheroidcylinder:: DS_3Dspheroidcylinder()
//---------------------------------------------------------------------------
  : DS_RigidBody()
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: DS_3Dspheroidcylinder" ) ;

}




//---------------------------------------------------------------------------
DS_3Dspheroidcylinder:: DS_3Dspheroidcylinder( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
  : DS_RigidBody( pgrb )
{
  MAC_LABEL( "DS_RigidBody:: DS_RigidBody" ) ;

}




//---------------------------------------------------------------------------
DS_3Dspheroidcylinder:: ~DS_3Dspheroidcylinder()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: ~DS_3Dspheroidcylinder" ) ;

}




//---------------------------------------------------------------------------
void DS_3Dspheroidcylinder:: update()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: update" ) ;

}




//---------------------------------------------------------------------------
void DS_3Dspheroidcylinder:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;

  out << space << "Geometric rigid body features" << endl;
  m_geometric_rigid_body->display( out, indent_width + 3 );
  out << space << "Direction splitting specific features" << endl;
  out << space << three << "None so far" << endl;

}




//---------------------------------------------------------------------------
void DS_3Dspheroidcylinder:: compute_rigid_body_halozone( double const& dx )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: compute_rigid_body_halozone" ) ;

  geomVector const* pgc = dynamic_cast<FS_3Dspheroidcylinder*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre();

  double r_equi = get_circumscribed_radius() + dx;

  geomVector delta(r_equi, r_equi, r_equi);

  m_halo_zone[0]->operator=(*pgc);
  m_halo_zone[1]->operator=(*pgc);

  m_halo_zone[0]->operator-=(delta);
  m_halo_zone[1]->operator+=(delta);

}




//---------------------------------------------------------------------------
bool DS_3Dspheroidcylinder:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: isIn(pt)" ) ;

  return ( m_geometric_rigid_body->isIn( pt ) );

}




//---------------------------------------------------------------------------
bool DS_3Dspheroidcylinder:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: isIn(x,y,z)" ) ;

  return ( m_geometric_rigid_body->isIn( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_3Dspheroidcylinder:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: level_set_value(pt)" ) ;

  return ( m_geometric_rigid_body->level_set_value( pt ) );

}




//---------------------------------------------------------------------------
double DS_3Dspheroidcylinder:: level_set_value( double const& x
                                     , double const& y
                                     , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: level_set_value(x,y,z)" ) ;

  return ( m_geometric_rigid_body->level_set_value( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_3Dspheroidcylinder:: get_distanceTo( geomVector const& source,
                                      geomVector const& rayDir,
                                      double const& delta ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: get_distanceTo" ) ;

  // return (m_geometric_rigid_body->distanceTo(source, rayDir, delta));
  return (m_geometric_rigid_body->analytical_distanceTo(source, rayDir));
}




//---------------------------------------------------------------------------
geomVector DS_3Dspheroidcylinder:: get_rigid_body_velocity( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: rigid_body_velocity(pt)" ) ;

  return (m_geometric_rigid_body->rigid_body_velocity(pt));

}




//---------------------------------------------------------------------------
geomVector DS_3Dspheroidcylinder:: get_rigid_body_angular_velocity( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: rigid_body_angular_velocity()" ) ;

  return (m_geometric_rigid_body->rigid_body_angular_velocity());

}




//---------------------------------------------------------------------------
std::tuple<double,double,double> DS_3Dspheroidcylinder:: get_mass_and_density_and_moi() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: get_mass_and_density()" ) ;

  return ( m_geometric_rigid_body->get_mass_and_density_and_moi() );

}




//---------------------------------------------------------------------------
double DS_3Dspheroidcylinder:: get_circumscribed_radius( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: get_circumscribed_radius()" ) ;

  return (m_geometric_rigid_body->get_circumscribed_radius());

}




//---------------------------------------------------------------------------
geomVector const* DS_3Dspheroidcylinder:: get_ptr_to_gravity_centre( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: get_ptr_to_gravity_centre( )" ) ;

  return (dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
                              ->get_ptr_to_gravity_centre());

}




//---------------------------------------------------------------------------
void DS_3Dspheroidcylinder:: update_RB_position_and_velocity(geomVector const& pos,
                                                    geomVector const& vel,
                                                    geomVector const& ang_vel,
                                   vector<geomVector> const& periodic_directions
                                   , double const& time_step)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: update_RB_position_and_velocity" ) ;

  return (m_geometric_rigid_body->update_RB_position_and_velocity(pos,vel
                                                                  ,ang_vel
                                                         ,periodic_directions
                                                         ,time_step));

}




//---------------------------------------------------------------------------
void DS_3Dspheroidcylinder:: update_additional_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: update_additional_parameters" ) ;

  m_geometric_rigid_body->update_additional_parameters();

}




//---------------------------------------------------------------------------
void DS_3Dspheroidcylinder:: compute_surface_points(  )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: compute_surface_points" ) ;

  // Pointers to location and additional parameters
  struct FS_3Dspheroidcylinder_Additional_Param const* pagp =
   dynamic_cast<FS_3Dspheroidcylinder*>(m_geometric_rigid_body)
      ->get_ptr_FS_3Dspheroidcylinder_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_3Dspheroidcylinder*>(m_geometric_rigid_body)
                            ->get_ptr_to_gravity_centre();

  // Reference: Becker and Becker, Computational Geometry 45 (2012) 275-283
  // Estimating the number of rings on the hemisphere, assuming Pmin=3
  // and aspect ratio(ar) as 1

  size_t kmax = (size_t) Nsphere_by_2;

  // Reference: Becker and Becker, Computational Geometry 45 (2012) 275-283
  double eta_temp = MAC::pi() / 2.;
  double k_temp = (double)kmax;
  double Ro_temp = MAC::sqrt(2);
  double Rn_temp = MAC::sqrt(2);
  size_t cntr = 0;

  // Estimating the number of rings on the hemisphere, assuming Pmin=3
  // and aspect ratio(ar) as 1
  while (k_temp > double(Pmin + 2))
  {
    eta_temp = eta_temp - (2. / ar) * MAC::sqrt(MAC::pi() / k_temp) * MAC::sin(eta_temp / 2.);
    Rn_temp = 2. * MAC::sin(eta_temp / 2.);
    k_temp = round(k_temp * MAC::pow(Rn_temp / Ro_temp, 2.));
    Ro_temp = Rn_temp;
    cntr++;
  }

  size_t Nrings = cntr + 1;

  // Summation of total discretized points
  // with increase in number of rings radially
  doubleVector k(Nrings, 0.);
  // Zenithal angle for the sphere
  doubleVector eta(Nrings, 0.);
  // Radius of the rings in lamber projection plane
  doubleVector Rring(Nrings, 0.);

  // Assigning the maximum number of discretized
  // points to the last element of the array
  k(Nrings - 1) = (double)kmax;
  // Zenithal angle for the last must be pi/2.
  eta(Nrings - 1) = MAC::pi() / 2.;
  // Radius of last ring in lamber projection plane
  Rring(Nrings - 1) = MAC::sqrt(2.);

  for (int i = int(Nrings) - 2; i >= 0; --i) {
    eta(i) = eta(i + 1) - (2. / ar) * MAC::sqrt(MAC::pi() / k(i + 1)) * MAC::sin(eta(i + 1) / 2.);
    Rring(i) = 2. * MAC::sin(eta(i) / 2.);
    k(i) = round(k(i + 1) * MAC::pow(Rring(i) / Rring(i + 1), 2.));
    if (i == 0)
      k(0) = (double)Pmin;
  }

  geomVector shift_Hby2(0., 0.5 * pagp->cylinder_height, 0.);

  size_t maxby2 = (size_t)k(Nrings - 1);
  // Calculation for all rings except at the pole
  for (int i = (int)Nrings - 1; i > 0; --i) {
    double Ri = Rring(i);
    Rring(i) = (Rring(i) + Rring(i - 1)) / 2.;
    eta(i) = (eta(i) + eta(i - 1)) / 2.;
    double d_theta = 2. * MAC::pi() / (k(i) - k(i - 1));
    // Initialize theta as 1% of the d_theta to avoid
    // point overlap with mesh gridlines
    double theta = 0.01 * d_theta;

    for (int j = (int)k(i - 1); j < k(i); j++) {
      theta = theta + d_theta;

      geomVector point(pagp->cylinder_radius * MAC::cos(theta) * MAC::sin(eta(i)),
                       pagp->cylinder_radius * MAC::cos(eta(i)),
                       pagp->cylinder_radius * MAC::sin(theta) * MAC::sin(eta(i)));

      geomVector mirror_point(pagp->cylinder_radius * MAC::cos(theta) * MAC::sin(eta(i)),
                     - pagp->cylinder_radius * MAC::cos(eta(i)),
                       pagp->cylinder_radius * MAC::sin(theta) * MAC::sin(eta(i)));

      m_surface_points[j]->operator=(point + shift_Hby2);
      m_surface_area[j]->operator()(0) = pagp->cylinder_radius * pagp->cylinder_radius *
                                         (0.5 * d_theta * (pow(Ri, 2.) - pow(Rring(i - 1), 2.)));
      // For second half of sphere
      m_surface_points[maxby2 + j]->operator=(mirror_point - shift_Hby2);
      m_surface_area[maxby2 + j] = m_surface_area[j];

      // Create surface normal vectors
      m_surface_normal[j]->operator=(point);
      m_surface_normal[maxby2 + j]->operator=(mirror_point);
    }
  }

  // Calculation at the ring on pole (i=0)
  double Ri = Rring(0);
  Rring(0) = Rring(0) / 2.;
  eta(0) = eta(0) / 2.;
  double d_theta = 2. * MAC::pi() / (k(0));

  // Initialize theta as 1% of the d_theta to avoid
  // point overlap with mesh gridlines
  double theta = 0.01 * d_theta;

  if (k(0) > 1) {
    for (int j = 0; j < k(0); j++) {
      theta = theta + d_theta;

      geomVector point(pagp->cylinder_radius * MAC::cos(theta) * MAC::sin(eta(0)),
                       pagp->cylinder_radius * MAC::cos(eta(0)),
                       pagp->cylinder_radius * MAC::sin(theta) * MAC::sin(eta(0)));

      geomVector mirror_point(pagp->cylinder_radius * MAC::cos(theta) * MAC::sin(eta(0)),
                     - pagp->cylinder_radius * MAC::cos(eta(0)),
                       pagp->cylinder_radius * MAC::sin(theta) * MAC::sin(eta(0)));

      m_surface_points[j]->operator=(point + shift_Hby2);
      m_surface_area[j]->operator()(0) =
          pagp->cylinder_radius * pagp->cylinder_radius * 0.5 * d_theta * pow(Ri, 2.);

      // For second half of sphere
      m_surface_points[maxby2 + j]->operator=(mirror_point - shift_Hby2);
      m_surface_area[maxby2 + j] = m_surface_area[j];

      // Create surface normal vectors
      m_surface_normal[j]->operator=(point);
      m_surface_normal[maxby2 + j]->operator=(mirror_point);
    }
  } else {
    geomVector point(0., pagp->cylinder_radius * 1., 0.);
    geomVector mirror_point(0., -pagp->cylinder_radius * 1., 0.);

    m_surface_points[0]->operator=(point + shift_Hby2);
    m_surface_area[0]->operator()(0) =
        pagp->cylinder_radius * pagp->cylinder_radius * 0.5 * d_theta * pow(Ri, 2.);

    //  For second half of sphere
    m_surface_points[maxby2]->operator=(mirror_point - shift_Hby2);
    m_surface_area[maxby2] = m_surface_area[0];

    // Create surface normal vectors
    m_surface_normal[0]->operator=(point);
    m_surface_normal[maxby2]->operator=(mirror_point);
  }

  // Generating one ring of points on cylindrical surface
  // Can be used to calculate stress on whole surface by
  // a constant shift of points

  int pts_1_ring = (int)(k(Nrings-1) - k(Nrings-2));
  int cyl_rings = (int)round((pagp->cylinder_height/pagp->cylinder_radius)/2.
                               /(1-MAC::sqrt(k(Nrings-2)/k(Nrings-1))));
  double cell_area = 2.*pagp->cylinder_radius*pagp->cylinder_height*MAC::pi()
                   / ((double) pts_1_ring*cyl_rings);

  d_theta = 2.*MAC::pi()/pts_1_ring;
  for (int j=0; j<cyl_rings; j++) {
     theta = 0.01*d_theta;
     for (int ij=0; ij<pts_1_ring; ij++) {
        theta = theta + d_theta;
        int n = 2*(int)k(Nrings-1) + j*pts_1_ring + ij;
        geomVector point (pagp->cylinder_radius*MAC::cos(theta)
                        , pagp->cylinder_height
                              *(-1./2. + (j+0.5)/(double(cyl_rings)))
                        , pagp->cylinder_radius*MAC::sin(theta) );
        geomVector normal(pagp->cylinder_radius * MAC::cos(theta)
                        , 0.
                        , pagp->cylinder_radius * MAC::sin(theta));
        m_surface_points[n]->operator=(point);
        m_surface_area[n]->operator()(0) = cell_area;
        m_surface_normal[n]->operator=(normal);
     }
  }

  // Translate and rotate
  for (size_t i = 0; i < m_surface_area.size(); i++) {
     m_geometric_rigid_body->rotate(m_surface_points[i]);
     m_geometric_rigid_body->rotate(m_surface_normal[i]);
     m_surface_points[i]->translate(*pgc);
  }

}




//---------------------------------------------------------------------------
void DS_3Dspheroidcylinder:: compute_number_of_surface_variables(
                                          double const& surface_cell_scale
                                        , double const& dx )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3Dspheroidcylinder:: compute_number_of_surface_variables" ) ;

  struct FS_3Dspheroidcylinder_Additional_Param const* pagp =
   dynamic_cast<FS_3Dspheroidcylinder*>(m_geometric_rigid_body)
      ->get_ptr_FS_3Dspheroidcylinder_Additional_Param();

  size_t temp = (size_t)((1. / surface_cell_scale) 
              * (4. * MAC::pi() * pagp->cylinder_radius * pagp->cylinder_radius) 
              / (dx * dx));

  // Getting the nearest even number for discrete points on a hemisphere
  Nsphere_by_2 = (size_t)(round((double)temp * 0.5) * 2.);

  double eta_temp = MAC::pi() / 2. - (2. / ar) * MAC::sqrt(MAC::pi() / Nsphere_by_2) * MAC::sin(MAC::pi() / 4.);
  double Rring_temp = 2. * MAC::sin(eta_temp / 2.);
  double Npm1 = round(Nsphere_by_2 * MAC::pow(Rring_temp / MAC::sqrt(2.), 2.));

  double dh = 1. - MAC::sqrt(Npm1 / Nsphere_by_2);
  double Nr = round((pagp->cylinder_height/pagp->cylinder_radius)/2./dh);
  Ntot = (size_t)(2. * Nsphere_by_2 + Nr * (Nsphere_by_2 - Npm1));

}
