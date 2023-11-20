#include <FS_3Dcylinder.hh>
#include <MAC.hh>
#include <math.h>
#define EPSILON 1.e-7
using std::endl;


//---------------------------------------------------------------------------
FS_3Dcylinder:: FS_3Dcylinder()
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_3Dcylinder:: FS_3Dcylinder" ) ;

  m_space_dimension = 3;
  m_shape_type = GEOM_3DCYLINDER;
  m_agp_3dcyl.BottomCenter.resize(3);
  m_agp_3dcyl.TopCenter.resize(3);
  m_agp_3dcyl.BottomToTopVec.resize(3);
  m_agp_3dcyl.RadialRefVec.resize(3);
  m_agp_3dcyl.cylinder_radius = 0.;
  m_agp_3dcyl.cylinder_height = 0.;
}




//---------------------------------------------------------------------------
FS_3Dcylinder:: FS_3Dcylinder( istream& in, size_t& id_ )
//---------------------------------------------------------------------------
  : FS_RigidBody()
{
  MAC_LABEL( "FS_3Dcylinder:: FS_3Dcylinder" ) ;

  // Default parameter
  m_space_dimension = 3;
  m_Id = id_;
  m_shape_type = GEOM_3DCYLINDER;

  // Resize parameters
  m_gravity_center.resize(3);
  m_translational_velocity.resize(3);
  m_angular_velocity.resize(3);
  m_hydro_force.resize(3);
  m_hydro_torque.resize(3);
  m_heat_flux.resize(3);
  m_agp_3dcyl.BottomCenter.resize(3);
  m_agp_3dcyl.TopCenter.resize(3);
  m_agp_3dcyl.BottomToTopVec.resize(3);
  m_agp_3dcyl.RadialRefVec.resize(3);

  // Set the rigid body features from the input stream
  set( in );

}




//---------------------------------------------------------------------------
FS_3Dcylinder:: ~FS_3Dcylinder()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: ~FS_3Dcylinder" ) ;

}




//---------------------------------------------------------------------------
void FS_3Dcylinder:: update( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: update" ) ;

  // Set the rigid body features from the input stream
  set( in );

}




//---------------------------------------------------------------------------
void FS_3Dcylinder:: set( istream& in )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: set" ) ;

  size_t ncorners, i, nfaces, nper ;
  geomVector vv(3);

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

  // Read the bottom and top disk centers
  // Second point is used to compute the radial reference vector
  // and the cylinder radius
  // The radial reference vector is a vector pointing from the axis to the
  // cylinder surface in the plane of the disks
  in >> m_agp_3dcyl.BottomCenter;
  in >> vv;
  m_agp_3dcyl.RadialRefVec = vv - m_agp_3dcyl.BottomCenter ;
  m_agp_3dcyl.cylinder_radius = m_agp_3dcyl.RadialRefVec.calcNorm() ;
  in >> m_agp_3dcyl.TopCenter;
  m_agp_3dcyl.BottomToTopVec = m_agp_3dcyl.TopCenter - m_agp_3dcyl.BottomCenter ;
  m_agp_3dcyl.cylinder_height = m_agp_3dcyl.BottomToTopVec.calcNorm() ;
  in >> nfaces;

  // Set volume
  m_volume = m_mass / m_density ;
}




//---------------------------------------------------------------------------
void FS_3Dcylinder:: display( ostream& out, size_t const& indent_width ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: display" ) ;

  string space( indent_width, ' ' ) ;
  string three( 3, ' ' ) ;
  out << space << "Shape type = " <<
  	FS_RigidBody::GEOMETRICSHAPE_name[m_shape_type] << endl;
  out << space << "Specific attributes" << endl;
  out << space << three << "Radius = " << m_agp_3dcyl.cylinder_radius << endl;
  out << space << three << "Height = " << m_agp_3dcyl.cylinder_height << endl;
  out << space << three << "Center of the bottom disk = "
  	                     << m_agp_3dcyl.BottomCenter(0)
               << three << m_agp_3dcyl.BottomCenter(1)
               << three << m_agp_3dcyl.BottomCenter(2) << endl;
  out << space << three << "Center of the top disk = "
  	                     << m_agp_3dcyl.TopCenter(0)
               << three << m_agp_3dcyl.TopCenter(1)
               << three << m_agp_3dcyl.TopCenter(2) << endl;
  out << space << three << "Axial vector from bottom to top disk center = "
  	                     << m_agp_3dcyl.BottomToTopVec(0)
               << three << m_agp_3dcyl.BottomToTopVec(1)
               << three << m_agp_3dcyl.BottomToTopVec(2) << endl;
  out << space << three << "Radial reference vector = "
  	                     << m_agp_3dcyl.RadialRefVec(0)
               << three << m_agp_3dcyl.RadialRefVec(1)
               << three << m_agp_3dcyl.RadialRefVec(2) << endl;
  out << space << "General attributes" << endl;
  display_general( out, indent_width + 3 );

}




//---------------------------------------------------------------------------
bool FS_3Dcylinder:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: isIn(pt)" ) ;

  bool b_isIn = false;

  geomVector BottomToPoint( pt - m_agp_3dcyl.BottomCenter );
  double dot = ( BottomToPoint , m_agp_3dcyl.BottomToTopVec )
	                        /  m_agp_3dcyl.cylinder_height ;

  double eps_local = EPSILON*m_agp_3dcyl.cylinder_height;

  if ( dot <= m_agp_3dcyl.cylinder_height && dot >= - eps_local )
    if ( BottomToPoint.calcNormSquare() - dot * dot  - eps_local <=
         m_agp_3dcyl.cylinder_radius * m_agp_3dcyl.cylinder_radius )
      b_isIn = true ;

  if (m_periodic_directions) {
     for (size_t i = 0; (i < m_periodic_directions->size()) && !b_isIn; ++i) {
        BottomToPoint = pt - (m_agp_3dcyl.BottomCenter + (*m_periodic_directions)[i]);
        dot = ( BottomToPoint , m_agp_3dcyl.BottomToTopVec  )
            /  m_agp_3dcyl.cylinder_height ;
        if ( dot <= m_agp_3dcyl.cylinder_height && dot >= - eps_local )
          if ( BottomToPoint.calcNormSquare() - dot * dot  - eps_local <=
               m_agp_3dcyl.cylinder_radius * m_agp_3dcyl.cylinder_radius )
            b_isIn = true ;
      }
   }

  return ( b_isIn );

}




//---------------------------------------------------------------------------
bool FS_3Dcylinder:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: isIn(x,y,z)" ) ;

  geomVector pt(x,y,z);

  return ( isIn(pt) );

}




//---------------------------------------------------------------------------
double FS_3Dcylinder:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: level_set_value(pt)" ) ;

  geomVector BottomToPoint( pt - m_agp_3dcyl.BottomCenter );
  double dot = ( BottomToPoint , m_agp_3dcyl.BottomToTopVec )
	                        /  m_agp_3dcyl.cylinder_height ;

  double rad_sq_dist = BottomToPoint.calcNormSquare() - dot * dot -
         m_agp_3dcyl.cylinder_radius * m_agp_3dcyl.cylinder_radius;

  double value = 0.;

  double eps_local = EPSILON*m_agp_3dcyl.cylinder_height;

  if ( dot < m_agp_3dcyl.cylinder_height &&
       dot > -eps_local &&
       rad_sq_dist < eps_local) {
     value = -1.*fabs(fabs(dot) - m_agp_3dcyl.cylinder_height)*fabs(rad_sq_dist);
  } else {
     value = MAC::sqrt(MAC::pow((fabs(dot)-m_agp_3dcyl.cylinder_height),2.)
                     + MAC::pow(fabs(rad_sq_dist),2.));
  }

  if (m_periodic_directions) {
     for (size_t i = 0; i < m_periodic_directions->size(); ++i) {
        BottomToPoint = pt - (m_agp_3dcyl.BottomCenter + (*m_periodic_directions)[i]);
        dot = ( BottomToPoint , m_agp_3dcyl.BottomToTopVec )
      	   /  m_agp_3dcyl.cylinder_height ;

        rad_sq_dist = BottomToPoint.calcNormSquare() - dot * dot -
               m_agp_3dcyl.cylinder_radius * m_agp_3dcyl.cylinder_radius;

        double temp = 0.;

        if ( dot < m_agp_3dcyl.cylinder_height &&
             dot > -eps_local &&
             rad_sq_dist < eps_local) {
           temp = -1.*fabs(fabs(dot) - m_agp_3dcyl.cylinder_height)*fabs(rad_sq_dist);
        } else {
           temp = MAC::sqrt(MAC::pow((fabs(dot)-m_agp_3dcyl.cylinder_height),2.)
                           + MAC::pow(fabs(rad_sq_dist),2.));
        }
        value = MAC::min(temp, value);
     }
  }

  return ( value );

}




//---------------------------------------------------------------------------
double FS_3Dcylinder:: level_set_value( double const& x
                                      , double const& y
                                      , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: level_set_value(x,y,z)" ) ;

  geomVector pt(x,y,z);

  return ( level_set_value(pt) );

}




//---------------------------------------------------------------------------
double FS_3Dcylinder::analytical_distanceTo(geomVector const &source,
                                            geomVector const &rayDir) const
//---------------------------------------------------------------------------
{
  MAC_LABEL("FS_3Dcylinder:: analytical_distanceTo");

  // Ref: https://hugi.scene.org/online/hugi24/coding%20graphics%20chris%20dragan%20raytracing%20shapes.htm

  geomVector D(rayDir);
  geomVector O(source);

  geomVector V(m_agp_3dcyl.BottomToTopVec);
  V = V * (1. / V.calcNorm());
  geomVector C(m_agp_3dcyl.BottomCenter);
  geomVector X(O - C);

  // Coefficients of quadratic equation
  double a = (D, D) - MAC::pow((D, V), 2.);
  double b = 2. * ((D, X) - (D, V) * (X, V));
  double c = (X, X) - MAC::pow((X, V), 2.) - MAC::pow(m_agp_3dcyl.cylinder_radius, 2.);

  double det = MAC::pow(b, 2) - 4. * a * c;

  if (det > 0)
  {
    double t1 = (-b + MAC::sqrt(det)) / (2. * a);
    double t2 = (-b - MAC::sqrt(det)) / (2. * a);
    double m1 = (D, V) * t1 + (X, V);
    double m2 = (D, V) * t2 + (X, V);
    // rayDir accounts for the direction, we only consider
    // positive t1 and t2
    if ((t1 > 0) && (t2 > 0))
    {
      if ((m1 >= 0. && m1 <= m_agp_3dcyl.cylinder_height) &&
          (m2 >= 0. && m2 <= m_agp_3dcyl.cylinder_height))
      {
        return (MAC::min(t1, t2));
      }
      else if (m1 >= 0. && m1 <= m_agp_3dcyl.cylinder_height)
      {
        // second point intersection with a plane
        geomVector Xt(O - C);
        t2 = ((D, V) != 0.) ? -(Xt, V) / (D, V) : 0.;
        geomVector P_t(O + t2 * D);
        if (P_t.calcDist(C) <= m_agp_3dcyl.cylinder_radius)
          if (t2 > 0 && t2 < t1 && t2 < m_agp_3dcyl.cylinder_height)
            return (t2);

        Xt = O - (C + V * m_agp_3dcyl.cylinder_height);
        t2 = ((D, V) != 0.) ? -(Xt, V) / (D, V) : 0.;
        P_t = O + t2 * D;
        if (P_t.calcDist(C + V * m_agp_3dcyl.cylinder_height) <= m_agp_3dcyl.cylinder_radius)
          if (t2 > 0 && t2 < t1 && t2 < m_agp_3dcyl.cylinder_height)
            return (t2);

        return (t1);
      }
      else if (m2 >= 0. && m2 <= m_agp_3dcyl.cylinder_height)
      {
        // first point intersection with a plane
        geomVector Xt(O - C);
        t1 = ((D, V) != 0.) ? -(Xt, V) / (D, V) : 0.;
        geomVector P_t(O + t1 * D);
        if (P_t.calcDist(C) <= m_agp_3dcyl.cylinder_radius)
          if (t1 > 0 && t1 < t2 && t1 < m_agp_3dcyl.cylinder_height)
            return (t1);

        Xt = O - (C + V * m_agp_3dcyl.cylinder_height);
        t1 = ((D, V) != 0.) ? -(Xt, V) / (D, V) : 0.;
        P_t = O + t1 * D;
        if (P_t.calcDist(C + V * m_agp_3dcyl.cylinder_height) <= m_agp_3dcyl.cylinder_radius)
          if (t1 > 0 && t1 < t2 && t1 < m_agp_3dcyl.cylinder_height)
            return (t1);

        return (t2);
      }
      else
      {
        double dist = 0.;
        geomVector Xt(O - C);
        t1 = ((D, V) != 0.) ? -(Xt, V) / (D, V) : 0.;
        geomVector P_t(O + t1 * D);
        if (P_t.calcDist(C) <= m_agp_3dcyl.cylinder_radius)
          if (t1 > 0 && t1 < m_agp_3dcyl.cylinder_height)
            dist = t1;

        Xt = O - (C + V * m_agp_3dcyl.cylinder_height);
        t2 = ((D, V) != 0.) ? -(Xt, V) / (D, V) : 0.;
        P_t = O + t2 * D;
        if (P_t.calcDist(C + V * m_agp_3dcyl.cylinder_height) <= m_agp_3dcyl.cylinder_radius)
          if (t2 > 0 && t2 < m_agp_3dcyl.cylinder_height)
            dist = (dist != 0) ? MAC::min(dist, t2) : t2;

        return (dist);
      }
      return (MAC::min(t1, t2));
    }
    else if (t1 > 0)
    {
      // second point intersection with a plane
      geomVector Xt(O - C);
      t2 = ((D, V) != 0.) ? -(Xt, V) / (D, V) : 0.;
      geomVector P_t(O + t2 * D);
      if (P_t.calcDist(C) <= m_agp_3dcyl.cylinder_radius)
        if (t2 > 0 && t2 < t1 && t2 < m_agp_3dcyl.cylinder_height)
          return (t2);

      Xt = O - (C + V * m_agp_3dcyl.cylinder_height);
      t2 = ((D, V) != 0.) ? -(Xt, V) / (D, V) : 0.;
      P_t = O + t2 * D;
      if (P_t.calcDist(C + V * m_agp_3dcyl.cylinder_height) <= m_agp_3dcyl.cylinder_radius)
        if (t2 > 0 && t2 < t1 && t2 < m_agp_3dcyl.cylinder_height)
          return (t2);

      return (t1);
    }
    else if (t2 > 0)
    {
      // first point intersection with a plane
      geomVector Xt(O - C);
      t1 = ((D, V) != 0.) ? -(Xt, V) / (D, V) : 0.;
      geomVector P_t(O + t1 * D);
      if (P_t.calcDist(C) <= m_agp_3dcyl.cylinder_radius)
        if (t1 > 0 && t1 < t2 && t1 < m_agp_3dcyl.cylinder_height)
          return (t1);

      Xt = O - (C + V * m_agp_3dcyl.cylinder_height);
      t1 = ((D, V) != 0.) ? -(Xt, V) / (D, V) : 0.;
      P_t = O + t1 * D;
      if (P_t.calcDist(C + V * m_agp_3dcyl.cylinder_height) <= m_agp_3dcyl.cylinder_radius)
        if (t1 > 0 && t1 < t2 && t1 < m_agp_3dcyl.cylinder_height)
          return (t1);

      return (t2);
    }
  }

  return (m_agp_3dcyl.cylinder_radius);
}




//---------------------------------------------------------------------------
struct FS_3Dcylinder_Additional_Param const* FS_3Dcylinder::
	get_ptr_FS_3Dcylinder_Additional_Param() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: isIn(x,y,z)" ) ;

  return( &m_agp_3dcyl );

}




//---------------------------------------------------------------------------
void FS_3Dcylinder::update_additional_parameters( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_3Dcylinder:: update_additional_parameters( )" ) ;



}
