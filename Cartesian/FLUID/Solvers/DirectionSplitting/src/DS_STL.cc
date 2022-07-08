#include <DS_STL.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <MAC.hh>
using std::endl;


//---------------------------------------------------------------------------
DS_STL:: DS_STL()
//---------------------------------------------------------------------------
  : m_geometric_rigid_body( NULL )
{
  MAC_LABEL( "DS_STL:: DS_STL" ) ;

}




//---------------------------------------------------------------------------
DS_STL:: DS_STL( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
  : m_geometric_rigid_body( pgrb )
{
  MAC_LABEL( "DS_STL:: DS_STL" ) ;



}




//---------------------------------------------------------------------------
DS_STL:: ~DS_STL()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: ~DS_STL" ) ;

  if ( !m_surface_points.empty() ) m_surface_points.clear();

}




//---------------------------------------------------------------------------
bool DS_STL:: isIn( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: isIn(pt)" ) ;

  // return ( m_geometric_rigid_body->isIn( pt ) );

}




//---------------------------------------------------------------------------
bool DS_STL:: isIn( double const& x, double const& y, double const& z )
	const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: isIn(x,y,z)" ) ;

  double x1, y1, z1, x2, y2, z2, x3, y3, z3;
  double xC, yC, zC;

  double Xmin = 0.;//UF->primary_grid()->get_main_domain_min_coordinate(0);
  double Xmax = 1.;//UF->primary_grid()->get_main_domain_max_coordinate(0);
  double Ymin = 0.;//UF->primary_grid()->get_main_domain_min_coordinate(1);
  double Ymax = 1.;//UF->primary_grid()->get_main_domain_max_coordinate(1);
  double Zmin = 0.;//UF->primary_grid()->get_main_domain_min_coordinate(2);
  double Zmax = 1.;//UF->primary_grid()->get_main_domain_max_coordinate(2);

  double trdx = (Xmax - Xmin) / Nopx;
  double trdz = (Zmax - Zmin) / Nopz;

  double eps = 0.01*trdx; double level_set;

  xC = x; yC = y; zC = z;

  int ii=-1; int kk=-1; 

  // Searching for the right 2D Halo
  for (int i=0; i<Nopx; i++)
  {
     double xmin = Xmin + i * trdx;     double xmax = Xmin + (i+1) * trdx;
	         	 
     for (int k=0; k<Nopz; k++)
     {
         double zmin = Zmin + k * trdz;     double zmax = Zmin + (k+1) * trdz;

         if ((xC > xmin-eps) && (xC < xmax+eps) && (zC > zmin-eps) && (zC < zmax+eps))
	 {
	     ii=i; 
	     kk=k;
         }
     }
  }

  int itrscount = 0;
  double tri1[3], tri2[3], tri3[3], q1[3], q2[3];
	     
  if (ii > -1 && kk > -1)
  {
     for (int m=0; m < tridx_xz[ii][kk]/3; m++)
     {
        // first point
	x1 = std::get<0>(ttrbox_xz[ii][kk][3*m  ]);
	y1 = std::get<1>(ttrbox_xz[ii][kk][3*m  ]);
	z1 = std::get<2>(ttrbox_xz[ii][kk][3*m  ]);

	tri1[0] = x1; tri1[1] = y1; tri1[2] = z1;
		
	// second point
	x2 = std::get<0>(ttrbox_xz[ii][kk][3*m + 1]);
	y2 = std::get<1>(ttrbox_xz[ii][kk][3*m + 1]);
	z2 = std::get<2>(ttrbox_xz[ii][kk][3*m + 1]);
	tri2[0] = x2; tri2[1] = y2; tri2[2] = z2;
	       	
	// third point
        x3 = std::get<0>(ttrbox_xz[ii][kk][3*m + 2]);
        y3 = std::get<1>(ttrbox_xz[ii][kk][3*m + 2]);
	z3 = std::get<2>(ttrbox_xz[ii][kk][3*m + 2]);
	tri3[0] = x3; tri3[1] = y3; tri3[2] = z3;

	q1[0] = xC; q1[1] = yC;        q1[2] = zC;
        q2[0] = xC; q2[1] = Ymin-1.0;  q2[2] = zC;

	if (intersect3d(q1,q2,tri1,tri2,tri3))
	{
	    itrscount = itrscount + 1;
	}
     }
  }

  if ( itrscount % 2 == 0 ) // even number -> outside  
     level_set = 1.0;
  else 
     level_set = -1.0;   

  return(level_set);

  // return ( m_geometric_rigid_body->isIn( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_STL:: level_set_value( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: level_set_value(pt)" ) ;

  // return ( m_geometric_rigid_body->level_set_value( pt ) );

}




//---------------------------------------------------------------------------
double DS_STL:: level_set_value( double const& x
                                     , double const& y
                                     , double const& z ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: level_set_value(x,y,z)" ) ;

  return isIn(x,y,z);

  // return ( m_geometric_rigid_body->level_set_value( x, y, z ) );

}




//---------------------------------------------------------------------------
double DS_STL:: get_distanceTo( geomVector const& source,
                                      geomVector const& rayDir,
                                      double const& delta ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: get_distanceTo" ) ;

  double tri1[3], tri2[3], tri3[3], q1[3], q2[3];
  double x1,y1,z1,x2,y2,z2,x3,y3,z3;

  double Xmin = 0.;//UF->primary_grid()->get_main_domain_min_coordinate(0) ;
  double Xmax = 1.;//UF->primary_grid()->get_main_domain_max_coordinate(0) ;
  double Ymin = 0.;//UF->primary_grid()->get_main_domain_min_coordinate(1) ;
  double Ymax = 1.;//UF->primary_grid()->get_main_domain_max_coordinate(1) ;
  double Zmin = 0.;//UF->primary_grid()->get_main_domain_min_coordinate(2) ;
  double Zmax = 1.;//UF->primary_grid()->get_main_domain_max_coordinate(2) ;

  double trdx = (Xmax - Xmin) / Nopx;
  double trdy = (Ymax - Ymin) / Nopy;
  double trdz = (Zmax - Zmin) / Nopz; 

  int ii, kk;

  double eps = 0.01*trdx;

  int dir = 1;

  double xvalue, yvalue, zvalue, xright2, xcenter2, xleft2;

  // yz
  if ( dir == 0 )
  {
 	for (int i=0; i<Nopy; i++)
  	{
	    double ymin = Ymin + i * trdy;     double ymax = Ymin + (i+1) * trdy;
	         	 
            for (int k=0; k<Nopz; k++)
	    {
                   double zmin = Zmin + k * trdz;     double zmax = Zmin + (k+1) * trdz;

                   if ((yvalue > ymin-eps) && (yvalue < ymax+eps) && (zvalue > zmin-eps) && (zvalue < zmax+eps))
	           {
		       ii=i; 
		       kk=k;
	           }
	    }
        } 

        for (int m=0; m<tridx_yz[ii][kk]/3; m++)
	{
		// first point
                x1 = std::get<0>(ttrbox_yz[ii][kk][3*m  ]);
		y1 = std::get<1>(ttrbox_yz[ii][kk][3*m  ]);
		z1 = std::get<2>(ttrbox_yz[ii][kk][3*m  ]);
		tri1[0] = x1; tri1[1] = y1; tri1[2] = z1;
		
		// second point
	        x2 = std::get<0>(ttrbox_yz[ii][kk][3*m + 1]);
		y2 = std::get<1>(ttrbox_yz[ii][kk][3*m + 1]);
		z2 = std::get<2>(ttrbox_yz[ii][kk][3*m + 1]);
		tri2[0] = x2; tri2[1] = y2; tri2[2] = z2;
	       	
	       	// third point
                x3 = std::get<0>(ttrbox_yz[ii][kk][3*m + 2]);
		y3 = std::get<1>(ttrbox_yz[ii][kk][3*m + 2]);
		z3 = std::get<2>(ttrbox_yz[ii][kk][3*m + 2]);
		tri3[0] = x3; tri3[1] = y3; tri3[2] = z3;

		q1[0] = xleft2;  q1[1] = yvalue;  q1[2] = zvalue;
                q2[0] = xright2; q2[1] = yvalue;  q2[2] = zvalue;

		if (intersect3d(q1,q2,tri1,tri2,tri3))
		{
		   // 974
		   double tmp1[3], tmp2[3], tmp3[3], tmp4[3];
	           diffProduct(tri2,tri1,tmp1);
		   diffProduct(tri3,tri1,tmp2);
                   double N[3];
		   crossProduct(tmp1,tmp2,N);
		   diffProduct(q1,tri1,tmp3);
		   diffProduct(q2,q1,tmp4);
		   double t = -dotProduct(tmp3,N)/dotProduct(tmp4,N);
		   xcenter2 = q1[0] + t * tmp4[0];
		}
	}
  }
  // xz
  if ( dir == 1 )
  {
 	for (int i=0; i<Nopx; i++)
  	{
	    double xmin = Xmin + i * trdx;     double xmax = Xmin + (i+1) * trdx;
	         	 
            for (int k=0; k<Nopz; k++)
	    {
                   double zmin = Zmin + k * trdz;     double zmax = Zmin + (k+1) * trdz;
                   if ((yvalue > xmin-eps) && (yvalue < xmax+eps) && (zvalue > zmin-eps) && (zvalue < zmax+eps))
	           {
		       ii=i; 
		       kk=k;
	           }
	    }
        } 

        for (int m=0; m<tridx_xz[ii][kk]/3; m++)
	{
		// first point
                x1 = std::get<0>(ttrbox_xz[ii][kk][3*m  ]);
		y1 = std::get<1>(ttrbox_xz[ii][kk][3*m  ]);
		z1 = std::get<2>(ttrbox_xz[ii][kk][3*m  ]);
		tri1[0] = x1; tri1[1] = y1; tri1[2] = z1;
		
		// second point
	        x2 = std::get<0>(ttrbox_xz[ii][kk][3*m + 1]);
		y2 = std::get<1>(ttrbox_xz[ii][kk][3*m + 1]);
		z2 = std::get<2>(ttrbox_xz[ii][kk][3*m + 1]);
		tri2[0] = x2; tri2[1] = y2; tri2[2] = z2;
	       	
	       	// third point
                x3 = std::get<0>(ttrbox_xz[ii][kk][3*m + 2]);
		y3 = std::get<1>(ttrbox_xz[ii][kk][3*m + 2]);
		z3 = std::get<2>(ttrbox_xz[ii][kk][3*m + 2]);
		tri3[0] = x3; tri3[1] = y3; tri3[2] = z3;

		q1[0] = yvalue; q1[1] = xleft2;   q1[2] = zvalue;
                q2[0] = yvalue; q2[1] = xright2;  q2[2] = zvalue;

		if (intersect3d(q1,q2,tri1,tri2,tri3))
		{
		   double tmp1[3], tmp2[3], tmp3[3], tmp4[3];
	           diffProduct(tri2,tri1,tmp1);
		   diffProduct(tri3,tri1,tmp2);
                   double N[3];
		   crossProduct(tmp1,tmp2,N);
		   diffProduct(q1,tri1,tmp3);
		   diffProduct(q2,q1,tmp4);
		   double t = -dotProduct(tmp3,N)/dotProduct(tmp4,N);
		   xcenter2 = q1[1] + t * tmp4[1];
		}
	}
  }
  // xy
  if ( dir == 2 )
  {
 	for (int i=0; i<Nopx; i++)
  	{
	    double xmin = Xmin + i * trdx;     double xmax = Xmin + (i+1) * trdx;
	         	 
            for (int k=0; k<Nopy; k++)
	    {
                   double ymin = Ymin + k * trdy;     double ymax = Ymin + (k+1) * trdy;

                   if ((yvalue > xmin-eps) && (yvalue < xmax+eps) && (zvalue > ymin-eps) && (zvalue < ymax+eps))
	           {
		       ii=i; 
		       kk=k;
	           }
	    }
        } 

        for (int m=0; m<tridx_xy[ii][kk]/3; m++)
	{
		// first point
                x1 = std::get<0>(ttrbox_xy[ii][kk][3*m  ]);
	        y1 = std::get<1>(ttrbox_xy[ii][kk][3*m  ]);
	        z1 = std::get<2>(ttrbox_xy[ii][kk][3*m  ]);
		tri1[0] = x1; tri1[1] = y1; tri1[2] = z1;
		
		// second point
	        x2 = std::get<0>(ttrbox_xy[ii][kk][3*m + 1]);
	        y2 = std::get<1>(ttrbox_xy[ii][kk][3*m + 1]);
	        z2 = std::get<2>(ttrbox_xy[ii][kk][3*m + 1]);
		tri2[0] = x2; tri2[1] = y2; tri2[2] = z2;
	       	
	       	// third point
                x3 = std::get<0>(ttrbox_xy[ii][kk][3*m + 2]);
	        y3 = std::get<1>(ttrbox_xy[ii][kk][3*m + 2]);
	        z3 = std::get<2>(ttrbox_xy[ii][kk][3*m + 2]);
		tri3[0] = x3; tri3[1] = y3; tri3[2] = z3;

		q1[0] = yvalue; q1[1] = zvalue;  q1[2] = xleft2;
                q2[0] = yvalue; q2[1] = zvalue;  q2[2] = xright2;

		if (intersect3d(q1,q2,tri1,tri2,tri3))
		{
		   double tmp1[3], tmp2[3], tmp3[3], tmp4[3];
	           diffProduct(tri2,tri1,tmp1);
		   diffProduct(tri3,tri1,tmp2);
                   double N[3];
		   crossProduct(tmp1,tmp2,N);
		   diffProduct(q1,tri1,tmp3);
		   diffProduct(q2,q1,tmp4);
		   double t = -dotProduct(tmp3,N)/dotProduct(tmp4,N);
		   xcenter2 = q1[2] + t * tmp4[2];		
		}
	}
  }

  // return (m_geometric_rigid_body->distanceTo(source, rayDir, delta));

}




//---------------------------------------------------------------------------
geomVector DS_STL:: get_rigid_body_velocity( geomVector const& pt ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: rigid_body_velocity(pt)" ) ;

  geomVector value(0.,0.,0.);

  return (value);

}




//---------------------------------------------------------------------------
geomVector DS_STL:: get_rigid_body_angular_velocity( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: rigid_body_angular_velocity()" ) ;

  geomVector value(0.,0.,0.);

  return (value);

}




//---------------------------------------------------------------------------
std::tuple<double,double> DS_STL:: get_mass_and_density() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: get_mass_and_density()" ) ;

  // return ( m_geometric_rigid_body->get_mass_and_density() );

}




//---------------------------------------------------------------------------
double DS_STL:: get_circumscribed_radius( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: get_circumscribed_radius()" ) ;

  //return (m_geometric_rigid_body->get_circumscribed_radius());

}




//---------------------------------------------------------------------------
geomVector const* DS_STL:: get_ptr_to_gravity_centre( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: get_ptr_to_gravity_centre( )" ) ;

  //return (dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
    //                          ->get_ptr_to_gravity_centre());

}




//---------------------------------------------------------------------------
void DS_STL:: initialize_surface_variables( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: initialize_surface_variables" ) ;

  if (m_surface_points.empty()) {
     m_surface_points.reserve( Ntot );
     m_surface_area.reserve( Ntot );
     m_surface_normal.reserve( Ntot );
     m_surface_Pforce.reserve( Ntot );
     m_surface_Vforce.reserve( Ntot );
     m_surface_Tgrad.reserve( Ntot );

     geomVector vvv(3);

     for (size_t i = 0; i < Ntot; ++i) {
        m_surface_points.push_back( new geomVector(3) );
        m_surface_area.push_back( new geomVector(1) );
        m_surface_normal.push_back( new geomVector(3) );
        m_surface_Pforce.push_back( vvv );
        m_surface_Vforce.push_back( vvv );
        m_surface_Tgrad.push_back( 0. );
     }
   }

}




//---------------------------------------------------------------------------
void DS_STL:: update_RB_position_and_velocity(geomVector const& pos,
                                                    geomVector const& vel,
                                                    geomVector const& ang_vel,
                                   vector<geomVector> const& periodic_directions)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: update_RB_position_and_velocity" ) ;

  //return (m_geometric_rigid_body->update_RB_position_and_velocity(pos,vel
    //                                                              ,ang_vel
      //                                                   ,periodic_directions));

}




//---------------------------------------------------------------------------
vector<geomVector*> DS_STL:: get_rigid_body_surface_points( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: get_rigid_body_surface_points" ) ;

  return (m_surface_points);

}




//---------------------------------------------------------------------------
vector<geomVector*> DS_STL:: get_rigid_body_surface_normals( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: get_rigid_body_surface_normals" ) ;

  return (m_surface_normal);

}




//---------------------------------------------------------------------------
vector<geomVector*> DS_STL:: get_rigid_body_surface_areas( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: get_rigid_body_surface_areas" ) ;

  return (m_surface_area);

}




//---------------------------------------------------------------------------
vector<geomVector*> DS_STL:: get_rigid_body_haloZone( ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: get_rigid_body_haloZone" ) ;

  return (m_halo_zone);

}



//---------------------------------------------------------------------------
void DS_STL:: update_Pforce_on_surface_point( size_t const& i
                                                  , geomVector const& value )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: update_Pforce_on_surface_point" ) ;

  m_surface_Pforce[i] = value;

}




//---------------------------------------------------------------------------
void DS_STL:: update_Vforce_on_surface_point( size_t const& i
                                                  , geomVector const& value )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: update_Vforce_on_surface_point" ) ;

  m_surface_Vforce[i] = value;

}




//---------------------------------------------------------------------------
void DS_STL:: update_Tgrad_on_surface_point( size_t const& i
                                                 , double const& value )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: update_Tgrad_on_surface_point" ) ;

  m_surface_Tgrad[i] = value;

}




//---------------------------------------------------------------------------
void DS_STL:: correct_surface_discretization( FV_Mesh const* MESH )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: correct_surface_discretization( )" ) ;

  // vector<geomVector> const* p_pbc =
  //                 dynamic_cast<FS_RigidBody*>(m_geometric_rigid_body)
  //                             ->get_ptr_to_periodic_directions();
  //
  // if (p_pbc) {
     boolVector const* periodic_comp = MESH->get_periodic_directions();
     size_t dim = MESH->nb_space_dimensions() ;

     for (size_t dir=0;dir < dim; dir++) {
        bool is_periodic = periodic_comp->operator()( dir );

        if (is_periodic) {
           double isize = MESH->get_main_domain_max_coordinate(dir)
                        - MESH->get_main_domain_min_coordinate(dir);
           double imin = MESH->get_main_domain_min_coordinate(dir);

           for (size_t i = 0; i < m_surface_area.size(); i++) {
              m_surface_points[i]->operator()(dir) =
                  m_surface_points[i]->operator()(dir)
                - MAC::floor((m_surface_points[i]->operator()(dir)-imin)/isize)
                  * isize;
           }
        }
     }
  // }
}




//---------------------------------------------------------------------------
void DS_STL:: write_surface_discretization( const std::string& file)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: write_surface_discretization" ) ;

  std::ofstream out;

  out.open(file.c_str());
  out << "x ,y ,z ,nx ,ny ,nz ,area ,Fpx ,Fpy ,Fpz ,Fvx ,Fvy ,Fvz, Tgrad" << endl;

  for (size_t i = 0; i < m_surface_area.size(); i++) {
     if ((m_surface_Pforce[i].calcNorm() != 0) ||
         (m_surface_Vforce[i].calcNorm() != 0) ||
         (m_surface_Tgrad[i] != 0))
        out << m_surface_points[i]->operator()(0) << " ,"
            << m_surface_points[i]->operator()(1) << " ,"
            << m_surface_points[i]->operator()(2) << " ,"
            << m_surface_normal[i]->operator()(0) << " ,"
            << m_surface_normal[i]->operator()(1) << " ,"
            << m_surface_normal[i]->operator()(2) << " ,"
            << m_surface_area[i]->operator()(0) << " ,"
            << m_surface_Pforce[i](0) << " ,"
            << m_surface_Pforce[i](1) << " ,"
            << m_surface_Pforce[i](2) << " ,"
            << m_surface_Vforce[i](0) << " ,"
            << m_surface_Vforce[i](1) << " ,"
            << m_surface_Vforce[i](2) << " ,"
            << m_surface_Tgrad[i] << endl;
  }

  out.close();

}

//---------------------------------------------------------------------------
double DS_STL:: dotProduct(double vect_A[], double vect_B[]) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: dotProduct" ) ;

  double product = 0;
 
  for (int i = 0; i < 3; i++)
     product = product + vect_A[i] * vect_B[i];

  return product;
}

//---------------------------------------------------------------------------
void DS_STL:: crossProduct(double vect_A[], double vect_B[] 
		                          , double cross_P[]) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: crossProduct" ) ;
 
  cross_P[0] = vect_A[1] * vect_B[2] - vect_A[2] * vect_B[1];
  cross_P[1] = vect_A[2] * vect_B[0] - vect_A[0] * vect_B[2];
  cross_P[2] = vect_A[0] * vect_B[1] - vect_A[1] * vect_B[0];
}

//---------------------------------------------------------------------------
void DS_STL:: diffProduct(double vect_A[], double vect_B[]
		                         , double diff_P[]) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: diffProduct" ) ;

  diff_P[0] = vect_A[0] - vect_B[0];
  diff_P[1] = vect_A[1] - vect_B[1];
  diff_P[2] = vect_A[2] - vect_B[2];
}

//---------------------------------------------------------------------------
int DS_STL:: orient3d(double vect_A[], double vect_B[], double vect_C[]
		                                 , double vect_D[]) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: orient3d" ) ;

  double tmp1[3], tmp2[3], tmp3[3], tmp4[3];

  diffProduct(vect_B,vect_A,tmp1);
  diffProduct(vect_C,vect_A,tmp2);
  diffProduct(vect_D,vect_A,tmp3);
  crossProduct(tmp1,tmp2,tmp4);
  if (dotProduct(tmp4,tmp3) > 0)
	  return 0;
  else
	  return 1;
}

//---------------------------------------------------------------------------
int DS_STL:: intersect3d(double q1[], double q2[], double tri1[]
		                    , double tri2[], double tri3[]) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: intersect3d" ) ;

  int s1 = orient3d(q1,tri1,tri2,tri3);
  int s2 = orient3d(q2,tri1,tri2,tri3);
  // Test whether the two extermities of the segment
  // are on the same side of the supporting plane of
  // the triangle
  if (s1 == s2) 
     return 0;

  // Now we know that the segment 'straddles' the supporing
  // plane. We need to test whether the three tetrahedra formed
  // by the segment and the three edges of the triangle have
  // the same orientation
  int s3 = orient3d(q1,q2,tri1,tri2);
  int s4 = orient3d(q1,q2,tri2,tri3);
  int s5 = orient3d(q1,q2,tri3,tri1);
  return (s3 == s4 && s4 == s5);
}

//---------------------------------------------------------------------------
void DS_STL:: readSTL()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_STL:: readSTL" ) ;

  double xp,yp,zp,Rp,vx,vy,vz,wx,wy,wz,Tp,off;
  double xn,yn,zn,x1,y1,z1,x2,y2,z2,x3,y3,z3;
  int kk;
  string notsodummy,dummy;
  clock_t start;
  double duration;
  int dim =3; // to be removed

  string solid_filename("q");
  	  
     kk = 0;	     
     start = clock(); // start time for reading the set of triangles
        
     // reading STL file format
     // filename
     std::ostringstream os2;
     os2 << "./InputFiles/" << solid_filename;
     std::string filename = os2.str();

     // check if the file is ASCII or binary
     std::ifstream inFilep(filename);
	
     int c;
     //if ( my_rank == is_master && field == 0 )
     {
        std::cout << endl;
	std::cout << "***" << endl;
	std::cout << "================= STL FILE =================" << endl;
	std::cout << "   Read STL file: *" << solid_filename << "*" << endl;
     }

     while( (c = inFilep.get()) != EOF && c <= 127);
	
     //if ( my_rank == is_master && field == 0 )
     {
        if( c == EOF ) 
	{
           std::cout << "   File is in ASCII format" << endl;
        }
	else 
	   std::cout << "   File is in Binary format" << endl;
     } 
     inFilep.close();
     // end reading STL file format

     // reading triangulation     	
     if ( c == EOF ) // ASCII
     {
	std::ifstream inFile(filename);	
        //inFile.open(filename.c_str());
     	string line;

	if ( dim == 3 )
	   getline(inFile,line);

     	while(getline(inFile,line))
     	{
	   std::istringstream iss1(line);
     	   iss1 >> notsodummy;
	   if ( dim == 3  && notsodummy.compare("endsolid") == 0 )
               break;
     	   iss1 >> dummy; 
     	   iss1 >> xn;
     	   iss1 >> yn;
     	   iss1 >> zn;
           getline(inFile,line);
           getline(inFile,line);
	   std::istringstream iss2(line);
           iss2 >> dummy;
     	   iss2 >> x1;
     	   iss2 >> y1;
     	   iss2 >> z1;
	   getline(inFile,line);
	   std::istringstream iss3(line);
           iss3 >> dummy;
     	   iss3 >> x2;
     	   iss3 >> y2;
      	   iss3 >> z2;
	   if ( dim == 3 )
	   {
	      getline(inFile,line);
	      std::istringstream iss4(line);
              iss4 >> dummy;
     	      iss4 >> x3;
     	      iss4 >> y3;
      	      iss4 >> z3;
	   }
	   getline(inFile,line);
	   getline(inFile,line);

	   // Vertices
           Llvls.push_back(std::make_tuple(x1, y1, z1));
	   Llvls.push_back(std::make_tuple(x2, y2, z2));
	   if ( dim == 3 )	   
	      Llvls.push_back(std::make_tuple(x3, y3, z3));
	   

	   // Normals
	   Llvns.push_back(std::make_tuple(xn, yn, zn));

	   if ( dim == 2 )
	      kk=kk+2;
	   if ( dim == 3 )
	      kk=kk+3;
	 
        }
         	
        inFile.close();	  
	Npls=kk;

     } // end ASCII
     else // binary
     {
	std::ifstream inFileb(filename, std::ifstream::binary);

        // rdbuf returns a streambuf object associated with the
	// input fstream object ifs.

	std::filebuf* pbuf = inFileb.rdbuf();

	// Calculate the file's size.

	auto size = pbuf->pubseekoff(0, inFileb.end);

	// Set the position pointer to the beginning of the file.

	pbuf->pubseekpos(0);

	// Allocate memory to contain file data.

	char* buffer = new char[(size_t)size];

	// Get file data. sgetn grabs all the characters from the streambuf 
	// object 'pbuf'. The return value of sgetn is the number of characters
	// obtained - ordinarily, this value should be checked for equality
	// against the number of characters requested.

	pbuf->sgetn(buffer, size);

	char * bufptr = buffer;

	bufptr += 80;  // Skip past the header.
	bufptr += 4;   // Skip past the number of triangles.

	double nx,ny,nz;

        while (bufptr < buffer + size)
       	{

		xn = *(float *)(bufptr);	
                yn = *(float *)(bufptr + 4);
		zn = *(float *)(bufptr + 8);
		bufptr += 12;

		Llvns.push_back(std::make_tuple(xn, yn, zn));

		x1 = *(float *)(bufptr);
		y1 = *(float *)(bufptr + 4);
		z1 = *(float *)(bufptr + 8);
		bufptr += 12;

		x2 = *(float *)(bufptr);
		y2 = *(float *)(bufptr + 4);
		z2 = *(float *)(bufptr + 8);
		bufptr += 12;

		x3 = *(float *)(bufptr);
		y3 = *(float *)(bufptr + 4);
		z3 = *(float *)(bufptr + 8);
		bufptr += 12;

		Llvls.push_back(std::make_tuple(x1, y1, z1));
	        Llvls.push_back(std::make_tuple(x2, y2, z2));
	        Llvls.push_back(std::make_tuple(x3, y3, z3));
	    
                kk=kk+3;
		bufptr += 2;
	}

        inFileb.close();
	Npls=kk;
        delete[] buffer;

     }
     // end reading triangulation
	
     //if ( my_rank == is_master && field == 0 )
     {        
         std::cout << "   Delaunay triangles: " << Npls/3 << endl;
	 duration = ( clock() - start ) / (double) CLOCKS_PER_SEC;
         std::cout << "   File read in " << duration << "s" << endl;
     }

     // Filtering operations

     int filter_trbox = 1; // 1D distances approach

     if ( filter_trbox )
     {

        double Xmin = 0.;//UF->primary_grid()->get_main_domain_min_coordinate(0) ;
        double Xmax = 1.;//UF->primary_grid()->get_main_domain_max_coordinate(0) ;
        double Ymin = 0.;//UF->primary_grid()->get_main_domain_min_coordinate(1) ;
        double Ymax = 1.;//UF->primary_grid()->get_main_domain_max_coordinate(1) ;
        double Zmin = 0.;//UF->primary_grid()->get_main_domain_min_coordinate(2) ;
        double Zmax = 1.;//UF->primary_grid()->get_main_domain_max_coordinate(2) ;

        //if ( field == 0 && my_rank == is_master )
        {
       	   std::cout << "   - Box dimensions" << endl;
       	   std::cout << "   x [" << Xmin  << " " << Xmax << "]" <<  endl;
	   std::cout << "   y [" << Ymin  << " " << Ymax << "]" <<  endl;
	   std::cout << "   z [" << Zmin  << " " << Zmax << "]" <<  endl << endl;
         }
     
        double trdx = (Xmax - Xmin) / Nopx;
        double trdy = (Ymax - Ymin) / Nopy;
        double trdz = (Zmax - Zmin) / Nopz;

        double cenh = 1.0; // if too small and triangles too large, 
                          // we may miss some triangles in the halo
       
        double enhx = trdx * cenh;
        double enhy = trdy * cenh;
        double enhz = trdz * cenh;

        // xz
        for (int i=0; i<Nopx; i++)
        {
	   double xmin = Xmin + i * trdx;     double xmax = Xmin + (i+1) * trdx;
           double xxmin = xmin - enhx;        double xxmax = xmax + enhx;
	         	 
           for (int k=0; k<Nopz; k++)
	   {
              double zmin = Zmin + k * trdz;     double zmax = Zmin + (k+1) * trdz;
              double zzmin = zmin - enhz;        double zzmax = zmax + enhz;

	      int m=0;
	      for (int l=0; l<Npls/3; l++)
              {
                 if  ( ((std::get<0>(Llvls[3*l])   > xxmin && std::get<0>(Llvls[3*l])   < xxmax)   ||
	                (std::get<0>(Llvls[3*l+1]) > xxmin && std::get<0>(Llvls[3*l+1]) < xxmax)   ||
	                (std::get<0>(Llvls[3*l+2]) > xxmin && std::get<0>(Llvls[3*l+2]) < xxmax))  &&
	               ((std::get<2>(Llvls[3*l])   > zzmin && std::get<2>(Llvls[3*l])   < zzmax)   || 
	                (std::get<2>(Llvls[3*l+1]) > zzmin && std::get<2>(Llvls[3*l+1]) < zzmax)   ||
	                (std::get<2>(Llvls[3*l+2]) > zzmin && std::get<2>(Llvls[3*l+2]) < zzmax)) )
	         {
                    x1 = std::get<0>(Llvls[3*l  ]);
	            x2 = std::get<0>(Llvls[3*l+1]); 
                    x3 = std::get<0>(Llvls[3*l+2]); 
	            y1 = std::get<1>(Llvls[3*l  ]);
	            y2 = std::get<1>(Llvls[3*l+1]); 
                    y3 = std::get<1>(Llvls[3*l+2]);
                    z1 = std::get<2>(Llvls[3*l  ]);
	            z2 = std::get<2>(Llvls[3*l+1]); 
                    z3 = std::get<2>(Llvls[3*l+2]);

		    ttrbox_xz[i][k].push_back(std::make_tuple(x1, y1, z1));
		    ttrbox_xz[i][k].push_back(std::make_tuple(x2, y2, z2));
		    ttrbox_xz[i][k].push_back(std::make_tuple(x3, y3, z3));

		    m=m+3;
	         }  
	     }
	     tridx_xz[i][k] = m;
           }
        }
        // xy
        for (int i=0; i<Nopx; i++)
        {
	   double xmin = Xmin + i * trdx;     double xmax = Xmin + (i+1) * trdx;
           double xxmin = xmin - enhx;        double xxmax = xmax + enhx;
	         	 
           for (int k=0; k<Nopy; k++)
	   {
              double ymin = Ymin + k * trdy;     double ymax = Ymin + (k+1) * trdy;
              double yymin = ymin - enhy;        double yymax = ymax + enhy;

	      int m=0;
	      for (int l=0; l<Npls/3; l++)
              {
                 if  ( ((std::get<0>(Llvls[3*l])   > xxmin && std::get<0>(Llvls[3*l])   < xxmax)   ||
	                (std::get<0>(Llvls[3*l+1]) > xxmin && std::get<0>(Llvls[3*l+1]) < xxmax)   ||
	                (std::get<0>(Llvls[3*l+2]) > xxmin && std::get<0>(Llvls[3*l+2]) < xxmax))  &&
		       ((std::get<1>(Llvls[3*l])   > yymin && std::get<1>(Llvls[3*l])   < yymax)   || 
	                (std::get<1>(Llvls[3*l+1]) > yymin && std::get<1>(Llvls[3*l+1]) < yymax)   ||
	                (std::get<1>(Llvls[3*l+2]) > yymin && std::get<1>(Llvls[3*l+2]) < yymax)) )
	         {

                    x1 = std::get<0>(Llvls[3*l  ]);
	            x2 = std::get<0>(Llvls[3*l+1]); 
                    x3 = std::get<0>(Llvls[3*l+2]); 
	            y1 = std::get<1>(Llvls[3*l  ]);
	            y2 = std::get<1>(Llvls[3*l+1]); 
                    y3 = std::get<1>(Llvls[3*l+2]);
		    z1 = std::get<2>(Llvls[3*l  ]);
	            z2 = std::get<2>(Llvls[3*l+1]); 
                    z3 = std::get<2>(Llvls[3*l+2]);

		    ttrbox_xy[i][k].push_back(std::make_tuple(x1, y1, z1));
		    ttrbox_xy[i][k].push_back(std::make_tuple(x2, y2, z2));
		    ttrbox_xy[i][k].push_back(std::make_tuple(x3, y3, z3));

                    m=m+3;
	         }
	      }    
              tridx_xy[i][k] = m;
           }   
        }
        // yz
        for (int i=0; i<Nopy; i++)
        {
	   double ymin = Ymin + i * trdy;     double ymax = Ymin + (i+1) * trdy;
           double yymin = ymin - enhy;        double yymax = ymax + enhy;
	         	 
           for (int k=0; k<Nopz; k++)
	   {
              double zmin = Zmin + k * trdz;     double zmax = Zmin + (k+1) * trdz;
              double zzmin = zmin - enhz;        double zzmax = zmax + enhz;

	      int m=0;
	      for (int l=0; l<Npls/3; l++)
              {
                 if  ( ((std::get<1>(Llvls[3*l])   > yymin && std::get<1>(Llvls[3*l])   < yymax)   ||
	                (std::get<1>(Llvls[3*l+1]) > yymin && std::get<1>(Llvls[3*l+1]) < yymax)   ||
	                (std::get<1>(Llvls[3*l+2]) > yymin && std::get<1>(Llvls[3*l+2]) < yymax))  &&
		       ((std::get<2>(Llvls[3*l])   > zzmin && std::get<2>(Llvls[3*l])   < zzmax)   || 
	                (std::get<2>(Llvls[3*l+1]) > zzmin && std::get<2>(Llvls[3*l+1]) < zzmax)   ||
	                (std::get<2>(Llvls[3*l+2]) > zzmin && std::get<2>(Llvls[3*l+2]) < zzmax)) )
	         {
	            x1 = std::get<0>(Llvls[3*l  ]);
	            x2 = std::get<0>(Llvls[3*l+1]); 
                    x3 = std::get<0>(Llvls[3*l+2]); 
	            y1 = std::get<1>(Llvls[3*l  ]);
	            y2 = std::get<1>(Llvls[3*l+1]); 
                    y3 = std::get<1>(Llvls[3*l+2]);
		    z1 = std::get<2>(Llvls[3*l  ]);
	            z2 = std::get<2>(Llvls[3*l+1]); 
                    z3 = std::get<2>(Llvls[3*l+2]);

		    ttrbox_yz[i][k].push_back(std::make_tuple(x1, y1, z1));
		    ttrbox_yz[i][k].push_back(std::make_tuple(x2, y2, z2));
		    ttrbox_yz[i][k].push_back(std::make_tuple(x3, y3, z3));

                    m=m+3;
	         }
	      }
              tridx_yz[i][k] = m;
           }
        }

        //if ( field == 0 && my_rank == is_master )
        {
           std::cout << "- Displaying Matrix of filtered triangles xz" << endl;
           for (int i=0; i<Nopx; i++)
           {	       	 
              for (int k=0; k<Nopz; k++)
	      {
                 std::cout << tridx_xz[i][k]/3 << "\t";
	      }
	      std::cout << endl;
           }
           std::cout << "- Displaying Matrix of filtered triangles xy" << endl;
           for (int i=0; i<Nopx; i++)
           {	       	 
              for (int k=0; k<Nopy; k++)
	      {
                 std::cout << tridx_xy[i][k]/3 << "\t";
	      }
	      std::cout << endl;
           }
           std::cout << "- Displaying Matrix of filtered triangles yz" << endl;
           for (int i=0; i<Nopy; i++)
           {	       	 
              for (int k=0; k<Nopz; k++)
	      {
                 std::cout << tridx_yz[i][k]/3 << "\t";
	      }
	      std::cout << endl;
           }

         }
      } // end filtertrbox 
}
     

