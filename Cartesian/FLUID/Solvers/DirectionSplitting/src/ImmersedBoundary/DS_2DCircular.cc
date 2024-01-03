#include <DS_2DCircular.hh>
#include <FS_RigidBody.hh>
#include <FS_Disc.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
using std::endl;


//---------------------------------------------------------------------------
DS_2DCircular:: DS_2DCircular()
//---------------------------------------------------------------------------
  : DS_ImmersedBoundary()
{
  MAC_LABEL( "DS_2DCircular:: DS_2DCircular" ) ;

}




//---------------------------------------------------------------------------
DS_2DCircular:: DS_2DCircular( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
  : DS_ImmersedBoundary( pgrb )
{
  MAC_LABEL( "DS_2DCircular:: DS_2DCircular" ) ;

}




//---------------------------------------------------------------------------
DS_2DCircular:: ~DS_2DCircular()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DCircular:: ~DS_2DCircular" ) ;

}




//---------------------------------------------------------------------------
void DS_2DCircular:: compute_number_of_surface_variables(
                                  double const& surface_cell_scale
                                , double const& dx)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DCircular:: compute_number_of_surface_variables" ) ;


  struct FS_Disc_Additional_Param const* pagp =
   dynamic_cast<FS_Disc*>(m_geometric_immersed_body)
      ->get_ptr_FS_Disc_Additional_Param();

  size_t temp = (size_t) ((1./surface_cell_scale)
               *(2.*MAC::pi()*pagp->radius)
               /(dx));

  // Getting the nearest even number
  Ntot = (size_t) (round((double)temp * 0.5) * 2.);

}




//---------------------------------------------------------------------------
void DS_2DCircular::compute_surface_points()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_2DCircular:: compute_surface_points");

  // Pointers to location and additional parameters
  struct FS_Disc_Additional_Param const* pagp =
   dynamic_cast<FS_Disc*>(m_geometric_immersed_body)
      ->get_ptr_FS_Disc_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_RigidBody*>(m_geometric_immersed_body)
                            ->get_ptr_to_gravity_centre();


  size_t Npoints = m_all_nodes.size();
  double d_theta = 2.*MAC::pi()/((double)Npoints);
  double theta = 0.5*d_theta;

  for (size_t i = 0; i < Npoints; i++) {
     theta = theta + d_theta;

     geomVector point( pagp->radius*MAC::cos(theta)
                     , pagp->radius*MAC::sin(theta)
                     , 0. );
      
     // Create surface point
     m_all_nodes[i].position[0] = point(0);
     m_all_nodes[i].position[1] = point(1);
     m_all_nodes[i].position[2] = point(2);
     //  m_surface_area[i]->operator()(0) = pagp->radius*d_theta;

     // Create surface normal vectors
     m_all_nodes[i].normal[0] = point(0);
     m_all_nodes[i].normal[1] = point(1);
     m_all_nodes[i].normal[2] = point(2);
  }

  // Translate and rotate
  for (size_t i = 0; i < m_all_nodes.size(); i++) {
//     m_geometric_rigid_body->rotate(m_surface_points[i]);
//     m_geometric_rigid_body->rotate(m_surface_normal[i]);
     m_all_nodes[i].position[0] += pgc->operator()(0);
     m_all_nodes[i].position[1] += pgc->operator()(1);
     m_all_nodes[i].position[2] += pgc->operator()(2);
     std::cout << m_all_nodes[i].position[0] << "," 
               << m_all_nodes[i].position[1] << "," 
               << m_all_nodes[i].position[2] << "," << endl; 
  }


}
