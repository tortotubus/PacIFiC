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
void DS_2DCircular::compute_surface_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_2DCircular:: compute_surface_parameters");

  // Pointers to location and additional parameters
  struct FS_Disc_Additional_Param const* pagp =
   dynamic_cast<FS_Disc*>(m_geometric_immersed_body)
      ->get_ptr_FS_Disc_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_RigidBody*>(m_geometric_immersed_body)
                            ->get_ptr_to_gravity_centre();


  size_t Npoints = m_all_nodes.size();
  size_t Nedges = Npoints;
  double d_theta = 2.*MAC::pi()/((double)Npoints);
  double theta = 0.5*d_theta;

  for (size_t i = 0; i < Npoints; i++) {
     theta = theta + d_theta;

     geomVector point( pagp->radius*MAC::cos(theta)
                     , pagp->radius*MAC::sin(theta)
                     , 0. );
      
     // Store ID
     m_all_nodes[i]->nodeID = i;

     // Create surface point
     m_all_nodes[i]->position = point;

     // Linking to neighbor nodes
     m_all_nodes[i]->neighbor[0] = (i != 0) ? m_all_nodes[i-1] 
                                            : m_all_nodes[Npoints - 1];
     m_all_nodes[i]->neighbor[1] = (i != Npoints - 1) ? m_all_nodes[i+1]
                                                      : m_all_nodes[0];

     // Create surface normal vectors
     m_all_nodes[i]->normal = point;
  }

  for (size_t i = 0; i < Nedges; i++) {
     m_all_edges[i]->edgeID = i;
     m_all_edges[i]->connecting_node[0] = m_all_nodes[i];
     m_all_edges[i]->connecting_node[1] = (i != Npoints - 1) ? m_all_nodes[i + 1]
                                                             : m_all_nodes[0];
     double dist = m_all_edges[i]->connecting_node[0]->position.calcDist(
                    m_all_edges[i]->connecting_node[1]->position);

     m_all_edges[i]->initial_length = dist;
     m_all_edges[i]->length = dist;
  }



  // Translate and rotate
  for (size_t i = 0; i < m_all_nodes.size(); i++) {
//     m_geometric_rigid_body->rotate(m_surface_points[i]);
//     m_geometric_rigid_body->rotate(m_surface_normal[i]);
     m_all_nodes[i]->position(0) += pgc->operator()(0);
     m_all_nodes[i]->position(1) += pgc->operator()(1);
     m_all_nodes[i]->position(2) += pgc->operator()(2);
  }


}
