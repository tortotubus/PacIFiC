#include <DS_2DBiconcave.hh>
#include <FS_RigidBody.hh>
#include <FS_Disc.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
using std::endl;


//---------------------------------------------------------------------------
DS_2DBiconcave:: DS_2DBiconcave()
//---------------------------------------------------------------------------
  : DS_ImmersedBoundary()
{
  MAC_LABEL( "DS_2DBiconcave:: DS_2DBiconcave" ) ;

}




//---------------------------------------------------------------------------
DS_2DBiconcave:: DS_2DBiconcave( FS_RigidBody* pgrb )
//---------------------------------------------------------------------------
  : DS_ImmersedBoundary( pgrb )
{
  MAC_LABEL( "DS_2DBiconcave:: DS_2DBiconcave" ) ;

}




//---------------------------------------------------------------------------
DS_2DBiconcave:: ~DS_2DBiconcave()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DBiconcave:: ~DS_2DBiconcave" ) ;

}




//---------------------------------------------------------------------------
void DS_2DBiconcave:: compute_number_of_surface_variables(
                                  double const& surface_cell_scale
                                , double const& dx)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DBiconcave:: compute_number_of_surface_variables" ) ;


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
void DS_2DBiconcave::compute_surface_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_2DBiconcave:: compute_surface_parameters");

  // Pointers to location and additional parameters
  struct FS_Disc_Additional_Param const* pagp =
   dynamic_cast<FS_Disc*>(m_geometric_immersed_body)
      ->get_ptr_FS_Disc_Additional_Param();

  geomVector const* pgc = dynamic_cast<FS_RigidBody*>(m_geometric_immersed_body)
                            ->get_ptr_to_gravity_centre();

  double c0 = 0.1035805;
  double c1 = 1.001279;
  double c2 = -0.561381;

  size_t Npoints = m_all_nodes.size();
  size_t Nedges = Npoints;
  double d_theta = 2.*MAC::pi()/((double)Npoints);
  double theta = 0.01*d_theta;

  for (size_t i = 0; i < Npoints; i++) {
     theta = theta + d_theta;

     double xi = pagp->radius*MAC::cos(theta);
     double yi = pagp->radius * MAC::sqrt(1. - MAC::pow(xi,2.) / MAC::pow(pagp->radius,2.))
               * (c0 + c1 * MAC::pow(xi,2.) / MAC::pow(pagp->radius,2.)
                     + c2 * MAC::pow(xi,4.) / MAC::pow(pagp->radius,4.));

     if (theta > MAC::pi()) yi *= -1.;

     geomVector point( xi
                     , yi
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
    m_geometric_immersed_body->rotate(&m_all_nodes[i]->position);
//     m_geometric_rigid_body->rotate(m_surface_normal[i]);
     m_all_nodes[i]->position(0) += pgc->operator()(0);
     m_all_nodes[i]->position(1) += pgc->operator()(1);
     m_all_nodes[i]->position(2) += pgc->operator()(2);
  }


}
