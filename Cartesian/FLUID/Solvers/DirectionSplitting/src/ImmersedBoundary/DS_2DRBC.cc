#include <DS_2DRBC.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
#include <cmath>
using std::endl;
using std::cout;
using std::cin;
using std::string;


//---------------------------------------------------------------------------
DS_2DRBC:: DS_2DRBC()
//---------------------------------------------------------------------------
  : DS_ImmersedBoundary()
{
  MAC_LABEL( "DS_2DRBC:: DS_2DRBC" ) ;

}




//---------------------------------------------------------------------------
DS_2DRBC:: ~DS_2DRBC()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: ~DS_2DRBC" ) ;

}




//---------------------------------------------------------------------------
void DS_2DRBC:: write_one_point_to_VTK( double const& time
                                   , size_t const& cyclenum )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: write_one_point_to_VTK()" ) ;

}




//---------------------------------------------------------------------------
void DS_2DRBC:: initialize_node_properties()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: initialize_node_properties()" ) ;

  m_all_nodes.reserve(shape_param.N_nodes);

  Node temp;
  temp.coordinates(2);
  temp.coordinates_pbc(2);
  temp.velocity(2);
  temp.angular_velocity(2);
  temp.sumforce(2);
  temp.sumforce_nm1(2);
  temp.spring_force(2);
  temp.bending_force(2);
  temp.viscous_force(2);
  temp.volume_force(2);
  temp.area_force(2);
  temp.unit_outwards_normal_vector(2);
  // temp.neighbors(2);
  temp.initial_angle = 0.;
  temp.angle_nm1 = 0.;
  temp.dangle_dt = 0.;
  temp.number = 0;

  for (size_t i = 0; i < shape_param.N_nodes; ++i) {
    m_all_nodes.push_back(temp);
  }

}




//---------------------------------------------------------------------------
void DS_2DRBC:: generate_membrane_mesh()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: generate_membrane_mesh" ) ;

  size_t num_nodes = shape_param.N_nodes;

  // Generate the node ID and coordinates
  for (size_t i=0; i<num_nodes; ++i)
  {
    // Node number or node ID
    m_all_nodes[i].number = i;

    // Coordinates
    m_all_nodes[i].coordinates(0) = shape_param.radius
                            * cos( 2. * M_PI * double(i) / double(num_nodes) ) ;
    m_all_nodes[i].coordinates(1) = shape_param.radius
                            * sin( 2. * M_PI * double(i) / double(num_nodes) ) ;
  }

  // Set the neighbor for each node


}
