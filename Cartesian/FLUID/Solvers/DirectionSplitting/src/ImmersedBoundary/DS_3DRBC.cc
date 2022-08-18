#include <DS_3DRBC.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
using std::endl;
using std::cout;
using std::cin;
using std::string;


//---------------------------------------------------------------------------
DS_3DRBC:: DS_3DRBC()
//---------------------------------------------------------------------------
  : DS_ImmersedBoundary()
{
  MAC_LABEL( "DS_3DRBC:: DS_3DRBC" ) ;

}




//---------------------------------------------------------------------------
DS_3DRBC:: ~DS_3DRBC()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: ~DS_3DRBC" ) ;

}




//---------------------------------------------------------------------------
void DS_3DRBC:: initialize_node_properties( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: initialize_node_properties()" ) ;

  m_all_nodes.reserve(shape_param.N_nodes);

  Node temp;
  temp.coordinates(3);
  temp.coordinates_pbc(3);
  temp.velocity(3);
  temp.angular_velocity(3);
  temp.sumforce(3);
  temp.sumforce_nm1(3);
  temp.spring_force(3);
  temp.bending_force(3);
  temp.viscous_force(3);
  temp.volume_force(3);
  temp.area_force(3);
  temp.unit_outwards_normal_vector(3);
  temp.initial_angle = 0.;
  temp.angle_nm1 = 0.;
  temp.dangle_dt = 0.;
  temp.number = 0;

  for (size_t i = 0; i < shape_param.N_nodes; i++) {
    m_all_nodes.push_back(temp);
  }

}




//---------------------------------------------------------------------------
void DS_3DRBC:: write_one_point_to_VTK( double const& time
                                   , size_t const& cyclenum )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: write_one_point_to_VTK()" ) ;

}




//---------------------------------------------------------------------------
void DS_3DRBC:: generate_membrane_mesh()
//---------------------------------------------------------------------------
{
  
}
