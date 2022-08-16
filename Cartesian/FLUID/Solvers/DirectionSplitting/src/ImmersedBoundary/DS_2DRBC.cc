#include <DS_2DRBC.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
using std::endl;


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
void DS_2DRBC:: initialize_node_properties( )
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
  temp.number = 0;

  for (size_t i = 0; i < shape_param.N_nodes; i++) {
    m_all_nodes.push_back(temp);
  }

}




//---------------------------------------------------------------------------
void DS_2DRBC:: create_RBC_structure( size_t const& nb_edges
                         , double const& radius
                         , double const& c0
                         , double const& c1
                         , double const& c2
                         , double const& rbc_orientation_angle)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: create_RBC_structure" ) ;

}
