#include <DS_ImmersedBoundary.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <MAC.hh>
using std::endl;


//---------------------------------------------------------------------------
DS_ImmersedBoundary:: DS_ImmersedBoundary()
//---------------------------------------------------------------------------
  : m_geometric_immersed_body( NULL )
{
  MAC_LABEL( "DS_ImmersedBoundary:: DS_ImmersedBoundary" ) ;



}




//---------------------------------------------------------------------------
DS_ImmersedBoundary::DS_ImmersedBoundary(FS_RigidBody *pgrb)
//---------------------------------------------------------------------------
  : m_geometric_immersed_body( pgrb )
{
  MAC_LABEL( "DS_ImmersedBoundary:: DS_ImmersedBoundary" ) ;



}




//---------------------------------------------------------------------------
DS_ImmersedBoundary:: ~DS_ImmersedBoundary()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_ImmersedBoundary:: ~DS_ImmersedBoundary" ) ;


}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: write_one_IB_to_VTK( double const& time
                                              , size_t const& cyclenum )
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: write_one_IB_to_VTK()");





}



//---------------------------------------------------------------------------
void DS_ImmersedBoundary::initialize_surface_variables()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: initialize_surface_variables");


  if (m_all_nodes.empty()) {
     m_all_nodes.reserve( Ntot );
     struct Node nnn;

     nnn.position.reserve(3);
     nnn.position.push_back(0.);
     nnn.position.push_back(0.);
     nnn.position.push_back(0.);

     nnn.velocity.reserve(3);
     nnn.velocity.push_back(0.);
     nnn.velocity.push_back(0.);
     nnn.velocity.push_back(0.);

     nnn.force.reserve(3);
     nnn.force.push_back(0.);
     nnn.force.push_back(0.);
     nnn.force.push_back(0.);

     nnn.normal.reserve(3);
     nnn.normal.push_back(0.);
     nnn.normal.push_back(0.);
     nnn.normal.push_back(0.);

     for (size_t i = 0; i < Ntot; ++i) m_all_nodes.push_back(nnn);
  }
}