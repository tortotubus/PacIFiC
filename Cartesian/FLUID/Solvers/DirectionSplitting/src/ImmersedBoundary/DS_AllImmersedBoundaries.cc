#include <DS_AllImmersedBoundaries.hh>
#include <FS_AllRigidBodies.hh>
#include <DS_ImmersedBoundary.hh>
#include <DS_ImmersedBoundary_BuilderFactory.hh>
#include <FV_TimeIterator.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <cmath>
using std::endl;


//---------------------------------------------------------------------------
DS_AllImmersedBoundaries:: DS_AllImmersedBoundaries()
//---------------------------------------------------------------------------
  : m_space_dimension( 2 )
  , m_nIB( 0 )
{
  MAC_LABEL( "DS_AllImmersedBoundaries:: DS_AllImmersedBoundaries" ) ;



}




//---------------------------------------------------------------------------
DS_AllImmersedBoundaries::DS_AllImmersedBoundaries(size_t &dimens
                                , istream &in
                                , bool const &b_IB_as_fixed_obstacles
                                , FV_DiscreteField const *arb_UF
                                , double const& arb_scs 
                                , MAC_Communicator const *arb_macCOMM)
//---------------------------------------------------------------------------
  : m_space_dimension(dimens)
  , UF(arb_UF)
  , MESH ( UF->primary_grid() )
  , surface_cell_scale (arb_scs)
  , m_macCOMM(arb_macCOMM)
{
  MAC_LABEL( "DS_AllImmersedBoundaries:: DS_AllImmersedBoundaries()" ) ;

  m_FSallrigidbodies = new FS_AllRigidBodies( m_space_dimension, in,
                        b_IB_as_fixed_obstacles );
  m_nIB = m_FSallrigidbodies->get_number_rigid_bodies();
  m_allDSimmersedboundaries.reserve(m_nIB);

  DS_ImmersedBoundary *dsib = NULL;

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundaries.push_back(dsib);
    m_allDSimmersedboundaries[i] = DS_ImmersedBoundary_BuilderFactory
                        ::create(m_FSallrigidbodies->get_ptr_rigid_body(i));
  }

  initialize_surface_variables_for_all_IB();
}




//---------------------------------------------------------------------------
DS_AllImmersedBoundaries:: ~DS_AllImmersedBoundaries()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundaries:: ~DS_AllImmersedBoundaries" ) ;

  for (size_t i = 0; i < m_nIB; ++i) 
      delete m_allDSimmersedboundaries[i];
  m_allDSimmersedboundaries.clear();

}




//---------------------------------------------------------------------------
size_t DS_AllImmersedBoundaries:: get_number_immersed_boundaries() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundaries:: get_number_immersed_boundaries" ) ;

  return ( m_nIB );

}




//---------------------------------------------------------------------------
DS_ImmersedBoundary* DS_AllImmersedBoundaries:: get_ptr_rigid_body( size_t i )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundaries:: get_ptr_rigid_body" ) ;

  return ( m_allDSimmersedboundaries[i] );

}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundaries::initialize_surface_variables_for_all_IB()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllImmersedBoundaries:: initialize_surface_variables_for_all_IB");

  double dx = MESH->get_smallest_grid_size();

  for (size_t i = 0; i < m_nIB; ++i) {
  // for (vector<size_t>::iterator it = local_RB_list.begin();
  //      it != local_RB_list.end(); ++it) {
    // size_t i = *it;
    m_allDSimmersedboundaries[i]->compute_number_of_surface_variables(
        surface_cell_scale, dx);
    m_allDSimmersedboundaries[i]->initialize_surface_variables();
    m_allDSimmersedboundaries[i]->compute_surface_points();
  }


}
