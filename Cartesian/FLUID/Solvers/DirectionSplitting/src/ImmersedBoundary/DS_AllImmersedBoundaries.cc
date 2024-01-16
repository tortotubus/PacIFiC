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
                                , FV_DiscreteField *arb_LF
                                , double const& arb_scs 
                                , MAC_Communicator const *arb_macCOMM)
//---------------------------------------------------------------------------
  : m_space_dimension(dimens)
  , UF(arb_UF)
  , LF(arb_LF)
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

  initialize_all_pvd();
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
    m_allDSimmersedboundaries[i]->compute_surface_parameters();
  }


}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundaries::initialize_all_pvd()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllImmersedBoundaries:: initialize_all_pvd");

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundaries[i]->initialize_pvd();
  }


}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundaries::finalize_all_pvd()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllImmersedBoundaries:: finalize_all_pvd");

  for (size_t i = 0; i < m_nIB; ++i) {
    string strID = m_allDSimmersedboundaries[i]->sizetToString(i);
    m_allDSimmersedboundaries[i]->finalize_pvd("./Res/saveIB" + strID + ".pvd");
  }


}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundaries::project_lagrangian_force_on_eulerian_grid()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllImmersedBoundaries:: project_lagrangian_force_on_eulerian_grid");

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundaries[i]->project_force_on_grid_from_oneIB(LF);
  }

  LF->synchronize(0);

}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundaries::reset_Lagrangian_and_Eulerian_Force_field()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllImmersedBoundaries:: reset_Lagrangian_and_Eulerian_Force_field");

  // Reset Eulerian force field
  for (size_t comp = 0; comp < LF->nb_components(); comp++) {
    LF->set_DOFs_value(comp, 0, 0.);
  }

  // Reset Lagrangian force field
  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundaries[i]->reset_Lagrangian_force_field();
  }

}





//---------------------------------------------------------------------------
void DS_AllImmersedBoundaries::write_all_IB_to_VTU(double const& time
                                              , size_t const& cyclenum)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllImmersedBoundaries:: write_all_IB_to_VTK");

  for (size_t i = 0; i < m_nIB; ++i) {
    string strID = m_allDSimmersedboundaries[i]->sizetToString(i);
    m_allDSimmersedboundaries[i]->write_one_IB_to_VTU("saveIB_" + strID,time,cyclenum);
  }

}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundaries::interpolate_velocity_on_all_IB()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllImmersedBoundaries:: interpolate_velocity_on_all_IB");

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundaries[i]->eulerian_velocity_on_lagrange_nodes(UF, m_macCOMM);
  }

}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundaries::advect_all_IB(double const& dt)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllImmersedBoundaries:: advect_all_IB");

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundaries[i]->advect_IB(dt);
    m_allDSimmersedboundaries[i]->check_and_update_periodic_clone(UF);
    m_allDSimmersedboundaries[i]->update_edge_length(UF);
  }

}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundaries::compute_force_on_all_lagrange_nodes(double const& Es)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_AllImmersedBoundaries:: compute_force_on_all_lagrange_nodes");

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundaries[i]->compute_elastic_force_on_lagrange_nodes(Es);
  }

}
