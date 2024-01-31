#include <DS_CommonMethods.hh>
#include <DS_DirectionSplitting.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DomainBuilder.hh>
#include <FV_DiscreteField.hh>
#include <FV_SystemNumbering.hh>
#include <PostProcessing.hh>
#include <FV_Mesh.hh>
#include <FV_TimeIterator.hh>
#include <MAC.hh>
#include <MAC_Root.hh>
#include <MAC_Error.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_Vector.hh>
#include <MAC_BoolArray2D.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <MAC_Application.hh>
#include <PAC_Misc.hh>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <FS_SolidPlugIn_BuilderFactory.hh>
#include <FS_SolidPlugIn.hh>
#include <FS_Grains3DPlugIn.hh>
#include <DS_AllRigidBodies.hh>
#include <DS_AllImmersedBoundaries.hh>

//---------------------------------------------------------------------------
DS_CommonMethods::DS_CommonMethods()
//---------------------------------------------------------------------------
  : dim(2)
{
    MAC_LABEL("DS_CommonMethods:: DS_CommonMethods");
}




//---------------------------------------------------------------------------
DS_CommonMethods::DS_CommonMethods(size_t &dimens
                                , MAC_Communicator const *arb_macCOMM
                                , FV_DiscreteField *arb_TF)
//---------------------------------------------------------------------------
  : dim( dimens )
  , allrigidbodies ( NULL )
  , macCOMM ( arb_macCOMM )
  , MESH (arb_TF->primary_grid())
  , TF (arb_TF)
{
    MAC_LABEL("DS_CommonMethods:: DS_CommonMethods()");

    // Periodic boundary condition check for velocity
    boolVector const* periodic_comp = MESH->get_periodic_directions();
    is_periodic[0] = periodic_comp->operator()(0);
    is_periodic[1] = periodic_comp->operator()(1);
    if (dim > 2) is_periodic[2] = periodic_comp->operator()(2);

}




//---------------------------------------------------------------------------
DS_CommonMethods::DS_CommonMethods(size_t &dimens
                                , MAC_Communicator const *arb_macCOMM
                                , FV_DiscreteField *arb_UF
                                , FV_DiscreteField *arb_PF)
//---------------------------------------------------------------------------
  : dim( dimens )
  , allrigidbodies ( NULL )
  , macCOMM ( arb_macCOMM )
  , MESH (arb_UF->primary_grid())
  , UF (arb_UF)
  , PF (arb_PF)
{
    MAC_LABEL("DS_CommonMethods:: DS_CommonMethods()");

    // Periodic boundary condition check for velocity
    boolVector const* periodic_comp = MESH->get_periodic_directions();
    is_periodic[0] = periodic_comp->operator()(0);
    is_periodic[1] = periodic_comp->operator()(1);
    if (dim > 2) is_periodic[2] = periodic_comp->operator()(2);

}




//---------------------------------------------------------------------------
DS_CommonMethods::DS_CommonMethods(size_t &dimens
                                , DS_AllRigidBodies *arb_allRB
                                , MAC_Communicator const *arb_macCOMM
                                , FV_DiscreteField *arb_TF)
//---------------------------------------------------------------------------
  : dim( dimens )
  , allrigidbodies ( arb_allRB )
  , macCOMM ( arb_macCOMM )
  , MESH (arb_TF->primary_grid())
  , TF (arb_TF)
{
    MAC_LABEL("DS_CommonMethods:: DS_CommonMethods()");

    // Periodic boundary condition check for velocity
    boolVector const* periodic_comp = MESH->get_periodic_directions();
    is_periodic[0] = periodic_comp->operator()(0);
    is_periodic[1] = periodic_comp->operator()(1);
    if (dim > 2) is_periodic[2] = periodic_comp->operator()(2);
}




//---------------------------------------------------------------------------
DS_CommonMethods::DS_CommonMethods(size_t &dimens
                                , DS_AllRigidBodies *arb_allRB
                                , MAC_Communicator const *arb_macCOMM
                                , FV_DiscreteField *arb_UF
                                , FV_DiscreteField *arb_PF)
//---------------------------------------------------------------------------
  : dim( dimens )
  , allrigidbodies ( arb_allRB )
  , macCOMM ( arb_macCOMM )
  , MESH (arb_UF->primary_grid())
  , UF (arb_UF)
  , PF (arb_PF)
{
    MAC_LABEL("DS_CommonMethods:: DS_CommonMethods()");

    // Periodic boundary condition check for velocity
    boolVector const* periodic_comp = MESH->get_periodic_directions();
    is_periodic[0] = periodic_comp->operator()(0);
    is_periodic[1] = periodic_comp->operator()(1);
    if (dim > 2) is_periodic[2] = periodic_comp->operator()(2);
}




//---------------------------------------------------------------------------
DS_CommonMethods::~DS_CommonMethods()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DS_CommonMethods:: ~DS_CommonMethods");


}




//---------------------------------------------------------------------------
void DS_CommonMethods::initialize_grid_nodes_on_rigidbody(FV_DiscreteField *FF
                                          , vector<size_t> const &list)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_CommonMethods::initialize_grid_nodes_on_rigidbody" ) ;

  size_t_vector min_unknown_index(3,0);
  size_t_vector max_unknown_index(3,0);

  // Vector for solid presence
  size_t_array2D* void_frac = allrigidbodies->get_void_fraction_on_grid(FF);
  size_t nb_comps = FF->nb_components();

  for (size_t comp = 0; comp < nb_comps; comp++) {
     // Get local min and max indices
     for (size_t dir = 0; dir < dim; dir++) {
        if (is_periodic[dir]) {
           min_unknown_index(dir) =
                     FF->get_min_index_unknown_handled_by_proc( comp, dir ) - 1;
           max_unknown_index(dir) =
                     FF->get_max_index_unknown_handled_by_proc( comp, dir ) + 1;
        } else {
           min_unknown_index(dir) =
                     FF->get_min_index_unknown_handled_by_proc( comp, dir );
           max_unknown_index(dir) =
                     FF->get_max_index_unknown_handled_by_proc( comp, dir );
        }
     }

     for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i) {
        double xC = FF->get_DOF_coordinate( i, comp, 0 ) ;
        for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j) {
           double yC = FF->get_DOF_coordinate( j, comp, 1 ) ;
           for (size_t k = min_unknown_index(2); k <= max_unknown_index(2); ++k) {
              double zC = (dim == 2) ? 0 : FF->get_DOF_coordinate( k, comp, 2 );
              geomVector pt(xC,yC,zC);
              size_t p = FF->DOF_local_number(i,j,k,comp);
              if (void_frac->operator()(p,0) != 0) {
                 size_t par_id = void_frac->operator()(p,0) - 1;
                 geomVector value = (FF == UF) ? allrigidbodies->rigid_body_velocity(par_id, pt)
                                               : allrigidbodies->rigid_body_temperature(par_id, pt);
                 for (size_t level : list)
                  FF->set_DOF_value(i, j, k, comp, level, value(comp));
              }
           }
        }
     }
  }

}