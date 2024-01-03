#include <DS_3DRBC.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
using std::endl;


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
void DS_3DRBC:: compute_number_of_surface_variables(
                                  double const& surface_cell_scale
                                , double const& dx)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_number_of_surface_variables" ) ;

  // struct FS_Disc_Additional_Param const* pagp =
  //  dynamic_cast<FS_Disc*>(m_geometric_rigid_body)
  //     ->get_ptr_FS_Disc_Additional_Param();

  // size_t temp = (size_t) ((1./surface_cell_scale)
  //              *(2.*MAC::pi()*pagp->radius)
  //              /(dx));

  // // Getting the nearest even number
  // Ntot = (size_t) (round((double)temp * 0.5) * 2.);

}




//---------------------------------------------------------------------------
void DS_3DRBC::compute_surface_points()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_3DRBC:: compute_surface_points");



}