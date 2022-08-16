#include <DS_ImmersedBoundary.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <MAC.hh>
using std::endl;


//---------------------------------------------------------------------------
DS_ImmersedBoundary:: DS_ImmersedBoundary()
//---------------------------------------------------------------------------
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
void DS_ImmersedBoundary:: create_RBC_structure()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_ImmersedBoundary:: create_RBC_structure" ) ;

  // Call the RBC2D or RBC3D create functions

}




//---------------------------------------------------------------------------
ShapeParameters* DS_ImmersedBoundary:: get_ptr_shape_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_ImmersedBoundary:: get_shape_parameters" ) ;

  return(&shape_param);

}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: display_parameters( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: display_parameters" ) ;

  std::cout << "Shape parameters" << endl;
  std::cout << shape_param.center(0) << " , "
            << shape_param.center(1) << " , "
            << shape_param.center(2) << " , "
            << shape_param.radius << " , "
            << shape_param.c0 << " , "
            << shape_param.c1 << " , "
            << shape_param.c2 << " , "
            << shape_param.N_levels << " , "
            << shape_param.N_nodes << endl;

}
