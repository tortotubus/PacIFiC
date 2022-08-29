#include <DS_ImmersedBoundary.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <MAC.hh>
#include <fstream>
#include <sstream>
using std::endl;
using std::cout;
using std::cin;
using std::string;
using std::ofstream;
using namespace std;


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
void DS_ImmersedBoundary:: display_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "FS_RigidBody:: display_parameters" ) ;

  std::cout << "Shape parameters" << endl;
  std::cout << shape_param.center(0) << "\t"
            << shape_param.center(1) << "\t"
            << shape_param.center(2) << "\t"
            << shape_param.xroll << "\t"
            << shape_param.ypitch << "\t"
            << shape_param.zyaw << "\t"
            << shape_param.radius << "\t"
            << shape_param.c0 << "\t"
            << shape_param.c1 << "\t"
            << shape_param.c2 << "\t"
            << shape_param.N_nodes << "\t"
            << shape_param.N_levels << "\t"
            << shape_param.node_spacing_with_dx << endl;

}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: position_membrane()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: position_membrane" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
  double x_center = shape_param.center(0);
  double y_center = shape_param.center(1);
  double z_center = shape_param.center(2);
  
  for (size_t i=0;i<num_nodes;++i)
    for (size_t dir=0;dir<3;++dir)
        m_all_nodes[i].coordinates(dir) += shape_param.center(dir);
}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: rotate_membrane()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: rotate_membrane" ) ;
  
  doubleArray2D rot_matrix(3,3,0);
  doubleVector coords(3), coords_rotated(3);
  
  size_t dim = 3;
  
  size_t num_nodes = shape_param.N_nodes;

  double roll_angle = shape_param.xroll;
  double pitch_angle = shape_param.ypitch;
  double yaw_angle = shape_param.zyaw;

  double degree_to_radians_conversion = MAC::pi() / 180.;

  double gamma = roll_angle * degree_to_radians_conversion;
  double beta = pitch_angle * degree_to_radians_conversion;
  double alpha = yaw_angle * degree_to_radians_conversion;
  
  // Refer https://en.wikipedia.org/wiki/Rotation_matrix 
  // with alpha = yaw, beta = pitch and gamma = roll
  rot_matrix(0, 0) =   MAC::cos(alpha) * MAC::cos(beta);
  rot_matrix(0, 1) =   MAC::cos(alpha) * MAC::sin(beta) 
                     * MAC::sin(gamma) - MAC::sin(alpha) * MAC::cos(gamma);
  rot_matrix(0, 2) =   MAC::cos(alpha) * MAC::sin(beta) 
                     * MAC::cos(gamma) + MAC::sin(alpha) * MAC::sin(gamma);
  rot_matrix(1, 0) =   MAC::sin(alpha) * MAC::cos(beta);
  rot_matrix(1, 1) =   MAC::sin(alpha) * MAC::sin(beta) 
                     * MAC::sin(gamma) + MAC::cos(alpha) * MAC::cos(gamma);
  rot_matrix(1, 2) =   MAC::sin(alpha) * MAC::sin(beta) 
                     * MAC::cos(gamma) - MAC::cos(alpha) * MAC::sin(gamma);
  rot_matrix(2, 0) = - MAC::sin(beta);
  rot_matrix(2, 1) =   MAC::cos(beta) * MAC::sin(gamma);
  rot_matrix(2, 2) =   MAC::cos(beta) * MAC::cos(gamma);

  for (size_t i=0;i<num_nodes;++i)
  {
    // Position membrane to (0, 0, 0)
    for (size_t dir=0;dir<dim;++dir)
        coords(dir) = m_all_nodes[i].coordinates(dir) 
                      - shape_param.center(dir);
                      
    // Rotate coordinates
    coords_rotated(0) = coords(0)*rot_matrix(0,0) 
                               + coords(1)*rot_matrix(0,1) 
                               + coords(2)*rot_matrix(0,2);
    coords_rotated(1) = coords(0)*rot_matrix(1,0) 
                               + coords(1)*rot_matrix(1,1) 
                               + coords(2)*rot_matrix(1,2);
    coords_rotated(2) = coords(0)*rot_matrix(2,0) 
                               + coords(1)*rot_matrix(2,1) 
                               + coords(2)*rot_matrix(2,2);
        
    // Re-position membrane back to (xcenter, ycenter, zcenter)
    for (size_t dir=0;dir<dim;++dir)
        m_all_nodes[i].coordinates(dir) = coords_rotated(dir) 
                                          + shape_param.center(dir);
  }
}




//---------------------------------------------------------------------------
string DS_ImmersedBoundary:: sizetToString( size_t const& figure ) const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: sizetToString" ) ;

    ostringstream oss;
    oss << figure; 

    return ( oss.str() ); 
}





//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: get_datatype_of_variable()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: get_datatype_of_variable" ) ;

  // cout << typeid(variable).name() << endl;
  // // cout << typeid(m_all_edges[i].ext_unit_normal).name() << endl;
}
