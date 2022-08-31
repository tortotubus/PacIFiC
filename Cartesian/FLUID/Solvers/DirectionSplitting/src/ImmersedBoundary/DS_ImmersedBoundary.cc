#include <DS_ImmersedBoundary.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <FV_TimeIterator.hh>
#include <MAC.hh>
#include <fstream>
#include <sstream>
#include <math.h>
#include <cmath>
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
MembraneParameters* DS_ImmersedBoundary:: get_ptr_membrane_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_ImmersedBoundary:: get_ptr_membrane_parameters" ) ;

  return(&membrane_param);

}




//---------------------------------------------------------------------------
IBMParameters* DS_ImmersedBoundary:: get_ptr_IBM_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_ImmersedBoundary:: get_ptr_IBM_parameters" ) ;

  return(&ibm_param);

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
            << shape_param.node_spacing_with_dx << "\t"
            << membrane_param.k_spring << "\t"
            << membrane_param.k_bending << "\t"
            << membrane_param.k_viscous << "\t"
            << membrane_param.k_area << "\t"
            << membrane_param.k_volume << endl;
}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: position_membrane()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: position_membrane" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
  
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
void DS_ImmersedBoundary:: do_one_inner_iteration
                           (FV_TimeIterator const* t_it
                          , FV_Mesh const* MESH
                          , size_t const& dim
                          , size_t const& periodic_dir)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: do_one_inner_iteration" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
  
  apply_periodic_boundary_conditions(MESH, dim, periodic_dir);
  eul_to_lag();
  /*
  doubleVector temp_lag_vel = copy_lag_velocity_to_vector(num_nodes);
  MAC_comm->reduce_vector(temp_lag_vel, 0);
  if(my_rank == is_master)
  {
    copy_vector_to_lag_vel(num_nodes);
    rbc_dynamics();
    copy_lag_position_and_force_to_vector(num_nodes);
  }
  copy_vector_to_lag_position_and_force(num_nodes);
  apply_periodic_boundary_conditions(MESH, dim, periodic_dir);
  lag_to_eul();
  */
}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: apply_periodic_boundary_conditions
                             ( FV_Mesh const* MESH
                             , size_t const& dim
                             , size_t const& periodic_dir)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: apply_periodic_boundary_conditions" ) ;
  
  size_t num_nodes = shape_param.N_nodes;

  // Assign current coordinates to coordinates_pbc variable
  for (size_t i=0;i<num_nodes;++i)
    for (size_t dir=0;dir<dim;++dir)
        m_all_nodes[i].coordinates_pbc(dir) = m_all_nodes[i].coordinates(dir);
  
  geomVector domain_length(dim);
  bool apply_periodic_bc = (periodic_dir >= 0) and (periodic_dir <= 2);
  if(apply_periodic_bc)
  {
      // Apply periodic boundary conditions
      for (size_t dir=0;dir<dim;++dir)
      {
        // Get min and max bounds of domain length
        double min = MESH->get_main_domain_min_coordinate(dir);
        double max = MESH->get_main_domain_max_coordinate(dir);
        domain_length(dir) = max - min;

        // Check if all nodes in immersed body moved 
        // out of periodic domain boundary?
        size_t num = 0;
        for (size_t i=0;i<num_nodes;++i)
        {
            m_all_nodes[i].coordinates_pbc(dir) -= 
                        MAC::floor(
                        (m_all_nodes[i].coordinates_pbc(dir) - min)
                        /domain_length(dir)) 
                        * domain_length(dir);
            
            // check if all RBC Lagrangian nodes are out of boundary
            if( (m_all_nodes[i].coordinates(dir) > max)
                 or 
                 (m_all_nodes[i].coordinates(dir) < min) )
                num += 1;
        }
        
        // apply periodic boundary conditions to membrane
        // coordinates array only when all Lagrangian 
        // nodes are out of the domain boundary
        if(num == num_nodes)
            for (size_t i=0;i<num_nodes;++i)
                m_all_nodes[i].coordinates(dir) -= 
                                 MAC::floor(
                                 (m_all_nodes[i].coordinates(dir) - min)
                                 /domain_length(dir)) 
                                 * domain_length(dir);
      }
  }
  else
  {
    // Do not apply periodic boundary conditions
    for (size_t i=0;i<num_nodes;++i)
      for (size_t dir=0;dir<dim;++dir)
        m_all_nodes[i].coordinates_pbc(dir) = m_all_nodes[i].coordinates(dir);
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




//---------------------------------------------------------------------------
double DS_ImmersedBoundary:: discrete_Dirac_delta ( double val, 
                                    string const& dirac_type, double dh, int n )
//---------------------------------------------------------------------------
{
    MAC_LABEL( "DS_ImmersedBoundary:: discrete_Dirac_delta" ) ;
    
    double abs_val = MAC::abs(val);
    
    if(dirac_type == "Balogh") // Same as Eggleton & Popel, Phy. of Fluids, 1998
    {
        if( ( abs_val >= 0.0 ) && ( abs_val <= 2.0 ) )
            return( 
            ( 1.0 / (4.0 * dh) ) 
            * ( 1.0 + MAC::cos( ( MAC::pi() / (2.0 * dh) ) * (val * dh) ) ) 
            ); // val*dh is done since val is computed as (xC-xp)*hx 
        else
            return 0.;
    }
    else if(dirac_type == "Roma")
    {
        if( abs_val <= 0.5 )
            return( ( 1.0 + MAC::sqrt( -3.0 * pow(val, 2) + 1 ) ) / 3.0 );
        else if( ( abs_val >= 0.5 ) && ( abs_val <= 1.5 ) )
            return(
            ( 5.0 - 3.0 * abs_val 
            - 
            MAC::sqrt( -3.0 * pow( 1.0 - abs_val, 2 ) + 1.0 ) ) / 6.0 
            );
        else
            return 0.;
    }
    else if(dirac_type == "Archer") // Uses more or less the equation 7 in 2011_Kruger_IBMBasics.pdf
    {
        if ( ( abs_val >= 0.0 ) && ( abs_val < 1.0 ) )
            return(
            0.125 
            * ( 3.0 - 2.0 * abs_val 
            + MAC::sqrt( MAC::abs( 1 + 4.0 * abs_val - 4.0 * pow(val, 2) ) ) ) 
            );
        else if ( ( abs_val >= 1.0 ) && ( abs_val < 2.0 ) )
            return(
            0.125 
            * ( 5.0 - 2.0 * abs_val - 
            MAC::sqrt( MAC::abs( -7.0 + 12.0 * abs_val - 4.0 * pow(val, 2) ) ) ) 
            );
        else
            return 0.;
    }
    else if(dirac_type == "Krueger") // Uses the equation 7 in 2011_Kruger_IBMBasics.pdf
    {
        if ( ( abs_val >= 0.0 ) && ( abs_val < 1.0 ) )
            return(
            0.125 
            * ( 3.0 - 2.0 * abs_val 
            + MAC::sqrt( 1 + 4.0 * abs_val - 4.0 * pow(val, 2) ) ) 
            );
        else if ( ( abs_val >= 1.0 ) && ( abs_val < 2.0 ) )
            return(
            0.125 
            * ( 5.0 - 2.0 * abs_val 
            - MAC::sqrt( -7.0 + 12.0 * abs_val - 4.0 * pow(val, 2) ) ) 
            );
        else
            return 0.;
    }
}




//---------------------------------------------------------------------------
bool DS_ImmersedBoundary::across_periodic(double p1, double p2, double length)
//---------------------------------------------------------------------------
{
    MAC_LABEL( "DS_ImmersedBoundary:: across_periodic" ) ;
    
    return(fabs(p1 - p2) > 0.5 * length);
}




//---------------------------------------------------------------------------
double DS_ImmersedBoundary::periodic_1D_distance(double p1, 
                            double p2, double length)
//---------------------------------------------------------------------------
{
    MAC_LABEL( "DS_ImmersedBoundary:: periodic_1D_distance" ) ;
    
    double dist = fabs(p1 - length - p2) > 0.5 * length 
                  ? 
                  p1 + length - p2 
                  : 
                  p1 - length - p2;
    return(dist);
}




//---------------------------------------------------------------------------
double DS_ImmersedBoundary::compute_dist_incl_pbc(double p1, 
                                                    double p2, double length)
//---------------------------------------------------------------------------
{
    MAC_LABEL( "DS_ImmersedBoundary:: compute_dist_incl_pbc" ) ;
    
    double dist = across_periodic(p1, p2, length) 
                  ? 
                  periodic_1D_distance(p1, p2, length) 
                  : 
                  p1 - p2;
    return(dist);
}







