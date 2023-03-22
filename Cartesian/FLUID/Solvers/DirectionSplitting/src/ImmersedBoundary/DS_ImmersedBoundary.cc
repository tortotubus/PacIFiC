#include <DS_ImmersedBoundary.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <FV_TimeIterator.hh>
#include <MAC.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <boolVector.hh>
#include <doubleVector.hh>
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
            << membrane_param.k_bending_visc << "\t"
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
void DS_ImmersedBoundary:: update_membrane_coordinates(size_t const& dim,
                                                  double const& dt)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: update_membrane_coordinates" ) ;
  
  size_t num_nodes = shape_param.N_nodes;

  for (size_t inode=0;inode<num_nodes;++inode)
    for (size_t j=0;j<dim;++j)
      m_all_nodes[inode].coordinates(j) += m_all_nodes[inode].velocity(j) * dt;
}
  
//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: do_one_inner_iteration
                           (FV_DiscreteField const* UF
                          , FV_DiscreteField* Eul_F
                          , FV_DiscreteField* F_Eul_tag
                          , FV_TimeIterator const* t_it
                          , FV_Mesh const* MESH
                          , size_t const& dim
                          , boolVector const* is_periodic
                          , string const& case_type
                          , bool const& Matlab_numbering
                          , bool const& combined_force_computation
                          , string const& model_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: do_one_inner_iteration" ) ;
  
  size_t num_nodes = shape_param.N_nodes;

  MAC_Communicator const* MAC_comm;
  MAC_comm = MAC_Exec::communicator();
  size_t my_rank = MAC_comm->rank();
  size_t nb_procs = MAC_comm->nb_ranks();
  size_t is_master = 0;
  
  // Allocate memory of temporary vectors
  size_t n1 = dim;
  size_t n2 = 2 * dim;
  doubleVector temp_lag_vel(n1 * num_nodes, 0.); // Lagrangian velocity
  doubleVector temp_lag_pos_and_force(n2 * num_nodes, 0.); // Lagrangian position & force
  
  // Apply periodic boundary conditions
  apply_periodic_boundary_conditions(MESH, dim, is_periodic);
  
  // Eulerian velocity to Lagrangian velocity interpolation
  size_t nb_comps = UF->nb_components();
  for (size_t comp=0;comp<nb_comps;comp++)
  {
    eul_to_lag(UF, dim, comp);
  }

  // MPI Reduce the Lagrangian velocity across all processors
  copy_lagrangian_velocity_to_vector(temp_lag_vel, dim);
  MAC_comm->reduce_vector(temp_lag_vel, 0);
  
  // Dynamics & deformation of the RBC
  if(my_rank == is_master)
  {
    // Copy MPI_Reduce'd Lagrangian velocity to master proc velocity variable
    copy_vector_to_lagrangian_velocity(temp_lag_vel, dim);
    
    // Solve for RBC dynamics using spring-dashpot model
    if(int(membrane_param.n_subtimesteps_RBC) > 0) // With sub time stepping
    {
      rbc_dynamics_solver(dim, t_it->iteration_number(), 
                          t_it->time_step(), case_type, 
                          Matlab_numbering, combined_force_computation, 
                          model_type);
    }
    else // no sub time stepping
    {
      update_membrane_coordinates(dim, t_it->time_step());
      
      rbc_dynamics_solver_no_sub_time_stepping(dim, t_it->iteration_number(),
                                               t_it->time_step(), case_type, 
                                               Matlab_numbering, 
                                               combined_force_computation, 
                                               model_type);
    }
    
    // Copy new Lagrangian positon & force into a doubleVector
    copy_lag_position_and_force_to_vector(temp_lag_pos_and_force, dim);
  }
  
  // Broadcast the Lagrangian position & force from master to all processors
  MAC_comm->broadcast(temp_lag_pos_and_force, 0);
  copy_vector_to_lag_position_and_force(temp_lag_pos_and_force, dim);
  
  // Apply periodic boundary conditions
  apply_periodic_boundary_conditions(MESH, dim, is_periodic);
  
  // Lagrangian to Eulerian force spreading
  // Initialising all Eulerian force cell tag value to zero
  for (size_t comp=0;comp<nb_comps;comp++)
  {
    Eul_F->set_DOFs_value(comp, 0, 0.0);
    F_Eul_tag->set_DOFs_value(comp, 0, 0.0);
  }
  for (size_t comp=0;comp<nb_comps;comp++) 
  {
    // // lag_to_eul(Eul_F, F_Eul_tag, dim, comp);
  }
}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: apply_periodic_boundary_conditions
                             ( FV_Mesh const* MESH
                             , size_t const& dim
                             , boolVector const* is_periodic)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: apply_periodic_boundary_conditions" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
  
  geomVector domain_length(dim, 0.);
  
  bool apply_periodic_bc = (is_periodic->operator()(0))
                           or
                           (is_periodic->operator()(1))
                           or
                           (is_periodic->operator()(2));
                           
  if(apply_periodic_bc)
  {
    // Assign current coordinates to coordinates_pbc variable
    for (size_t i=0;i<num_nodes;++i)
      for (size_t dir=0;dir<dim;++dir)
          m_all_nodes[i].coordinates_pbc(dir) = m_all_nodes[i].coordinates(dir);
    
    // Apply periodic boundary conditions
    for (size_t dir=0;dir<dim;++dir)
    {
      if(is_periodic->operator()(dir))
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
        {
          for (size_t i=0;i<num_nodes;++i)
          {
            m_all_nodes[i].coordinates(dir) -= 
                             MAC::floor(
                             (m_all_nodes[i].coordinates(dir) - min)
                             /domain_length(dir)) 
                             * domain_length(dir);
          }
        }
      }
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
void DS_ImmersedBoundary:: copy_lagrangian_velocity_to_vector
                                      (doubleVector& lag_vel, size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_ImmersedBoundary:: copy_lagrangian_velocity_to_vector" ) ;

  size_t num_nodes = shape_param.N_nodes;

  size_t j = 0;
  for(size_t inode=0; inode<num_nodes; ++inode)
  {
    for(size_t dir=0; dir<dim; ++dir)
    {
      lag_vel(j) = m_all_nodes[inode].velocity(dir);
      j++;
    }
  }
}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: copy_vector_to_lagrangian_velocity
                                      (doubleVector& lag_vel, size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_ImmersedBoundary:: copy_vector_to_lag_vel" ) ;

  size_t num_nodes = shape_param.N_nodes;
  
  size_t j = 0;
  for(size_t inode=0; inode<num_nodes; ++inode)
  {
    for(size_t dir=0; dir<dim; ++dir)
    {
        m_all_nodes[inode].velocity(dir) = lag_vel(j);
        j++;
    }
  }
}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: copy_lag_position_and_force_to_vector
                            (doubleVector& lag_pos_and_force, size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_ImmersedBoundary:: copy_lag_position_and_force_to_vector" ) ;

  size_t num_nodes = shape_param.N_nodes;
  
  size_t j = 0;
  for(size_t inode=0; inode<num_nodes; ++inode)
  {
    // storing coordinates
    for(size_t dir=0; dir<dim; ++dir)
    {
        lag_pos_and_force(j) = m_all_nodes[inode].coordinates(dir);
        j++;
    }
    // storing Lagrangian force
    for(size_t dir=0; dir<dim; ++dir)
    {
        lag_pos_and_force(j) = m_all_nodes[inode].sumforce(dir);
        j++;
    }
  }
}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: copy_vector_to_lag_position_and_force
                            (doubleVector& lag_pos_and_force, size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_ImmersedBoundary:: copy_vector_to_lag_position_and_force" ) ;

  size_t num_nodes = shape_param.N_nodes;
  
  size_t j = 0;
  for(size_t inode=0; inode<num_nodes; ++inode)
  {
    // storing coordinates
    for(size_t dir=0; dir<dim; ++dir)
    {
        m_all_nodes[inode].coordinates(dir) = lag_pos_and_force(j);
        j++;
    }
    // storing Lagrangian force
    for(size_t dir=0; dir<dim; ++dir)
    {
        m_all_nodes[inode].sumforce(dir) = lag_pos_and_force(j);
        j++;
    }
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
    
  if(dirac_type == "Balogh")
  {
      // Same as Eggleton & Popel, Phy. of Fluids, 1998
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
  else if(dirac_type == "Archer")
  {
      // Uses more or less the equation 7 in 2011_Kruger_IBMBasics.pdf
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
  
  double distance = fabs(p1 - length - p2) > 0.5 * length 
                ? 
                p1 + length - p2 
                : 
                p1 - length - p2;

  return(distance);
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




//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: write_rbc_dot_pvd_file(size_t IB_number)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_ImmersedBoundary:: write_rbc_dot_pvd_file" ) ;
  
  string filename = "rbc" + sizetToString( IB_number ) + ".pvd";
  string rootname = "Res/";
  string file_to_write = rootname + filename;
  ofstream rbc_pvd_file( file_to_write, ios::out );
  
  // initialize pvd
  m_rbc_pvd << "<?xml version=\"1.0\"?>" << endl;
  m_rbc_pvd << "<VTKFile type=\"Collection\" version=\"0.1\""
            << " byte_order=\"LittleEndian\"";
  m_rbc_pvd << ">" << endl;
  m_rbc_pvd << "<Collection>" << endl;
  
  // write time steps
  m_rbc_pvd << m_vtk_to_pvd.str();
            
  // finalize pvd
  m_rbc_pvd << "</Collection>" << endl;
  m_rbc_pvd << "</VTKFile>" << endl;

  // Write to rbc.pvd file inside Res/ directory
  rbc_pvd_file << m_rbc_pvd.str();
  
  // clear the instance so that next iteration of rbc.pvd
  // writing is free of previous iteration data
  m_rbc_pvd.str(""); // resets string to empty
  m_rbc_pvd.clear(); // clear any error flags

  rbc_pvd_file.close();
}






