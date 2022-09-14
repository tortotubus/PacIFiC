#include <DS_3DRBC.hh>
#include <FV_Mesh.hh>
#include <FV_DiscreteField.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <doubleArray2D.hh>
#include <math.h>
#include <cmath>
#include <typeinfo>
#include <iomanip>      // std::setprecision
using std::endl;
using std::cout;
using std::cin;
using std::string;
using std::max;
using namespace std;


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
void DS_3DRBC:: initialize_node_properties( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: initialize_node_properties()" ) ;

  m_all_nodes.reserve(shape_param.N_nodes);

  Node temp;
  temp.coordinates(3);
  temp.coordinates_pbc(3);
  temp.velocity(3);
  temp.angular_velocity(3);
  temp.sumforce(3);
  temp.sumforce_nm1(3);
  temp.spring_force(3);
  temp.bending_force(3);
  temp.viscous_force(3);
  temp.volume_force(3);
  temp.area_force(3);
  temp.unit_outwards_normal_vector(3);
  // // temp.neighbors.resize(3);
  // // temp.neighbors.push_back(0);
  // // temp.neighbors.push_back(0);
  // // temp.neighbors.push_back(0);
  temp.initial_angle = 0.;
  temp.angle_nm1 = 0.;
  temp.dangledt = 0.;
  temp.number = 0;

  for (size_t i = 0; i < shape_param.N_nodes; i++) {
    m_all_nodes.push_back(temp);
  }

}




//---------------------------------------------------------------------------
void DS_3DRBC:: write_one_point_to_VTK( double const& time
                                   , size_t const& cyclenum )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: write_one_point_to_VTK()" ) ;

}




//---------------------------------------------------------------------------
void DS_3DRBC:: set_all_nodes()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: set_all_nodes" ) ;

  
}




//---------------------------------------------------------------------------
void DS_3DRBC:: set_all_trielements()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: set_all_trielements" ) ;

  
}




//---------------------------------------------------------------------------
void DS_3DRBC:: initialize_edge_properties( )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: initialize_edge_properties" ) ;

}



//---------------------------------------------------------------------------
void DS_3DRBC:: set_all_edges()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: set_all_edges" ) ;

  
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_spring_lengths(bool init)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_spring_lengths" ) ;

  
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_edge_normals()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_edge_normals" ) ;


}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_edge_angle(bool init)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_edge_angle" ) ;


}




//---------------------------------------------------------------------------
void DS_3DRBC:: project_membrane_shape()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: project_membrane_shape" ) ;

  double c0 = shape_param.c0;
  double c1 = shape_param.c1;
  double c2 = shape_param.c2;
  size_t num_nodes = shape_param.N_nodes;
  double radius = shape_param.radius;
  
  // 3D RBC shape from Eq. 14 in Li et al, Biophysical Journal, 2005
  for (size_t i=0;i<num_nodes;++i)
  {
      double x = m_all_nodes[i].coordinates(0);
      double y = m_all_nodes[i].coordinates(1);
      double ww = ( pow( x, 2.) + pow( y, 2. ) ) / pow( radius, 2. );
      double z = ( m_all_nodes[i].coordinates(2) > 0. ? 1. : -1. ) * radius 
                 * sqrt( max( 1. - ww, 0. ) ) 
                 * ( c0 + c1 * ww + c2 * pow( ww, 2.) ) ;
      m_all_nodes[i].coordinates(2) = z;
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC:: preprocess_membrane_parameters(string const& case_type,
                                             double const& mu,
                                             size_t const& num_subtimesteps_RBC)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: preprocess_membrane_parameters" ) ;
    
    
    
}




//---------------------------------------------------------------------------
void DS_3DRBC:: write_mesh_to_vtk_file( size_t IB_number, double const& time,
                                        size_t const& cyclenum )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: write_mesh_to_vtk_file()" ) ;

  size_t num_nodes = shape_param.N_nodes;

  /*
  // File name
  ofstream fileOUT;
  string filename = "rbc" + sizetToString( IB_number ) + "_T" 
                    + sizetToString( cyclenum ) + ".vtu";
  string directory = "Res/";
  string file_to_write = directory + filename.c_str();
  fileOUT.open( file_to_write, ios::out );

  // Add a line to pvd oss
  m_vtk_to_pvd << "<DataSet timestep=\"" << time
          	   << "\" " << "group=\"\" part=\"0\" file=\"" 
          	   << filename << "\"/>\n";

  // Header
  fileOUT << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" " 
          << "byte_order=\"LittleEndian\">" << endl;
  fileOUT << "<UnstructuredGrid>" << endl;

  // Number of nodes and triangular elements         	  
  fileOUT << "<Piece NumberOfPoints=\"" 
          << num_nodes 
          << "\" NumberOfCells=\"" 
          << m_nTriangles 
          << "\">" 
          << endl; 

  // Write node coordinates
  fileOUT << "<Points>" << endl;
  fileOUT << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" " 
          << "format=\"ascii\">" << endl;
  for (size_t i=0;i<num_nodes;++i)
  {
      fileOUT << std::scientific << std::setprecision(12) 
              << m_all_nodes[i].coordinates(0) 
              << " " 
              << m_all_nodes[i].coordinates(1) 
              << " " 
              << m_all_nodes[i].coordinates(2) << endl;
  }
  fileOUT << "</DataArray>" << endl;
  fileOUT << "</Points>" << endl;

  // Write triangular elements 
  fileOUT << "<Cells>" << endl;
  fileOUT << "<DataArray type=\"Int32\" Name=\"connectivity\" "
  << "format=\"ascii\">" << endl;
  for (size_t i=0;i<m_nTriangles;++i)
      fileOUT << m_all_trielements[i].pnodes(0)->number 
              << " " 
              << m_all_trielements[i].pnodes(1)->number 
              << " " 
              << m_all_trielements[i].pnodes(2)->number 
              << " ";
  fileOUT << endl;  
  fileOUT << "</DataArray>" << endl;
  fileOUT << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" 
          << endl;
  size_t offset = 3;
  for (size_t i=0;i<m_nTriangles;++i)
  {
      fileOUT << offset << " ";
      offset += 3;
  }
  fileOUT << endl;   
  fileOUT << "</DataArray>" << endl;
  fileOUT << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" 
  << endl;
  for (size_t i=0;i<m_nTriangles;++i) 
  {
      fileOUT << "5 ";
  }
  fileOUT << endl;   
  fileOUT << "</DataArray>" << endl;
  fileOUT << "</Cells>" << endl;
  fileOUT << "</Piece>" << endl;
  fileOUT << "</UnstructuredGrid>" << endl;
  fileOUT << "</VTKFile>" << endl;

  fileOUT.close();   
  */ 
}



    
//---------------------------------------------------------------------------
// Interpolates Eulerian to Lagrangian velocity
//---------------------------------------------------------------------------
// get imin, imax, jmin, jmax which include ghost cells
// get istart, iend, jstart, jend
// do all_RBCs
//  do all_nodes
//      rbc->velocity(comp) = 0.
//      get xp, yp
//      if((xp, yp) belongs to current proc)
//          compute ipi, ipf, jpi, jpf - bounds of stencil
//          do cells_in_stencil(ii, jj)
//              get xC, yC, dxC, dyC
//              r1 = (xC - xp) * hxC
//              r1 = dirac(r1, dx, Nx)
//              p1 = (yC - yp) * hyC
//              p1 = dirac(p1, dy, Ny)
//              weight = r1 * p1
//              rbc->velocity(comp) += eul_velocity * weight * dxC * dyC
//          end do
//      end if
//  end do
// end do
//---------------------------------------------------------------------------
void DS_3DRBC:: eul_to_lag(FV_DiscreteField const* FF
                         , size_t const& dim
                         , size_t const& comp)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: eul_to_lag" ) ;

  double xC, yC, zC;
  double dxC, dyC, dzC;
  double hxC, hyC, hzC;
  int Nx, Ny, Nz;
  double r1, r2, p1, p2, q1, q2, delt1, delt2, delt;
  size_t istart, iend, jstart, jend, kstart, kend;

  MAC_Communicator const* PAC_comm;
  PAC_comm = MAC_Exec::communicator();
  size_t my_rank = PAC_comm->rank();
  size_t nb_procs = PAC_comm->nb_ranks();
  size_t is_master = 0;
  
  FV_Mesh const* fvm = FF->primary_grid();
  
  size_t_vector min_unknown_index_with_halozone(dim, 0);
  size_t_vector max_unknown_index_with_halozone(dim, 0);
  size_t_vector min_unknown_index_without_halozone(dim, 0);
  size_t_vector max_unknown_index_without_halozone(dim, 0);
  doubleVector Dmin(dim), Dmax(dim);
  doubleVector domain_min(dim), domain_max(dim);
  doubleVector domain_length(dim);
  
  // Get proc's "LOCAL" indices
  for (size_t l=0;l<dim;++l) 
  {
      // getting boundary indices of procs WITH ghost cells/halozone cells
      min_unknown_index_with_halozone(l) = 
                                   FF->get_min_index_unknown_on_proc( comp, l );
      max_unknown_index_with_halozone(l) = 
                                   FF->get_max_index_unknown_on_proc( comp, l );

      // getting boundary indices of procs WITHOUT ghost cells/halozone cells
      min_unknown_index_without_halozone(l) = 
                           FF->get_min_index_unknown_handled_by_proc( comp, l );
      max_unknown_index_without_halozone(l) = 
                           FF->get_max_index_unknown_handled_by_proc( comp, l );

      // getting coordinates of domain bounds
      Dmin(l) = fvm->get_min_coordinate_on_current_processor(l);
      Dmax(l) = fvm->get_max_coordinate_on_current_processor(l);
      domain_min(l) = fvm->get_main_domain_min_coordinate(l);
      domain_max(l) = fvm->get_main_domain_max_coordinate(l);
      domain_length(l) = domain_max(l) - domain_min(l);
  }
  
  // x-direction start & end proc indices
  istart = min_unknown_index_with_halozone(0);
  iend = max_unknown_index_with_halozone(0);

  // y-direction start & end proc indices
  jstart = min_unknown_index_with_halozone(1);
  jend = max_unknown_index_with_halozone(1);
  
  kstart = min_unknown_index_with_halozone(2);
  kend = max_unknown_index_with_halozone(2);
  
  // Number of cells along each Cartesian direction
  Nx = iend - istart + 1;
  Ny = jend - jstart + 1;
  Nz = kend - kstart + 1;

  // Mesh spacing
  dxC = FF->get_cell_size( istart, 0, 0 ) ;
  dyC = FF->get_cell_size( jstart, 1, 1 ) ;
  dzC = FF->get_cell_size( kstart, 2, 2 ) ;
  hxC = 1.0 / dxC ;
  hyC = 1.0 / dyC ;
  hzC = 1.0 / dzC ;

  // Number of nodes on each RBC
  size_t num_nodes = shape_param.N_nodes;
  
  // loop over all Lagrangian marker in iRBC'th particle
  for(size_t inode=0; inode<num_nodes; ++inode)
  {
    // Initialising Lagrangian velocity to 0.0
    m_all_nodes[inode].velocity(comp) = 0.0;

    // Get coordinates of inode's Lagrangian marker
    double xp = m_all_nodes[inode].coordinates_pbc(0);
    double yp = m_all_nodes[inode].coordinates_pbc(1);
    double zp = m_all_nodes[inode].coordinates_pbc(2);
    
    // Condition to check if a Lagrangian node belongs to a proc
    bool lag_cell_within_processor = 
                             fvm->is_in_domain_on_current_processor(xp, yp, zp);
    // // // bool eul_cell_in_domain = (xp > Dmin(0)) && 
    // // //                           (xp <= Dmax(0)) && 
    // // //                           (yp > Dmin(1)) && 
    // // //                           (yp <= Dmax(1));
    
    // if a node belongs to a processor, then compute it's
    // Lagrangian velocity from neighbouring Eulerian nodes
    if(lag_cell_within_processor)
    {
      double sum_dirac_delta = 0.0;
      
      for (size_t ii=min_unknown_index_with_halozone(0);
                  ii<=max_unknown_index_with_halozone(0);
                  ++ii) 
      {
        for (size_t jj=min_unknown_index_with_halozone(1);
                    jj<=max_unknown_index_with_halozone(1);
                    ++jj) 
        {
          for (size_t kk=min_unknown_index_with_halozone(2);
                      kk<=max_unknown_index_with_halozone(2);
                      ++kk)
          {
            xC = FF->get_DOF_coordinate( ii, comp, 0 ) ;
            yC = FF->get_DOF_coordinate( jj, comp, 1 ) ;
            zC = FF->get_DOF_coordinate( kk, comp, 2 ) ;
            
            // Check if Eulerian cell is within Dirac delta 2x2 stencil
            double dist_x = 
                compute_dist_incl_pbc(xC, xp, domain_length(0)) * hxC;
            dist_x = (xC - xp) * hxC;
            double dist_y = 
                compute_dist_incl_pbc(yC, yp, domain_length(1)) * hyC;
            dist_y = (yC - yp) * hyC;
            double dist_z = 
                compute_dist_incl_pbc(zC, zp, domain_length(2)) * hzC;
            dist_z = (zC - zp) * hzC;
            bool eul_cell_within_Dirac_delta_stencil = 
                                                 (fabs(dist_x) <= 2.) 
                                             and (fabs(dist_y) <= 2.) 
                                             and (fabs(dist_z) <= 2.);
            
            if( eul_cell_within_Dirac_delta_stencil )
            {
                r1 = dist_x;
                r1 = discrete_Dirac_delta(r1, ibm_param.dirac_type, dxC, Nx);
                p1 = dist_y;
                p1 = discrete_Dirac_delta(p1, ibm_param.dirac_type, dyC, Ny);
                q1 = dist_z;
                q1 = discrete_Dirac_delta(q1, ibm_param.dirac_type, dzC, Nz);

                // Dirac delta function value
                delt1 = r1 * p1 * q1;

                // Numerical integration of Dirac delta function value
                sum_dirac_delta += delt1 * dxC * dyC * dzC;

                // Computing Lagrangian velocity
                m_all_nodes[inode].velocity(comp) += 
                                  FF->DOF_value( ii, jj, kk, comp, 0 ) 
                                  * delt1 * dxC * dyC * dzC;
            }
          }
        }
      }
    }
  }
}





//---------------------------------------------------------------------------
// Lagrangian to Eulerian force spreading
//---------------------------------------------------------------------------
// do all_RBCs
//  do all_nodes
//      get imin, imax, jmin, jmax which include ghost cells
//      compute ipi, ipf, jpi, jpf - bounds of stencil
//      get xp, yp
//      do cells_in_stencil(ii, jj)
//          get xC, yC
//          if(cell belongs to current proc)
//              get xp, yp
//              get xC, yC, dxC, dyC
//              r1 = (xC - xp) * hxC
//              r1 = dirac(r1, dx, Nx)
//              p1 = (yC - yp) * hyC
//              p1 = dirac(p1, dy, Ny)
//              weight = r1 * p1
//              euler_force += rbc->lag_force * weight * dxC * dyC
//          end if
//      end do
//  end do
// end do
//---------------------------------------------------------------------------
void DS_3DRBC:: lag_to_eul(FV_DiscreteField* FF, FV_DiscreteField* FF_tag,
                           size_t const& dim, size_t const& comp)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: lag_to_eul" ) ;
    
  double xC, yC, zC;
  double dxC, dyC, dzC;
  double hxC, hyC, hzC; // Reciprocal of dxC, dyC, dzC
  int Nx, Ny, Nz;
  double r1, r2, p1, p2, q1, q2, delt1, delt2, delt; // Dirac delta variables
  size_t istart, iend, jstart, jend, kstart, kend;
  double euler_force; // temporary Eulerian force summation variable
  double euler_force_tag; // temporary Eulerian force tag for each cell for debugging

  FV_Mesh const* fvm = FF->primary_grid() ;

  size_t_vector min_unknown_index_with_halozone(dim,0);
  size_t_vector max_unknown_index_with_halozone(dim,0);
  size_t_vector min_unknown_index_without_halozone(dim, 0);
  size_t_vector max_unknown_index_without_halozone(dim, 0);
  doubleVector Dmin(dim), Dmax(dim);
  doubleVector domain_min(dim), domain_max(dim);
  doubleVector domain_length(dim);
  
  // Get proc's "LOCAL" indices
  for (size_t l=0;l<dim;++l) 
  {
    // getting boundary indices of procs WITH ghost cells/halozone cells
    min_unknown_index_with_halozone(l) = 
                                 FF->get_min_index_unknown_on_proc( comp, l );
    max_unknown_index_with_halozone(l) = 
                                 FF->get_max_index_unknown_on_proc( comp, l );

    // getting boundary indices of procs WITHOUT ghost cells/halozone cells
    min_unknown_index_without_halozone(l) = 
                         FF->get_min_index_unknown_handled_by_proc( comp, l );
    max_unknown_index_without_halozone(l) = 
                         FF->get_max_index_unknown_handled_by_proc( comp, l );

    // getting coordinates of domain bounds
    Dmin(l) = fvm->get_min_coordinate_on_current_processor(l);
    Dmax(l) = fvm->get_max_coordinate_on_current_processor(l);
    domain_min(l) = fvm->get_main_domain_min_coordinate(l);
    domain_max(l) = fvm->get_main_domain_max_coordinate(l);
    domain_length(l) = domain_max(l) - domain_min(l);
  }
  
  // x-direction start & end proc indices
  istart = min_unknown_index_with_halozone(0);
  iend = max_unknown_index_with_halozone(0);

  // y-direction start & end proc indices
  jstart = min_unknown_index_with_halozone(1);
  jend = max_unknown_index_with_halozone(1);
  
  kstart = min_unknown_index_with_halozone(2);
  kend = max_unknown_index_with_halozone(2);
  
  // Number of cells along each Cartesian direction
  Nx = iend - istart + 1;
  Ny = jend - jstart + 1;
  Nz = kend - kstart + 1;

  // Mesh spacing
  dxC = FF->get_cell_size( istart, 0, 0 ) ;
  dyC = FF->get_cell_size( jstart, 1, 1 ) ;
  dzC = FF->get_cell_size( kstart, 2, 2 ) ;
  hxC = 1.0 / dxC ;
  hyC = 1.0 / dyC ;
  hzC = 1.0 / dzC ;

  size_t num_nodes = shape_param.N_nodes;

  for(size_t inode=0; inode<num_nodes; ++inode)
  {
    // Get coordinates of inode's Lagrangian marker
    double xp = m_all_nodes[inode].coordinates_pbc(0);
    double yp = m_all_nodes[inode].coordinates_pbc(1);
    double zp = m_all_nodes[inode].coordinates_pbc(2);
    
    double sum_dirac_delta = 0.;
    double sum_euler_force = 0., euler_force_tag = 0.;
    
    for (size_t ii=min_unknown_index_with_halozone(0);
                ii<=max_unknown_index_with_halozone(0);
                ++ii) 
    {
      for (size_t jj=min_unknown_index_with_halozone(1);
                  jj<=max_unknown_index_with_halozone(1);
                  ++jj) 
      {
        for (size_t kk=min_unknown_index_with_halozone(2);
                    kk<=max_unknown_index_with_halozone(2);
                    ++kk) 
        {
          xC = FF->get_DOF_coordinate( ii, comp, 0 ) ;
          yC = FF->get_DOF_coordinate( jj, comp, 1 ) ;
          zC = FF->get_DOF_coordinate( kk, comp, 2 ) ;
          
          // Check if Eulerian cell is within processor domain
          bool eul_cell_within_proc_domain = 
                     fvm->is_in_domain_on_current_processor(xC, yC, zC);

          // Check if Eulerian cell is within Dirac delta 2x2 stencil
          double dist_x = 
                  compute_dist_incl_pbc(xC, xp, domain_length(0)) * hxC;
          dist_x = (xC - xp) * hxC;
          double dist_y = 
                  compute_dist_incl_pbc(yC, yp, domain_length(1)) * hyC;
          dist_y = (yC - yp) * hyC;
          double dist_z = 
                  compute_dist_incl_pbc(zC, zp, domain_length(2)) * hzC;
          dist_z = (zC - zp) * hzC;
          bool eul_cell_within_Dirac_delta_stencil = 
                                                   (fabs(dist_x) <= 2.) 
                                                   and 
                                                   (fabs(dist_y) <= 2.) 
                                                   and 
                                                   (fabs(dist_z) <= 2.);
          
          // // if( eul_cell_within_proc_domain 
          // //     and 
          // //     eul_cell_within_Dirac_delta_stencil )
          if( eul_cell_within_Dirac_delta_stencil )
          {
            r1 = dist_x;
            r1 = discrete_Dirac_delta(r1, ibm_param.dirac_type, dxC, Nx);
            p1 = dist_y;
            p1 = discrete_Dirac_delta(p1, ibm_param.dirac_type, dyC, Ny);
            q1 = dist_z;
            q1 = discrete_Dirac_delta(q1, ibm_param.dirac_type, dzC, Nz);

            // Dirac delta function value
            delt1 = r1 * p1 * q1;

            // Numerical integration of Dirac delta function value
            sum_dirac_delta += delt1 * dxC * dyC * dzC;
            
            // Computing Eulerian force
            euler_force = FF->DOF_value(ii, jj, kk, comp, 0) 
                          + 
                          m_all_nodes[inode].sumforce(comp) 
                          * delt1 * dxC * dyC * dzC;
            FF->set_DOF_value( ii, jj, kk, comp, 0, euler_force );
            sum_euler_force += euler_force;
            
            // Assigning Eulerian force tag for each cell for debugging
            euler_force_tag = FF_tag->DOF_value(ii, jj, kk, comp, 0) 
                              + 1.0;
            FF_tag->set_DOF_value(ii, jj, kk, comp, 0, 
                                  euler_force_tag);
          }
        }
      }
    }
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_spring_force(size_t const& dim, 
                                     double const& spring_constant)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_spring_force" ) ;
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_linear_spring_force(size_t const& dim, 
                                     double const& spring_constant)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_linear_spring_force" ) ;
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_bending_resistance( size_t const& dim, 
                                        double const& bending_spring_constant,
                                        double const& bending_viscous_constant, 
                                        double const& dt )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_bending_resistance" ) ;
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_viscous_drag_force( size_t const& dim,
                                           double const& viscous_drag_constant )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_viscous_drag_force" ) ;
    
}




//---------------------------------------------------------------------------
void DS_3DRBC:: rbc_dynamics_solver(size_t const& dim, 
                                    double const& dt_fluid, 
                                    string const& case_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: rbc_dynamics_solver" ) ;

  size_t num_nodes = shape_param.N_nodes;
  
  

}




//-----------------------------------------------------------
void DS_3DRBC:: compute_centroid(size_t const& dim)
//-----------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_centroid" ) ;
}




//---------------------------------------------------------------------------
double DS_3DRBC:: compute_axial_diameter()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_axial_diameter" ) ;
}




//---------------------------------------------------------------------------
double DS_3DRBC:: compute_transverse_diameter()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_transverse_diameter" ) ;
}




//-----------------------------------------------------------
void DS_3DRBC:: compute_tdp_orientation_angle()
//-----------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_tdp_orientation_angle" ) ;
}




//---------------------------------------------------------------------------
double DS_3DRBC:: compute_avg_tangential_velocity()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_avg_tangential_velocity" ) ;
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_stats(string const& directory, string const& filename, 
                              size_t const& dim, double const& time, 
                              size_t const& cyclenum)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_stats" ) ;

  m_rootdir = "Res";
  m_rootname = "rbc";
  m_video_rootname = "video_rbc";
  m_kinetic_energy_rootname = "kinetic_energy.res";
  m_morphology_rootname = "membrane_morphology.datnew";
  m_force_stats_rootname = "force_stats.txt";
  m_diameter_stats_rootname = "diameter_and_morphology.txt";
  m_rbc_one_point_rootname = "rbc_one_point";
  m_triangle_unit_normals_rootname = "triangle_unit_normals";
  m_node_unit_normals_rootname = "node_unit_normals";


}




//---------------------------------------------------------------------------
double DS_3DRBC:: norm( double const* v )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: 3D_norm" ) ;
    
  return ( pow( v[0]*v[0] + v[1]*v[1] + v[2]*v[2], 0.5 ) );
}



    
//---------------------------------------------------------------------------
double DS_3DRBC:: scalar( double const* v0, double const* v1 )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: scalar" ) ;
    
  return ( v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2] ); 
}



    
//---------------------------------------------------------------------------
void DS_3DRBC:: cross_3D( double const* v0, double const* v1, double* res )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: cross" ) ;
    
  res[0] = v0[1] * v1[2] - v0[2] * v1[1];
  res[1] = v0[2] * v1[0] - v0[0] * v1[2];
  res[2] = v0[0] * v1[1] - v0[1] * v1[0];
} 




//---------------------------------------------------------------------------
double DS_3DRBC:: perimeter()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: perimeter" ) ;

  double perimeter = 0.;
  for (size_t i=0;i<m_nEdges;++i)
      perimeter += m_all_edges[i].length;

  return ( perimeter );  
}
