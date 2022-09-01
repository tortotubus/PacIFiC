#include <DS_2DRBC.hh>
#include <FV_Mesh.hh>
#include <FV_DiscreteField.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <doubleArray2D.hh>
#include <math.h>
#include <cmath>
#include <typeinfo>
using std::endl;
using std::cout;
using std::cin;
using std::string;
using std::max;
using namespace std;


//---------------------------------------------------------------------------
DS_2DRBC:: DS_2DRBC()
//---------------------------------------------------------------------------
  : DS_ImmersedBoundary()
{
  MAC_LABEL( "DS_2DRBC:: DS_2DRBC" ) ;

}




//---------------------------------------------------------------------------
DS_2DRBC:: ~DS_2DRBC()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: ~DS_2DRBC" ) ;

}




//---------------------------------------------------------------------------
void DS_2DRBC:: initialize_node_properties()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: initialize_node_properties()" ) ;

  m_all_nodes.reserve(shape_param.N_nodes);

  Node temp;
  temp.coordinates(2);
  temp.coordinates_pbc(2);
  temp.velocity(2);
  temp.angular_velocity(2);
  temp.sumforce(2);
  temp.sumforce_nm1(2);
  temp.spring_force(2);
  temp.bending_force(2);
  temp.viscous_force(2);
  temp.volume_force(2);
  temp.area_force(2);
  temp.unit_outwards_normal_vector(2);
  temp.initial_angle = 0.;
  temp.angle_nm1 = 0.;
  temp.dangle_dt = 0.;
  temp.number = 0;

  for (size_t i = 0; i < shape_param.N_nodes; ++i) {
    m_all_nodes.push_back(temp);
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: set_all_nodes()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: set_all_nodes" ) ;

  size_t num_nodes = shape_param.N_nodes;
  
  // Generate the node ID and coordinates
  for (size_t i=0; i<num_nodes; ++i)
  {
    // Node number or node ID
    m_all_nodes[i].number = i;

    // Coordinates
    m_all_nodes[i].coordinates(0) = shape_param.radius
                            * cos( 2. * MAC::pi() * double(i) / double(num_nodes) ) ;
    m_all_nodes[i].coordinates(1) = shape_param.radius
                            * sin( 2. * MAC::pi() * double(i) / double(num_nodes) ) ;

    // Node neighbors
    m_all_nodes[i].neighbors[0] = 
                                &(m_all_nodes[i == 0 ? num_nodes - 1 : i - 1 ]);
    m_all_nodes[i].neighbors[1] = 
                                &(m_all_nodes[i == num_nodes - 1 ? 0 : i + 1 ]);
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: project_membrane_shape()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: project_membrane_shape" ) ;

  double c0 = shape_param.c0;
  double c1 = shape_param.c1;
  double c2 = shape_param.c2;
  size_t num_nodes = shape_param.N_nodes;
  double radius = shape_param.radius;
  
  // 2D RBC shape from modifying Eq. 14 in Li et al, Bioph. Journal, 2005
  for (size_t i=0;i<num_nodes;++i)
  {
      double x = m_all_nodes[i].coordinates(0);
      double y = m_all_nodes[i].coordinates(1);
      double ww = pow( x, 2.) / pow( radius, 2. );
      m_all_nodes[i].coordinates(1) = ( y > 0. ? 1. : -1. ) * radius
                                      * sqrt( max( 1. - ww, 0. ) ) 
                                      * ( c0 + c1 * ww + c2 * pow( ww, 2.) );
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: set_all_trielements()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: set_all_trielements" ) ;
}




//---------------------------------------------------------------------------
void DS_2DRBC:: initialize_edge_properties()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: initialize_edge_properties" ) ;
  
  m_nEdges = shape_param.N_nodes;
  
  // Allocate memory for m_all_edges object
  m_all_edges.reserve(m_nEdges);
  Edge temp;
  temp.ext_unit_normal(2);
  temp.initial_length = 0.;
  temp.length = 0.;
  temp.sintheta0 = 0.;
  temp.costheta0 = 0.;
  temp.l0 = 0.;
  temp.l = 0.;
  temp.lmax = 0.;
  temp.k = 0.;
  temp.kp = 0.;
  temp.initial_angle = 0.;
  temp.angle = 0.;
  temp.angle_nm1 = 0.;
  temp.dangledt = 0.;

  for (size_t i=0;i<m_nEdges;++i) m_all_edges.push_back(temp);
}




//---------------------------------------------------------------------------
void DS_2DRBC:: set_all_edges()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: set_all_edges" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
  
  m_nEdges = shape_param.N_nodes;
  
  for (size_t i=0;i<m_nEdges;++i)
  {
    // Edge ID/number
    m_all_edges[i].number = i;
    
    // Nodes forming the edge
    m_all_edges[i].n2 = &(m_all_nodes[i]);
    m_all_edges[i].n3 = &(m_all_nodes[i == m_nEdges - 1 ? 0 : i + 1 ]);

  }

  for (size_t i=0; i<num_nodes; ++i)
  {
    // Edge of neighbors
    m_all_nodes[i].edge_of_neighbors[0] = 
                                 &(m_all_edges[i == 0 ? m_nEdges - 1 : i - 1 ]);
    m_all_nodes[i].edge_of_neighbors[1] = &(m_all_edges[i]);
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: compute_spring_lengths(bool init)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_spring_lengths" ) ;
  
  double diff[2];
  
  for (size_t i=0;i<m_nEdges;++i)
  {
    for (size_t j = 0; j < 2; ++j)
    {
        diff[j] = m_all_edges[i].n3->coordinates(j) 
                  - m_all_edges[i].n2->coordinates(j);
    }
    m_all_edges[i].length = norm( diff );
    if(init) m_all_edges[i].initial_length = m_all_edges[i].length;
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: compute_edge_normals()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_edge_normals" ) ;
  
  double diff[2];
  
  for (size_t i=0;i<m_nEdges;++i)
  {
    for (size_t j = 0; j < 2; ++j)
    {
        diff[j] = m_all_edges[i].n3->coordinates(j) 
                  - m_all_edges[i].n2->coordinates(j);
    }

    m_all_edges[i].length = norm( diff );

    m_all_edges[i].ext_unit_normal(0) =   diff[1]
                                          / m_all_edges[i].length;
    m_all_edges[i].ext_unit_normal(1) = - diff[0]
                                          / m_all_edges[i].length;
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: compute_edge_angle(bool init)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_edge_angle" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
  geomVector n1(2), n2(2);
  
  for (size_t i=0;i<num_nodes;++i)
  {
    n1.operator=(m_all_nodes[i].edge_of_neighbors[0]->ext_unit_normal);
    n2.operator=(m_all_nodes[i].edge_of_neighbors[1]->ext_unit_normal);
    double scalar_prod = n1.operator,(n2);
    double vector_prod = cross_2D(n1, n2);
    
    m_all_nodes[i].angle = (vector_prod > 0. ? 1. : -1.) * acos( scalar_prod ); 
    m_all_nodes[i].angle_nm1 = m_all_nodes[i].angle;
    if(init) m_all_nodes[i].initial_angle = m_all_nodes[i].angle;
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: eul_to_lag(FV_DiscreteField const* FF
                         , size_t const& dim
                         , size_t const& comp)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: eul_to_lag" ) ;

  double xC, yC, zC;
  double dxC, dyC, dzC;
  double hxC, hyC, hzC; // Reciprocal of dxC, dyC, dzC
  int Nx, Ny, Nz;
  double r1, r2, p1, p2, q1, q2, delt1, delt2, delt; // Dirac delta variables
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
  
  size_t kk = 0;
          
  // Number of cells along each Cartesian direction
  Nx = iend - istart + 1;
  Ny = jend - jstart + 1;
  
  // Mesh spacing
  dxC = FF->get_cell_size( istart, 0, 0 ) ;
  dyC = FF->get_cell_size( jstart, 1, 1 ) ;
  hxC = 1.0 / dxC ;
  hyC = 1.0 / dyC ;

  // Number of nodes on each RBC
  size_t num_nodes = shape_param.N_nodes;
  
  // Generate the node ID and coordinates
  for (size_t inode=0; inode<num_nodes; ++inode)
  {
    // Initialising Lagrangian velocity to 0.0
    m_all_nodes[inode].velocity(comp) = 0.;
     
    // Get coordinates of inode's Lagrangian marker
    double xp = m_all_nodes[inode].coordinates_pbc(0);
    double yp = m_all_nodes[inode].coordinates_pbc(1);
      
    // Condition to check if a Lagrangian node belongs to a proc
    bool lag_cell_within_processor = 
                               fvm->is_in_domain_on_current_processor(xp, yp);
      
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
          xC = FF->get_DOF_coordinate( ii, comp, 0 ) ;
          yC = FF->get_DOF_coordinate( jj, comp, 1 ) ;
                  
          // Check if Eulerian cell is within Dirac delta 2x2 stencil
          double dist_x = 
                compute_dist_incl_pbc(xC, xp, domain_length(0)) * hxC;
          double dist_y = 
                compute_dist_incl_pbc(yC, yp, domain_length(1)) * hyC;
          bool eul_cell_within_Dirac_delta_stencil = 
                        (fabs(dist_x) <= 2.) and (fabs(dist_y) <= 2.);
                  
          if( eul_cell_within_Dirac_delta_stencil )
          {
            r1 = dist_x;
            r1 = discrete_Dirac_delta(r1, 
                                      ibm_param.dirac_type, 
                                      dxC, 
                                      Nx);
            p1 = dist_y;
            p1 = discrete_Dirac_delta(p1, 
                                      ibm_param.dirac_type, 
                                      dyC, 
                                      Ny);

            // Dirac delta function value
            delt1 = r1 * p1;

            // Numerical integration of Dirac delta function value
            sum_dirac_delta += delt1 * dxC * dyC;

            kk = 0;

            // Computing Lagrangian velocity
            m_all_nodes[inode].velocity(comp) += 
                                FF->DOF_value( ii, jj, kk, comp, 0 ) 
                                * delt1 * dxC * dyC;
          }
        }
      }
    }
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: copy_lagrangian_velocity_to_vector()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: write_mesh_to_vtk_file()" ) ;

  size_t num_nodes = shape_param.N_nodes;

}




//---------------------------------------------------------------------------
void DS_2DRBC:: write_mesh_to_vtk_file( size_t IB_number, double const& time,
                                        size_t const& cyclenum )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: write_mesh_to_vtk_file()" ) ;

  size_t num_nodes = shape_param.N_nodes;

  // File name
  ofstream fileOUT;
  string filename = "rbc" + sizetToString( IB_number ) + "_T0" + ".vtu";
  string directory = "Res/";
  string file_to_write = directory + filename.c_str();
  fileOUT.open( file_to_write, ios::out );
  
  /*
  // Add a line to pvd oss
  m_vtk_to_pvd << "<DataSet timestep=\"" << time
          	   << "\" " << "group=\"\" part=\"0\" file=\"" 
          	   << filename << "\"/>\n";
  */
  
  
  // Header
  fileOUT << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
  	<< "byte_order=\"LittleEndian\">" << endl;
  fileOUT << "<UnstructuredGrid>" << endl;
  
  // Number of vertices 
  size_t ncell = 1 + num_nodes;           	  
  fileOUT << "<Piece NumberOfPoints=\"" << num_nodes 
  	<< "\" NumberOfCells=\"" << ncell << "\">" << endl; 

  // Write vertex coordinates
  fileOUT << "<Points>" << endl;
  fileOUT << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
  	<< "format=\"ascii\">" << endl;
  for (size_t i=0;i<num_nodes;++i)
    fileOUT << m_all_nodes[i].coordinates(0) << " " 
    	<< m_all_nodes[i].coordinates(1) << " 0." << endl;  
  fileOUT << "</DataArray>" << endl;
  fileOUT << "</Points>" << endl;

  // Write cells = 1 polyline + num_nodes VTK_vertices
  fileOUT << "<Cells>" << endl;
  fileOUT << "<DataArray type=\"Int32\" Name=\"connectivity\" "
  	<< "format=\"ascii\">" << endl;
  for (size_t i=0;i<num_nodes;++i)
    fileOUT << i << " ";          
  fileOUT << "0 "; 
  for (size_t i=0;i<num_nodes;++i)
    fileOUT << i << " ";    
  fileOUT << endl;
  fileOUT << "</DataArray>" << endl;
  fileOUT << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
  	<< endl;
  size_t offset = num_nodes + 1;
  fileOUT << offset;
  offset++;
  for (size_t i=0;i<num_nodes;++i,offset++)
    fileOUT << " " << offset;
  fileOUT << endl;          
  fileOUT << "</DataArray>" << endl;
  fileOUT << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">" 
  	<< endl;
  fileOUT << "4 "; 
  for (size_t i=0;i<num_nodes;++i)
    fileOUT << "1 ";
  fileOUT << endl;          
  fileOUT << "</DataArray>" << endl;
  fileOUT << "</Cells>" << endl;
  fileOUT << "</Piece>" << endl;
  fileOUT << "</UnstructuredGrid>" << endl;
  fileOUT << "</VTKFile>" << endl;
  
  fileOUT.close();   

}




//---------------------------------------------------------------------------
void DS_2DRBC:: preprocess_membrane_parameters(string const& case_type
                                           , size_t const& num_subtimesteps_RBC)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: preprocess_membrane_parameters" ) ;
  
  // Membrane mass to Node mass
  membrane_param.node_mass = membrane_param.mass / double(shape_param.N_nodes);
  
  // Number of subtimesteps for RBC dynamics iterations
  membrane_param.n_subtimesteps_RBC = num_subtimesteps_RBC;
  
  // Build the Direction Splitting immersed boundary
  if(case_type.compare("Breyannis2000case") != 0)
  {
    // Membrane constants to edge & node based constants
    membrane_param.membrane_spring_constant = membrane_param.k_spring;
    membrane_param.k_spring *= double(m_nEdges); // edge spring constant
    membrane_param.k_bending *= double(m_nEdges); // node bending constant
    membrane_param.k_viscous /= double(m_nEdges); // node viscous constant
      
    // Timescales
    membrane_param.membrane_mass_spring_timescale = 2. * MAC::pi() 
                               * pow( membrane_param.mass / 
                               membrane_param.membrane_spring_constant, 0.5 );
    membrane_param.edge_mass_spring_timescale = 2. * MAC::pi() 
                               * pow( membrane_param.node_mass / 
                               membrane_param.k_spring, 0.5 );
    membrane_param.node_bending_mass_spring_timescale = 2. * MAC::pi() 
              * ( 2. * MAC::pi() * shape_param.radius / double(m_nEdges) ) 
              * pow( membrane_param.node_mass / membrane_param.k_bending, 0.5 );
    
    membrane_param.ntimescales = 20.;
    membrane_param.dt = min( membrane_param.edge_mass_spring_timescale, 
                        membrane_param.node_bending_mass_spring_timescale ) 
                        / 50. ;
    membrane_param.tmax = double(membrane_param.ntimescales) 
                          * membrane_param.membrane_mass_spring_timescale;
    membrane_param.ntimesteps = size_t(membrane_param.tmax / membrane_param.dt);
   }
}




//---------------------------------------------------------------------------
double DS_2DRBC:: norm( double const* v )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: norm" ) ;
    
  return ( pow( v[0]*v[0] + v[1]*v[1], 0.5 ) );
}




//---------------------------------------------------------------------------
double DS_2DRBC:: scalar( double const* v0, double const* v1 )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "RBC2D:: scalar" ) ;
    
  return ( v0[0] * v1[0] + v0[1] * v1[1] ); 
}




//---------------------------------------------------------------------------
double DS_2DRBC::cross_2D( geomVector const v0, geomVector const v1 )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "RBC2D:: cross_2D" ) ;

  return ( v0(0) * v1(1) - v0(1) * v1(0) );
} 

