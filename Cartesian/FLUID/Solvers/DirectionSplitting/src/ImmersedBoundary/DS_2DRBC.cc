#include <DS_2DRBC.hh>
#include <FV_Mesh.hh>
#include <FV_DiscreteField.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <doubleArray2D.hh>
#include <math.h>
#include <cmath>
#include <typeinfo>
#include <iomanip>      // std::setprecision
#include <fstream>
#include <sstream>
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
void DS_2DRBC:: initialize_node_properties(string const& mesh_filename,
                                           size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: initialize_node_properties()" ) ;

  m_all_nodes.reserve(shape_param.N_nodes);

  Node temp;
  temp.coordinates(dim);
  temp.coordinates_pbc(dim);
  temp.velocity(dim);
  temp.angular_velocity(dim);
  temp.sumforce(dim);
  temp.sumforce_nm1(dim);
  temp.spring_force(dim);
  temp.bending_force(dim);
  temp.viscous_force(dim);
  temp.volume_force(dim);
  temp.area_force(dim);
  temp.unit_outwards_normal_vector(dim);
  temp.initial_angle = 0.;
  temp.angle_nm1 = 0.;
  temp.dangledt = 0.;
  temp.number = 0;

  for (size_t i = 0; i < shape_param.N_nodes; ++i) {
    m_all_nodes.push_back(temp);
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: set_all_nodes(istream& fileIN, size_t const& dim)
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
void DS_2DRBC:: initialize_triangle_properties(size_t const& num_triangles,
                                               size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: initialize_triangle_properties()" ) ;
}




//---------------------------------------------------------------------------
void DS_2DRBC:: set_all_trielements(istream& fileIN, size_t const& dim,
                                    bool const& MatlabNumb)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: set_all_trielements" ) ;
}




//---------------------------------------------------------------------------
void DS_2DRBC:: compute_triangle_area_normals_centre_of_mass(bool init,
                                                    size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_triangle_area_normals_centre_of_mass" ) ;
}




//---------------------------------------------------------------------------
void DS_2DRBC:: set_all_node_neighbors()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: set_all_node_neighbors" ) ;
}




//---------------------------------------------------------------------------
void DS_2DRBC:: initialize_edge_properties(size_t const& dim)
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
void DS_2DRBC:: set_all_edges(size_t const& dim)
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

  // Edge of neighbors
  for (size_t i=0; i<num_nodes; ++i)
  {
    m_all_nodes[i].edge_of_neighbors[0] = 
                                 &(m_all_edges[i == 0 ? m_nEdges - 1 : i - 1 ]);
    m_all_nodes[i].edge_of_neighbors[1] = &(m_all_edges[i]);
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: compute_spring_lengths(bool init, size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_spring_lengths" ) ;
  
  geomVector diff(2);
  
  for (size_t i=0;i<m_nEdges;++i)
  {
    for (size_t j=0;j<dim;++j)
    {
        diff(j) = m_all_edges[i].n3->coordinates(j) 
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
  
  geomVector diff(2);
  
  for (size_t i=0;i<m_nEdges;++i)
  {
    for (size_t j=0;j<2;++j)
    {
        diff(j) = m_all_edges[i].n3->coordinates(j) 
                  - m_all_edges[i].n2->coordinates(j);
    }

    m_all_edges[i].length = norm( diff );

    m_all_edges[i].ext_unit_normal(0) =   diff(1)
                                          / m_all_edges[i].length;
    m_all_edges[i].ext_unit_normal(1) = - diff(0)
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
void DS_2DRBC:: preprocess_membrane_parameters(string const& model_type
                                           , string const& case_type
                                           , double const& mu
                                           , size_t const& num_subtimesteps_RBC
                                           , size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: preprocess_membrane_parameters" ) ;
  
  if(model_type.compare("NumericalMembraneModel") == 0) // if it is detailed numerical membrane model
  {
    // Nothing to do since NMM is valid for 3D RBC
  }
  else
  {
    // Membrane mass to Node mass
    membrane_param.node_mass = membrane_param.mass / double(shape_param.N_nodes);
    
    // Number of subtimesteps for RBC dynamics iterations
    membrane_param.n_subtimesteps_RBC = num_subtimesteps_RBC;
    
    // Membrane constants to edge & node based constants
    membrane_param.membrane_spring_constant = membrane_param.k_spring;

    // Set the spring constant values
    if(case_type.compare("Breyannis2000case") != 0)
    {
      membrane_param.k_spring *= double(m_nEdges); // edge spring constant
      membrane_param.k_bending *= double(m_nEdges); // node bending constant
      membrane_param.k_bending_visc /= double(m_nEdges);
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
    else
    {
      membrane_param.k_spring = mu * membrane_param.ShearRate * shape_param.radius 
                                / membrane_param.CapillaryNumber;
      membrane_param.membrane_spring_constant = membrane_param.k_spring;
    }
  }
}




//---------------------------------------------------------------------------
// Interpolates Eulerian to Lagrangian velocity
//---------------------------------------------------------------------------
// get imin, imax, jmin, jmax which include ghost cells
// get istart, iend, jstart, jend
// do all_RBCs
//  do all_nodes
//      rbc->velocity[comp] = 0.
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
//              rbc->velocity[comp] += eul_velocity * weight * dxC * dyC
//          end do
//      end if
//  end do
// end do
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
  double r1, p1, delt1; // Dirac delta variables
  size_t istart, iend, jstart, jend;

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
  
  // Eulerian to Lagrangian velocity interpolation for all nodes
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
          // // dist_x = (xC - xp) * hxC;
          double dist_y = 
                compute_dist_incl_pbc(yC, yp, domain_length(1)) * hyC;
          // // dist_y = (yC - yp) * hyC;
          bool eul_cell_within_Dirac_delta_stencil = (fabs(dist_x) <= 2.) 
                                                     and 
                                                     (fabs(dist_y) <= 2.);
                  
          if( eul_cell_within_Dirac_delta_stencil )
          {
            r1 = dist_x;
            r1 = discrete_Dirac_delta(r1, ibm_param.dirac_type, dxC, Nx);
            p1 = dist_y;
            p1 = discrete_Dirac_delta(p1, ibm_param.dirac_type, dyC, Ny);

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
void DS_2DRBC:: lag_to_eul(FV_DiscreteField* FF, FV_DiscreteField* FF_tag,
                           size_t const& dim, size_t const& comp)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: lag_to_eul" ) ;

  double xC, yC, zC;
  double dxC, dyC, dzC;
  double hxC, hyC, hzC; // Reciprocal of dxC, dyC, dzC
  int Nx, Ny, Nz;
  double r1, p1, delt1; // Dirac delta variables
  size_t istart, iend, jstart, jend;

  FV_Mesh const* fvm = FF->primary_grid() ;
  
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

  size_t num_nodes = shape_param.N_nodes;

  // Lagrangian to Eulerian force spreading for all nodes
  for(size_t inode=0; inode<num_nodes; ++inode)
  {
    // Get coordinates of inode's Lagrangian marker
    double xp = m_all_nodes[inode].coordinates_pbc(0);
    double yp = m_all_nodes[inode].coordinates_pbc(1);
    
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
        xC = FF->get_DOF_coordinate( ii, comp, 0 ) ;
        yC = FF->get_DOF_coordinate( jj, comp, 1 ) ;
        
        // Check if Eulerian cell is within processor domain
        bool eul_cell_within_proc_domain = 
                           fvm->is_in_domain_on_current_processor(xC, yC);

        // Check if Eulerian cell is within Dirac delta 2x2 stencil
        double dist_x = compute_dist_incl_pbc(xC, xp, domain_length(0)) * hxC;
        // // dist_x = (xC - xp) * hxC;
        double dist_y = compute_dist_incl_pbc(yC, yp, domain_length(1)) * hyC;
        // // dist_y = (yC - yp) * hyC;
        bool eul_cell_within_Dirac_delta_stencil = (fabs(dist_x) <= 2.) 
                                                   and 
                                                   (fabs(dist_y) <= 2.);
        
        // if( eul_cell_within_proc_domain and eul_cell_within_Dirac_delta_stencil )
        if( eul_cell_within_Dirac_delta_stencil )
        {
          r1 = dist_x;
          r1 = discrete_Dirac_delta(r1, ibm_param.dirac_type, dxC, Nx);
          p1 = dist_y;
          p1 = discrete_Dirac_delta(p1, ibm_param.dirac_type, dyC, Ny);

          // Dirac delta function value
          delt1 = r1 * p1;

          // Numerical integration of Dirac delta function value
          sum_dirac_delta += delt1 * dxC * dyC;
          
          kk = 0;

          // Computing Eulerian force
          double euler_force = FF->DOF_value(ii, jj, kk, comp, 0) 
                               + 
                               m_all_nodes[inode].sumforce(comp)
                               * delt1 * dxC * dyC;
          FF->set_DOF_value( ii, jj, kk, comp, 0, euler_force );
          sum_euler_force += euler_force;
          
          // Assigning Eulerian force tag for each cell for debugging
          double euler_force_tag = FF_tag->DOF_value(ii, jj, kk, comp, 0) + 1.0;
          FF_tag->set_DOF_value(ii, jj, kk, comp, 0, euler_force_tag);
        }
      }
    }
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: compute_spring_force(size_t const& dim, 
                                     double const& spring_constant,
                                     string const& model_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_spring_force" ) ;

  size_t num_nodes = shape_param.N_nodes;
  geomVector Iij(dim);
  double length = 0.;
  
  for (size_t inode=0;inode<num_nodes;++inode)
  {
    // Loop over neighboring nodes
    // k = 0 => forward connected node
    // k = 1 => backward connected node
    for (size_t k=0;k<dim;++k)
    {
      // spring vector
      for (size_t j=0;j<dim;++j)
        Iij(j) = m_all_nodes[inode].neighbors[k]->coordinates(j) 
                 - m_all_nodes[inode].coordinates(j);
          
      // spring length
      length = norm(Iij);
      
      // normalization of unit spring vector for components of spring force
      for (size_t j=0;j<dim;++j)
        Iij(j) /= length;

      // spring force computation = k(l - l_0)
      double initial_spring_length = 
                        m_all_nodes[inode].edge_of_neighbors[k]->initial_length;
      for (size_t j=0;j<dim;++j)
        m_all_nodes[inode].sumforce(j) += spring_constant *
                                    ( length - initial_spring_length ) * Iij(j);
    }
  }
} 




//---------------------------------------------------------------------------
void DS_2DRBC:: compute_linear_spring_force( size_t const& dim, 
                                             double const& spring_constant )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_linear_spring_force" ) ;

  size_t num_nodes = shape_param.N_nodes;
  geomVector Iij(dim);
  double length = 0.;
  
  for (size_t inode=0;inode<num_nodes;++inode)
  {
    // Loop over neighboring nodes
    // k = 0 => forward connected node
    // k = 1 => backward connected node
    for (size_t k=0;k<dim;++k)
    {
      // spring vector
      for (size_t j=0;j<dim;++j)
        Iij(j) = m_all_nodes[inode].neighbors[k]->coordinates(j) 
                 - m_all_nodes[inode].coordinates(j);
          
      // spring length
      length = m_all_nodes[inode].edge_of_neighbors[k]->length;
      length = norm(Iij);
      
      // normalization of unit spring vector for components of spring force
      for (size_t j=0;j<dim;++j)
        Iij(j) /= length;

      // initial spring length
      double initial_spring_length = 
                        m_all_nodes[inode].edge_of_neighbors[k]->initial_length;
      
      // tension along the spring
      double tension = spring_constant 
                       * ( (length / initial_spring_length) - 1. );
      
      // compute spring force
      for (size_t j=0;j<dim;++j)
        m_all_nodes[inode].sumforce(j) += tension * Iij(j);
    }
  }
} 




//---------------------------------------------------------------------------
void DS_2DRBC:: compute_bending_resistance( size_t const& dim, 
                                         double const& bending_spring_constant,
                                         double const& bending_viscous_constant, 
                                         double const& dt,
                                         string const& model_type )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_bending_resistance" ) ;

  size_t num_nodes = shape_param.N_nodes;
  geomVector n1(2), n2(2);
  double f = 0.;
  double top = 0.;
  
  for (size_t inode=0;inode<num_nodes;++inode)
  {
    // Compute the angle, pay attention to the sign determined by 
    // the sign of the cross product
    n1.operator=(m_all_nodes[inode].edge_of_neighbors[0]->ext_unit_normal);
    n2.operator=(m_all_nodes[inode].edge_of_neighbors[1]->ext_unit_normal);
    double scalar_prod = n1.operator,(n2);
    double vector_prod = cross_2D(n1, n2);
    double theta = (vector_prod > 0. ? 1. : -1.) * acos( scalar_prod);

    // Bending moment
    double moment_of_inertia = bending_spring_constant 
                               * ( theta - m_all_nodes[inode].initial_angle );
    moment_of_inertia += bending_viscous_constant * m_all_nodes[inode].dangledt ;

    // Bending spring forces
    // Neighbor 0
    for (size_t j=0;j<dim;++j)
    {
      f = (moment_of_inertia / m_all_nodes[inode].edge_of_neighbors[0]->length) 
          * n1(j);
      m_all_nodes[inode].neighbors[0]->sumforce(j) += f ;
      m_all_nodes[inode].sumforce(j) -= f ;
    }  

    // Neighbor 1
    for (size_t j=0;j<dim;++j)
    {
      f = (moment_of_inertia / m_all_nodes[inode].edge_of_neighbors[1]->length) 
          * n2(j);
      m_all_nodes[inode].neighbors[1]->sumforce(j) += f ;
      m_all_nodes[inode].sumforce(j) -= f ;
      top += f;
    }
    
    // Update previous angle and d(angle)/dt
    m_all_nodes[inode].dangledt = ( theta - m_all_nodes[inode].angle_nm1 ) / dt;
    m_all_nodes[inode].angle_nm1 = theta;
  }
}




//---------------------------------------------------------------------------
// Computes viscous drag force
void DS_2DRBC:: compute_viscous_drag_force( size_t const& dim,
                                           double const& viscous_drag_constant,
                                           double const& dt,
                                           string const& model_type )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_viscous_drag_force" ) ;
    
  size_t num_nodes = shape_param.N_nodes;
  
  // viscous force = - c.v where c is damping constant and v is node velocity
  for (size_t inode=0;inode<num_nodes;++inode)
    for (size_t j=0;j<dim;++j)
        m_all_nodes[inode].sumforce(j) -= viscous_drag_constant 
                                         * m_all_nodes[inode].velocity(j);
} 




//---------------------------------------------------------------------------
void DS_2DRBC:: compute_volume_conservation_force(size_t const& dim,
                                                  string const& model_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_volume_conservation_force" ) ;

}




//---------------------------------------------------------------------------
void DS_2DRBC:: compute_area_conservation_force(size_t const& dim,
                                                bool const& Matlab_numbering,
                                                string const& model_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_area_conservation_force" ) ;

}




//---------------------------------------------------------------------------
void DS_2DRBC:: rbc_dynamics_solver(size_t const& dim, 
                                    size_t const& fluid_iter_num,
                                    double const& dt_fluid, 
                                    string const& case_type,
                                    bool const& Matlab_numbering,
                                    string const& model_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: rbc_dynamics_solver" ) ;

  size_t num_nodes = shape_param.N_nodes;
  double node_mass = membrane_param.node_mass;
  double spring_constant = membrane_param.k_spring;
  double bending_spring_constant = membrane_param.k_bending;
  double bending_viscous_constant = membrane_param.k_bending_visc;
  double viscous_drag_constant = membrane_param.k_viscous;
  size_t n_sub_timesteps = membrane_param.n_subtimesteps_RBC;
  
  double time = 0.;
  double total_kinetic_energy = 0.;
  double dt = dt_fluid / n_sub_timesteps;
  
  // Compute edge normals
  compute_edge_normals();
  
  // Initial perimeter
  membrane_param.initial_perimeter = perimeter();
  

  // Time loop
  for (size_t iter_num=0;iter_num<n_sub_timesteps;++iter_num)
  {
    // Time
    time += dt;
    
    // Initialize forces on all nodes
    for (size_t inode=0;inode<num_nodes;++inode)
      for (size_t j=0;j<dim;++j)
        m_all_nodes[inode].sumforce(j) = 0.0;
            
    // Compute external unit normals to edges and edge length
    compute_edge_normals();
    
    // Compute forces on all nodes
    // Spring force
    if(case_type.compare("Breyannis2000case") != 0)
      compute_spring_force( dim, spring_constant, model_type );
    else
      compute_linear_spring_force( dim, spring_constant );
        
    // Bending resistance
    if(case_type.compare("Breyannis2000case") != 0)
      compute_bending_resistance( dim,
                                  bending_spring_constant, 
                                  bending_viscous_constant, 
                                  dt,
                                  model_type );

    // Viscous drag force
    compute_viscous_drag_force( dim, viscous_drag_constant, dt, model_type );


    for (size_t inode=0;inode<num_nodes;++inode)
    {
      // Solve momentum conservation
      // 1st order Euler explicit method
      if ( !iter_num ) // ==> same as if(iter_num == 0)
      {
        for (size_t j=0;j<dim;++j)
        {
          m_all_nodes[inode].velocity(j) += ( dt /  node_mass ) 
                                            * m_all_nodes[inode].sumforce(j);
        }
      }
      else // 2nd order Adams-Bashforth from the 2nd time
      {
        for (size_t j=0;j<dim;++j)
        {
          m_all_nodes[inode].velocity(j) += ( dt /  node_mass ) 
                                   * ( 1.5 * m_all_nodes[inode].sumforce(j) 
                                   - 0.5 * m_all_nodes[inode].sumforce_nm1(j) );
          m_all_nodes[inode].sumforce_nm1(j) = m_all_nodes[inode].sumforce(j);
        }
      }

      // Update position/coordinates with 
      // 2nd order Taylor series expansion
      // x = x0 + u*dt + (1/2) * a * (dt)^2
      for (size_t j=0;j<dim;++j)
        m_all_nodes[inode].coordinates(j) += dt * m_all_nodes[inode].velocity(j) 
                                        + 0.5 * ( m_all_nodes[inode].sumforce(j) 
                                                 / node_mass ) * pow( dt, 2.0 );

      // Compute total kinetic energy
      for (size_t j=0;j<dim;++j)
        total_kinetic_energy += 0.5 * node_mass * 
                               pow( m_all_nodes[inode].velocity(j), 2.0 );
    }
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: rbc_dynamics_solver_no_sub_time_stepping(size_t const& dim, 
                                                   size_t const& fluid_iter_num,
                                                   double const& dt_fluid, 
                                                   string const& case_type,
                                                   bool const& Matlab_numbering,
                                                   string const& model_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: rbc_dynamics_solver_no_sub_time_stepping" ) ;

}




//---------------------------------------------------------------------------
void DS_2DRBC:: write_mesh_to_vtk_file( size_t IB_number, double const& time,
                                        size_t const& cyclenum )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: write_mesh_to_vtk_file" ) ;

  size_t num_nodes = shape_param.N_nodes;

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




//-----------------------------------------------------------
void DS_2DRBC:: compute_centroid(size_t const& dim)
//-----------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_centroid" ) ;
    
  size_t num_nodes = shape_param.N_nodes;
  
  // centroid or center of mass
  for (size_t j=0;j<dim;++j)
  {
    // re-initialisation to avoid wrong 
    // centroid computation from previous 
    // time step
    membrane_param.centroid_coordinates(j) = 0.0;
    
    for (size_t inode=0;inode<num_nodes;++inode)
      membrane_param.centroid_coordinates(j) += m_all_nodes[inode].coordinates(j);
    
    membrane_param.centroid_coordinates(j) /= num_nodes;
  }
}




//---------------------------------------------------------------------------
double DS_2DRBC:: compute_axial_diameter()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_axial_diameter" ) ;
    
  size_t num_nodes = shape_param.N_nodes;
  
  double dist_axial, euclid_dist_axial = -HUGE_VAL;
  double x_node, y_node;
  size_t node_index;
  double axial_radius, axial_dia;

  double x_centroid = 0.; // membrane_param.centroid_coordinates(0);
  double y_centroid = 0.; // membrane_param.centroid_coordinates(1);
  
  for (size_t inode=0;inode<num_nodes;++inode)
  {
    double x = m_all_nodes[inode].coordinates(0);
    double y = m_all_nodes[inode].coordinates(1);
    
    dist_axial = MAC::sqrt( pow(x - x_centroid, 2.) 
                          + pow(y - y_centroid, 2.) );
    
    if(dist_axial > euclid_dist_axial)
    {
      euclid_dist_axial = dist_axial;
      node_index = inode;
      x_node = x;
      y_node = y;
    }
  }
  
  axial_radius = MAC::sqrt( pow(x_node - x_centroid, 2.) 
                          + pow(y_node - y_centroid, 2.) );
  axial_dia = 2. * axial_radius;
      
  return(axial_dia);
}




//---------------------------------------------------------------------------
double DS_2DRBC:: compute_transverse_diameter()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_transverse_diameter" ) ;
    
  size_t num_nodes = shape_param.N_nodes;
  
  double dist_transverse, euclid_dist_transverse = HUGE_VAL;
  double x_node, y_node;
  size_t node_index;
  double transverse_radius, transverse_dia;

  double x_centroid = 0.; // membrane_param.centroid_coordinates(0);
  double y_centroid = 0.; // membrane_param.centroid_coordinates(1);
  
  for (size_t inode=0;inode<num_nodes;++inode)
  {
    double x = m_all_nodes[inode].coordinates(0);
    double y = m_all_nodes[inode].coordinates(1);
    
    dist_transverse = MAC::sqrt( pow(x - x_centroid, 2.) 
                               + pow(y - y_centroid, 2.) );
    
    if(dist_transverse < euclid_dist_transverse)
    {
      euclid_dist_transverse = dist_transverse;
      node_index = inode;
      x_node = x;
      y_node = y;
    }
  }
  
  transverse_radius = MAC::sqrt( pow(x_node - x_centroid, 2.) 
                               + pow(y_node - y_centroid, 2.) );
  transverse_dia = 2. * transverse_radius;
      
  return(transverse_dia);
}




//-----------------------------------------------------------
void DS_2DRBC:: compute_tdp_orientation_angle()
//-----------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_tdp_orientation_angle" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
    
  double rmax = -HUGE_VAL, rmin = HUGE_VAL, theta = 0.;
  double pi = MAC::pi();
  double radians_to_angle_conversion = 180. / pi;
  
  for (size_t inode=0;inode<num_nodes;++inode)
  {
    double x = m_all_nodes[inode].coordinates(0);
    double y = m_all_nodes[inode].coordinates(1);
    
    double r = MAC::sqrt( pow(x, 2.) + pow(y, 2.) );
    if(r > rmax)
    {
      rmax = r;
      theta = (fabs(x) < 1.e-14) ? (y > 0. ? pi/2. : 3.*pi/2.) : atan2(y, x);
      theta = (theta >= 0.) ? ((theta > pi) ? theta - pi : theta) : theta + 2.*pi;
      theta = (theta > pi) ? theta - pi : theta; // to keep the angle between 0 and 90 degrees
      theta *= radians_to_angle_conversion;
    }
    
    if(r < rmin)
      rmin = r;
  }
  
  membrane_param.axial_diameter = 2. * rmax;
  membrane_param.transverse_diameter = 2. * rmin;
  membrane_param.taylor_deformation_parameter = (rmax - rmin) / (rmax + rmin);
  membrane_param.orientation_angle = theta;
}




//---------------------------------------------------------------------------
double DS_2DRBC:: compute_avg_tangential_velocity()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_avg_tangential_velocity" ) ;

  size_t num_nodes = shape_param.N_nodes;
  
  double avg_tangential_velocity = 0.;
  double u, v;
  
  for (size_t inode=0;inode<num_nodes;++inode)
  {
      u = m_all_nodes[inode].velocity(0);
      v = m_all_nodes[inode].velocity(1);
      
      avg_tangential_velocity += MAC::sqrt( pow(u, 2.) + pow(v, 2.) );
  }
  
  avg_tangential_velocity /= num_nodes;
  
  return(avg_tangential_velocity);
}




//---------------------------------------------------------------------------
void DS_2DRBC:: compute_stats(string const& directory, string const& filename, 
                              size_t const& dim, double const& time, 
                              double const& final_time, size_t const& cyclenum,
                              string const& case_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: compute_stats" ) ;
  
  m_rootdir = "Res";
  m_rootname = "rbc";
  m_video_rootname = "video_rbc";
  m_kinetic_energy_rootname = "kinetic_energy.txt";
  m_morphology_rootname = "membrane_morphology.datnew";
  m_force_stats_rootname = "force_stats.txt";
  m_diameter_stats_rootname = "diameter_and_morphology.txt";
  m_rbc_one_point_rootname = "rbc_one_point";
  m_triangle_unit_normals_rootname = "triangle_unit_normals";
  m_node_unit_normals_rootname = "node_unit_normals";


  ofstream rbc_stats_file;
  string file_to_write = directory + "/" + filename;
    
  // Opening the file
  if(cyclenum == 1)
  {
    rbc_stats_file.open( file_to_write, ios::out );
    rbc_stats_file << "Iteration_number" << "\t"
                    << "Time" << "\t"
                    << "Axial_diameter" << "\t"
                    << "Transverse_diameter" << "\t"
                    << "Taylor_deformation_parameter" << "\t"
                    << "Orientation_angle" << "\t"
                    << "Average_tangential_velocity" << "\t"
                    << "Initial_perimeter" << "\t"
                    << "Final_perimeter" << "\t"
                    << "Final_perimeter-Initial_perimeter" << "\t"
                    << "Final_perimeter/Initial_perimeter" << "\t"
                    << "Initial area" << "\t"
                    << "Final_area" << "\t"
                    << "Centroid_x" << "\t"
                    << "Centroid_y" << "\t"
                    << "Total_kinetic_energy"
                    << endl;
                    
    // Initial perimeter
    compute_edge_normals();
    membrane_param.initial_perimeter = perimeter();
    membrane_param.initial_area = pow(membrane_param.initial_perimeter, 2) 
                                  / (4.0 * MAC::pi());
  }
  else
  {
    rbc_stats_file.open( file_to_write, ios::app );
  }
  
  // Centroid of membrane
  compute_centroid(dim);

  // Axial diameter
  membrane_param.axial_diameter = compute_axial_diameter();
  
  // Transverse diameter
  membrane_param.transverse_diameter = compute_transverse_diameter();
  
  // Taylor deformation parameter and orientation angle computation
  compute_tdp_orientation_angle();
  
  // Final perimeter
  compute_edge_normals();
  double final_perimeter = perimeter();
  
  // Final area
  membrane_param.final_area = pow(final_perimeter, 2) / (4.0 * MAC::pi());
  
  // Average tangential velocity for tank treading = avg magnitude of velocity over all nodes
  membrane_param.avg_tangential_velocity = compute_avg_tangential_velocity();
  
  // Writing to file
  rbc_stats_file << std::scientific // << setprecision(12)
    << cyclenum << "\t"
    << time << "\t"
    << membrane_param.axial_diameter << "\t"
    << membrane_param.transverse_diameter << "\t"
    << membrane_param.taylor_deformation_parameter << "\t"
    << membrane_param.orientation_angle << "\t"
    << membrane_param.avg_tangential_velocity << "\t"
    << membrane_param.initial_perimeter << "\t"
    << membrane_param.final_perimeter << "\t"
    << membrane_param.final_perimeter - membrane_param.initial_perimeter << "\t"
    << membrane_param.final_perimeter / membrane_param.initial_perimeter << "\t"
    << membrane_param.initial_area << "\t"
    << membrane_param.final_area << "\t"
    << membrane_param.centroid_coordinates(0) << "\t"
    << membrane_param.centroid_coordinates(1) << "\t"
    << membrane_param.total_kinetic_energy
    << endl;

  // Closing the file
  rbc_stats_file.close();
}




//---------------------------------------------------------------------------
double DS_2DRBC:: norm( geomVector& v )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: norm" ) ;
    
  return ( pow( v(0)*v(0) + v(1)*v(1), 0.5 ) );
}




//---------------------------------------------------------------------------
double DS_2DRBC:: scalar( geomVector const& v0, geomVector const& v1 )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: scalar" ) ;
    
  return ( v0(0) * v1(0) + v0(1) * v1(1) ); 
}




//---------------------------------------------------------------------------
double DS_2DRBC:: cross_2D( geomVector const v0, geomVector const v1 )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: cross_2D" ) ;

  return ( v0(0) * v1(1) - v0(1) * v1(0) );
} 




//---------------------------------------------------------------------------
double DS_2DRBC:: perimeter()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: perimeter" ) ;

  double perimeter = 0.;
  for (size_t i=0;i<m_nEdges;++i)
      perimeter += m_all_edges[i].length;

  return ( perimeter );  
}
