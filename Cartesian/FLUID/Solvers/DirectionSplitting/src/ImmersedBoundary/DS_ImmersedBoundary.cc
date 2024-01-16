#include <DS_ImmersedBoundary.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <MAC.hh>
#include <math.h>
using std::endl;
using std::ofstream;

//---------------------------------------------------------------------------
DS_ImmersedBoundary::DS_ImmersedBoundary()
    //---------------------------------------------------------------------------
    : m_geometric_immersed_body(NULL)
{
  MAC_LABEL( "DS_ImmersedBoundary:: DS_ImmersedBoundary" ) ;



}




//---------------------------------------------------------------------------
DS_ImmersedBoundary::DS_ImmersedBoundary(FS_RigidBody *pgrb)
//---------------------------------------------------------------------------
  : m_geometric_immersed_body( pgrb )
{
  MAC_LABEL( "DS_ImmersedBoundary:: DS_ImmersedBoundary" ) ;



}




//---------------------------------------------------------------------------
DS_ImmersedBoundary:: ~DS_ImmersedBoundary()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_ImmersedBoundary:: ~DS_ImmersedBoundary" ) ;


  if (!m_all_nodes.empty()) {
    for (size_t i = 0; i < m_all_nodes.size(); ++i) {
       delete m_all_nodes[i];
    }
  }

}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary:: write_one_IB_to_VTU( string const& rootname
                                              , double const& time
                                              , size_t const& cyclenum )
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: write_one_IB_to_VTU()");

  // File name
  ofstream fileOUT;
  string filename = "./Res/" + rootname + "_T" + sizetToString(cyclenum) + ".vtu";
  fileOUT.open(filename.c_str(), std::ios::out);

  string filepvd = rootname + "_T" + sizetToString(cyclenum) + ".vtu";
  // Add a line to pvd oss
  m_IB_pvd << "<DataSet timestep=\"" << time
           << "\" "
           << "group=\"\" part=\"0\" file=\"" << filepvd << "\"/>\n";

  // Header
  fileOUT << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
          << "byte_order=\"LittleEndian\">" << endl;
  fileOUT << "<UnstructuredGrid>" << endl;

  // Number of vertices
  size_t ncell = 1 + Ntot;
  fileOUT << "<Piece NumberOfPoints=\"" << Ntot
          << "\" NumberOfCells=\"" << ncell << "\">" << endl;

  // Write vertex coordinates
  fileOUT << "<Points>" << endl;
  fileOUT << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
          << "format=\"ascii\">" << endl;
  for (size_t i = 0; i < Ntot; ++i)
    fileOUT << m_all_nodes[i]->position(0) << " "
            << m_all_nodes[i]->position(1) << " "
            << m_all_nodes[i]->position(2) << " "
            << endl;
  fileOUT << "</DataArray>" << endl;
  fileOUT << "</Points>" << endl;

  fileOUT << "<PointData Vectors=\"U,F\">" << endl;

  // Write velocity (vector field)
  fileOUT << "<DataArray type=\"Float32\" Name=\"U\" NumberOfComponents=\"3\" "
          << "format=\"ascii\">" << endl;
  for (size_t i = 0; i < Ntot; ++i)
    fileOUT << m_all_nodes[i]->velocity(0) << " "
            << m_all_nodes[i]->velocity(1) << " "
            << m_all_nodes[i]->velocity(2) << " "
            << endl;
  fileOUT << "</DataArray>" << endl;
  // Write force (vector field)
  fileOUT << "<DataArray type=\"Float32\" Name=\"F\" NumberOfComponents=\"3\" "
          << "format=\"ascii\">" << endl;
  for (size_t i = 0; i < Ntot; ++i)
    fileOUT << m_all_nodes[i]->force(0) << " "
            << m_all_nodes[i]->force(1) << " "
            << m_all_nodes[i]->force(2) << " "
            << endl;
  fileOUT << "</DataArray>" << endl;

  fileOUT << "</PointData>" << endl;

  // // Write scalar field
  // fileOUT << "<PointData Scalars=\"procID\">" << endl;
  // // Write the local proc ID's
  // fileOUT << "<DataArray type=\"Float32\" Name=\"procID\" "
  //         << "format=\"ascii\">" << endl;
  // for (size_t i = 0; i < Ntot; ++i)
  //   fileOUT << m_all_nodes[i]->local_procID << " "
  //           << endl;
  // fileOUT << "</DataArray>" << endl;
  // fileOUT << "</PointData>" << endl;

  // Write cells = 1 polyline + m_number_of_nodes VTK_bertices
  fileOUT << "<Cells>" << endl;
  fileOUT << "<DataArray type=\"Int32\" Name=\"connectivity\" "
          << "format=\"ascii\">" << endl;
  for (size_t i = 0; i < Ntot; ++i)
    fileOUT << i << " ";
  fileOUT << "0 ";
  for (size_t i = 0; i < Ntot; ++i)
    fileOUT << i << " ";
  fileOUT << endl;
  fileOUT << "</DataArray>" << endl;
  fileOUT << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
          << endl;
  size_t offset = Ntot + 1;
  fileOUT << offset;
  offset++;
  for (size_t i = 0; i < Ntot; ++i, offset++)
    fileOUT << " " << offset;
  fileOUT << endl;
  fileOUT << "</DataArray>" << endl;
  fileOUT << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">"
          << endl;
  fileOUT << "4 ";
  for (size_t i = 0; i < Ntot; ++i)
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
void DS_ImmersedBoundary::initialize_pvd()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: initialize_pvd");

  m_IB_pvd << "<?xml version=\"1.0\"?>" << endl;
  m_IB_pvd <<
	"<VTKFile type=\"Collection\" version=\"0.1\""
	<< " byte_order=\"LittleEndian\"";
  m_IB_pvd << ">" << endl;
  m_IB_pvd << "<Collection>" << endl;

}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary::finalize_pvd(string const& filename)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: initialize_pvd");

  m_IB_pvd << "</Collection>" << endl;
  m_IB_pvd << "</VTKFile>" << endl;
  ofstream f(filename, std::ios::out);
  f << m_IB_pvd.str();
  f.close();
}




//---------------------------------------------------------------------------
string DS_ImmersedBoundary::sizetToString( size_t const &value )
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: sizetToString");

  ostringstream oss;
  oss << value; 
  
  return ( oss.str() );

}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary::advect_IB(double const& dt)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: advect_IB");

  for (size_t i = 0; i < Ntot; ++i) {
    m_all_nodes[i]->position += m_all_nodes[i]->velocity * dt;
  }

}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary::project_force_on_grid_from_oneIB(FV_DiscreteField *LF)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary::project_force_on_grid_from_oneIB");

  geomVector delta(3);

  for (size_t i = 0; i < Ntot; ++i) {
    for (size_t comp = 0; comp < LF->nb_components(); comp++) {
      size_t_vector i0(3,0);
      for (size_t dir = 0; dir < LF->primary_grid()->nb_space_dimensions(); dir++) {
          size_t i_temp = 0;
          FV_Mesh::between(LF->get_DOF_coordinates_vector(comp, dir),
                                        m_all_nodes[i]->position(dir), i_temp);
          i0(dir) = i_temp;
      }

      size_t ix_min = LF->get_min_index_unknown_handled_by_proc(comp, 0);
      size_t ix_max = LF->get_max_index_unknown_handled_by_proc(comp, 0);
      size_t iy_min = LF->get_min_index_unknown_handled_by_proc(comp, 1);
      size_t iy_max = LF->get_max_index_unknown_handled_by_proc(comp, 1);

      for (size_t ix = i0(0)-2; ix <= i0(0)+2; ix++) {
        if (ix >= ix_min && ix <= ix_max) {
          for (size_t iy = i0(1)-2; iy <= i0(1)+2; iy++) {
            if (iy >= iy_min && iy <= iy_max) {
              double dx = m_all_nodes[i]->position(0) 
                        - LF->get_DOF_coordinate(ix, comp, 0);
              double dy = m_all_nodes[i]->position(1) 
                        - LF->get_DOF_coordinate(iy, comp, 1);

              double deltax = LF->get_cell_size(ix, comp, 0);
              double deltay = LF->get_cell_size(iy, comp, 1);

              double weight = (1. + MAC::cos(0.5 * MAC::pi() * dx / deltax))
                            * (1. + MAC::cos(0.5 * MAC::pi() * dy / deltay))
                            / (4. * deltax) / (4. * deltay);

              if ((MAC::abs(dx) > 2. * deltax) || (MAC::abs(dy) > 2. * deltay)) 
                  weight = 0.;
              
              double value = m_all_nodes[i]->force(comp) * weight;
              LF->add_value_to_DOF(ix, iy, 0, comp, 0, value);
            }
          }
        }
      }
    }
  }
}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary::eulerian_velocity_on_lagrange_nodes
                                      ( FV_DiscreteField const *UF
                                      , MAC_Communicator const* macCOMM)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: eulerian_velocity_on_lagrange_nodes");

  geomVector delta(3);
  doubleArray2D data_for_MPI(Ntot,3,0.);

  for (size_t i = 0; i < Ntot; ++i) {

    bool is_on_local_proc = UF->primary_grid()
                          ->is_in_domain_on_current_processor(m_all_nodes[i]->position(0)
                                                            , m_all_nodes[i]->position(1));

    if (is_on_local_proc) {
      for (size_t comp = 0; comp < UF->nb_components(); comp++) {
        m_all_nodes[i]->velocity(comp) = 0.;
        size_t_vector i0(3,0);
        for (size_t dir = 0; dir < UF->primary_grid()->nb_space_dimensions(); dir++) {
            size_t i_temp = 0;
            FV_Mesh::between(UF->get_DOF_coordinates_vector(comp, dir),
                                          m_all_nodes[i]->position(dir), i_temp);
            i0(dir) = i_temp;
        }

        for (size_t ix = i0(0)-2; ix <= i0(0)+2; ix++) {
          for (size_t iy = i0(1)-2; iy <= i0(1)+2; iy++) {
            double dx = m_all_nodes[i]->position(0) 
                      - UF->get_DOF_coordinate(ix, comp, 0);
            double dy = m_all_nodes[i]->position(1) 
                      - UF->get_DOF_coordinate(iy, comp, 1);

            double deltax = UF->get_cell_size(ix, comp, 0);
            double deltay = UF->get_cell_size(iy, comp, 1);

            double weight = (1. + MAC::cos(0.5 * MAC::pi() * dx / deltax))
                          * (1. + MAC::cos(0.5 * MAC::pi() * dy / deltay))
                          / (4. * deltax) / (4. * deltay);

            if ((MAC::abs(dx) > 2. * deltax) || (MAC::abs(dy) > 2. * deltay)) 
                weight = 0.;
            
            double value = UF->DOF_value(ix, iy, 0, comp, 0) * weight * deltax * deltay;
            m_all_nodes[i]->velocity(comp) += value;
          }
        }
      }
      // Store values on an doubleArray2D for MPI communications
      data_for_MPI(i, 0) = m_all_nodes[i]->velocity(0);
      data_for_MPI(i, 1) = m_all_nodes[i]->velocity(1);
      data_for_MPI(i, 2) = m_all_nodes[i]->velocity(2);
    }
  }

  // Transfer and get velocity on all nodes
  macCOMM->sum_array(data_for_MPI);

  // Update the values to local nodal data structure
  for (size_t i = 0; i < m_all_nodes.size(); i++) {
    m_all_nodes[i]->velocity(0) = data_for_MPI(i, 0);
    m_all_nodes[i]->velocity(1) = data_for_MPI(i, 1);
    m_all_nodes[i]->velocity(2) = data_for_MPI(i, 2);
  }

}

//---------------------------------------------------------------------------
void DS_ImmersedBoundary::update_edge_length(FV_DiscreteField const *UF)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: update_edge_length");

  boolVector const *periodic_comp = UF->primary_grid()->get_periodic_directions();
  size_t dim = UF->primary_grid()->nb_space_dimensions();

  for (size_t i = 0; i < m_all_edges.size(); i++) {
     geomVector di = m_all_edges[i]->connecting_node[0]->position
                   - m_all_edges[i]->connecting_node[1]->position;
     for (size_t dir = 0; dir < dim; dir++) {
        bool is_periodic = periodic_comp->operator()(dir);

        if (is_periodic) {
          double isize = UF->primary_grid()->get_main_domain_max_coordinate(dir) 
                      - UF->primary_grid()->get_main_domain_min_coordinate(dir);
          di(dir) = di(dir) - round(di(dir)/isize) * isize;
        }
     }
     m_all_edges[i]->length = MAC::sqrt(di(0)*di(0) + di(1)*di(1) + di(2)*di(2));
  }

}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary::reset_Lagrangian_force_field()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: reset_Lagrangian_force_field");

  for (size_t i = 0; i < m_all_nodes.size(); i++) {
     m_all_nodes[i]->force.set(0.); 
  }

}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary::compute_elastic_force_on_lagrange_nodes
                                      ( double const& Es)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: compute_elastic_force_on_lagrange_nodes");

  for (size_t i = 0; i < m_all_edges.size(); ++i) {
      double lambda = m_all_edges[i]->length/m_all_edges[i]->initial_length;
      double fmag = (lambda > 1.e-14) ?
                  Es * (MAC::pow(lambda,3) - 1.) / MAC::pow(lambda,1.5) : 0.;

      geomVector nhat = m_all_edges[i]->connecting_node[1]->position
                      - m_all_edges[i]->connecting_node[0]->position;

      nhat = nhat * (1./nhat.calcNorm());

      // Add force to connecting nodes
      m_all_edges[i]->connecting_node[0]->force += fmag * nhat;
      m_all_edges[i]->connecting_node[1]->force += - fmag * nhat;
  }

}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary::check_and_update_periodic_clone(
                                                  FV_DiscreteField const *UF)
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: check_and_update_periodic_clone");

  boolVector const *periodic_comp = UF->primary_grid()->get_periodic_directions();
  size_t dim = UF->primary_grid()->nb_space_dimensions();

  for (size_t dir = 0; dir < dim; dir++) {
    bool is_periodic = periodic_comp->operator()(dir);

    if (is_periodic) {
      double isize = UF->primary_grid()->get_main_domain_max_coordinate(dir) 
                   - UF->primary_grid()->get_main_domain_min_coordinate(dir);
      double imin = UF->primary_grid()->get_main_domain_min_coordinate(dir);

      for (size_t i = 0; i < Ntot; i++) {
        m_all_nodes[i]->position(dir) += - MAC::floor((m_all_nodes[i]->position(dir) 
                                                        - imin) / isize) * isize;
      }
    }
  }
}




//---------------------------------------------------------------------------
void DS_ImmersedBoundary::initialize_surface_variables()
//---------------------------------------------------------------------------
{
  MAC_LABEL("DS_ImmersedBoundary:: initialize_surface_variables");

  geomVector vvv(3);

  // Initialize all nodes
  if (m_all_nodes.empty()) {
     m_all_nodes.reserve( Ntot );

     // Initialize structure of node
     for (size_t i = 0; i < Ntot; ++i) {
        struct Node *nnn = NULL;

        nnn = new struct Node();

        nnn->nodeID = 0;

        nnn->position = vvv;
        nnn->velocity = vvv;
        nnn->force = vvv;
        nnn->normal = vvv;

        nnn->neighbor.reserve(3);
        nnn->neighbor.push_back(new struct Node());
        nnn->neighbor.push_back(new struct Node());
        nnn->neighbor.push_back(new struct Node());

        m_all_nodes.push_back(nnn);
     }
  }


  // Initialize all edges
  if (m_all_edges.empty()) {
     m_all_edges.reserve( Ntot );

     // Initialize structure of edge
     for (size_t i = 0; i < Ntot; ++i) {
        struct Edge *eee = NULL;

        eee = new struct Edge();

        eee->edgeID = 0;
        eee->initial_length = 0.;
        eee->length = 0.;

        eee->connecting_node.reserve(2);
        eee->connecting_node.push_back(new struct Node());
        eee->connecting_node.push_back(new struct Node());

        m_all_edges.push_back(eee);
     }
  }

}