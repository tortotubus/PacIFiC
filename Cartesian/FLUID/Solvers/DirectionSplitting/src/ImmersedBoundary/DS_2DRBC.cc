#include <DS_2DRBC.hh>
#include <doubleArray2D.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
#include <cmath>
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
  temp.neighbors.resize(2);
  temp.neighbors.push_back(0);
  temp.neighbors.push_back(0);
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
                            * cos( 2. * M_PI * double(i) / double(num_nodes) ) ;
    m_all_nodes[i].coordinates(1) = shape_param.radius
                            * sin( 2. * M_PI * double(i) / double(num_nodes) ) ;

    // Set the neighbor for each node
    m_all_nodes[i].neighbors[0] = i == 0 ? num_nodes - 1 : i - 1;
    m_all_nodes[i].neighbors[1] = i == num_nodes - 1 ? 0 : i + 1;
  }
}




//---------------------------------------------------------------------------
void DS_2DRBC:: set_all_trielements()
//---------------------------------------------------------------------------
{
  
}




//---------------------------------------------------------------------------
void DS_2DRBC:: set_all_edges()
//---------------------------------------------------------------------------
{
  
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
void DS_2DRBC:: eul_to_lag()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_2DRBC:: eul_to_lag_2D" ) ;

  size_t num_nodes = shape_param.N_nodes;
  // FV_Mesh const* fvm = FF->primary_grid();


  // Generate the node ID and coordinates
  for (size_t i=0; i<num_nodes; ++i)
  {
  }
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





