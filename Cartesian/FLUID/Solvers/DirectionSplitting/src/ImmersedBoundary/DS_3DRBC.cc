#include <DS_3DRBC.hh>
#include <FV_Mesh.hh>
#include <FV_DiscreteField.hh>
#include <MAC_Communicator.hh>
#include <MAC_Exec.hh>
#include <doubleArray2D.hh>
#include <doubleVector.hh>
#include <math.h>
#include <cmath>
#include <typeinfo>
#include <iomanip>      // std::setprecision
#include <set>
#include <fstream>
#include <sstream>
#include <FV_Mesh.hh> // for using nb_space_dimensions() to get dim value & FV_Mesh::between function
using std::endl;
using std::cout;
using std::cin;
using std::string;
using std::max;
using namespace std;

class doubleVector;

#define TOLERANCE 0.05
#define SMALL     0.001
#define SMALLER   0.00001


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
void DS_3DRBC:: initialize_node_properties(string const& mesh_filename,
                                           size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: initialize_node_properties" ) ;
  
  // Read number of nodes for each RBC3D into the variable
  ifstream fileIN( mesh_filename.c_str(), ios::in );
  if(fileIN.fail())
  {
    //File does not exist
    cout << "Mesh file does not exist -- Exiting the program!!" << endl;
    exit(3);
  }
  fileIN >> shape_param.N_nodes;
  fileIN.close();

  m_all_nodes.reserve(shape_param.N_nodes);

  Node temp;
  temp.coordinates(dim);
  temp.coordinates_pbc(dim);
  temp.velocity(dim);
  temp.angular_velocity(dim);
  temp.sumforce(dim);
  temp.sumforce_nm1(dim);
  temp.WLC_force(dim);
  temp.POW_force(dim);
  temp.spring_force(dim);
  temp.bending_force(dim);
  temp.area_force(dim);
  temp.volume_force(dim);
  temp.viscous_force(dim);
  temp.unit_outwards_normal_vector(dim);
  set<Node const*> ww;
  vector< set<Node const*> > neighbors_3D( shape_param.N_nodes, ww );
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
  MAC_LABEL( "DS_3DRBC:: write_one_point_to_VTK" ) ;

}




//---------------------------------------------------------------------------
void DS_3DRBC:: set_all_nodes(istream& fileIN, size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: set_all_nodes" ) ;
  
  fileIN >> shape_param.N_nodes;
  
  size_t num_nodes = shape_param.N_nodes;
  
  // Read node number/ID and their coordinates
  for (size_t i=0;i<num_nodes;++i)
  {
    m_all_nodes[i].number = i;
    fileIN >> m_all_nodes[i].coordinates(0) 
           >> m_all_nodes[i].coordinates(1) 
           >> m_all_nodes[i].coordinates(2);
  }
  
  set_radius(num_nodes, dim);
  
  scale_all_node_coordinates(num_nodes, dim);
}




//---------------------------------------------------------------------------
void DS_3DRBC:: set_radius(size_t const& num_nodes, size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: set_radius" ) ;

  // Compute radius of membrane ASSUMING membrane centroid is (0, 0, 0)
  shape_param.radius = 0.;
  geomVector coords(dim);
  for (size_t i=0;i<num_nodes;++i)
  {
    for (size_t j=0;j<dim;++j) coords(j) = m_all_nodes[i].coordinates(j);

    shape_param.radius = max( shape_param.radius, norm(coords) );
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC:: scale_all_node_coordinates(size_t const& num_nodes, 
                                           size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: scale_all_node_coordinates" ) ;
  
  // Scaling of node coordinates
  for (size_t i=0;i<num_nodes;++i)
    for (size_t j=0;j<dim;++j)
      m_all_nodes[i].coordinates(j) *= shape_param.scaling_factor;
  
  // Scaling the radius based on scaling factor
  shape_param.radius *= shape_param.scaling_factor;
  
  // Order of magnitude of radius of 3D RBC
  shape_param.order_of_magnitude_of_radius 
                                   = pow(10., floor(log10(shape_param.radius)));
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
    double z = ( m_all_nodes[i].coordinates(2) > 0. ? 1. : -1. ) 
               * radius 
               * sqrt( max( 1. - ww, 0. ) ) 
               * ( c0 + c1 * ww + c2 * pow( ww, 2.) ) ;
    m_all_nodes[i].coordinates(2) = z;
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC:: initialize_triangle_properties(size_t const& num_triangles,
                                               size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: initialize_triangle_properties" ) ;
  
  m_all_trielements.reserve(num_triangles);
  
  TriElement temp;
  temp.number = 0;
  temp.twice_area_outwards_normal_vector(dim);
  temp.center_of_mass(dim);
  temp.tri_area = 0.;
  temp.tri_initial_area = 0.;
  temp.tri_volume = 0.;
  temp.tri_initial_volume = 0.;

  for (size_t i = 0; i < num_triangles; i++) {
    m_all_trielements.push_back(temp);
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC:: set_all_trielements(istream& fileIN, size_t const& dim,
                                    bool const& MatlabNumb)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: set_all_trielements" ) ;
  
  size_t n1, n2, n3;
  
  // Reads the number of nodes
  fileIN >> m_nTriangles;
  
  // Allocate memory for members of "TriElement" structure
  initialize_triangle_properties(m_nTriangles, dim);
  
  // Reads triangles
  for (size_t i=0;i<m_nTriangles;++i) 
  {
    fileIN >> n1 >> n2 >> n3;
    if ( MatlabNumb )
    {
        n1--;
        n2--;
        n3--;
    }
    m_all_trielements[i].number = i;
    m_all_trielements[i].pnodes[0] = &(m_all_nodes[n1]);
    m_all_trielements[i].pnodes[1] = &(m_all_nodes[n2]);
    m_all_trielements[i].pnodes[2] = &(m_all_nodes[n3]);
  }
  
  // Define the 2 vectors to compute the outwards oriented twice area normal 
  // vector
  Node const* pnode = NULL;
  pair<Node const*,Node const*> pp(pnode,pnode);
  geomVector twav(dim);
  for (size_t i=0;i<m_nTriangles;++i)
  {
    // Compute center of mass
    for (size_t j=0;j<dim;++j)
    {
        m_all_trielements[i].center_of_mass(j) = ( 
        m_all_trielements[i].pnodes[0]->coordinates(j) +
        m_all_trielements[i].pnodes[1]->coordinates(j) +
        m_all_trielements[i].pnodes[2]->coordinates(j) ) / 3.;
    }

    m_all_trielements[i].varea.reserve(2);
    for (size_t j=0;j<2;++j) m_all_trielements[i].varea.push_back(pp);

    // Set first node in each pair
    m_all_trielements[i].varea[0].first = m_all_trielements[i].pnodes[0];
    m_all_trielements[i].varea[1].first = m_all_trielements[i].pnodes[0];

    // Set arbitrarily second node in each pair
    m_all_trielements[i].varea[0].second = m_all_trielements[i].pnodes[1];
    m_all_trielements[i].varea[1].second = m_all_trielements[i].pnodes[2];

    // Compute twice area vector 
    compute_twice_area_vector( m_all_trielements[i].varea, twav );

    // Compute scalar product with center of the triangle as
    // an approximation of a radial vector at the center of the triangle
    double prod = scalar( m_all_trielements[i].center_of_mass, twav );

    // If hav is negative, swap second node in each pair
    if ( prod < 0. )
    {
        m_all_trielements[i].varea[0].second = m_all_trielements[i].pnodes[2];
        m_all_trielements[i].varea[1].second = m_all_trielements[i].pnodes[1];    
    }

    // Compute twice area normal vector
    compute_twice_area_vector( m_all_trielements[i].varea, 
                       m_all_trielements[i].twice_area_outwards_normal_vector );
    if(isinf(m_all_trielements[i].twice_area_outwards_normal_vector(0)) or isinf(m_all_trielements[i].twice_area_outwards_normal_vector(1)) or isinf(m_all_trielements[i].twice_area_outwards_normal_vector(2)))
    {
      cout << "normals are infinity\n";
      exit(3);
    }
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC::compute_triangle_area_normals_centre_of_mass(bool init,
                                                    size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_triangle_area_normals_centre_of_mass" ) ;

  for (size_t i=0;i<m_nTriangles;++i)
  {
    // Compute normals of triangle
    compute_twice_area_vector( m_all_trielements[i].varea, 
                       m_all_trielements[i].twice_area_outwards_normal_vector );

    m_all_trielements[i].tri_area = 0.5 
                 * norm(m_all_trielements[i].twice_area_outwards_normal_vector);
    
    if(init) 
      m_all_trielements[i].tri_initial_area = m_all_trielements[i].tri_area;

    // Compute center of mass
    for (size_t j=0;j<dim;++j)
    {
      m_all_trielements[i].center_of_mass(j) = 
      ( m_all_trielements[i].pnodes[0]->coordinates(j) +
        m_all_trielements[i].pnodes[1]->coordinates(j) +
        m_all_trielements[i].pnodes[2]->coordinates(j) ) / 3.;
    }
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC::set_all_node_neighbors()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: set_all_node_neighbors" ) ;
  
  size_t num_nodes = shape_param.N_nodes;

  set<Node const*> ww;
  vector< set<Node const*> > neighbors( num_nodes, ww );
  
  // Searches the neighbors
  for (vector<TriElement>::const_iterator il=m_all_trielements.begin();
       il!=m_all_trielements.end();
       il++)
  {
    size_t n0 = ((il->pnodes)[0])->number;
    size_t n1 = ((il->pnodes)[1])->number;
    size_t n2 = ((il->pnodes)[2])->number;

    neighbors[n0].insert( &(m_all_nodes[n1]) );
    neighbors[n0].insert( &(m_all_nodes[n2]) ); 
    neighbors[n1].insert( &(m_all_nodes[n0]) );
    neighbors[n1].insert( &(m_all_nodes[n2]) );     
    neighbors[n2].insert( &(m_all_nodes[n0]) );
    neighbors[n2].insert( &(m_all_nodes[n1]) );
  }                       


  // Transfer sets to vectors and compute initial length
  size_t nn = 0;
  Node const* pp = NULL;
  double zero = 0.;
  for (size_t i=0;i<num_nodes;++i)
  {
    // Allocate
    nn = neighbors[i].size();
    m_all_nodes[i].neighbors_3D.reserve(nn);
    for (size_t j=0;j<nn;++j) m_all_nodes[i].neighbors_3D.push_back(pp);

    m_all_nodes[i].initial_spring_length.reserve(nn); 
    for (size_t j=0;j<nn;++j) 
      m_all_nodes[i].initial_spring_length.push_back(zero);


    // Sets pointers to neighbor and initial length 
    size_t j = 0;
    for (set<Node const*>::iterator is=neighbors[i].begin();
         is!=neighbors[i].end();
         is++,++j)
    {
      m_all_nodes[i].neighbors_3D[j] = *is;
    }
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC:: initialize_edge_properties(size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: initialize_edge_properties" ) ;

  /*
  // Allocate memory for m_all_edges object
  m_all_edges.reserve(m_nEdges);
  
  Edge temp;
  temp.ext_unit_normal(dim);
  temp.initial_length = 0.;
  temp.t1v1.first = &(m_all_trielements[0]);
  temp.t1v1.second = m_all_trielements[0].pnodes[0];
  temp.t2v4.first = &(m_all_trielements[0]);
  temp.t2v4.second = m_all_trielements[0].pnodes[0];
  // // temp.n2 = m_all_trielements[0].pnodes[0];
  // // temp.n3 = m_all_trielements[9].pnodes[1];
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

  for (size_t i=0;i<m_nEdges;++i) 
    m_all_edges.push_back(temp);
  */
}



//---------------------------------------------------------------------------
void DS_3DRBC:: set_all_edges(size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: set_all_edges" ) ;

  Edge ed;
  ed.n2 = NULL;
  ed.n3 = NULL;
  ed.t1v1.first = NULL;
  ed.t1v1.second = NULL;
  ed.t2v4.first = NULL;
  ed.t2v4.second = NULL;
  ed.sintheta0 = ed.costheta0 = 0.;
  Edge* ped = NULL;

  // Create edges
  ed.number = -1;
  for (size_t i=0;i<m_nTriangles;++i)
  {
    for (size_t j=0;j<dim;++j)
    {
      ped = does_edge_exist( m_all_trielements[i].pnodes[j], m_all_trielements[i].pnodes[j == 2 ? 0 : j+1] );

      if ( ped ) // If the edge already exists
      {
        ped->t2v4.first = &(m_all_trielements[i]);
        ped->t2v4.second = m_all_trielements[i].pnodes[ j == 0 ? 2 : j == 1 ? 0 : 1 ];
      }
      else // If the edge does not already exist, add it
      {
        ed.n2 = m_all_trielements[i].pnodes[j];
        ed.n3 = m_all_trielements[i].pnodes[j == 2 ? 0 : j+1];
        ed.t1v1.first = &(m_all_trielements[i]);
        ed.t1v1.second = m_all_trielements[i].pnodes[ j == 0 ? 2 : j == 1 ? 0 : 1 ];
        ed.number++;
        m_all_edges.push_back( ed );
      }
    }
  }
  
  double order_of_magnitude_of_radius = shape_param.order_of_magnitude_of_radius;
  
  // Check if n2-n3 are correctly ordered else swap them
  // i.e., check if normals of triangles on either side 
  // of edge match original normals from set_all_trielements()
  Node* buffer = NULL;
  for (vector<Edge>::iterator il=m_all_edges.begin();il!=m_all_edges.end();il++)
  {
    geomVector normal_t1(dim);
    for (size_t j=0;j<dim;++j) 
      normal_t1(j) = il->t1v1.first->twice_area_outwards_normal_vector(j) / pow(order_of_magnitude_of_radius, 2.);
    
    geomVector normal_t2(dim);
    for (size_t j=0;j<dim;++j) 
      normal_t2(j) = il->t2v4.first->twice_area_outwards_normal_vector(j) / pow(order_of_magnitude_of_radius, 2.);
    
    // a21 = a1 - a2 = vector from node 2 to node 1
    geomVector a21(dim);
    for (size_t j=0;j<dim;++j) 
      a21(j) = (il->t1v1.second->coordinates(j) - il->n2->coordinates(j))/order_of_magnitude_of_radius;
    
    // a31 = a1 - a3 = vector from node 3 to node 1
    geomVector a31(dim);
    for (size_t j=0;j<dim;++j) 
      a31(j) = (il->t1v1.second->coordinates(j) - il->n3->coordinates(j))/order_of_magnitude_of_radius;
    
    // Chi = a21 x a31
    geomVector Chi(dim);
    cross_3D(a21, a31, Chi);
    
    if (!compare( Chi, normal_t1, 1.e-10 ))
    //if(scalar(Chi, normal_t1) < 0.)
    {
      buffer = il->n2;
      il->n2 = il->n3;
      il->n3 = buffer;
      
      // a21 = a1 - a2 = vector from node 2 to node 1
      for (size_t j=0;j<dim;++j) a21(j) = (il->t1v1.second->coordinates(j) - il->n2->coordinates(j))/order_of_magnitude_of_radius;
      
      // a31 = a1 - a3 = vector from node 3 to node 1
      for (size_t j=0;j<dim;++j) a31(j) = (il->t1v1.second->coordinates(j) - il->n3->coordinates(j))/order_of_magnitude_of_radius;
      
      // Chi = a21 x a31
      cross_3D(a21, a31, Chi);
    }
    
    
    
    // a34 = a4 - a3 = vector from node 3 to node 4
    geomVector a34(dim);
    for (size_t j=0;j<dim;++j) a34(j) = (il->t2v4.second->coordinates(j) - il->n3->coordinates(j))/order_of_magnitude_of_radius;
    
    // a24 = a4 - a2 = vector from node 2 to node 4
    geomVector a24(dim);
    for (size_t j=0;j<dim;++j) a24(j) = (il->t2v4.second->coordinates(j) - il->n2->coordinates(j))/order_of_magnitude_of_radius;
    
    // Gamma = a34 x a24
    geomVector Gamma(dim);
    cross_3D(a34, a24, Gamma);
    
    
    if ( (!compare(Chi, normal_t1, 1.e-10)) or (!compare(Gamma, normal_t2, 1.e-10)) )
    //if( (scalar(Chi, normal_t1) < 0.) || (scalar(Gamma, normal_t2) < 0.) )
    {
      cout << "Numbering problem in edge connected by nodes " << il->n2->number << "-" << il->n3->number << endl;
      exit(3);
    }


    double A1 = 0.5 * norm(Chi);
    double A2 = 0.5 * norm(Gamma);
    
    // Compute chi - gamma
    geomVector Chi_minus_Gamma(dim);
    for (size_t j=0;j<dim;++j) Chi_minus_Gamma(j) = Chi(j) - Gamma(j);
    

    // Compute tc1 - tc2
    geomVector tc1(dim), tc2(dim);
    for (size_t j=0;j<dim;++j) tc1(j) = (il->t1v1.second->coordinates(j) + il->n2->coordinates(j) + il->n3->coordinates(j))/3.;
    for (size_t j=0;j<dim;++j) tc2(j) = (il->n2->coordinates(j) + il->n3->coordinates(j) + il->t2v4.second->coordinates(j))/3.;
    geomVector tc1_minus_tc2(dim);
    for (size_t j=0;j<dim;++j) tc1_minus_tc2(j) = (tc1(j) - tc2(j))/order_of_magnitude_of_radius;
    double dot_product = scalar(Chi_minus_Gamma, tc1_minus_tc2);


    // Compute costheta
    il->costheta0 = scalar(Chi, Gamma)/A1/A2/4.0;
    if(il->costheta0 >  1.) il->costheta0 =  1.;
    if(il->costheta0 < -1.) il->costheta0 = -1.;
    il->initial_angle = acos(il->costheta0);
    if(dot_product < 0.) il->initial_angle *= -1.;
    

    // Compute sintheta (with a clipping for values above 0.5 or so)
    il->sintheta0 = sin(il->initial_angle);
    if(fabs(il->sintheta0) < SMALLER) il->sintheta0 = SMALLER;
    
    il->angle_nm1 = il->initial_angle;
    il->dangledt = 0.;
  }
}




//---------------------------------------------------------------------------
Edge* DS_3DRBC::does_edge_exist( Node const* n2_, Node const* n3_ )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: does_edge_exist" ) ;
  
  Edge* ped = NULL;
  bool found = false;
  for (vector<Edge>::iterator il=m_all_edges.begin();
       il!=m_all_edges.end() && !found;
       il++)
  {
    if ( ( il->n2 == n2_ && il->n3 == n3_ ) 
    || ( il->n2 == n3_ && il->n3 == n2_ ) )
    {
        ped = &(*il);
        found = true;
    }
  }

  return ( ped );
} 




//---------------------------------------------------------------------------
bool DS_3DRBC:: compare( geomVector const& v0, 
                         geomVector const& v1, 
                         double const& eps )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compare" ) ;
  
  return ( fabs( v0(0) - v1(0) ) < eps
        && fabs( v0(1) - v1(1) ) < eps
        && fabs( v0(2) - v1(2) ) < eps );
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_spring_lengths(bool init, size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_spring_lengths" );
  
  size_t num_nodes = shape_param.N_nodes;
  
  size_t nn = 0;
  for (size_t i=0;i<num_nodes;++i)
  {
    nn = m_all_nodes[i].initial_spring_length.size();
    for (size_t j=0;j<nn;++j)
    {
      double length = 0.;
      for (size_t k=0;k<dim;++k)
      {
        length += pow(m_all_nodes[i].coordinates(k)
                    - m_all_nodes[i].neighbors_3D[j]->coordinates(k), 2.);
      }
      m_all_nodes[i].initial_spring_length[j] = pow( length, 0.5 );
    }
  }
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
void DS_3DRBC:: init_membrane_parameters_in_physical_units()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: init_membrane_parameters_in_physical_units" );
  
  // membrane_param.mu0_P = 6.3e-6; // 29.e-6; // 6.3e-6; // Shear modulus in N/m
  membrane_param.Y_P = 18.9e-6; // 13.3437e-6; // 18.9e-6; // Young's modulus in N/m
  // membrane_param.x0 = 1./2.2; // 1./1.8; // malaria cell = 1./1.8 // healthy cell = 1./2.2; // Maximum allowable spring extension --> x0=l/lmax
  
  cout << membrane_param.x0 << "\t" << membrane_param.mu0_P << endl; exit(3);
  membrane_param.D0_P = 2. * shape_param.radius; // 7.82e-6; // Diameter of RBC micro-metre
  membrane_param.kc_P = 2.4e-19; // 4.8e-19; // 2.4e-19; // bending rigidity - Joules
  membrane_param.kbending_P = (2./sqrt(3.)) * membrane_param.kc_P; // bending spring constant
  membrane_param.eta_P = 0.022; // 0.022 * 10.; // 6.e-3; // Viscosity of membrane Newton-second/metre^2 = Pa.s = 6 centipoise as given in Pozrikidis book Pg. 124
  membrane_param.eta_plasma_P = 1.2e-3; // Viscosity of plasma fluid surrounding RBC Newton-second/metre^2 = Pa.s = 1.2 centipoise as given in Pozrikidis book Pg. 124
  membrane_param.eta_cytoplasm_P = 5.e-3; // Viscosity of fluid inside RBC Pa.s // internal fluid inside the RBC
  membrane_param.rho_plasma_P = 1023.9; // Density of fluid surrounding RBC in kg/m^3 as given in Pozrikidis book Pg. 253 // fluid viscosity
  membrane_param.rho_P = membrane_param.rho_plasma_P; // Density of RBC
  membrane_param.m = 2.;
  membrane_param.alpha = 1.;
}




//---------------------------------------------------------------------------
void DS_3DRBC:: init_membrane_parameters_in_model_units()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: init_membrane_parameters_model_units" );
  
  membrane_param.node_mass_M = 1.; // TO BE VERIFIED
  membrane_param.mu0_M = 100.; // model unit's shear modulus
  membrane_param.k_area = 4900.; // ka
  membrane_param.k_area_local = 100.; // kd
  membrane_param.k_volume = 5000.; // kv
  membrane_param.D0_M = 7.82; // 8.25; // given by "user"
  membrane_param.eta_M = 6.6; // 1.8; // 6.6 * 10.; // in model units // I think this is the value we choose just like choosing membrane_param.D0_M
  membrane_param.eta_plasma_M = 10.; // in model units -- need better value prediction mathematically or analytically
}




//---------------------------------------------------------------------------
void DS_3DRBC:: scaling_membrane_params_from_physical_to_model_units()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: scaling_membrane_params_from_physical_to_model_units" ) ;
  
  double mu0_M = membrane_param.mu0_M;
  double ka = membrane_param.k_area;
  double kd = membrane_param.k_area_local;
  double Y_P = membrane_param.Y_P;
  double D0_P = membrane_param.D0_P;
  double D0_M = membrane_param.D0_M;
  double eta_P = membrane_param.eta_P;
  double eta_M = membrane_param.eta_M;
  double alpha = membrane_param.alpha;
  
  // Compression/Volumetric/Bulk modulus of membrane in model units
  double K_M = 2. * mu0_M + ka + kd;
                               
  // Young's modulus of membrane in model units
  double Y_M = 4. * K_M * mu0_M / (K_M + mu0_M);
  
  
  // // time scale
  membrane_param.t = pow( ((D0_P * eta_P/Y_P) / (D0_M * eta_M/Y_M)), alpha );
  
  // // membrane_param.t = pow( ((D0_P * eta_P/mu0_P) / (D0_M * eta_M/mu0_M)), alpha );
  // // membrane_param.t = pow( ((D0_P * eta_plasma_P/mu0_P) / (D0_M * eta_plasma_M/mu0_M)), alpha );
  
  
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_spring_constant_values_in_model_units(size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_spring_constant_values" ) ;
  
  double x0 = membrane_param.x0;
  double alpha_1 = 1./(4.*pow(1-x0, 2.));
  double alpha_2 = 1./4.;
  double alpha_3 = x0;
  double alpha = alpha_1 - alpha_2 + alpha_3;
  
  double beta_1 = x0/(2.*pow(1-x0, 3.));
  double beta_2 = 1./(4.*pow(1-x0, 2.));
  double beta_3 = 1./4.;
  double beta = beta_1 - beta_2 + beta_3;
  
  double denominator = sqrt(3.) * (beta + 3. * alpha);
  
  double O_R = shape_param.order_of_magnitude_of_radius;
  double mu0_M = membrane_param.mu0_M;
  
  geomVector length_component(dim);
  double edge_length;
  
  for (vector<Edge>::iterator il=m_all_edges.begin();il!=m_all_edges.end();il++)
  {
    for (size_t j=0;j<dim;++j) 
      length_component(j) = il->n2->coordinates(j) - il->n3->coordinates(j);

    edge_length = norm(length_component);
    
    // Initial length of each spring
    il->l0 = edge_length;
    
    // Max allowed length of each spring
    il->lmax = il->l0 / x0;
    
    // Spring length in model units with O_R being scale of reference
    double l0_M = il->l0 / O_R;
    
    // Spring constant values
    il->k =  ( 4. * mu0_M * l0_M ) / denominator; // kBT/p in WLC spring model
    il->kp = ( 4. * mu0_M * alpha * pow(l0_M, 3.) ) / denominator; // = kp in POW spring model
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_bending_constant_values_in_model_units()
//---------------------------------------------------------------------------
{
  double kbending_P = membrane_param.kbending_P;
  double D0_P = membrane_param.D0_P;
  double mu0_P = membrane_param.mu0_P;
  double D0_M = membrane_param.D0_M;
  double mu0_M = membrane_param.mu0_M;
  
  // Computing non-dimensional quantity
  double gamma = kbending_P / (pow(D0_P, 2.) * mu0_P);
  
  // Computing bending constant in model units
  membrane_param.kbending_M = gamma * pow(D0_M, 2.) * mu0_M;
}




//---------------------------------------------------------------------------
void DS_3DRBC:: preprocess_membrane_parameters(string const& model_type,
                                             string const& case_type,
                                             double const& mu,
                                             size_t const& num_subtimesteps_RBC,
                                             size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: preprocess_membrane_parameters" ) ;
  
  if(model_type.compare("NumericalMembraneModel") == 0) // if it is detailed numerical membrane model (NMM)
  {
    // Set the spring constant values
    // Initialize membrane material properties in physical units
    init_membrane_parameters_in_physical_units();
    
    // Initialize membrane material properties in model units
    init_membrane_parameters_in_model_units();
    
    // If it's shear flow tank treading case - recompute the mu0_P
    if(case_type.compare("Breyannis2000case") == 0) // if it's not Breyainnis case
    {
      membrane_param.mu0_P = mu * membrane_param.ShearRate * shape_param.radius 
                                / membrane_param.CapillaryNumber;
      membrane_param.D0_P = 2. * shape_param.radius;
      membrane_param.Y_P = membrane_param.FopplVonKarmanNumber 
                            * membrane_param.kc_P/pow(shape_param.radius, 2.);
      // for incompressible sheet, K >> mu0 ==> in Y = 4Kmu0/(K+mu0), we've Y = 4mu0
      membrane_param.mu0_P = membrane_param.Y_P / 4.;
      // // membrane_param.eta_P = 0.;
      // // membrane_param.eta_M = 0.;
    }

    // Scaling membrane material properties from physical units to model units
    scaling_membrane_params_from_physical_to_model_units();
    
    // Compute spring constant values in model units
    compute_spring_constant_values_in_model_units(dim);
    
    // Compute bending constant value in model units
    compute_bending_constant_values_in_model_units();
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
    if(case_type.compare("Breyannis2000case") != 0) // if it's not Breyainnis case
    {
      /*
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
      */
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
void DS_3DRBC::compute_twice_area_vector
      (vector< pair<Node const*,Node const*> > const& ppnodes, geomVector& res )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_twice_area_vector" ) ;
  
  geomVector v0(3);
  geomVector v1(3);

  for (size_t j=0;j<3;++j)
  {
    v0(j) = ppnodes[0].second->coordinates(j) 
            - ppnodes[0].first->coordinates(j);
    v1(j) = ppnodes[1].second->coordinates(j) 
            - ppnodes[1].first->coordinates(j);
  }

  cross_3D( v0, v1, res );
}




//---------------------------------------------------------------------------
void DS_3DRBC:: write_mesh_to_vtk_file( size_t IB_number, double const& time,
                                        size_t const& cyclenum )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: write_mesh_to_vtk_file()" ) ;

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
      fileOUT << m_all_trielements[i].pnodes[0]->number 
              << " " 
              << m_all_trielements[i].pnodes[1]->number 
              << " " 
              << m_all_trielements[i].pnodes[2]->number 
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
  double hxC, hyC, hzC; // Reciprocal of dxC, dyC, dzC
  int Nx, Ny, Nz;
  double dirac_x, dirac_y, dirac_z, dirac_delta; // Dirac delta variables
  size_t istart, iend, jstart, jend, kstart, kend;
  size_t iLag, jLag, kLag;
  size_t ipi, ipf, jpi, jpf, kpi, kpf;
  double sum_dirac_delta;
  double dist_x, dist_y, dist_z;
  bool eul_cell_within_Dirac_delta_stencil;

  // Mesh pertaining to the Eulerian velocity component
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
  
  // z-direction start & end proc indices
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
    
    // if a node belongs to a processor, then compute it's
    // Lagrangian velocity from neighbouring Eulerian nodes
    if(lag_cell_within_processor)
    {
      sum_dirac_delta = 0.0;
      
      // Find the Eulerian cell to which Lagrangian node belong to
      bool found_i = FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,0), xp, iLag);
      bool found_j = FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,1), yp, jLag);
      bool found_k = FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,2), zp, kLag);
      
      // Compute the stencil support bounds
      size_t ipi = iLag - 3;
      size_t ipf = iLag + 3;
      size_t jpi = jLag - 3;
      size_t jpf = jLag + 3;
      size_t kpi = kLag - 3;
      size_t kpf = kLag + 3;

      // Loop over Dirac stencil support
      for (size_t ii=ipi;ii<=ipf;++ii)
      {
        for (size_t jj=jpi;jj<=jpf;++jj)
        {
          for (size_t kk=kpi;kk<=kpf;++kk)
          {
            bool euler_cell_in_stencil_inside_domain = FF->DOF_in_domain(ii, jj, kk, comp);
            
            if(euler_cell_in_stencil_inside_domain)
            {
              xC = FF->get_DOF_coordinate( ii, comp, 0 ) ;
              yC = FF->get_DOF_coordinate( jj, comp, 1 ) ;
              zC = FF->get_DOF_coordinate( kk, comp, 2 ) ;
              
              // Check if Eulerian cell is within Dirac delta 2x2 stencil
              dist_x = compute_dist_incl_pbc(xC, xp, domain_length(0)) * hxC;
              // // dist_x = (xC - xp) * hxC;
              dist_y = compute_dist_incl_pbc(yC, yp, domain_length(1)) * hyC;
              // // dist_y = (yC - yp) * hyC;
              dist_z = compute_dist_incl_pbc(zC, zp, domain_length(2)) * hzC;
              // // dist_z = (zC - zp) * hzC;
              eul_cell_within_Dirac_delta_stencil = (fabs(dist_x) <= 2.) 
                                                    and 
                                                    (fabs(dist_y) <= 2.) 
                                                    and 
                                                    (fabs(dist_z) <= 2.);
              
              // Compute the Dirac delta value for each direction
              dirac_x = discrete_Dirac_delta(dist_x, ibm_param.dirac_type, dxC, Nx);
              dirac_y = discrete_Dirac_delta(dist_y, ibm_param.dirac_type, dyC, Ny);
              dirac_z = discrete_Dirac_delta(dist_z, ibm_param.dirac_type, dzC, Nz);

              // Dirac delta function value for cell (ii, jj, kk)
              dirac_delta = dirac_x * dirac_y * dirac_z;

              // Numerical integration of Dirac delta function value
              sum_dirac_delta += dirac_delta * dxC * dyC * dzC;

              // Computing Lagrangian velocity
              m_all_nodes[inode].velocity(comp) += 
                                FF->DOF_value( ii, jj, kk, comp, 0 ) 
                                * dirac_delta * dxC * dyC * dzC;
            }
          }
        }
      }
      
      // Normalizing the Lagrangian velocity to the support size (specifically for boundary nodes)
      // total_lag_velocity = \sum\limits_{ii, jj, kk} dirac_delta(ii, jj, kk) * eulerian_velocity(ii, jj, kk)
      // total_dirac_delta = \sum\limits_{ii, jj, kk} dirac_delta(ii, jj, kk)
      // lagrangian_velocity = total_lag_velocity / total_dirac_delta
      // // // m_all_nodes[inode].velocity(comp) /= sum_dirac_delta;
      
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
  double dirac_x, dirac_y, dirac_z, dirac_delta; // Dirac delta variables
  size_t istart, iend, jstart, jend, kstart, kend;
  size_t iLag, jLag, kLag;
  double dist_x, dist_y, dist_z;
  double sum_dirac_delta = 0.;
  double sum_euler_force = 0., euler_force_tag = 0.;
  bool eul_cell_within_Dirac_delta_stencil;

  // Mesh pertaining to the Eulerian force component
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
  
  // z-direction start & end proc indices
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
    sum_dirac_delta = 0.;
    sum_euler_force = 0.;
    euler_force_tag = 0.;
    
    // Get coordinates of inode's Lagrangian marker
    double xp = m_all_nodes[inode].coordinates_pbc(0);
    double yp = m_all_nodes[inode].coordinates_pbc(1);
    double zp = m_all_nodes[inode].coordinates_pbc(2);
    
    // Find the Eulerian cell to which Lagrangian node belong to
    bool found_i = FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,0), xp, iLag);
    bool found_j = FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,1), yp, jLag);
    bool found_k = FV_Mesh::between(FF->get_DOF_coordinates_vector(comp,2), zp, kLag);
    
    // Compute the stencil support bounds
    size_t ipi = iLag - 3;
    size_t ipf = iLag + 3;
    size_t jpi = jLag - 3;
    size_t jpf = jLag + 3;
    size_t kpi = kLag - 3;
    size_t kpf = kLag + 3;
    
    // Loop over Dirac stencil support
    for (size_t ii=ipi;ii<=ipf;++ii)
    {
      for (size_t jj=jpi;jj<=jpf;++jj)
      {
        for (size_t kk=kpi;kk<=kpf;++kk)
        {
          // Is the Eulerian cell (ii, jj, kk) within the simulation domain 
          // AND 
          // within the current processor's domain?
          bool euler_cell_within_x_proc_and_domain_bounds 
                                 = (ii >= min_unknown_index_without_halozone(0) 
                                   and 
                                   ii <= max_unknown_index_without_halozone(0));

          bool euler_cell_within_y_proc_and_domain_bounds 
                                 = (jj >= min_unknown_index_without_halozone(1) 
                                   and 
                                   jj <= max_unknown_index_without_halozone(1));

          bool euler_cell_within_z_proc_and_domain_bounds 
                                 = (kk >= min_unknown_index_without_halozone(2) 
                                   and 
                                   kk <= max_unknown_index_without_halozone(2));

          if( euler_cell_within_x_proc_and_domain_bounds
              and
              euler_cell_within_y_proc_and_domain_bounds
              and
              euler_cell_within_z_proc_and_domain_bounds )
          {
            bool euler_cell_is_unknown = FF->DOF_is_unknown(ii, jj, kk, comp);
            
            if(euler_cell_is_unknown)
            {
              xC = FF->get_DOF_coordinate( ii, comp, 0 ) ;
              yC = FF->get_DOF_coordinate( jj, comp, 1 ) ;
              zC = FF->get_DOF_coordinate( kk, comp, 2 ) ;
            
              // Check if Eulerian cell is within Dirac delta 2x2 stencil
              dist_x = compute_dist_incl_pbc(xC, xp, domain_length(0)) * hxC;
              // // dist_x = (xC - xp) * hxC;
              dist_y = compute_dist_incl_pbc(yC, yp, domain_length(1)) * hyC;
              // // dist_y = (yC - yp) * hyC;
              dist_z = compute_dist_incl_pbc(zC, zp, domain_length(2)) * hzC;
              // // dist_z = (zC - zp) * hzC;
              eul_cell_within_Dirac_delta_stencil = (fabs(dist_x) <= 2.) 
                                                    and 
                                                    (fabs(dist_y) <= 2.) 
                                                    and 
                                                    (fabs(dist_z) <= 2.);
              
              // Compute the Dirac delta value for each direction
              dirac_x = discrete_Dirac_delta(dist_x, ibm_param.dirac_type, dxC, Nx);
              dirac_y = discrete_Dirac_delta(dist_y, ibm_param.dirac_type, dyC, Ny);
              dirac_z = discrete_Dirac_delta(dist_z, ibm_param.dirac_type, dzC, Nz);

              // Dirac delta function value
              dirac_delta = dirac_x * dirac_y * dirac_z;

              // Numerical integration of Dirac delta function value
              sum_dirac_delta += dirac_delta * dxC * dyC * dzC;
              
              // Computing Eulerian force
              double euler_force = FF->DOF_value(ii, jj, kk, comp, 0) 
                                   + 
                                   m_all_nodes[inode].sumforce(comp) 
                                   * dirac_delta * dxC * dyC * dzC;
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
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC::compute_total_surface_area_total_volume( bool init,
                                                        size_t const& dim )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_total_surface_area_total_volume" ) ;
  
  // Total surface area
  membrane_param.total_area = 0.;
  for (size_t i=0;i<m_nTriangles;++i)
    membrane_param.total_area += m_all_trielements[i].tri_area;
  if ( init ) membrane_param.initial_area = membrane_param.total_area;
    
  /*
  // Compute center of mass of each triangle
  for (size_t i=0;i<m_nTriangles;++i)
  {
    for (size_t j=0;j<dim;++j)
    {
      m_all_trielements[i].center_of_mass(j) = 
      ( m_all_trielements[i].pnodes[0]->coordinates(j) +
        m_all_trielements[i].pnodes[1]->coordinates(j) +
        m_all_trielements[i].pnodes[2]->coordinates(j) ) / 3.;
    }
  }
  */


  // Center of mass of the membrane
  for (size_t j=0;j<dim;++j) 
    membrane_param.centroid_coordinates(j) = 0.;
  for (size_t i=0;i<m_nTriangles;++i)
  {
    for (size_t j=0;j<dim;++j)
    {
      membrane_param.centroid_coordinates(j) += m_all_trielements[i].tri_area 
                                       * m_all_trielements[i].center_of_mass(j);
    }
  }
  for (size_t j=0;j<dim;++j) 
    membrane_param.centroid_coordinates(j) /= membrane_param.total_area;


  // Total volume
  membrane_param.total_volume = 0.;
  geomVector relpos(dim);
  double oneoversix = 1./6.;
  for (size_t i=0;i<m_nTriangles;++i)
  {
    for (size_t j=0;j<dim;++j) relpos(j) = m_all_trielements[i].center_of_mass(j) - membrane_param.centroid_coordinates(j);
    membrane_param.total_volume += oneoversix * scalar( m_all_trielements[i].twice_area_outwards_normal_vector, relpos );
  }
  if ( init ) membrane_param.initial_volume = membrane_param.total_volume; 
}




//---------------------------------------------------------------------------
void DS_3DRBC::compute_membrane_area_centroid_volume( size_t const& dim )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_membrane_area_centroid_volume" ) ;
  
  // Membrane surface area
  membrane_param.total_area = 0.;
  for (size_t i=0;i<m_nTriangles;++i)
    membrane_param.total_area += m_all_trielements[i].tri_area;
    
  // Membrane centroid
  for (size_t j=0;j<dim;++j) 
    membrane_param.centroid_coordinates(j) = 0.;
  for (size_t i=0;i<m_nTriangles;++i)
  {
    for (size_t j=0;j<dim;++j)
    {
      membrane_param.centroid_coordinates(j) += m_all_trielements[i].tri_area 
                                       * m_all_trielements[i].center_of_mass(j);
    }
  }
  for (size_t j=0;j<dim;++j) 
    membrane_param.centroid_coordinates(j) /= membrane_param.total_area;


  // Membrane volume
  membrane_param.total_volume = 0.;
  geomVector relpos(dim);
  double oneoversix = 1./6.;
  for (size_t i=0;i<m_nTriangles;++i)
  {
    for (size_t j=0;j<dim;++j) relpos(j) = m_all_trielements[i].center_of_mass(j) - membrane_param.centroid_coordinates(j);
    membrane_param.total_volume += oneoversix * scalar( m_all_trielements[i].twice_area_outwards_normal_vector, relpos );
  }
}



    
//---------------------------------------------------------------------------
void DS_3DRBC::compute_spring_bending_viscous_forces( size_t const& dim, 
                                            double const& spring_constant, 
                                            double const& bending_spring_constant,
                                            double const& bending_viscous_constant,
                                            double const& viscous_drag_constant, 
                                            double const& dt, 
                                            string const& model_type )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_spring_bending_viscous_forces" ) ;
  
  size_t num_nodes = shape_param.N_nodes;

  if(model_type.compare("NumericalMembraneModel") == 0)
  {
    double order_of_magnitude_of_radius = shape_param.order_of_magnitude_of_radius;
    
    geomVector Iij(dim);
    double length = 0.;
    geomVector f_WLC(dim);
    geomVector f_POW(dim);

    geomVector a21(dim), a31(dim), a24(dim), a34(dim);
    geomVector Chi(dim), Gamma(dim), crossp(dim), n2n3(dim);
    geomVector Chi_minus_Gamma(dim), tc1_minus_tc2(dim), tc1(dim), tc2(dim);
    
    membrane_param.mean_WLC_spring_force_magnitude = 0.;
    membrane_param.mean_POW_spring_force_magnitude = 0.;
    membrane_param.mean_bending_force_magnitude = 0.;
    membrane_param.mean_viscous_force_magnitude = 0.;
    
    for (size_t i=0;i<num_nodes;++i)
    {
        for (size_t j=0;j<dim;++j)
        {
            m_all_nodes[i].WLC_force(j) = 0.;
            m_all_nodes[i].POW_force(j) = 0.;
            m_all_nodes[i].bending_force(j) = 0.;
            m_all_nodes[i].viscous_force(j) = 0.;
        }
    }
        
    for (vector<Edge>::iterator il=m_all_edges.begin();il!=m_all_edges.end();il++)
    {
        //----------------------//
        /**** Spring force   ****/
        //----------------------//
        // Computing the spring vector from node 3 to node 2, i.e., node j to i
        for (size_t j=0;j<dim;++j)
          Iij(j) = il->n2->coordinates(j) - il->n3->coordinates(j); 

        // Spring length
        length = norm ( Iij );

        // Normalizing the spring vector to unit vector
        for (size_t j=0;j<dim;++j)
          Iij(j) /= length;
        
        // WLC spring prefactor
        double x = length / il->lmax;
        double alpha_1 = 1. / (4. * pow(1.-x, 2.));
        double alpha_2 = 1./4.;
        double alpha_3 = x;
        double alpha = alpha_1 - alpha_2 + alpha_3;

        // POW spring prefactor
        double l_M = length / order_of_magnitude_of_radius;
        double beta = 1. / pow(l_M, 2.);

        for (size_t j=0;j<dim;++j)
        {
            // WLC force
            f_WLC(j) = - il->k * alpha * Iij(j);
            il->n2->sumforce(j)  +=  f_WLC(j);
            il->n3->sumforce(j)  += -f_WLC(j);
            il->n2->WLC_force(j) +=  f_WLC(j);
            il->n3->WLC_force(j) += -f_WLC(j);

            // POW force
            f_POW(j) = il->kp * beta * Iij(j);
            il->n2->sumforce(j)  +=  f_POW(j);
            il->n3->sumforce(j)  += -f_POW(j);
            il->n2->WLC_force(j) +=  f_POW(j);
            il->n3->WLC_force(j) += -f_POW(j);
        }
        
        
        
        //----------------------//
        /**** Bending force  ****/
        //----------------------//
        // Steps for computing bending force
        
        //  1. Compute a21, a31
        // a21 = a1 - a2 = vector from node 2 to node 1
        for (size_t j=0;j<dim;++j) 
          a21(j) = (-il->t1v1.second->coordinates(j) + il->n2->coordinates(j))/order_of_magnitude_of_radius;
        
        // a31 = a1 - a3 = vector from node 3 to node 1
        for (size_t j=0;j<dim;++j) 
          a31(j) = (-il->t1v1.second->coordinates(j) + il->n3->coordinates(j))/order_of_magnitude_of_radius;
        
        
        //  2. Compute chi
        cross_3D(a21, a31, Chi);
        
        
        //  3. Compute a34, a24
        // a34 = a4 - a3 = vector from node 3 to node 4
        for (size_t j=0;j<dim;++j) 
          a34(j) = (-il->t2v4.second->coordinates(j) + il->n3->coordinates(j))/order_of_magnitude_of_radius;
        
        // a24 = a4 - a2 = vector from node 2 to node 4
        for (size_t j=0;j<dim;++j) 
          a24(j) = (-il->t2v4.second->coordinates(j) + il->n2->coordinates(j))/order_of_magnitude_of_radius;
        
        
        //  4. Compute gamma
        cross_3D(a34, a24, Gamma);
            
            
        double A1 = 0.5 * norm(Chi);
        double A2 = 0.5 * norm(Gamma);
            
        //  6. Compute chi - gamma
        for (size_t j=0;j<dim;++j) Chi_minus_Gamma(j) = Chi(j) - Gamma(j);
            
        
        //  7. Compute tc1 - tc2
        for (size_t j=0;j<dim;++j) tc1(j) = (il->t1v1.second->coordinates(j) + il->n2->coordinates(j) + il->n3->coordinates(j))/3.;
        for (size_t j=0;j<dim;++j) tc2(j) = (il->n2->coordinates(j) + il->n3->coordinates(j) + il->t2v4.second->coordinates(j))/3.;
        for (size_t j=0;j<dim;++j) tc1_minus_tc2(j) = (tc1(j) - tc2(j))/order_of_magnitude_of_radius;
        double dot_product = scalar(Chi_minus_Gamma, tc1_minus_tc2);


        //  5. Compute costheta
        double costheta = scalar(Chi, Gamma)/A1/A2/4.0;
            
        
        //  6. Compute sintheta (with a clipping for values above 0.5 or so)
        if(costheta >  1.) costheta =  1.;
        if(costheta < -1.) costheta = -1.;
        double theta = acos(costheta);
        if(dot_product < 0.) theta *= -1.;
        double sintheta = sin(theta);
        if(fabs(sintheta) < SMALLER) sintheta = SMALLER;
        double siinv = 1. / sintheta;
            



            
        //  7. Compute a12, a13, a43, a42, a32, a23
        // a12 = a2 - a1 = vector from node 1 to node 2
        geomVector a12(dim);
        for (size_t j=0;j<dim;++j) a12(j) = (-il->n2->coordinates(j) + il->t1v1.second->coordinates(j))/order_of_magnitude_of_radius;
        
        // a13 = a3 - a1 = vector from node node 1 to node 3
        geomVector a13(dim);
        for (size_t j=0;j<dim;++j) a13(j) = (-il->n3->coordinates(j) + il->t1v1.second->coordinates(j))/order_of_magnitude_of_radius;
        
        // a43 = a3 - a4 = vector from node 4 to node 3
        geomVector a43(dim);
        for (size_t j=0;j<dim;++j) a43(j) = (-il->n3->coordinates(j) + il->t2v4.second->coordinates(j))/order_of_magnitude_of_radius;
        
        // a42 = a2 - a4 = vector from node 4 to node 2
        geomVector a42(dim);
        for (size_t j=0;j<dim;++j) a42(j) = (-il->n2->coordinates(j) + il->t2v4.second->coordinates(j))/order_of_magnitude_of_radius;
        
        // a32 = a2 - a3 = vector from node 3 to node 2
        geomVector a32(dim);
        for (size_t j=0;j<dim;++j) a32(j) = (-il->n2->coordinates(j) + il->n3->coordinates(j))/order_of_magnitude_of_radius;
        
        // a23 = a3 - a2 = vector from node to node 3
        geomVector a23(dim);
        for (size_t j=0;j<dim;++j) a23(j) = (-il->n3->coordinates(j) + il->n2->coordinates(j))/order_of_magnitude_of_radius;
            
            
        //  8. Compute beta_b
        double beta_b = membrane_param.kbending_M * (sintheta * il->costheta0 - costheta * il->sintheta0) * siinv;
        
        //  9. Compute b11, b12, b22
        double b11 = - beta_b*costheta/A1/A1/4.0;
        double b12 =   beta_b/A1/A2/4.0;
        double b22 = - beta_b*costheta/A2/A2/4.0;
            

        //  10. Compute nodal forces
        // Node 1
        geomVector Chi_cross_a32(dim);
        cross_3D(Chi, a32, Chi_cross_a32);
        geomVector Gamma_cross_a32(dim);
        cross_3D(Gamma, a32, Gamma_cross_a32);
        for (size_t j=0;j<dim;++j) il->t1v1.second->sumforce(j) += b11 * Chi_cross_a32(j) + b12 * Gamma_cross_a32(j);
        for (size_t j=0;j<dim;++j) il->t1v1.second->bending_force(j) += b11 * Chi_cross_a32(j) + b12 * Gamma_cross_a32(j);
        
        // Node 2
        geomVector Chi_cross_a13(dim);
        cross_3D(Chi, a13, Chi_cross_a13);
        geomVector Chi_cross_a34(dim);
        cross_3D(Chi, a34, Chi_cross_a34);
        geomVector Gamma_cross_a13(dim);
        cross_3D(Gamma, a13, Gamma_cross_a13);
        geomVector Gamma_cross_a34(dim);
        cross_3D(Gamma, a34, Gamma_cross_a34);
        for (size_t j=0;j<dim;++j) il->n2->sumforce(j) += b11 * Chi_cross_a13(j) + b12 * (Chi_cross_a34(j) + Gamma_cross_a13(j)) + b22 * Gamma_cross_a34(j);
        for (size_t j=0;j<dim;++j) il->n2->bending_force(j) += b11 * Chi_cross_a13(j) + b12 * (Chi_cross_a34(j) + Gamma_cross_a13(j)) + b22 * Gamma_cross_a34(j);
        
        // Node 3
        geomVector Chi_cross_a21(dim);
        cross_3D(Chi, a21, Chi_cross_a21);
        geomVector Chi_cross_a42(dim);
        cross_3D(Chi, a42, Chi_cross_a42);
        geomVector Gamma_cross_a21(dim);
        cross_3D(Gamma, a21, Gamma_cross_a21);
        geomVector Gamma_cross_a42(dim);
        cross_3D(Gamma, a42, Gamma_cross_a42);
        for (size_t j=0;j<dim;++j) il->n3->sumforce(j) += b11 * Chi_cross_a21(j) + b12 * (Chi_cross_a42(j) + Gamma_cross_a21(j)) + b22 * Gamma_cross_a42(j);
        for (size_t j=0;j<dim;++j) il->n3->bending_force(j) += b11 * Chi_cross_a21(j) + b12 * (Chi_cross_a42(j) + Gamma_cross_a21(j)) + b22 * Gamma_cross_a42(j);
        
        // Node 4
        geomVector Chi_cross_a23(dim);
        cross_3D(Chi, a23, Chi_cross_a23);
        geomVector Gamma_cross_a23(dim);
        cross_3D(Gamma, a23, Gamma_cross_a23);
        for (size_t j=0;j<dim;++j) il->t2v4.second->sumforce(j) += b12 * Chi_cross_a23(j) + b22 * Gamma_cross_a23(j);
        for (size_t j=0;j<dim;++j) il->t2v4.second->bending_force(j) += b12 * Chi_cross_a23(j) + b22 * Gamma_cross_a23(j);
        
        
        
        //----------------------//
        /**** Viscous force ****/
        //----------------------//
        double eta_M = membrane_param.eta_M;

        geomVector eij(dim);
        geomVector vij(dim);
        geomVector vij_M(dim);

        // (x) ------------------- (x)
        //  i                       j
        // Nodes i and j of a spring
        // Computing viscous force as two lines of code for nodes i and j
        geomVector vjmvi(dim);
        geomVector vimvj(dim);

        // Vector of each spring
        for (size_t j=0;j<dim;++j)
            eij(j) = il->n3->coordinates(j) - il->n2->coordinates(j);
        
        // Length of each spring
        double length = norm ( eij );
        
        // Unit vector of each spring
        for (size_t j=0;j<dim;++j) 
            eij(j) /= length;

        // Relative velocity between nodes i and j
        for (size_t j=0;j<dim;++j)
        {
            vimvj(j) = (il->n2->velocity(j) - il->n3->velocity(j)) * (dt/order_of_magnitude_of_radius);
            vjmvi(j) = (il->n3->velocity(j) - il->n2->velocity(j)) * (dt/order_of_magnitude_of_radius);
        }
    
        // Projection of relative velocity on the spring unit vector
        double vimvj_dot_eij = scalar(vimvj, eij);
        double vjmvi_dot_eij = scalar(vjmvi, eij);
        
        // Viscous force computation
        double constant_1 = 12.0 / (13.0 * sqrt(3.0));
        double constant_2 =  4.0 / (13.0 * sqrt(3.0)); // gamma_C = gamma_T/3
        double gamma_T = constant_1 * eta_M; // T = translational component
        double gamma_C = constant_2 * eta_M; // C = central component
        for (size_t j=0;j<dim;++j)
        {
            // Viscous force on node i
            double f_i = - (gamma_T * vimvj(j) + gamma_C * vimvj_dot_eij * eij(j));
            il->n2->sumforce(j) += f_i;
            il->n2->viscous_force(j) += f_i;

            // Viscous force on node j
            double f_j = - (gamma_T * vjmvi(j) + gamma_C * vjmvi_dot_eij * eij(j));
            il->n3->sumforce(j) += f_j;
            il->n3->viscous_force(j) += f_j;
        }
    }
    
    // Mean spring, bending & viscous forces
    for (size_t i=0;i<num_nodes;++i)
    {
        membrane_param.mean_WLC_spring_force_magnitude += norm(m_all_nodes[i].WLC_force);
        membrane_param.mean_POW_spring_force_magnitude += norm(m_all_nodes[i].POW_force);
        membrane_param.mean_bending_force_magnitude += norm(m_all_nodes[i].bending_force);
        membrane_param.mean_viscous_force_magnitude += norm(m_all_nodes[i].viscous_force);
    }
    membrane_param.mean_WLC_spring_force_magnitude /= num_nodes;
    membrane_param.mean_POW_spring_force_magnitude /= num_nodes;
    membrane_param.mean_bending_force_magnitude /= num_nodes;
    membrane_param.mean_viscous_force_magnitude /= num_nodes;
  }
}



    
//---------------------------------------------------------------------------
void DS_3DRBC::compute_volume_area_conservation_forces( size_t const& dim, 
                                              bool const& Matlab_numbering, 
                                              string const& model_type )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_volume_area_conservation_forces" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
  double order_of_magnitude_of_radius = shape_param.order_of_magnitude_of_radius;
  double kv = membrane_param.k_volume;
  double ka = membrane_param.k_area;
  double kd = membrane_param.k_area_local;

  // Force reversal coefficient
  double force_reversal_coefficient = (Matlab_numbering) ? -1. : 1.;

  geomVector a21(dim), a31(dim), a13(dim), a32(dim), xi(dim);
  geomVector node_1_comp_volume(dim), node_2_comp_volume(dim), node_3_comp_volume(dim);
  geomVector node_1_comp_area(dim), node_2_comp_area(dim), node_3_comp_area(dim);
  double f1, f2, f3, f4, f5, f6;
  geomVector triangle_centroid(dim);
  double beta_v;
  geomVector a21_M(dim), a13_M(dim), a32_M(dim);
  geomVector xi_M(dim);
  geomVector center_of_mass_M(dim);
  
  // Initialising volume & area forces to 0.0
  for (size_t i=0;i<num_nodes;++i)
  {
    for (size_t j=0;j<dim;++j)
    {
      m_all_nodes[i].volume_force(j) = 0.;
      m_all_nodes[i].area_force(j) = 0.;
    }
  }
  
  // Nullifying mean volume & area force magnitudes
  membrane_param.mean_volume_force_magnitude = 0.;
  membrane_param.mean_area_force_magnitude = 0.;
          
  // Computing area forces based on looping over triangles
  for (size_t i=0;i<m_nTriangles;++i) 
  {
      size_t n1_num = m_all_trielements[i].pnodes[0]->number; // Node 1's number
      size_t n2_num = m_all_trielements[i].pnodes[1]->number; // Node 2's number
      size_t n3_num = m_all_trielements[i].pnodes[2]->number; // Node 3's number
      
      // Get each triangle's center of mass vector pointing from membrane centroid to triangle center of mass
      for (size_t j=0;j<dim;++j) 
        triangle_centroid(j) = (m_all_trielements[i].center_of_mass(j) 
                                - 
                                membrane_param.centroid_coordinates(j)) 
                                / order_of_magnitude_of_radius;
      
      // a21
      for (size_t j=0;j<dim;++j)
          a21(j) = (m_all_trielements[i].pnodes[0]->coordinates(j) 
                    - 
                    m_all_trielements[i].pnodes[1]->coordinates(j))
                    / order_of_magnitude_of_radius;
          
      // a31
      for (size_t j=0;j<dim;++j)
          a13(j) = (m_all_trielements[i].pnodes[2]->coordinates(j) 
                    - 
                    m_all_trielements[i].pnodes[0]->coordinates(j))
                    / order_of_magnitude_of_radius;
      
      // a32
      for (size_t j=0;j<dim;++j)
          a32(j) = (m_all_trielements[i].pnodes[1]->coordinates(j) 
                    - 
                    m_all_trielements[i].pnodes[2]->coordinates(j))
                    / order_of_magnitude_of_radius;
          

      // outward pointing normal to the triangle
      for (size_t j=0;j<dim;++j) 
        xi(j) = (m_all_trielements[i].twice_area_outwards_normal_vector(j)) 
                / pow(order_of_magnitude_of_radius, 2.);

      
      // Converting quantities to non-dimensional units
      for (size_t j=0;j<dim;++j)
      {
          xi_M(j) = xi(j); // / pow(order_of_magnitude_of_radius, 2);
          a21_M(j) = a21(j); // / order_of_magnitude_of_radius;
          a13_M(j) = a13(j); // / order_of_magnitude_of_radius;
          a32_M(j) = a32(j); // / order_of_magnitude_of_radius;
          center_of_mass_M(j) = triangle_centroid(j); // - membrane_param.centroid_coordinates(j)) / order_of_magnitude_of_radius;
      }
      
      // t_c x a32, tc x a13, tc x a21
      cross_3D(center_of_mass_M, a32_M, node_1_comp_volume);
      cross_3D(center_of_mass_M, a13_M, node_2_comp_volume);
      cross_3D(center_of_mass_M, a21_M, node_3_comp_volume);
      
      //----------------------------------//
      /**** Volume conservation force  ****/
      //----------------------------------//
      force_reversal_coefficient = 1.;
      
      // Nodal volume conservation force
      beta_v = - kv * ((membrane_param.total_volume 
                        - 
                        membrane_param.initial_volume)
                        /membrane_param.initial_volume);
      // Node 1
      for (size_t j=0;j<dim;++j)
      {
          f1 = (beta_v/6.) * ((1./3.)*xi_M(j) + node_1_comp_volume(j) * force_reversal_coefficient);
          m_all_trielements[i].pnodes[0]->sumforce(j) += f1;
          m_all_trielements[i].pnodes[0]->volume_force(j) += f1;
      }
      // Node 2
      for (size_t j=0;j<dim;++j)
      {
          f2 = (beta_v/6.) * ((1./3.)*xi_M(j) + node_2_comp_volume(j) * force_reversal_coefficient);
          m_all_trielements[i].pnodes[1]->sumforce(j) += f2;
          m_all_trielements[i].pnodes[1]->volume_force(j) += f2;
      }
      // Node 3
      for (size_t j=0;j<dim;++j)
      {
          f3 = (beta_v/6.) * ((1./3.)*xi_M(j) + node_3_comp_volume(j) * force_reversal_coefficient);
          m_all_trielements[i].pnodes[2]->sumforce(j) += f3;
          m_all_trielements[i].pnodes[2]->volume_force(j) += f3;
      }
        
        
        
      //--------------------------------//
      /**** Area conservation force ****/
      //--------------------------------//
      // Global area force contribution
      double ratio_global = (membrane_param.total_area - membrane_param.initial_area) 
                            / membrane_param.initial_area;
      double beta_global = - ka * ratio_global;
      double A_k = m_all_trielements[i].tri_area / pow(order_of_magnitude_of_radius, 2.);
      double global = beta_global / ( 4. * A_k );
      
      // Local area force contribution
      double ratio_local = (m_all_trielements[i].tri_area - m_all_trielements[i].tri_initial_area) 
                           / m_all_trielements[i].tri_initial_area;
      double beta_local = - kd * ratio_local;
      double local = beta_local / ( 4. * A_k );
      
      for (size_t j=0;j<dim;++j)
      {
          xi_M(j) = xi(j); // / pow(order_of_magnitude_of_radius, 2.);
          a21_M(j) = a21(j); // / order_of_magnitude_of_radius;
          a13_M(j) = a13(j); // / order_of_magnitude_of_radius;
          a32_M(j) = a32(j); // / order_of_magnitude_of_radius;
      }

      // xi x a32, xi x a13, xi x a21
      cross_3D( xi_M, a32_M, node_1_comp_area );
      cross_3D( xi_M, a13_M, node_2_comp_area );
      cross_3D( xi_M, a21_M, node_3_comp_area );
      
      force_reversal_coefficient = 1.;
      
      // Nodal area conservation force
      // Node 1
      for (size_t j=0;j<dim;++j)
      {
          f4 = (global + local) * node_1_comp_area(j) * force_reversal_coefficient;
          m_all_trielements[i].pnodes[0]->sumforce(j) += f4; // node 1
          m_all_trielements[i].pnodes[0]->area_force(j) += f4; // node 1
      }
      // Node 2
      for (size_t j=0;j<dim;++j)
      {
          f5 = (global + local) * node_2_comp_area(j) * force_reversal_coefficient;
          m_all_trielements[i].pnodes[1]->sumforce(j) += f5; // node 2
          m_all_trielements[i].pnodes[1]->area_force(j) += f5; // node 2
      }
      // Node 3
      for (size_t j=0;j<dim;++j)
      {
          f6 = (global + local) * node_3_comp_area(j) * force_reversal_coefficient;
          m_all_trielements[i].pnodes[2]->sumforce(j) += f6; // node 3
          m_all_trielements[i].pnodes[2]->area_force(j) += f6; // node 3
      }
  }
  
  
  // Mean force magnitude
  for (size_t i=0;i<num_nodes;++i)
  {
      membrane_param.mean_volume_force_magnitude += norm(m_all_nodes[i].volume_force);
      membrane_param.mean_area_force_magnitude += norm(m_all_nodes[i].area_force);
  }
  membrane_param.mean_volume_force_magnitude /= num_nodes;
  membrane_param.mean_area_force_magnitude /= num_nodes;
  
}



    
//---------------------------------------------------------------------------
void DS_3DRBC::compute_spring_force( size_t const& dim,
                                     double const& spring_constant,
                                     string const& model_type )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_spring_force" ) ;
  
  size_t num_nodes = shape_param.N_nodes;

  geomVector Iij(dim);
  double length = 0.;
  
  if(model_type.compare("Simplified") == 0)
  {
    for (size_t i=0;i<num_nodes;++i)
    {
      size_t nn = m_all_nodes[i].neighbors_3D.size();
      for (size_t j=0;j<nn;++j)
      {
        for (size_t k=0;k<dim;++k) 
          Iij(k) = m_all_nodes[i].neighbors_3D[j]->coordinates(k) 
                   - m_all_nodes[i].coordinates(k);
                   
        length = norm( Iij );
        
        for (size_t k=0;k<dim;++k) Iij(k) /= length;
        
        double initial_spring_length = m_all_nodes[i].initial_spring_length[j];
        
        for (size_t k=0;k<dim;++k)
          m_all_nodes[i].sumforce(k) += spring_constant 
                * ( length - initial_spring_length ) * Iij(k);
      }
    }
  }
  else if(model_type.compare("NumericalMembraneModel") == 0)
  {
    double order_of_magnitude_of_radius = shape_param.order_of_magnitude_of_radius;
    
    geomVector f_WLC(dim);
    geomVector f_POW(dim);
    
    membrane_param.mean_WLC_spring_force_magnitude = 0.;
    membrane_param.mean_POW_spring_force_magnitude = 0.;
    
    for (size_t i=0;i<num_nodes;++i)
    {
        for (size_t j=0;j<dim;++j)
        {
            m_all_nodes[i].WLC_force(j) = 0.;
            m_all_nodes[i].POW_force(j) = 0.;
        }
    }
        
    for (vector<Edge>::iterator il=m_all_edges.begin();il!=m_all_edges.end();il++)
    {
        // Computing the spring vector from node 3 to node 2, i.e., node j to i
        for (size_t j=0;j<dim;++j)
          Iij(j) = il->n2->coordinates(j) - il->n3->coordinates(j); 

        // Spring length
        length = norm ( Iij );

        // Normalizing the spring vector to unit vector
        for (size_t j=0;j<dim;++j)
          Iij(j) /= length;
        
        // WLC spring prefactor
        double x = length / il->lmax;
        double alpha_1 = 1. / (4. * pow(1.-x, 2.));
        double alpha_2 = 1./4.;
        double alpha_3 = x;
        double alpha = alpha_1 - alpha_2 + alpha_3;

        // POW spring prefactor
        double l_M = length / order_of_magnitude_of_radius;
        double beta = 1. / pow(l_M, 2.);

        for (size_t j=0;j<dim;++j)
        {
            // WLC force
            f_WLC(j) = - il->k * alpha * Iij(j);
            il->n2->sumforce(j)  +=  f_WLC(j);
            il->n3->sumforce(j)  += -f_WLC(j);
            il->n2->WLC_force(j) +=  f_WLC(j);
            il->n3->WLC_force(j) += -f_WLC(j);

            // POW force
            f_POW(j) = il->kp * beta * Iij(j);
            il->n2->sumforce(j)  +=  f_POW(j);
            il->n3->sumforce(j)  += -f_POW(j);
            il->n2->WLC_force(j) +=  f_POW(j);
            il->n3->WLC_force(j) += -f_POW(j);
        }
    }
    
    for (size_t i=0;i<num_nodes;++i)
    {
        membrane_param.mean_WLC_spring_force_magnitude += norm(m_all_nodes[i].WLC_force);
        membrane_param.mean_POW_spring_force_magnitude += norm(m_all_nodes[i].POW_force);
    }
    
    // Mean spring force
    membrane_param.mean_WLC_spring_force_magnitude /= num_nodes;
    membrane_param.mean_POW_spring_force_magnitude /= num_nodes;
  }
}



    
//---------------------------------------------------------------------------
void DS_3DRBC:: compute_linear_spring_force( size_t const& dim, 
                                             double const& spring_constant )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_linear_spring_force" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
  geomVector Iij(dim);
  double length = 0.;
  
  for (size_t inode=0;inode<num_nodes;++inode)
  {
    size_t nn = m_all_nodes[inode].neighbors_3D.size();
    
    for (size_t j=0;j<nn;++j)
    {
      // spring force vector
      for (size_t k=0;k<dim;++k) 
        Iij(k) = m_all_nodes[inode].neighbors_3D[j]->coordinates(k) 
                 - m_all_nodes[inode].coordinates(k);
                 
      // spring length
      length = norm(Iij);
      
      // normalised spring force vector
      for (size_t k=0;k<dim;++k)
        Iij(k) /= length;
      
      // initial spring length
      double initial_spring_length = m_all_nodes[inode].initial_spring_length[j];
      
      // tension along the spring
      double tension = spring_constant 
                       * ( (length / initial_spring_length) - 1. );
      
      // Compute spring force
      for (size_t k=0;k<dim;++k)
        m_all_nodes[inode].sumforce(k) += tension * Iij(k);
    }
  }
} 




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_bending_resistance( size_t const& dim, 
                                        double const& bending_spring_constant,
                                        double const& bending_viscous_constant, 
                                        double const& dt,
                                        string const& model_type )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_bending_resistance" ) ;
  
  geomVector a21(dim), a31(dim), a24(dim), a34(dim);
  geomVector Chi(dim), Gamma(dim), crossp(dim), n2n3(dim);
  geomVector Chi_minus_Gamma(dim), tc1_minus_tc2(dim), tc1(dim), tc2(dim);
  double angle = 0., dot = 0., mi = 0.;
  double lever1 = 0., lever4 = 0., nxi = 0., nzeta = 0.;
  
  size_t num_nodes = shape_param.N_nodes;


  if(model_type.compare("Simplified") == 0)
  {
    for (vector<Edge>::iterator il=m_all_edges.begin();il!=m_all_edges.end();il++)
    {
      // a21
      for (size_t j=0;j<dim;++j)
        a21(j) =  il->t1v1.second->coordinates(j) - il->n2->coordinates(j);

      // a31
      for (size_t j=0;j<dim;++j)
        a31(j) = il->t1v1.second->coordinates(j) - il->n3->coordinates(j);

      // Compute cross product Chi = a21 x a31
      cross_3D( a21, a31, Chi );
      nxi = norm(Chi);
      for (size_t j=0;j<dim;++j) Chi(j) /= nxi;

      // a34
      for (size_t j=0;j<dim;++j)
        a34(j) =  il->t2v4.second->coordinates(j) - il->n3->coordinates(j);

      // a24
      for (size_t j=0;j<dim;++j)
        a24(j) = il->t2v4.second->coordinates(j) - il->n2->coordinates(j);

      // Compute cross product Gamma = a34 x a24
      cross_3D( a34, a24, Gamma );
      nzeta = norm(Gamma);
      for (size_t j=0;j<dim;++j) Gamma(j) /= nzeta;        


      double A1 = 0.5 * norm(Chi);
      double A2 = 0.5 * norm(Gamma);
      
      // Compute chi - gamma
      geomVector Chi_minus_Gamma(dim);
      for (size_t j=0;j<dim;++j) Chi_minus_Gamma(j) = Chi(j) - Gamma(j);
      

      // Compute tc1 - tc2
      geomVector tc1(dim), tc2(dim);
      for (size_t j=0;j<dim;++j) tc1(j) = (il->t1v1.second->coordinates(j) + il->n2->coordinates(j) + il->n3->coordinates(j))/3.;
      for (size_t j=0;j<dim;++j) tc2(j) = (il->n2->coordinates(j) + il->n3->coordinates(j) + il->t2v4.second->coordinates(j))/3.;
      geomVector tc1_minus_tc2(dim);
      for (size_t j=0;j<dim;++j) tc1_minus_tc2(j) = tc1(j) - tc2(j); // /order_of_magnitude_of_radius;
      double dot_product = scalar(Chi_minus_Gamma, tc1_minus_tc2);


      // Compute costheta
      double costheta = scalar(Chi, Gamma)/A1/A2/4.0;
      if(costheta >  1.) costheta =  1.;
      if(costheta < -1.) costheta = -1.;
      angle = acos(costheta);
      if(dot_product < 0.) angle *= -1.;

      // Compute sintheta (with a clipping for values above 0.5 or so)
      double sintheta = sin(angle);
      if(fabs(sintheta) < SMALLER) sintheta = SMALLER;
      
      il->angle_nm1 = angle;
      il->dangledt = 0.;
      
      /*
      // Compute angle
      dot = scalar( Chi, Gamma ) ;
      cross_3D( Chi, Gamma, crossp );
      for (size_t j=0;j<dim;++j)
        n2n3(j) =  il->n3->coordinates(j) - il->n2->coordinates(j);
      double vnxizeta = scalar( crossp, n2n3 );    
      angle = ( vnxizeta > 0. ? 1. : -1. ) 
              * acos( max( - 1., min( il->costheta0, 1. ) ) );
      */

      // Bending moment
      mi = bending_spring_constant * ( angle - il->initial_angle );    

      // Bending spring forces
      // Triangle 1
      lever1 = 0.5 * ( norm( a21 ) + norm( a31 ) );
      for (size_t j=0;j<dim;++j)
      {
        il->t1v1.second->sumforce(j) += ( mi / lever1 ) * Chi(j);
        il->n2->sumforce(j) -= 0.5 * ( mi / lever1 ) * Chi(j);
        il->n3->sumforce(j) -= 0.5 * ( mi / lever1 ) * Chi(j);
      }

      // Triangle 2
      lever4 = 0.5 * ( norm( a34 ) + norm( a24 ) );
      for (size_t j=0;j<dim;++j)
      {
        il->t2v4.second->sumforce(j) += ( mi / lever4 ) * Gamma(j);
        il->n2->sumforce(j) -= 0.5 * ( mi / lever4 ) * Gamma(j);
        il->n3->sumforce(j) -= 0.5 * ( mi / lever4 ) * Gamma(j);
      }
    }
  }
  else if(model_type.compare("NumericalMembraneModel") == 0)
  {
    double order_of_magnitude_of_radius = shape_param.order_of_magnitude_of_radius;
    
    // Setting bending_force variable to 0.0
    for (size_t inode=0;inode<num_nodes;++inode)
        for (size_t j=0;j<dim;++j)
            m_all_nodes[inode].bending_force(j) = 0.;
            
    
    for (vector<Edge>::iterator il=m_all_edges.begin();il!=m_all_edges.end();il++)
    {
        // Steps for computing bending force
        
        //  1. Compute a21, a31
        // a21 = a1 - a2 = vector from node 2 to node 1
        for (size_t j=0;j<dim;++j) 
          a21(j) = (-il->t1v1.second->coordinates(j) + il->n2->coordinates(j))/order_of_magnitude_of_radius;
        
        // a31 = a1 - a3 = vector from node 3 to node 1
        for (size_t j=0;j<dim;++j) 
          a31(j) = (-il->t1v1.second->coordinates(j) + il->n3->coordinates(j))/order_of_magnitude_of_radius;
        
        
        //  2. Compute chi
        cross_3D(a21, a31, Chi);
        
        
        //  3. Compute a34, a24
        // a34 = a4 - a3 = vector from node 3 to node 4
        for (size_t j=0;j<dim;++j) 
          a34(j) = (-il->t2v4.second->coordinates(j) + il->n3->coordinates(j))/order_of_magnitude_of_radius;
        
        // a24 = a4 - a2 = vector from node 2 to node 4
        for (size_t j=0;j<dim;++j) 
          a24(j) = (-il->t2v4.second->coordinates(j) + il->n2->coordinates(j))/order_of_magnitude_of_radius;
        
        
        //  4. Compute gamma
        cross_3D(a34, a24, Gamma);
            
            
        double A1 = 0.5 * norm(Chi);
        double A2 = 0.5 * norm(Gamma);
            
        //  6. Compute chi - gamma
        for (size_t j=0;j<dim;++j) Chi_minus_Gamma(j) = Chi(j) - Gamma(j);
            
        
        //  7. Compute tc1 - tc2
        for (size_t j=0;j<dim;++j) tc1(j) = (il->t1v1.second->coordinates(j) + il->n2->coordinates(j) + il->n3->coordinates(j))/3.;
        for (size_t j=0;j<dim;++j) tc2(j) = (il->n2->coordinates(j) + il->n3->coordinates(j) + il->t2v4.second->coordinates(j))/3.;
        for (size_t j=0;j<dim;++j) tc1_minus_tc2(j) = (tc1(j) - tc2(j))/order_of_magnitude_of_radius;
        double dot_product = scalar(Chi_minus_Gamma, tc1_minus_tc2);


        //  5. Compute costheta
        double costheta = scalar(Chi, Gamma)/A1/A2/4.0;
            
        
        //  6. Compute sintheta (with a clipping for values above 0.5 or so)
        if(costheta >  1.) costheta =  1.;
        if(costheta < -1.) costheta = -1.;
        double theta = acos(costheta);
        if(dot_product < 0.) theta *= -1.;
        double sintheta = sin(theta);
        if(fabs(sintheta) < SMALLER) sintheta = SMALLER;
        double siinv = 1. / sintheta;
            



            
        //  7. Compute a12, a13, a43, a42, a32, a23
        // a12 = a2 - a1 = vector from node 1 to node 2
        geomVector a12(dim);
        for (size_t j=0;j<dim;++j) a12(j) = (-il->n2->coordinates(j) + il->t1v1.second->coordinates(j))/order_of_magnitude_of_radius;
        
        // a13 = a3 - a1 = vector from node node 1 to node 3
        geomVector a13(dim);
        for (size_t j=0;j<dim;++j) a13(j) = (-il->n3->coordinates(j) + il->t1v1.second->coordinates(j))/order_of_magnitude_of_radius;
        
        // a43 = a3 - a4 = vector from node 4 to node 3
        geomVector a43(dim);
        for (size_t j=0;j<dim;++j) a43(j) = (-il->n3->coordinates(j) + il->t2v4.second->coordinates(j))/order_of_magnitude_of_radius;
        
        // a42 = a2 - a4 = vector from node 4 to node 2
        geomVector a42(dim);
        for (size_t j=0;j<dim;++j) a42(j) = (-il->n2->coordinates(j) + il->t2v4.second->coordinates(j))/order_of_magnitude_of_radius;
        
        // a32 = a2 - a3 = vector from node 3 to node 2
        geomVector a32(dim);
        for (size_t j=0;j<dim;++j) a32(j) = (-il->n2->coordinates(j) + il->n3->coordinates(j))/order_of_magnitude_of_radius;
        
        // a23 = a3 - a2 = vector from node to node 3
        geomVector a23(dim);
        for (size_t j=0;j<dim;++j) a23(j) = (-il->n3->coordinates(j) + il->n2->coordinates(j))/order_of_magnitude_of_radius;
            
            
        //  8. Compute beta_b
        double beta_b = membrane_param.kbending_M * (sintheta * il->costheta0 - costheta * il->sintheta0) * siinv;
        
        //  9. Compute b11, b12, b22
        double b11 = - beta_b*costheta/A1/A1/4.0;
        double b12 =   beta_b/A1/A2/4.0;
        double b22 = - beta_b*costheta/A2/A2/4.0;
            

        //  10. Compute nodal forces
        // Node 1
        geomVector Chi_cross_a32(dim);
        cross_3D(Chi, a32, Chi_cross_a32);
        geomVector Gamma_cross_a32(dim);
        cross_3D(Gamma, a32, Gamma_cross_a32);
        for (size_t j=0;j<dim;++j) il->t1v1.second->sumforce(j) += b11 * Chi_cross_a32(j) + b12 * Gamma_cross_a32(j);
        for (size_t j=0;j<dim;++j) il->t1v1.second->bending_force(j) += b11 * Chi_cross_a32(j) + b12 * Gamma_cross_a32(j);
        
        // Node 2
        geomVector Chi_cross_a13(dim);
        cross_3D(Chi, a13, Chi_cross_a13);
        geomVector Chi_cross_a34(dim);
        cross_3D(Chi, a34, Chi_cross_a34);
        geomVector Gamma_cross_a13(dim);
        cross_3D(Gamma, a13, Gamma_cross_a13);
        geomVector Gamma_cross_a34(dim);
        cross_3D(Gamma, a34, Gamma_cross_a34);
        for (size_t j=0;j<dim;++j) il->n2->sumforce(j) += b11 * Chi_cross_a13(j) + b12 * (Chi_cross_a34(j) + Gamma_cross_a13(j)) + b22 * Gamma_cross_a34(j);
        for (size_t j=0;j<dim;++j) il->n2->bending_force(j) += b11 * Chi_cross_a13(j) + b12 * (Chi_cross_a34(j) + Gamma_cross_a13(j)) + b22 * Gamma_cross_a34(j);
        
        // Node 3
        geomVector Chi_cross_a21(dim);
        cross_3D(Chi, a21, Chi_cross_a21);
        geomVector Chi_cross_a42(dim);
        cross_3D(Chi, a42, Chi_cross_a42);
        geomVector Gamma_cross_a21(dim);
        cross_3D(Gamma, a21, Gamma_cross_a21);
        geomVector Gamma_cross_a42(dim);
        cross_3D(Gamma, a42, Gamma_cross_a42);
        for (size_t j=0;j<dim;++j) il->n3->sumforce(j) += b11 * Chi_cross_a21(j) + b12 * (Chi_cross_a42(j) + Gamma_cross_a21(j)) + b22 * Gamma_cross_a42(j);
        for (size_t j=0;j<dim;++j) il->n3->bending_force(j) += b11 * Chi_cross_a21(j) + b12 * (Chi_cross_a42(j) + Gamma_cross_a21(j)) + b22 * Gamma_cross_a42(j);
        
        // Node 4
        geomVector Chi_cross_a23(dim);
        cross_3D(Chi, a23, Chi_cross_a23);
        geomVector Gamma_cross_a23(dim);
        cross_3D(Gamma, a23, Gamma_cross_a23);
        for (size_t j=0;j<dim;++j) il->t2v4.second->sumforce(j) += b12 * Chi_cross_a23(j) + b22 * Gamma_cross_a23(j);
        for (size_t j=0;j<dim;++j) il->t2v4.second->bending_force(j) += b12 * Chi_cross_a23(j) + b22 * Gamma_cross_a23(j);
    }
    
    // Compute mean bending force magnitude
    membrane_param.mean_bending_force_magnitude = 0.;
    for (size_t inode=0;inode<num_nodes;++inode)
        membrane_param.mean_bending_force_magnitude += norm(m_all_nodes[inode].bending_force);
    membrane_param.mean_bending_force_magnitude /= num_nodes;
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_viscous_drag_force( size_t const& dim,
                                            double const& viscous_drag_constant,
                                            double const& dt,
                                            string const& model_type )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_viscous_drag_force" ) ;
    
  if(model_type.compare("Simplified") == 0)
  {
    // This is not Fedosov's model but a simplified model of isotropic 
    // compression/extension of the surface area using the vectors
    // connection the nodes to the center of mass of each triangle

    double coef = 0., norm_unit_relpos = 0.;
    geomVector unit_relpos(dim);
    for (size_t i=0;i<m_nTriangles;++i)
    {
      coef = membrane_param.k_area 
             * ( m_all_trielements[i].tri_area - m_all_trielements[i].tri_initial_area ) 
             / m_all_trielements[i].tri_initial_area ;
             
      for (size_t k=0;k<dim;++k)
      {
        for (size_t j=0;j<dim;++j)
          unit_relpos(j) = m_all_trielements[i].center_of_mass(j) 
                           - m_all_trielements[i].pnodes[k]->coordinates(j);
          
        norm_unit_relpos = norm( unit_relpos );
        
        for (size_t j=0;j<dim;++j) 
          unit_relpos(j) /= norm_unit_relpos;

        for (size_t j=0;j<dim;++j)
          m_all_trielements[i].pnodes[k]->sumforce(j) += coef * unit_relpos(j);
      }
    }
  }
  else if(model_type.compare("NumericalMembraneModel") == 0)
  {
    size_t num_nodes = shape_param.N_nodes;
    double order_of_magnitude_of_radius = shape_param.order_of_magnitude_of_radius;
    double eta_M = membrane_param.eta_M;

    geomVector eij(dim);
    geomVector vij(dim);
    geomVector vij_M(dim);

    membrane_param.mean_viscous_force_magnitude = 0.;
    
    // (x) ------------------- (x)
    //  i                       j
    // Nodes i and j of a spring
    // Computing viscous force as two lines of code for nodes i and j
    geomVector vjmvi(dim);
    geomVector vimvj(dim);

    // Making all the viscous force to 0.0
    for (size_t i=0;i<num_nodes;++i)
        for (size_t j=0;j<dim;++j) 
            m_all_nodes[i].viscous_force(j) = 0.;
        
    for (vector<Edge>::iterator il=m_all_edges.begin();il!=m_all_edges.end();il++)
    {
        // Vector of each spring
        for (size_t j=0;j<dim;++j)
            eij(j) = il->n3->coordinates(j) - il->n2->coordinates(j);
        
        // Length of each spring
        double length = norm ( eij );
        
        // Unit vector of each spring
        for (size_t j=0;j<dim;++j) 
            eij(j) /= length;

        // Relative velocity between nodes i and j
        for (size_t j=0;j<dim;++j)
        {
            vimvj(j) = (il->n2->velocity(j) - il->n3->velocity(j)) * (dt/order_of_magnitude_of_radius);
            vjmvi(j) = (il->n3->velocity(j) - il->n2->velocity(j)) * (dt/order_of_magnitude_of_radius);
        }
    
        // Projection of relative velocity on the spring unit vector
        double vimvj_dot_eij = scalar(vimvj, eij);
        double vjmvi_dot_eij = scalar(vjmvi, eij);
        
        // Viscous force computation
        double constant_1 = 12.0 / (13.0 * sqrt(3.0));
        double constant_2 =  4.0 / (13.0 * sqrt(3.0)); // gamma_C = gamma_T/3
        double gamma_T = constant_1 * eta_M; // T = translational component
        double gamma_C = constant_2 * eta_M; // C = central component
        for (size_t j=0;j<dim;++j)
        {
            // Viscous force on node i
            double f_i = - (gamma_T * vimvj(j) + gamma_C * vimvj_dot_eij * eij(j));
            il->n2->sumforce(j) += f_i;
            il->n2->viscous_force(j) += f_i;

            // Viscous force on node j
            double f_j = - (gamma_T * vjmvi(j) + gamma_C * vjmvi_dot_eij * eij(j));
            il->n3->sumforce(j) += f_j;
            il->n3->viscous_force(j) += f_j;
        }
    }

    // Computing mean viscous force magnitude
    for (size_t i=0;i<num_nodes;++i)
        membrane_param.mean_viscous_force_magnitude += norm(m_all_nodes[i].viscous_force);
    
    // Mean viscous force
    membrane_param.mean_viscous_force_magnitude /= num_nodes;
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_volume_conservation_force(size_t const& dim,
                                                  string const& model_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_volume_conservation_force" ) ;

  geomVector a21(dim), a31(dim), a13(dim), a32(dim), xi(dim);
  geomVector node_1_comp(dim), node_2_comp(dim), node_3_comp(dim);
  double f1, f2, f3;
  geomVector triangle_centroid(dim);
  double beta_v;
  
  size_t num_nodes = shape_param.N_nodes;
  double order_of_magnitude_of_radius = shape_param.order_of_magnitude_of_radius;
  double kv = membrane_param.k_volume;
  
  // Initialising all area forces to 0.0
  for (size_t i=0;i<num_nodes;++i)
      for (size_t j=0;j<dim;++j)
          m_all_nodes[i].volume_force(j) = 0.;
          
  // Computing area forces based on looping over triangles
  for (size_t i=0;i<m_nTriangles;++i) 
  {
      size_t n1_num = m_all_trielements[i].pnodes[0]->number; // Node 1's number
      size_t n2_num = m_all_trielements[i].pnodes[1]->number; // Node 2's number
      size_t n3_num = m_all_trielements[i].pnodes[2]->number; // Node 3's number
      
      // Get each triangle's center of mass vector pointing from membrane centroid to triangle center of mass
      for (size_t j=0;j<dim;++j) triangle_centroid(j) = (m_all_trielements[i].center_of_mass(j) - membrane_param.centroid_coordinates(j)) / order_of_magnitude_of_radius;
      
      
      // a21
      for (size_t j=0;j<dim;++j)
          a21(j) = (m_all_trielements[i].pnodes[0]->coordinates(j) - m_all_trielements[i].pnodes[1]->coordinates(j))/order_of_magnitude_of_radius;
          
      // a31
      for (size_t j=0;j<dim;++j)
          a13(j) = (m_all_trielements[i].pnodes[2]->coordinates(j) - m_all_trielements[i].pnodes[0]->coordinates(j))/order_of_magnitude_of_radius;
      
      // a32
      for (size_t j=0;j<dim;++j)
          a32(j) = (m_all_trielements[i].pnodes[1]->coordinates(j) - m_all_trielements[i].pnodes[2]->coordinates(j))/order_of_magnitude_of_radius;
      
      
      // outward pointing normal to the triangle
      for (size_t j=0;j<dim;++j) xi(j) = (m_all_trielements[i].twice_area_outwards_normal_vector(j)) / pow(order_of_magnitude_of_radius, 2.);

      
      // Converting quantities to non-dimensional units
      geomVector a21_M(dim), a13_M(dim), a32_M(dim);
      geomVector xi_M(dim);
      geomVector center_of_mass_M(dim);
      for (size_t j=0;j<dim;++j)
      {
          xi_M(j) = xi(j); // / pow(order_of_magnitude_of_radius, 2);
          a21_M(j) = a21(j); // / order_of_magnitude_of_radius;
          a13_M(j) = a13(j); // / order_of_magnitude_of_radius;
          a32_M(j) = a32(j); // / order_of_magnitude_of_radius;
          center_of_mass_M(j) = triangle_centroid(j); // - membrane_param.centroid_coordinates(j)) / order_of_magnitude_of_radius;
      }

      
      // t_c x a32, tc x a13, tc x a21
      cross_3D(center_of_mass_M, a32_M, node_1_comp);
      cross_3D(center_of_mass_M, a13_M, node_2_comp);
      cross_3D(center_of_mass_M, a21_M, node_3_comp);
      
      double force_reversal_coefficient = 1.;
      
      // Nodal force
      beta_v = - kv * ((membrane_param.total_volume - membrane_param.initial_volume)/membrane_param.initial_volume);
      // Node 1
      for (size_t j=0;j<dim;++j)
      {
          f1 = (beta_v/6.) * ((1./3.)*xi_M(j) + node_1_comp(j) * force_reversal_coefficient);
          m_all_trielements[i].pnodes[0]->sumforce(j) += f1;
          m_all_trielements[i].pnodes[0]->volume_force(j) += f1;
      }
      // Node 2
      for (size_t j=0;j<dim;++j)
      {
          f2 = (beta_v/6.) * ((1./3.)*xi_M(j) + node_2_comp(j) * force_reversal_coefficient);
          m_all_trielements[i].pnodes[1]->sumforce(j) += f2;
          m_all_trielements[i].pnodes[1]->volume_force(j) += f2;
      }
      // Node 3
      for (size_t j=0;j<dim;++j)
      {
          f3 = (beta_v/6.) * ((1./3.)*xi_M(j) + node_3_comp(j) * force_reversal_coefficient);
          m_all_trielements[i].pnodes[2]->sumforce(j) += f3;
          m_all_trielements[i].pnodes[2]->volume_force(j) += f3;
      }
  }
  
  // magnitude of volume force on node "i"
  membrane_param.mean_volume_force_magnitude = 0.;
  for (size_t i=0;i<num_nodes;++i)
      membrane_param.mean_volume_force_magnitude += norm(m_all_nodes[i].volume_force);
  
  // Mean area force
  membrane_param.mean_volume_force_magnitude /= num_nodes;
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_area_conservation_force(size_t const& dim,
                                                bool const& Matlab_numbering,
                                                string const& model_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_area_conservation_force" ) ;

  size_t num_nodes = shape_param.N_nodes;
  double order_of_magnitude_of_radius = shape_param.order_of_magnitude_of_radius;
  double ka = membrane_param.k_area;
  double kd = membrane_param.k_area_local;

  geomVector a21(dim), a31(dim), a13(dim), a32(dim), xi(dim);
  geomVector node_1_comp(dim), node_2_comp(dim), node_3_comp(dim);
  double f1, f2, f3;
  

  // Force reversal coefficient
  double force_reversal_coefficient = (Matlab_numbering) ? -1. : 1.;

  
  // Initialising all area forces to 0.0
  for (size_t i=0;i<num_nodes;++i)
      for (size_t j=0;j<dim;++j) 
          m_all_nodes[i].area_force(j) = 0.;
  
  // Computing area forces based on looping over triangles
  for (size_t i=0;i<m_nTriangles;++i) 
  {
      size_t n1_num = m_all_trielements[i].pnodes[0]->number; // Node 1's number/index
      size_t n2_num = m_all_trielements[i].pnodes[1]->number; // Node 2's number/index
      size_t n3_num = m_all_trielements[i].pnodes[2]->number; // Node 3's number/index
      
      // a21
      for (size_t j=0;j<dim;++j)
          a21(j) = (m_all_trielements[i].pnodes[0]->coordinates(j) 
                    - 
                    m_all_trielements[i].pnodes[1]->coordinates(j))
                    / order_of_magnitude_of_radius;
          
      // a31
      for (size_t j=0;j<dim;++j)
          a13(j) = (m_all_trielements[i].pnodes[2]->coordinates(j) 
                    - m_all_trielements[i].pnodes[0]->coordinates(j))
                    / order_of_magnitude_of_radius;
      
      // a32
      for (size_t j=0;j<dim;++j)
          a32(j) = (m_all_trielements[i].pnodes[1]->coordinates(j) 
                    - 
                    m_all_trielements[i].pnodes[2]->coordinates(j))
                    / order_of_magnitude_of_radius;
          
      // outward pointing normal to the triangle
      for (size_t j=0;j<dim;++j) 
        xi(j) = (m_all_trielements[i].twice_area_outwards_normal_vector(j)) 
                / pow(order_of_magnitude_of_radius, 2.);
      
      
      // Global area force contribution
      double ratio_global = (membrane_param.total_area - membrane_param.initial_area) 
                            / membrane_param.initial_area;
      double beta_global = - ka * ratio_global;
      double A_k = m_all_trielements[i].tri_area / pow(order_of_magnitude_of_radius, 2.); // scaling the area of the triangle to model length units with 1e-6 as model length
      double global = beta_global / ( 4. * A_k );
      
      // Local area force contribution
      double ratio_local = (m_all_trielements[i].tri_area - m_all_trielements[i].tri_initial_area) 
                           / m_all_trielements[i].tri_initial_area;
      double beta_local = - kd * ratio_local;
      double local = beta_local / ( 4. * A_k );
      
      // Converting quantities to non-dimensional units
      geomVector a21_M(dim), a13_M(dim), a32_M(dim);
      geomVector xi_M(dim);
      for (size_t j=0;j<dim;++j)
      {
          xi_M(j) = xi(j); // / pow(order_of_magnitude_of_radius, 2.);
          a21_M(j) = a21(j); // / order_of_magnitude_of_radius;
          a13_M(j) = a13(j); // / order_of_magnitude_of_radius;
          a32_M(j) = a32(j); // / order_of_magnitude_of_radius;
      }
      
      // xi x a32, xi x a13, xi x a21
      cross_3D( xi_M, a32_M, node_1_comp );
      cross_3D( xi_M, a13_M, node_2_comp );
      cross_3D( xi_M, a21_M, node_3_comp );
      
      double force_reversal_coefficient = 1.;
      
      // Nodal force
      // Node 1
      for (size_t j=0;j<dim;++j)
      {
          f1 = (global + local) * node_1_comp(j) * force_reversal_coefficient;
          m_all_trielements[i].pnodes[0]->sumforce(j) += f1; // node 1
          m_all_trielements[i].pnodes[0]->area_force(j) += f1; // node 1
      }
      // Node 2
      for (size_t j=0;j<dim;++j)
      {
          f2 = (global + local) * node_2_comp(j) * force_reversal_coefficient;
          m_all_trielements[i].pnodes[1]->sumforce(j) += f2; // node 2
          m_all_trielements[i].pnodes[1]->area_force(j) += f2; // node 2
      }
      // Node 3
      for (size_t j=0;j<dim;++j)
      {
          f3 = (global + local) * node_3_comp(j) * force_reversal_coefficient;
          m_all_trielements[i].pnodes[2]->sumforce(j) += f3; // node 3
          m_all_trielements[i].pnodes[2]->area_force(j) += f3; // node 3
      }
  }
  
  // magnitude of area force on node "i"
  membrane_param.mean_area_force_magnitude = 0.;
  for (size_t i=0;i<num_nodes;++i)
  {
      // magnitude of spring force on node "i"
      membrane_param.mean_area_force_magnitude += norm(m_all_nodes[i].area_force);
  }
  
  // Mean area force
  membrane_param.mean_area_force_magnitude /= num_nodes;
}




//---------------------------------------------------------------------------
// Converts the external force from model units to physical units
double DS_3DRBC:: convert_model_to_physical_units(double f_M)
//---------------------------------------------------------------------------
{
  double mu0_M = membrane_param.mu0_M;
  double D0_M = membrane_param.D0_M;
  double mu0_P = membrane_param.mu0_P;
  double D0_P = membrane_param.D0_P;
  
  // kappa = f_M / (mu0_M * D0_M) = f_P / (mu0_P * D0_P) = const.
  
  // non-dimensional quantity
  double kappa =  f_M / (mu0_M * D0_M);
  
  // external force in physical units
  double f_P = kappa * mu0_P * D0_P;
  
  return ( f_P );
}




//---------------------------------------------------------------------------
void DS_3DRBC:: rbc_dynamics_solver(size_t const& dim, 
                                    size_t const& fluid_iter_num,
                                    double const& dt_fluid, 
                                    string const& case_type,
                                    bool const& Matlab_numbering,
                                    bool const& combined_force_computation,
                                    string const& model_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: rbc_dynamics_solver" ) ;

  size_t num_nodes = shape_param.N_nodes;
  double node_mass = membrane_param.node_mass_M; // TO BE VERIFIED OF ITS VALUE
  double spring_constant = membrane_param.k_spring;
  double bending_spring_constant = membrane_param.k_bending;
  double bending_viscous_constant = membrane_param.k_bending_visc;
  double viscous_drag_constant = membrane_param.k_viscous;
  double area_force_constant = membrane_param.k_area;
  double volume_force_constant = membrane_param.k_volume;
  size_t n_sub_timesteps = membrane_param.n_subtimesteps_RBC;
  
  double time = 0.;
  double dt = dt_fluid / n_sub_timesteps;
  membrane_param.total_kinetic_energy = 0.;
  
  MAC_Communicator const* MAC_comm;
  MAC_comm = MAC_Exec::communicator();
  size_t my_rank = MAC_comm->rank();
  size_t nb_procs = MAC_comm->nb_ranks();
  size_t is_master = 0;
  
  bool compute_init_values; 
  compute_init_values = (int(fluid_iter_num) == 1) ? true : false;
  
  // Compute initial area & volume of membrane
  compute_total_surface_area_total_volume(true, dim);
  
  // Time loop
  for (size_t iter_num=0;iter_num<n_sub_timesteps;++iter_num)
  {
    // Time
    time += dt;
    
    // Initialize forces on all nodes
    for (size_t inode=0;inode<num_nodes;++inode)
      for (size_t j=0;j<dim;++j)
        m_all_nodes[inode].sumforce(j) = 0.0;
        
    // Compute normals, areas and center of mass of triangles
    compute_triangle_area_normals_centre_of_mass(false, dim);
    
    // Compute total surface area and total volume and write to file
    compute_total_surface_area_total_volume(false, dim);
    
    // Compute forces on all nodes
    // Spring force
    compute_spring_force( dim, spring_constant, model_type );

    // Bending resistance
    compute_bending_resistance( dim, bending_spring_constant, 
                                bending_viscous_constant, dt, model_type );

    // Viscous drag force
    compute_viscous_drag_force( dim, viscous_drag_constant, dt, model_type );

    // Volume conservation force
    compute_volume_conservation_force( dim, model_type );

    // Triangle surface area conservation force
    compute_area_conservation_force( dim, Matlab_numbering, model_type );


    membrane_param.total_kinetic_energy = 0.;

    // Compute new velocity and position
    for (size_t i=0;i<num_nodes;++i)
    {
      // Solve momentum conservation
      if ( !iter_num ) // First order explicit at the 1st time loop iteration
      {
        for (size_t j=0;j<dim;++j)
          m_all_nodes[i].velocity(j) += ( dt /  node_mass ) 
                                        * m_all_nodes[i].sumforce(j);
      }
      else // 2nd order Adams-Bashforth from the 2nd time loop iteration
      {
        for (size_t j=0;j<dim;++j)
        {
          m_all_nodes[i].velocity(j) += ( dt /  node_mass ) 
                                        * 
                                        ( 1.5 * m_all_nodes[i].sumforce(j) 
                                        - 0.5 * m_all_nodes[i].sumforce_nm1(j) );
                                        
          m_all_nodes[i].sumforce_nm1(j) = m_all_nodes[i].sumforce(j);
        } 			       
      }
      
      // Update position with 2nd order Taylor series expansion
      for (size_t j=0;j<dim;++j)
      {
        m_all_nodes[i].coordinates(j) += dt 
                                         * m_all_nodes[i].velocity(j) 
                                         + 0.5 
                                         * ( m_all_nodes[i].sumforce(j) / node_mass ) 
                                         * pow( dt, 2. );
      }

      // Compute total kinetic energy
      for (size_t j=0;j<dim;++j)
        membrane_param.total_kinetic_energy += 0.5 * node_mass * 
                                     pow( m_all_nodes[i].velocity(j), 2.0 );

      // Convert force from model to physical units
      // // // cout << "Node = " << i << "\t";
      for (size_t j=0;j<dim;++j)
      {
        m_all_nodes[i].sumforce(j) = convert_model_to_physical_units(m_all_nodes[i].sumforce(j));
        // // // cout << m_all_nodes[i].sumforce(j) << "\t";
      }
      // // // cout << endl;
    }
    // // // cin >> pot;
  }
}




//---------------------------------------------------------------------------
void DS_3DRBC:: rbc_dynamics_solver_no_sub_time_stepping(size_t const& dim, 
                                                   size_t const& fluid_iter_num,
                                                   double const& dt_fluid, 
                                                   string const& case_type,
                                                   bool const& Matlab_numbering,
                                                   bool const& combined_force_computation, 
                                                   string const& model_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: rbc_dynamics_solver_no_sub_time_stepping" ) ;

  size_t num_nodes = shape_param.N_nodes;
  double node_mass = membrane_param.node_mass_M; // TO BE VERIFIED OF ITS VALUE
  double spring_constant = membrane_param.k_spring;
  double bending_spring_constant = membrane_param.k_bending;
  double bending_viscous_constant = membrane_param.k_bending_visc;
  double viscous_drag_constant = membrane_param.k_viscous;
  double area_force_constant = membrane_param.k_area;
  double volume_force_constant = membrane_param.k_volume;
  
  double dt = dt_fluid;
  
  bool compute_init_values; 
  compute_init_values = (int(fluid_iter_num) == 1) ? true : false;
  
  // Compute initial area & volume of membrane
  if(int(fluid_iter_num) == 1) 
    compute_total_surface_area_total_volume(compute_init_values, dim);
  
  // Initialize forces on all nodes
  for (size_t inode=0;inode<num_nodes;++inode)
    for (size_t j=0;j<dim;++j)
      m_all_nodes[inode].sumforce(j) = 0.0;
      
  // Compute normals, areas and center of mass of triangles
  compute_triangle_area_normals_centre_of_mass(compute_init_values, dim);
  
  // Compute total surface area and total volume and write to file
  compute_total_surface_area_total_volume(compute_init_values, dim);
  
  // Compute forces on all nodes
  if(combined_force_computation)
  {
    // Spring + Bending + Viscous forces
    compute_spring_bending_viscous_forces( dim, spring_constant, 
                                           bending_spring_constant,
                                           bending_viscous_constant,
                                           viscous_drag_constant, dt, 
                                           model_type );
                                           
    // Volume + Area force computation
    compute_volume_area_conservation_forces( dim, Matlab_numbering, model_type );
  }
  else
  {
    // Spring force
    compute_spring_force( dim, spring_constant, model_type );

    // Bending resistance
    compute_bending_resistance( dim, bending_spring_constant, 
                                bending_viscous_constant, dt, model_type );

    // Viscous drag force
    compute_viscous_drag_force( dim, viscous_drag_constant, dt, model_type );

    // Volume conservation force
    compute_volume_conservation_force( dim, model_type );

    // Triangle surface area conservation force
    compute_area_conservation_force( dim, Matlab_numbering, model_type );
  }

  for (size_t inode=0;inode<num_nodes;++inode)
    for (size_t j=0;j<dim;++j)
      m_all_nodes[inode].sumforce(j) = convert_model_to_physical_units(m_all_nodes[inode].sumforce(j));
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

  size_t num_nodes = shape_param.N_nodes;
  
  double axial_radius = 0.;
  
  double x_centroid = 0.; // membrane_param.centroid_coordinates(0);
  double y_centroid = 0.; // membrane_param.centroid_coordinates(1);
  double z_centroid = 0.; // membrane_param.centroid_coordinates(2);
  
  for (size_t inode=0;inode<num_nodes;++inode)
  {
    double x = m_all_nodes[inode].coordinates(0);
    double y = m_all_nodes[inode].coordinates(1);
    double z = m_all_nodes[inode].coordinates(2);
    
    double dist_axial = MAC::sqrt( pow(x - x_centroid, 2.)
                                 + pow(y - y_centroid, 2.)
                                 + pow(z - z_centroid, 2.) );
    
    axial_radius = max(axial_radius, dist_axial);
  }

  double axial_diameter = 2. * axial_radius;
  
  return(axial_diameter);
}




//---------------------------------------------------------------------------
double DS_3DRBC:: compute_transverse_diameter()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_transverse_diameter" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
  
  double euclid_dist_transverse = -HUGE_VAL;
  double transverse_radius = 0.;
  double y_node, z_node;
  size_t node_index;
  
  double y_centroid = 0.; // membrane_param.centroid_coordinates(1);
  double z_centroid = 0.; // membrane_param.centroid_coordinates(2);
  
  for (size_t inode=0;inode<num_nodes;++inode)
  {
    double y = m_all_nodes[inode].coordinates(1);
    double z = m_all_nodes[inode].coordinates(2);
    
    double dist_transverse = sqrt( pow(y - y_centroid, 2.)
                                 + pow(z - z_centroid, 2.) );
  
    transverse_radius = max(transverse_radius, dist_transverse);
  }
  
  double transverse_diameter = 2. * transverse_radius;
  
  return(transverse_diameter);
}




//-----------------------------------------------------------
void DS_3DRBC:: compute_tdp_orientation_angle()
//-----------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_tdp_orientation_angle" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
  
  // Taylor Deformation Parameter (TDP)
  double ad = membrane_param.axial_diameter;
  double td = membrane_param.transverse_diameter;
  membrane_param.taylor_deformation_parameter = (ad - td) / (ad + td);


  // Orientation angle
  double x_centroid = 0.; // membrane_param.centroid_coordinates(0);
  double y_centroid = 0.; // membrane_param.centroid_coordinates(1);
  double z_centroid = 0.; // membrane_param.centroid_coordinates(2);
  
  double axial_radius = 0.;
  double pi = MAC::pi();
  double radians_to_angle_conversion = 180. / pi;
  double theta;
  
  for (size_t inode=0;inode<num_nodes;++inode)
  {
    double x = m_all_nodes[inode].coordinates(0);
    double y = m_all_nodes[inode].coordinates(1);
    double z = m_all_nodes[inode].coordinates(2);
    
    double dist_axial = MAC::sqrt( pow(x - x_centroid, 2.)
                                 + pow(y - y_centroid, 2.)
                                 + pow(z - z_centroid, 2.) );
    
    if(dist_axial > axial_radius)
    {
      // Orientation angle with respect to x-axis
      double roll = x / MAC::sqrt(pow(x, 2.) + pow(y, 2.) + pow(z, 2.));
      theta = acos(roll);
      theta = (theta > pi) ? theta - pi : theta; // to keep the angle between 0 and 90 degrees
      membrane_param.orientation_roll_angle = theta * radians_to_angle_conversion;
      
      // Orientation angle with respect to y-axis
      double pitch = y / MAC::sqrt(pow(x, 2.) + pow(y, 2.) + pow(z, 2.));
      theta = acos(pitch);
      theta = (theta > pi) ? theta - pi : theta; // to keep the angle between 0 and 90 degrees
      membrane_param.orientation_pitch_angle = theta * radians_to_angle_conversion;
      
      // Orientation angle with respect to z-axis
      double yaw = z / MAC::sqrt(pow(x, 2.) + pow(y, 2.) + pow(z, 2.));
      theta = acos(yaw);
      theta = (theta > pi) ? theta - pi : theta; // to keep the angle between 0 and 90 degrees
      membrane_param.orientation_yaw_angle = theta * radians_to_angle_conversion;
    }
  }
}




//---------------------------------------------------------------------------
double DS_3DRBC:: compute_avg_tangential_velocity()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_avg_tangential_velocity" ) ;
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_gyration_tensor(size_t const& cyclenum,
                                        double const& time,
                                        string const& directory, 
                                        string const& filename)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_gyration_tensor" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
  
  double x_centroid = membrane_param.centroid_coordinates(0);
  double y_centroid = membrane_param.centroid_coordinates(1);
  double z_centroid = membrane_param.centroid_coordinates(2);
  
  for (size_t j=0;j<6;++j)
    membrane_param.gyration_tensor[j] = 0.;
  
  for (size_t inode=0;inode<num_nodes;++inode)
  {
    double x = m_all_nodes[inode].coordinates(0);
    double y = m_all_nodes[inode].coordinates(1);
    double z = m_all_nodes[inode].coordinates(2);

    // Gxx
    membrane_param.gyration_tensor[0] += (x - x_centroid) * (x - x_centroid);
    
    // Gyy
    membrane_param.gyration_tensor[1] += (y - y_centroid) * (y - y_centroid);
    
    // Gzz
    membrane_param.gyration_tensor[2] += (z - z_centroid) * (z - z_centroid);
    
    // Gxy
    membrane_param.gyration_tensor[3] += (x - x_centroid) * (y - y_centroid);
    
    // Gxz
    membrane_param.gyration_tensor[4] += (x - x_centroid) * (z - z_centroid);
    
    // Gyz
    membrane_param.gyration_tensor[5] += (y - y_centroid) * (z - z_centroid);
  }
  
  for (size_t j=0;j<6;++j)
    membrane_param.gyration_tensor[j] /= num_nodes;
    
  // Write gyration tensor to a text file
  ofstream gyration_tensor_file;
  string file_to_write = directory + "/" + filename;
  if(int(cyclenum) == 1)
  {
    gyration_tensor_file.open( file_to_write, ios::out );
    gyration_tensor_file << "Cyclenum" << "\t"
                         << "Time" << "\t"
                         << "Gxx" << "\t"
                         << "Gyy" << "\t"
                         << "Gzz" << "\t"
                         << "Gxy" << "\t"
                         << "Gxz" << "\t"
                         << "Gyz"
                         << endl;
  }
  else
  {
    gyration_tensor_file.open( file_to_write, ios::app );
  }
  
  gyration_tensor_file << cyclenum << "\t"
                       << time << "\t"
                       << membrane_param.gyration_tensor[0] << "\t"
                       << membrane_param.gyration_tensor[1] << "\t"
                       << membrane_param.gyration_tensor[2] << "\t"
                       << membrane_param.gyration_tensor[3] << "\t"
                       << membrane_param.gyration_tensor[4] << "\t"
                       << membrane_param.gyration_tensor[5]
                       << endl;
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_stats(string const& directory, string const& filename, 
                              size_t const& dim, double const& time, 
                              double const& final_time, size_t const& cyclenum,
                              string const& case_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_stats" ) ;

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
  m_gyration_tensor_rootname = "gyration_tensor_init_and_final.txt";


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
                    << "Orientation_roll_angle" << "\t"
                    << "Orientation_pitch_angle" << "\t"
                    << "Orientation_yaw_angle" << "\t"
                    // << "Average_tangential_velocity" << "\t"
                    << "Initial_area" << "\t"
                    << "Final_area" << "\t"
                    << "Initial_volume" << "\t"
                    << "Final_volume" << "\t"
                    << "Centroid_x" << "\t"
                    << "Centroid_y" << "\t"
                    << "Centroid_z" << "\t"
                    // << "Total_kinetic_energy"
                    << endl;
                    
  }
  else
  {
    rbc_stats_file.open( file_to_write, ios::app );
  }
  
  // Membrane area, volume & centroid
  // // compute_membrane_area_centroid_volume(dim);

  // Axial diameter
  membrane_param.axial_diameter = compute_axial_diameter();
  
  // Transverse diameter
  membrane_param.transverse_diameter = compute_transverse_diameter();
  
  // Taylor deformation parameter and orientation angle computation
  compute_tdp_orientation_angle();
  
  // Compute smallest eigenvalue of Gyration tensor
  if( (case_type.compare("Parabolic_flow") == 0) )
      // // and
      // // ((cyclenum == 1) or (abs(time - final_time) <= 1.e-8)) )
    compute_gyration_tensor(cyclenum, time, directory, m_gyration_tensor_rootname);
  
  /*
  // Average tangential velocity for tank treading = avg magnitude of velocity over all nodes
  membrane_param.avg_tangential_velocity = compute_avg_tangential_velocity();
  */
  
  // Writing to file
  rbc_stats_file << std::scientific // << setprecision(12)
                 << cyclenum << "\t"
                 << time << "\t"
                 << membrane_param.axial_diameter << "\t"
                 << membrane_param.transverse_diameter << "\t"
                 << membrane_param.taylor_deformation_parameter << "\t"
                 << membrane_param.orientation_roll_angle << "\t"
                 << membrane_param.orientation_pitch_angle << "\t"
                 << membrane_param.orientation_yaw_angle << "\t"
                 // << membrane_param.avg_tangential_velocity << "\t"
                 << membrane_param.initial_area << "\t"
                 << membrane_param.total_area << "\t"
                 << membrane_param.initial_volume << "\t"
                 << membrane_param.total_volume << "\t"
                 << membrane_param.centroid_coordinates(0) << "\t"
                 << membrane_param.centroid_coordinates(1) << "\t"
                 << membrane_param.centroid_coordinates(2) << "\t"
                 // << membrane_param.total_kinetic_energy
                 << endl;

  // Closing the file
  rbc_stats_file.close();
}




//---------------------------------------------------------------------------
double DS_3DRBC:: norm( geomVector& v )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: 3D_norm" ) ;
    
  return ( pow( v(0)*v(0) + v(1)*v(1) + v(2)*v(2), 0.5 ) );
}



    
//---------------------------------------------------------------------------
double DS_3DRBC:: scalar( geomVector const& v0, geomVector const& v1 )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: scalar" ) ;
    
  return ( v0(0) * v1(0) + v0(1) * v1(1) + v0(2) * v1(2) ); 
}



    
//---------------------------------------------------------------------------
void DS_3DRBC:: cross_3D( geomVector const& v0, 
                          geomVector const& v1, 
                          geomVector& res )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: cross" ) ;
    
  res(0) = v0(1) * v1(2) - v0(2) * v1(1);
  res(1) = v0(2) * v1(0) - v0(0) * v1(2);
  res(2) = v0(0) * v1(1) - v0(1) * v1(0);
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
