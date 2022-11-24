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
#include <set>
#include <fstream>
#include <sstream>
using std::endl;
using std::cout;
using std::cin;
using std::string;
using std::max;
using namespace std;

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
  temp.spring_force(dim);
  temp.bending_force(dim);
  temp.viscous_force(dim);
  temp.volume_force(dim);
  temp.area_force(dim);
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
  cout << "Number of triangles = " << m_nTriangles << endl;
  
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

    if(isinf(m_all_trielements[i].twice_area_outwards_normal_vector(0)) or isinf(m_all_trielements[i].twice_area_outwards_normal_vector(1)) or isinf(m_all_trielements[i].twice_area_outwards_normal_vector(2)))
    {
      cout << "normals are infinity\n";
      exit(3);
    }

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
  
  // Check if n2-n3 are correctly ordered else swap them
  // i.e., check if normals of triangles on either side 
  // of edge match original normals from set_all_trielements()
  Node* buffer = NULL;
  for (vector<Edge>::iterator il=m_all_edges.begin();il!=m_all_edges.end();il++)
  {
    geomVector normal_t1(dim);
    for (size_t j=0;j<dim;++j) 
      normal_t1(j) = il->t1v1.first->twice_area_outwards_normal_vector(j);
    
    geomVector normal_t2(dim);
    for (size_t j=0;j<dim;++j) 
      normal_t2(j) = il->t2v4.first->twice_area_outwards_normal_vector(j);
    
    // a21 = a1 - a2 = vector from node 2 to node 1
    geomVector a21(dim);
    for (size_t j=0;j<dim;++j) 
      a21(j) = il->t1v1.second->coordinates(j) - il->n2->coordinates(j); // /order_of_magnitude_of_radius;
    
    // a31 = a1 - a3 = vector from node 3 to node 1
    geomVector a31(dim);
    for (size_t j=0;j<dim;++j) a31(j) = il->t1v1.second->coordinates(j) - il->n3->coordinates(j); // /order_of_magnitude_of_radius;
    
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
      for (size_t j=0;j<dim;++j) a21(j) = il->t1v1.second->coordinates(j) - il->n2->coordinates(j); // /order_of_magnitude_of_radius;
      
      // a31 = a1 - a3 = vector from node 3 to node 1
      for (size_t j=0;j<dim;++j) a31(j) = il->t1v1.second->coordinates(j) - il->n3->coordinates(j); // /order_of_magnitude_of_radius;
      
      // Chi = a21 x a31
      cross_3D(a21, a31, Chi);
    }
    
    
    
    // a34 = a4 - a3 = vector from node 3 to node 4
    geomVector a34(dim);
    for (size_t j=0;j<dim;++j) a34(j) = il->t2v4.second->coordinates(j) - il->n3->coordinates(j); // /order_of_magnitude_of_radius;
    
    // a24 = a4 - a2 = vector from node 2 to node 4
    geomVector a24(dim);
    for (size_t j=0;j<dim;++j) a24(j) = il->t2v4.second->coordinates(j) - il->n2->coordinates(j); // /order_of_magnitude_of_radius;
    
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
    for (size_t j=0;j<dim;++j) tc1_minus_tc2(j) = tc1(j) - tc2(j); // /order_of_magnitude_of_radius;
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
  
  membrane_param.mu0_P = 6.3e-6; // Shear modulus in N/m
  membrane_param.Y_P = 18.9e-6; // Young's modulus in N/m
  membrane_param.x0 = 1./2.2; // Maximum allowable spring extension --> x0=l/lmax
  membrane_param.D0_P = 7.82e-6; // Diameter of RBC micro-metre
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
  MAC_LABEL( "DS_3DRBC:: preprocess_membrane_parameters" ) ;
  
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
void DS_3DRBC:: compute_spring_constant_values(size_t const& dim)
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
void DS_3DRBC:: preprocess_membrane_parameters(string const& model_type,
                                             string const& case_type,
                                             double const& mu,
                                             size_t const& num_subtimesteps_RBC,
                                             size_t const& dim)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: preprocess_membrane_parameters" ) ;
  
  if(model_type.compare("NumericalMembraneModel") == 0) // if it is detailed numerical membrane model
  {
    // Initialize membrane material properties in physical units
    init_membrane_parameters_in_physical_units();
    
    // Initialize membrane material properties in model units
    init_membrane_parameters_in_model_units();
    
    // Scaling membrane material properties from physical units to model units
    scaling_membrane_params_from_physical_to_model_units();
    
    // Compute WLC and POW spring constant values for each spring using mu0, x0, l0, lmax
    compute_spring_constant_values(dim);
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
  double r1, p1, q1, delt1; // Dirac delta variables
  size_t istart, iend, jstart, jend, kstart, kend;

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
            // // dist_x = (xC - xp) * hxC;
            double dist_y = 
                compute_dist_incl_pbc(yC, yp, domain_length(1)) * hyC;
            // // dist_y = (yC - yp) * hyC;
            double dist_z = 
                compute_dist_incl_pbc(zC, zp, domain_length(2)) * hzC;
            // // dist_z = (zC - zp) * hzC;
            bool eul_cell_within_Dirac_delta_stencil = 
                                                 (fabs(dist_x) <= 2.) 
                                                 and 
                                                 (fabs(dist_y) <= 2.) 
                                                 and 
                                                 (fabs(dist_z) <= 2.);
            
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
  double r1, p1, q1, delt1; // Dirac delta variables
  size_t istart, iend, jstart, jend, kstart, kend;

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
          // // dist_x = (xC - xp) * hxC;
          double dist_y = 
                  compute_dist_incl_pbc(yC, yp, domain_length(1)) * hyC;
          // // dist_y = (yC - yp) * hyC;
          double dist_z = 
                  compute_dist_incl_pbc(zC, zp, domain_length(2)) * hzC;
          // // dist_z = (zC - zp) * hzC;
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
            double euler_force = FF->DOF_value(ii, jj, kk, comp, 0) 
                                 + 
                                 m_all_nodes[inode].sumforce(comp) 
                                 * delt1 * dxC * dyC * dzC;
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
  if ( init ) membrane_param.initial_area = membrane_param.total_area;


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
void DS_3DRBC::compute_spring_force( size_t const& dim,
                                     double const& spring_constant,
                                     string const& force_type )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_spring_force" ) ;
  
  size_t num_nodes = shape_param.N_nodes;
  
  if(force_type.compare("Anthony") == 0)
  {
    geomVector Iij(dim);
    double length = 0.;
    
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
                                        string const& force_type )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_bending_resistance" ) ;
  
  geomVector a21(dim), a31(dim), a24(dim), a34(dim);
  geomVector Chi(dim), Gamma(dim), crossp(dim), n2n3(dim);
  double angle = 0., dot = 0., mi = 0.;
  double lever1 = 0., lever4 = 0., nxi = 0., nzeta = 0.;

  if(force_type.compare("Anthony") == 0)
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
}




//---------------------------------------------------------------------------
void DS_3DRBC:: compute_viscous_drag_force( size_t const& dim,
                                            double const& viscous_drag_constant,
                                            string const& force_type )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: compute_viscous_drag_force" ) ;
    
  if(force_type.compare("Anthony") == 0)
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
}




//---------------------------------------------------------------------------
void DS_3DRBC:: rbc_dynamics_solver(size_t const& dim, 
                                    double const& dt_fluid, 
                                    string const& case_type)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: rbc_dynamics_solver" ) ;

  size_t num_nodes = shape_param.N_nodes;
  double node_mass = membrane_param.node_mass;
  double spring_constant = membrane_param.k_spring;
  double bending_spring_constant = membrane_param.k_bending;
  double bending_viscous_constant = membrane_param.k_bending_visc;
  double viscous_drag_constant = membrane_param.k_viscous;
  double area_force_constant = membrane_param.k_area;
  double volume_force_constant = membrane_param.k_volume;
  size_t n_sub_timesteps = membrane_param.n_subtimesteps_RBC;
  
  double time = 0.;
  double total_kinetic_energy = 0.;
  double dt = dt_fluid / n_sub_timesteps;
  
  MAC_Communicator const* MAC_comm;
  MAC_comm = MAC_Exec::communicator();
  size_t my_rank = MAC_comm->rank();
  size_t nb_procs = MAC_comm->nb_ranks();
  size_t is_master = 0;
  
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
    
    // cout << "Proc id = " << my_rank << "\ttime = " << time << "\tMembrane centroid = " << membrane_param.centroid_coordinates(0) << "\t" << membrane_param.centroid_coordinates(1) << "\t" << membrane_param.centroid_coordinates(2) << endl;
    if(isnan(membrane_param.centroid_coordinates(0)) or isnan(membrane_param.centroid_coordinates(1)) or isnan(membrane_param.centroid_coordinates(2)))
    {
      cout << "Membrane centroid NaN\n";
      exit(3);
    }

    // Compute forces on all nodes
    // Spring force
    if(case_type.compare("Breyannis2000case") != 0)
      compute_spring_force( dim, spring_constant, "Anthony" );
    else
      compute_linear_spring_force( dim, spring_constant );
        

    // Bending resistance
    if(case_type.compare("Breyannis2000case") != 0)
    {
      compute_bending_resistance( dim, bending_spring_constant, 
                                  bending_viscous_constant, dt, "none" );
    }

    /*
    // Viscous drag force
    compute_viscous_drag_force( viscous_drag_constant );

    // Volume conservation force
    compute_volume_conservation_force( volume_conservation_constant );

    // Triangle surface area conservation force
    compute_triangle_surface_area_conservation_force( triangle_surface_area_conservation_constant );
    */

    // Compute new velocity and position
    for (size_t i=0;i<num_nodes;++i)
    {
      // Solve momentum conservation
      if ( !iter_num ) // First order explicit at the 1st time
      {
        for (size_t j=0;j<dim;++j)
          m_all_nodes[i].velocity(j) += ( dt /  node_mass ) 
                                        * m_all_nodes[i].sumforce(j);
      }
      else // 2nd order Adams-Bashforth from the 2nd time
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
      
      // cout << "time = " << time << "\tinode = " << i << "\tforce = " << m_all_nodes[i].sumforce(0) << "\t" << m_all_nodes[i].sumforce(1) << "\t" << m_all_nodes[i].sumforce(2) << endl;
      // cout << "time = " << time << "\tinode = " << i << "\tvelocity = " << m_all_nodes[i].velocity(0) << "\t" << m_all_nodes[i].velocity(1) << "\t" << m_all_nodes[i].velocity(2) << endl;

      // Update position with 2nd order Taylor series expansion
      for (size_t j=0;j<dim;++j)
      {
        m_all_nodes[i].coordinates(j) += dt 
                                         * m_all_nodes[i].velocity(j) 
                                         + 0.5 
                                         * ( m_all_nodes[i].sumforce(j) / node_mass ) 
                                         * pow( dt, 2. );
      }
    }
  }
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
void DS_3DRBC:: compute_stats(string const& directory, string const& filename, 
                              size_t const& dim, double const& time, 
                              size_t const& cyclenum)
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
