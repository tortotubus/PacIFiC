#include <DS_3DRBC.hh>
#include <FV_DiscreteField.hh>
#include <math.h>
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
  temp.dangle_dt = 0.;
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
void DS_3DRBC:: preprocess_membrane_parameters(string const& case_type
                                           , size_t const& num_subtimesteps_RBC)
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
}



    
//---------------------------------------------------------------------------
void DS_3DRBC:: apply_periodic_boundary_conditions()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: apply_periodic_boundary_conditions" ) ;
    
}




//---------------------------------------------------------------------------
void DS_3DRBC:: eul_to_lag()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_3DRBC:: eul_to_lag_3D()" ) ;

}





//---------------------------------------------------------------------------
double DS_3DRBC::norm( double const* v )
//---------------------------------------------------------------------------
{
    MAC_LABEL( "DS_3DRBC:: 3D_norm" ) ;
    
    return ( pow( v[0]*v[0] + v[1]*v[1] + v[2]*v[2], 0.5 ) );
}



    
//---------------------------------------------------------------------------
double DS_3DRBC::scalar( double const* v0, double const* v1 )
//---------------------------------------------------------------------------
{
    MAC_LABEL( "DS_3DRBC:: scalar" ) ;
    
    return ( v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2] ); 
}



    
//---------------------------------------------------------------------------
void DS_3DRBC::cross_3D( double const* v0, double const* v1, double* res )
//---------------------------------------------------------------------------
{
    MAC_LABEL( "DS_3DRBC:: cross" ) ;
    
    res[0] = v0[1] * v1[2] - v0[2] * v1[1];
    res[1] = v0[2] * v1[0] - v0[0] * v1[2];
    res[2] = v0[0] * v1[1] - v0[1] * v1[0];
} 

