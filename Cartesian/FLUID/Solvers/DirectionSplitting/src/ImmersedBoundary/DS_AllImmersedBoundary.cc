#include <DS_AllImmersedBoundary.hh>
#include <DS_ImmersedBoundary.hh>
#include <DS_ImmersedBoundary_BuilderFactory.hh>
#include <FV_Mesh.hh>
#include <FV_TimeIterator.hh>
#include <FV_DiscreteField.hh>
#include <cmath>
using std::endl;
using std::cout;
using std::cin;
using std::string;


//-----------------------------------------------------------------------------
DS_AllImmersedBoundary:: DS_AllImmersedBoundary(size_t const& space_dimension
                                               , string const& IB_file
                                               , size_t const& N_IB
                                               , string const& case_type
                                               , FV_DiscreteField const* arb_UF
                                               , FV_DiscreteField* arb_EulF
                                               , FV_DiscreteField* arb_EulF_tag
                                               , size_t const& nRBC_subtimesteps
                                               , string const& dirac_type
                                               , size_t const& periodic_dir)
//-----------------------------------------------------------------------------
: m_space_dimension ( space_dimension )
, m_IB_file ( IB_file )
, m_nIB ( N_IB )
, UF ( arb_UF )
, Eul_F ( arb_EulF )
, F_Eul_tag ( arb_EulF_tag )
, MESH ( UF->primary_grid() )
, m_IB_case_type ( case_type )
, m_subtimesteps_RBC ( nRBC_subtimesteps )
, m_dirac_type ( dirac_type )
, m_periodic_dir ( periodic_dir )
{
  MAC_LABEL( "DS_AllImmersedBoundary:: DS_AllImmersedBoundary" ) ;

  // // m_nIB = get_num_lines_in_IB_file();

  m_allDSimmersedboundary.reserve( m_nIB );

  DS_ImmersedBoundary* dsib = NULL;

  for (size_t i = 0; i < m_nIB; ++i)
  {
     m_allDSimmersedboundary.push_back( dsib );
     m_allDSimmersedboundary[i] = DS_ImmersedBoundary_BuilderFactory::
                                                   create(m_space_dimension);
  }

  read_shape_and_membrane_parameters();

  generate_immersed_body_mesh();
  
  // write_immersed_body_mesh_to_vtk_file();
  
  preprocess_immersed_body_parameters(m_IB_case_type, m_subtimesteps_RBC);
  
  set_IBM_parameters(m_dirac_type, m_periodic_dir);
}




//---------------------------------------------------------------------------
DS_AllImmersedBoundary:: ~DS_AllImmersedBoundary()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: ~DS_AllImmersedBoundary" ) ;

  for (size_t i = 0; i < m_nIB; ++i) delete m_allDSimmersedboundary[i];
  m_allDSimmersedboundary.clear();

}




//---------------------------------------------------------------------------
size_t DS_AllImmersedBoundary:: get_number_of_immersed_boundaries() const
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: get_number_of_immersed_boundaries" ) ;

  return ( m_nIB );

}




//---------------------------------------------------------------------------
DS_ImmersedBoundary* DS_AllImmersedBoundary:: get_ptr_immersed_body( size_t i )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: get_ptr_immersed_body" ) ;

  return ( m_allDSimmersedboundary[i] );

}




//---------------------------------------------------------------------------
size_t DS_AllImmersedBoundary:: get_num_lines_in_IB_file()
//---------------------------------------------------------------------------
{
  std::ifstream inFile;
  std::ostringstream os2;
  os2 << "./InputFiles/" << m_IB_file;
  std::string filename = os2.str();

  inFile.open(filename.c_str());
  string line;
  size_t number_of_lines = 0;

  if(inFile.is_open())
  {
    while(inFile.peek()!=EOF)
    {
      getline(inFile, line);
      number_of_lines++;
    }
    inFile.close();
    cout<<"Number of lines in the file are: " << number_of_lines << endl;
  }
  else
    cout << "Couldn't open the file" << endl;

  return ( number_of_lines - 1 );
}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundary:: read_shape_and_membrane_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: read_shape_and_membrane_parameters" ) ;

  double xp, yp, zp;
  double Rp;
  double x_roll_angle, y_pitch_angle, z_yaw_angle;
  double c0, c1, c2;
  size_t N_nodes, N_levels;
  size_t node_spacing_with_dx;
  double k_spring, k_bending, k_bending_visc, k_viscous, k_area, k_volume;
  double membrane_mass;
  double dx = 1.; // CHANGE THIS

  std::ifstream inFile;
  std::ostringstream os2;
  os2 << "./InputFiles/" << m_IB_file;
  std::string filename = os2.str();

  inFile.open(filename.c_str());
  string line;
  getline(inFile,line); // read header line of input data file
  for (size_t i = 0; i < m_nIB; ++i) {
     inFile >> xp >> yp >> zp >> x_roll_angle >> y_pitch_angle
            >> z_yaw_angle >> Rp >> c0 >> c1 >> c2 >> N_nodes >> N_levels
            >> node_spacing_with_dx >> k_spring >> k_bending >> k_bending_visc
            >> k_viscous >> k_area >> k_volume >> membrane_mass;
            
     ShapeParameters* p_shape_param =  m_allDSimmersedboundary[i]
                                                ->get_ptr_shape_parameters();
     
     geomVector position(xp,yp,zp);
     p_shape_param->center = position;
     p_shape_param->xroll = x_roll_angle;
     p_shape_param->ypitch = y_pitch_angle;
     p_shape_param->zyaw = z_yaw_angle;
     p_shape_param->radius = Rp;
     p_shape_param->c0 = c0;
     p_shape_param->c1 = c1;
     p_shape_param->c2 = c2;
     p_shape_param->N_levels = N_levels;
     p_shape_param->node_spacing_with_dx = node_spacing_with_dx;
     p_shape_param->N_nodes = (node_spacing_with_dx != 0) 
                              ? round ( 2.0 * MAC::pi() * Rp 
                                / (double(node_spacing_with_dx) * dx) )
                              : N_nodes;
     
     MembraneParameters* p_membrane_param = m_allDSimmersedboundary[i]
                                                ->get_ptr_membrane_parameters();
     
     p_membrane_param->k_spring = k_spring;
     p_membrane_param->k_bending = k_bending;
     p_membrane_param->k_bending_visc = k_bending_visc;
     p_membrane_param->k_viscous = k_viscous;
     p_membrane_param->k_area = k_area;
     p_membrane_param->k_volume = k_volume;
     p_membrane_param->mass = membrane_mass;
     
     // m_allDSimmersedboundary[i]->display_parameters();
  }
  inFile.close();
}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundary:: initialize_variables()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: initialize_variables" ) ;

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundary[i]->initialize_node_properties();
  }

}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundary:: generate_immersed_body_mesh()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: generate_immersed_body_mesh" ) ;

  for (size_t i = 0; i < m_nIB; ++i) {
    // Node properties
    m_allDSimmersedboundary[i]->initialize_node_properties();
    m_allDSimmersedboundary[i]->set_all_nodes();
    m_allDSimmersedboundary[i]->project_membrane_shape();
    m_allDSimmersedboundary[i]->position_membrane();
    m_allDSimmersedboundary[i]->rotate_membrane();
    
    // Triangle properties
    m_allDSimmersedboundary[i]->set_all_trielements();
    
    // Edge properties
    m_allDSimmersedboundary[i]->initialize_edge_properties();
    m_allDSimmersedboundary[i]->set_all_edges();
    m_allDSimmersedboundary[i]->compute_spring_lengths(true);
    m_allDSimmersedboundary[i]->compute_edge_normals();
    m_allDSimmersedboundary[i]->compute_edge_angle(true);
  }
}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundary:: project_shape_of_immersed_body()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: project_shape_of_immersed_body" ) ;
  
  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundary[i]->project_membrane_shape();
  }
}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundary:: position_immersed_body()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: position_immersed_body" ) ;
  
  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundary[i]->position_membrane();
  }
}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundary:: rotate_immersed_body()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: rotate_immersed_body" ) ;

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundary[i]->rotate_membrane();
  }
}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundary:: set_IBM_parameters(string const& dirac_type
                                               , size_t const& periodic_dir)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: set_IBM_parameters" ) ;

  for (size_t i = 0; i < m_nIB; ++i) {
     IBMParameters* p_ibm_param = m_allDSimmersedboundary[i]
                                                ->get_ptr_IBM_parameters();
     p_ibm_param->dirac_type = dirac_type;
     p_ibm_param->periodic_dir = periodic_dir;
  }
}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundary:: write_immersed_body_mesh_to_vtk_file()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: write_immersed_body_mesh_to_vtk_file" ) ;

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundary[i]->write_mesh_to_vtk_file(i, 0., 0);
  }

}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundary:: preprocess_immersed_body_parameters
                   (string const& case_type, size_t const& num_subtimesteps_RBC)
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: preprocess_immersed_body_parameters" ) ;
  
  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundary[i]->preprocess_membrane_parameters(case_type
                                                        , num_subtimesteps_RBC);
  }
}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundary:: do_one_inner_iteration
                             ( FV_TimeIterator const* t_it )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: do_one_inner_iteration" ) ;

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundary[i]->do_one_inner_iteration(UF
                                                     , Eul_F
                                                     , F_Eul_tag
                                                     , t_it
                                                     , MESH
                                                     , m_space_dimension
                                                     , m_periodic_dir
                                                     , m_IB_case_type);
  }
}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundary:: do_additional_savings
                              ( FV_TimeIterator const* t_it,
                                size_t const& cycleNumber )
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: do_additional_savings" ) ;

  for (size_t i = 0; i < m_nIB; ++i) {
    // Writing Res/rbc_T*.vtu files for every iterations
    m_allDSimmersedboundary[i]->write_mesh_to_vtk_file(i, t_it->time(), 
                                                       cycleNumber);
    // Writing Res/rbc.pvd file
    m_allDSimmersedboundary[i]->write_rbc_dot_pvd_file();
    

    /*
    // Writing Res/rbc_one_point_T*.vtu files for every sampled iterations
    m_allDSimmersedboundary[i]->write_one_point_to_VTK(t_it->time(), 
                                                       cycleNumber);
    // Writing one point to another rbc_one_point.pvd file for tank treading
    m_allDSimmersedboundary[i]->write_one_point_of_rbc_mesh_to_pvd_file(
                                                     t_it->time(), cycleNumber);
    

    // Writing statistics of RBC membrane
    m_allDSimmersedboundary[i]->compute_stats(t_it->time(), cycleNumber);


    // Writing triangle unit normals to .vtu file
    m_allDSimmersedboundary[i]->write_triangle_normals_to_VTK
                                                           ("triangle_normals");
    // writing node unit normals to .vtu file
    m_allDSimmersedboundary[i]->write_node_normals_to_VTK("node_normals");
    */
  }
}










































