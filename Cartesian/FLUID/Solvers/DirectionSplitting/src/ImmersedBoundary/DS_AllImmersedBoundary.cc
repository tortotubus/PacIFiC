#include <DS_AllImmersedBoundary.hh>
#include <DS_ImmersedBoundary.hh>
#include <DS_ImmersedBoundary_BuilderFactory.hh>
#include <FV_TimeIterator.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <cmath>
using std::endl;
using std::cout;
using std::cin;
using std::string;


//---------------------------------------------------------------------------
DS_AllImmersedBoundary:: DS_AllImmersedBoundary(size_t const& space_dimension
                                              , string const& IB_file
                                              , size_t const& N_IB)
//---------------------------------------------------------------------------
: m_space_dimension ( space_dimension )
, m_IB_file ( IB_file )
, m_nIB ( N_IB )
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

  read_shape_parameters();

  initialize_variables();

  generate_immersed_body_mesh();
  
  write_immersed_body_mesh_to_vtk_file();

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
void DS_AllImmersedBoundary:: read_shape_parameters()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: read_shape_parameters" ) ;

  double xp, yp, zp;
  double Rp;
  double x_roll_angle, y_pitch_angle, z_yaw_angle;
  double c0, c1, c2;
  size_t N_nodes, N_levels;

  std::ifstream inFile;
  std::ostringstream os2;
  os2 << "./InputFiles/" << m_IB_file;
  std::string filename = os2.str();

  inFile.open(filename.c_str());
  string line;
  getline(inFile,line); // read header line of input data file
  for (size_t i = 0; i < m_nIB; ++i) {
     inFile >> xp >> yp >> zp >> x_roll_angle >> y_pitch_angle
            >> z_yaw_angle >> Rp >> c0 >> c1 >> c2 >> N_nodes >> N_levels;
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
     p_shape_param->N_nodes = N_nodes;
     p_shape_param->N_levels = N_levels;

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
  MAC_LABEL( "DS_AllImmersedBoundary:: initialize_variables" ) ;

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundary[i]->set_all_nodes();
    m_allDSimmersedboundary[i]->project_membrane_shape();
    m_allDSimmersedboundary[i]->position_membrane();
    m_allDSimmersedboundary[i]->rotate_membrane();
    m_allDSimmersedboundary[i]->set_all_trielements();
    m_allDSimmersedboundary[i]->set_all_edges();
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
void DS_AllImmersedBoundary:: write_immersed_body_mesh_to_vtk_file()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: write_immersed_body_mesh_to_vtk_file" ) ;

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundary[i]->write_mesh_to_vtk_file(i, 0., 0);
  }

}




//---------------------------------------------------------------------------
void DS_AllImmersedBoundary:: eul_to_lag_velocity_interpolate()
//---------------------------------------------------------------------------
{
  MAC_LABEL( "DS_AllImmersedBoundary:: eul_to_lag_velocity_interpolate" ) ;

  for (size_t i = 0; i < m_nIB; ++i) {
    m_allDSimmersedboundary[i]->eul_to_lag();
  }

}





































