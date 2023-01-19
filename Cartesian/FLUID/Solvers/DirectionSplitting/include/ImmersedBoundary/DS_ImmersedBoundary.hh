#ifndef _DS_IMMERSEDBOUNDARY__
#define _DS_IMMERSEDBOUNDARY__

#include <geomVector.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <doubleArray2D.hh>
#include <boolVector.hh>
#include <doubleVector.hh>
#include <FV_Mesh.hh>
#include <FV_TimeIterator.hh>
#include <MAC_assertions.hh>
#include <MAC_Communicator.hh>
#include <vector>
#include <iostream>
#include <map>
#include <fstream>
#include <sstream>
using std::ostream;
using std::vector;
using std::tuple;
using std::string;
using namespace std;


class FV_DiscreteField;
class FV_Mesh;
struct Edge;

struct ShapeParameters 
{
  geomVector center;
  double radius;
  double xroll, ypitch, zyaw;
  double c0, c1, c2;
  size_t N_nodes;
  size_t N_levels;
  double node_spacing_with_dx; // deciding N_nodes based on this value
  double scaling_factor; // for scaling the coordinates & radius of immersed body
  double order_of_magnitude_of_radius;
};


struct MembraneParameters
{
  // Physical parameters
  double mass, node_mass;
  double k_spring, k_bending, k_bending_visc, k_viscous, k_area, k_area_local, k_volume;
  double membrane_spring_constant;
  double membrane_mass_spring_timescale, edge_mass_spring_timescale;
  double node_bending_mass_spring_timescale;

  // Flow physics parameters
  double ReynoldsNumber, CapillaryNumber, ShearRate;
  
  // Detailed "numerical membrane model (NMM)" parameters in 'physical' units
  double mu0_P; // Shear modulus in N/m 
  double Y_P; // Young's modulus in N/m
  double x0; // Maximum allowable extension of spring length --> x0 = l/lmax
  double D0_P; // Diameter of immersed body in micro-metre
  double kc_P; // Bending rigidity of immersed body in Joules
  double kbending_P; // Bending constant of immersed body
  double eta_P; // membrane viscosity
  double eta_plasma_P; // Dynamic viscosity of fluid surrounding immersed body
  double eta_cytoplasm_P; // Dynamic viscosity of fluid inside immersed body
  double rho_plasma_P; // Density of fluid surrounding immersed body in in kg/m^3
  double rho_P; // Density of immersed body
  double m; // exponent in POW spring force model
  double alpha; // exponent in time scale calculation
  double t; // timescale of the immersed body in seconds
  
  // // NMM parameters in 'model' units
  double node_mass_M;
  double mu0_M;
  double D0_M;
  double kc_M;
  double kbending_M;
  double eta_M;
  double eta_plasma_M;
  double eta_cytoplasm_M;
  double rho_plasma_M;
  double rho_M;

  // Temporal parameters
  double tmax, dt;
  size_t ntimescales;
  size_t ntimesteps;
  size_t n_subtimesteps_RBC;
  
  // Morphology parameters
  double axial_diameter, transverse_diameter, taylor_deformation_parameter;
  double orientation_angle;
  double orientation_roll_angle, orientation_pitch_angle, orientation_yaw_angle;
  double avg_tangential_velocity;
  double initial_perimeter, final_perimeter;
  double initial_area, final_area, total_area;
  double initial_volume, final_volume, total_volume;
  geomVector centroid_coordinates;
  vector<double> gyration_tensor;
  double smallest_eigen_value; // smallest eigenvalue of gyration tensor
  double init_smallest_eigen_value; // used in computing shifted eigenvalue
  
  // Statistics variables
  double mean_WLC_spring_force_magnitude;
  double mean_POW_spring_force_magnitude;
  double mean_bending_force_magnitude;
  double mean_viscous_force_magnitude;
  double mean_volume_force_magnitude;
  double mean_area_force_magnitude;
  double total_kinetic_energy;
};


struct IBMParameters
{
  string dirac_type;
  size_t dim;
};


struct Node
{
  size_t number;
  geomVector coordinates;
  geomVector coordinates_pbc;
  geomVector sumforce;
  geomVector sumforce_nm1;
  geomVector WLC_force;
  geomVector POW_force;
  geomVector spring_force;
  geomVector bending_force;
  geomVector area_force;
  geomVector volume_force;
  geomVector viscous_force;
  geomVector velocity;
  geomVector angular_velocity;
  double angle, initial_angle, angle_nm1, dangledt;
  geomVector unit_outwards_normal_vector;
  Node* neighbors[2]; // FOR 2D
  Edge const* edge_of_neighbors[2]; // FOR 2D
  vector<Node const*> neighbors_3D; // FOR 3D
  vector<double> initial_spring_length;
  
  double mass;
};


struct TriElement
{
  Node* pnodes[3]; // FOR 3D
  size_t number;
  vector<size_t> node_number; // contain the node numbers forming the edge
  vector< pair<Node const*,Node const*> > varea;
  geomVector twice_area_outwards_normal_vector;
  geomVector center_of_mass;

  double tri_area, tri_initial_area;
  double tri_volume, tri_initial_volume;
};


struct Edge
{
  size_t number;
  Node* n2;
  Node* n3;
  geomVector ext_unit_normal; // only needed for 2D membranes
  double initial_length;
  double length;
  
  double sintheta0, costheta0;
  pair<TriElement const*,Node*> t1v1;
  pair<TriElement const*,Node*> t2v4;

  double l0; // initial length of the spring
  double l; // instantaneous length of the spring
  double lmax; // maximum allowed length of the spring (often = 2.2 * l0 = l0/x0)
  double k; // = kBT/p for the WLC spring force for each edge
  double kp; // for the POW spring force each edge
  double initial_angle;
  double angle;
  double angle_nm1;
  double dangledt;
};


/** @brief The class DS_ImmersedBoundary.

A moving or stationary rigid body in the Direction Splitting solver.

@author A. Goyal - Pacific project 2022 */

class DS_ImmersedBoundary
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_ImmersedBoundary();

      /** @brief Destructor */
      virtual ~DS_ImmersedBoundary();
      //@}


   //-- Set methods

      /**@name Set methods */
      //@{

      //@}

   //-- Get Methods
      /**@name Get methods */
      //@{
      /** @brief Returns the shape parameters of IB */
      ShapeParameters* get_ptr_shape_parameters();
      
      /** @brief Returns the object of MembraneParameters structure */
      MembraneParameters* get_ptr_membrane_parameters();
      
      /** @brief Returns the object of IBMParameters structure */
      IBMParameters* get_ptr_IBM_parameters();

      //@}

   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Creating the surface points on a RBC **/
      void create_RBC_structure( );

      /** @brief Outputs the shape and membrane parameters to screen */
      void display_parameters();

      /** @brief Allocates memory & initializes node properties & attributes */
      virtual void initialize_node_properties(string const& mesh_filename,
                                              size_t const& dim) = 0;

      /** @brief Generates the mesh for 2D and 3D immersed bodies **/
      virtual void set_all_nodes(istream& fileIN, size_t const& dim) = 0;
      
      /** @brief Allocates memory & initializes all variables of 'Struct TriElements' */
      virtual void initialize_triangle_properties(size_t const& num_triangles,
                                                  size_t const& dim) = 0;
      
      /** @brief Sets triangular elements for 2D and 3D immersed bodies */
      virtual void set_all_trielements(istream& fileIN, size_t const& dim,
                                       bool const& MatlabNumb) = 0;
      
      /** @brief Sets the neighbors of each node from triangle nodes */
      virtual void set_all_node_neighbors() = 0;
      
      /** @brief Computes the area, normals and centre of mass of each triangle */
      virtual void compute_triangle_area_normals_centre_of_mass(bool init,
                                                        size_t const& dim) = 0;
      
      /** @brief Allocates memory & initializes edge properties & attributes */
      virtual void initialize_edge_properties(size_t const& dim) = 0;
      
      /** @brief Sets the edges and nodes for 2D and 3D immersed bodies */
      virtual void set_all_edges(size_t const& dim) = 0;
      
      /** @brief Projects shape of the membrane to sphere or biconcave **/
      virtual void project_membrane_shape() = 0;
      
      /** @brief Positions the membrane from centers of mass in RBC_data.csv */
      void position_membrane();
      
      /** @brief Rotates the membrane based on 3D rotation matrix 
      2D: rotation is achieved using ONLY "yaw" angle
      3D: rotation is achieved using roll, pitch and yaw angles */
      void rotate_membrane();
      
      /** @brief Computes current (& initial) spring lengths */
      virtual void compute_spring_lengths(bool init, size_t const& dim) = 0;
      
      /** @brief Computes the unit normal of each edge/spring */
      virtual void compute_edge_normals() = 0;
      
      /** @brief Computes angle between edges */
      virtual void compute_edge_angle(bool init) = 0;

      /** @brief Computes norm of a vector variable */
      virtual double norm(geomVector& v) = 0;
      
      /** @brief Scalar product of two vectors */
      virtual double scalar( geomVector const& v0, geomVector const& v1 ) = 0;
      
      /** @brief Checks if point p1 and point p2 has 
      connection across periodic boundary? **/
      bool across_periodic(double p1, double p2, double length);

      /** @brief Computes the periodic distance between two points */
      double periodic_1D_distance(double p1, double p2, double length);

      /** Computes distance between two points based 
      on periodic boundary condition **/
      double compute_dist_incl_pbc(double p1, double p2, double length);
      
      /** @brief Computes node based spring and bending constants **/
      virtual void preprocess_membrane_parameters(string const& model_type,
                                        string const& case_type,
                                        double const& mu,
                                        size_t const& num_subtimesteps_RBC,
                                        size_t const& dim) = 0;
      
      /** @brief Function which calls RBC and IBM functions along with
      periodic boundary conditions and parallelisation temporary variables */
      void do_one_inner_iteration( FV_DiscreteField const* UF,
                                   FV_DiscreteField* Eul_F,
                                   FV_DiscreteField* F_Eul_tag,
                                   FV_TimeIterator const* t_it,
                                   FV_Mesh const* MESH,
                                   size_t const& dim,
                                   boolVector const* is_periodic,
                                   string const& case_type,
                                   bool const& Matlab_numbering,
                                   string const& model_type = "NumericalMembraneModel" );
        
      /** @brief Discretised Dirac delta function
      @param val -> the value which is to be converted using Dirac delta
      @param Dirac_type -> Balogh, Roma, Archer **/
      double discrete_Dirac_delta( double val, string const& Dirac_type, 
                                   double dh, int n ) ;

      /** @brief Applies periodic boundary condition to each immersed body */
      void apply_periodic_boundary_conditions(FV_Mesh const* MESH,
                                              size_t const& dim,
                                              boolVector const* is_periodic);
      
      /** @brief IBM: Eulerian velocity to Lagrangian velocity interpolation **/
      virtual void eul_to_lag(FV_DiscreteField const* FF, size_t const& dim, 
                              size_t const& periodic_dir) = 0;

      /** @brief Copies the Lagrangain velocity to a doubleVector */
      void copy_lagrangian_velocity_to_vector(doubleVector& lag_vel, 
                                              size_t const& dim);
      
      /** @brief Copies the Lagrangain position & force to a doubleVector */
      void copy_vector_to_lagrangian_velocity (doubleVector& lag_vel, 
                                               size_t const& dim);
                            
      /** @brief Computes RBC deformation using spring-dashpot model */
      virtual void rbc_dynamics_solver(size_t const& dim, 
                       double const& dt_fluid,
                       string const& case_type,
                       bool const& Matlab_numbering,
                       string const& model_type = "NumericalMembraneModel" ) = 0;
      
      /** @brief Computes RBC deformation using spring-dashpot model */
      virtual void rbc_dynamics_solver_no_sub_time_stepping(size_t const& dim, 
                      double const& dt_fluid,
                      string const& case_type,
                      bool const& Matlab_numbering,
                      string const& model_type = "NumericalMembraneModel" ) = 0;
      
      /** @brief Computes spring force */
      virtual void compute_spring_force( size_t const& dim,
                      double const& spring_constant,
                      string const& model_type = "NumericalMembraneModel" ) = 0;
      
      /** @brief Computes linear spring force Breyannis2000 and Bagchi2007 
      1. Compute tension T = Es * (l/l0 -1) for each edge of a node
      2. Compute spring force F = T_i * e_i - T_j * e_j 
      where i and j are neighbouring edge indices of a node */
      virtual void compute_linear_spring_force( size_t const& dim,
                                            double const& spring_constant ) = 0;
      
      /** @brief Computes bending force */
      virtual void compute_bending_resistance (size_t const& dim, 
                      double const& bending_spring_constant,
                      double const& bending_viscous_constant, 
                      double const& dt,
                      string const& model_type = "NumericalMembraneModel" ) = 0;

      /** @brief Computes viscous force */
      virtual void compute_viscous_drag_force(size_t const& dim, 
                      double const& viscous_drag_constant,
                      double const& dt,
                      string const& model_type = "NumericalMembraneModel" ) = 0;

      /** @brief Computes volume conservation force */
      virtual void compute_volume_conservation_force(size_t const& dim,
                      string const& model_type = "NumericalMembraneModel" ) = 0;
      
      /** @brief Computes area conservation force */
      virtual void compute_area_conservation_force(size_t const& dim,
                      bool const& Matlab_numbering,
                      string const& model_type = "NumericalMembraneModel" ) = 0;
      
      /** @brief Computes the statistics of each immersed boundary */
      virtual void compute_stats(string const& directory, 
                                 string const& filename, 
                                 size_t const& dim,
                                 double const& time, 
                                 double const& final_time,
                                 size_t const& cyclenum) = 0;
      
      /** @brief Computes the centroid of each immersed boundary */
      virtual void compute_centroid(size_t const& dim) = 0;
      
      /** @brief Copies the Lagrangain position & force to a doubleVector */
      void copy_lag_position_and_force_to_vector
                            (doubleVector& lag_pos_and_force, size_t const& dim);
                            
      /** @brief Copies vector to Lagrangian position & force */
      void copy_vector_to_lag_position_and_force
                           (doubleVector& lag_pos_and_force, size_t const& dim);

      /** @brief Lagrangian to Eulerian force spreading **/
      virtual void lag_to_eul(FV_DiscreteField* FF, 
                              FV_DiscreteField* FF_tag,
                              size_t const& dim, 
                              size_t const& comp) = 0;
      
      /** @brief Computes axial diameter of each immersed boundary */
      virtual double compute_axial_diameter() = 0;
      
      /** @brief Computes transverse diameter of each immersed boundary */
      virtual double compute_transverse_diameter() = 0;
      
      /** @brief Computes Taylor Deformation Parameter & Orientation angle for
      each immersed boundary */
      virtual void compute_tdp_orientation_angle() = 0;
      
      /** @brief Computes the average tangential velocity of immersed boundary*/
      virtual double compute_avg_tangential_velocity() = 0;
      
      /** @brief Computes the perimeter = sum of all edge lengths */
      virtual double perimeter() = 0;
      
      /** @brief Writes one point of RBC mesh to .vtu file **/
      virtual void write_mesh_to_vtk_file( size_t IB_number, double const& time,
                                           size_t const& cyclenum ) = 0;
                                           
      /** @brief Writes the rbc.pvd file for each RBC membrane */
      void write_rbc_dot_pvd_file(size_t IB_number);
      
      /** @brief Converts size_t to string datatype **/
      string sizetToString( size_t const& figure ) const;
      
      /** @brief Outputs the datatype of a variable **/
      void get_datatype_of_variable();

      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      // Shape parameters
      ShapeParameters shape_param;
      
      // Membrane physical property parameters
      MembraneParameters membrane_param;
      
      // IBM parameters
      IBMParameters ibm_param;

      // Vector containing all nodes properties
      vector<Node> m_all_nodes;
      vector<Edge> m_all_edges;
      vector<TriElement> m_all_trielements;

      size_t m_nEdges;
      size_t m_nTriangles;
      
      // File handlers for writing RBC coordinates files
      ostringstream m_rbc_pvd;
      ostringstream m_vtk_to_pvd;
      ostringstream m_rbc_one_point_pvd;
      ostringstream m_one_point_vtk_to_pvd;

      // Filename variables for statistics computation & writing
      string m_rootdir;
      string m_rootname;
      string m_video_rootname;
      string m_kinetic_energy_rootname;
      string m_morphology_rootname;
      string m_force_stats_rootname;
      string m_diameter_stats_rootname;
      string m_rbc_one_point_rootname;
      string m_triangle_unit_normals_rootname;
      string m_node_unit_normals_rootname;
      string m_gyration_tensor_rootname;
      //@}


   //-- Methods

      /**@name Methods */
      //@{

      //@}


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Copy constructor
      @param copy copied DS_ImmersedBoundary object */
      DS_ImmersedBoundary( DS_ImmersedBoundary const& copy );
      
      
      //@}
};

#endif
