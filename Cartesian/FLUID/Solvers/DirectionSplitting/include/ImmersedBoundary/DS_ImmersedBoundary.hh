#ifndef _DS_IMMERSEDBOUNDARY__
#define _DS_IMMERSEDBOUNDARY__

#include <geomVector.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <doubleArray2D.hh>
#include <boolVector.hh>
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
  size_t node_spacing_with_dx; // deciding N_nodes based on this value
};


struct Node
{
  size_t number;
  geomVector coordinates;
  geomVector coordinates_pbc;
  geomVector sumforce;
  geomVector sumforce_nm1;
  geomVector velocity;
  geomVector angular_velocity;
  double angle, initial_angle, angle_nm1, dangle_dt;
  geomVector spring_force, bending_force;
  geomVector area_force, volume_force, viscous_force;
  geomVector unit_outwards_normal_vector;
  Node* neighbors[2]; // CHANGE FOR 3D
  Edge const* edge_of_neighbors[2]; // CHANGE FOR 3D
};


struct TriElement
{
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

      void display_parameters();

      //@}

   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Creating the surface points on a RBC **/
      void create_RBC_structure( );

      /** @brief Allocates memory & initializes node properties & attributes */
      virtual void initialize_node_properties() = 0;

      /** @brief Generates the mesh for 2D and 3D immersed bodies **/
      virtual void set_all_nodes() = 0;
      
      /** @brief Sets triangular elements for 2D and 3D immersed bodies */
      virtual void set_all_trielements() = 0;
      
      /** @brief Allocates memory & initializes edge properties & attributes */
      virtual void initialize_edge_properties() = 0;
      
      /** @brief Sets the edges and nodes for 2D and 3D immersed bodies */
      virtual void set_all_edges() = 0;
      
      /** @brief Projects shape of the membrane to sphere or biconcave **/
      virtual void project_membrane_shape() = 0;
      
      /** @brief Positions the membrane from centers of mass in RBC_data.csv */
      void position_membrane();
      
      /** @brief Rotates the membrane based on 3D rotation matrix 
      2D: rotation is achieved using ONLY "yaw" angle
      3D: rotation is achieved using roll, pitch and yaw angles */
      void rotate_membrane();
      
      /** @brief Computes current (& initial) spring lengths */
      virtual void compute_spring_lengths(bool init) = 0;
      
      /** @brief Computes the unit normal of each edge/spring */
      virtual void compute_edge_normals() = 0;
      
      /** @brief Computes angle between edges */
      virtual void compute_edge_angle(bool init) = 0;

      /** @brief Computes norm of a vector variable */
      virtual double norm(double const* v) = 0;
      
      /** @brief Scalar product of two vectors */
      virtual double scalar( double const* v0, double const* v1 ) = 0;
      
      /** @brief IBM: Eulerian velocity to Lagrangian velocity interpolation **/
      virtual void eul_to_lag() = 0;

      /** @brief writes one point of RBC mesh to .vtu file **/
      virtual void write_mesh_to_vtk_file( size_t IB_number, double const& time,
                                           size_t const& cyclenum ) = 0;
                                           
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

      // Vector containing all nodes properties
      vector<Node> m_all_nodes;
      vector<Edge> m_all_edges;

      size_t m_nEdges;
      size_t m_nTriangles;
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
