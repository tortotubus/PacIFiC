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

struct ShapeParameters 
{
  geomVector center;
  double radius;
  double xroll, ypitch, zyaw;
  double c0, c1, c2;
  size_t N_nodes;
  size_t N_levels;
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
  double initial_angle, angle_nm1, dangle_dt;
  geomVector spring_force, bending_force, area_force, volume_force, viscous_force;
  geomVector unit_outwards_normal_vector;
  vector<size_t> neighbors;
  // size_t_vector neighbors;
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

      /** @brief Initialize node properties **/
      virtual void initialize_node_properties() = 0;

      /** @brief Generates the mesh for 2D and 3D immersed bodies **/
      virtual void set_all_nodes() = 0;
      
      /** @brief Sets triangular elements for 2D and 3D immersed bodies */
      virtual void set_all_trielements() = 0;
      
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
      
      /** @brief IBM: Eulerian velocity to Lagrangian velocity interpolation **/
      virtual void eul_to_lag() = 0;

      /** @brief writes one point of RBC mesh to .vtu file **/
      virtual void write_mesh_to_vtk_file( size_t IB_number, double const& time,
                                           size_t const& cyclenum ) = 0;
                                           
      /** @brief Converts size_t to string datatype **/
      string sizetToString( size_t const& figure ) const;

      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      // Shape parameters
      ShapeParameters shape_param;

      // Vector containing all nodes properties
      vector<Node> m_all_nodes;

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
