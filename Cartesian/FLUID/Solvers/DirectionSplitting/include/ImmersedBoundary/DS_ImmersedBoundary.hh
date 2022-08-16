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
#include <iostream>
#include <fstream>
using std::ostream;
using std::vector;
using std::tuple;
class FV_DiscreteField;
class FV_Mesh;

struct ShapeParameters {
  double c0;
  double c1;
  double c2;
  size_t N_nodes;
  size_t N_levels;
  // Gravity centre and radius
  geomVector center;
  double radius;
};

struct Node {
  size_t number;
  geomVector coordinates;
  geomVector coordinates_pbc;
  geomVector sumforce;
  geomVector velocity;
  geomVector angular_velocity;
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
      ShapeParameters* get_ptr_shape_parameters( );

      void display_parameters();

      //@}

   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Creating the surface points on a RBC **/
      void create_RBC_structure( );

      /** @brief Initialize node properties **/
      virtual void initialize_node_properties( ) = 0;

      /** @brief writes one point of RBC mesh to .vtu file **/
      virtual void write_one_point_to_VTK( double const& time
                                         , size_t const& cyclenum ) = 0;

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
