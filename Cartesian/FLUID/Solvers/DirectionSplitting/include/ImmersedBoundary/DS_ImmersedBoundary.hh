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
class FS_RigidBody;
class FV_DiscreteField;
class FV_Mesh;

struct Node
{
   vector<double> position;
   vector<double> velocity;
   vector<double> normal;
   vector<double> force;
};



/** @brief The class DS_ImmersedBoundary.

A moving or stationary rigid body in the Direction Splitting solver.

@author A. Goyal - Pacific project 2024 */

class DS_ImmersedBoundary
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_ImmersedBoundary();

      /** @brief Constructor with arguments
      @param pgrb pointer to the corresponding geometric rigid body */
      DS_ImmersedBoundary(FS_RigidBody *pgrb);

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


      //@}

   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Initializa the surface variables for a immersed boundary
      such as points, normals, and area */
      void initialize_surface_variables();

      /** @brief writes one point of RBC mesh to .vtu file **/
      void write_one_IB_to_VTK( double const& time
                              , size_t const& cyclenum );

      /** @brief Compute the surface points by discretizing the immersed boundary
      surface in approximately equal areas (if possible) */
      virtual void compute_surface_points() = 0;

      /** @brief Compute number of discrete points on IBM
      @param surface_cell_scale scale of surface cell compared with the grid
      @param dx grid size */
      virtual void compute_number_of_surface_variables(
               double const& surface_cell_scale, double const& dx ) = 0;

      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      FS_RigidBody *m_geometric_immersed_body; /**< Pointer to the corresponding
            geometric immersed body */
      vector<struct Node> m_all_nodes;
      size_t Ntot; /** < Stores the total number of surface points 
      on a Immersed boundary */
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
      DS_ImmersedBoundary(DS_ImmersedBoundary const &copy);
      //@}
};

#endif
