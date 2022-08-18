#ifndef _DS_2DRBC__
#define _DS_2DRBC__

#include <DS_ImmersedBoundary.hh>
#include <string>
using std::string;


/** @brief The class DS_2DRBC.

A moving or stationary rigid 2D RBC of axisymmetric cross-section in the
Direction Splitting solver.

@author A. Goyal - Pacific project 2022 */

class DS_2DRBC: public DS_ImmersedBoundary
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_2DRBC();

      /** @brief Destructor */
      ~DS_2DRBC();
      //@}


   //-- Set methods

      /**@name Set methods */
      //@{
      /** @brief writes one point of RBC mesh to .vtu file **/
      void write_one_point_to_VTK( double const& time
                                         , size_t const& cyclenum );
      //@}


   //-- Methods

      /**@name Methods */
      //@{
      /** @brief writes one point of RBC mesh to .vtu file
      @param nb_edges asgasdg
      @param nb_edges asgasdg
      @param nb_edges asgasdg
      @param nb_edges asgasdg
       **/
      void create_RBC_structure( size_t const& nb_edges
                               , double const& radius
                               , double const& c0
                               , double const& c1
                               , double const& c2
                               , double const& rbc_orientation_angle);

      /** @brief Initializes all variables of 'Struct Node' */
      void initialize_node_properties();

      /** @brief Generates the mesh for the 3D immersed body */
      void generate_membrane_mesh();

      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{

      //@}


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Copy constructor
      @param copy copied DS_2DRBC object */
      DS_2DRBC( DS_2DRBC const& copy );
      //@}
};

#endif
