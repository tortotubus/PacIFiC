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
      void write_mesh_to_vtk_file( size_t IB_number, double const& time,
                                   size_t const& cyclenum );
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

      /** @brief Sets all nodes - number, coordinates, neighbours */
      void set_all_nodes();
      
      /** @brief Sets all trielements - nodes, centers of mass, normals */
      void set_all_trielements();
      
      /** @brief Sets all edges and its nodes and normals */
      void set_all_edges();
      
      /** @brief Projects shape of the membrane to sphere or biconcave **/
      void project_membrane_shape();
      
      /** @brief IBM: Eulerian velocity to Lagrangian velocity interpolation */
      void eul_to_lag();

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
