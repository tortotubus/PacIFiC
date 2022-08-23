#ifndef _DS_3DRBC__
#define _DS_3DRBC__

#include <DS_ImmersedBoundary.hh>
#include <string>
using std::string;


/** @brief The class DS_3DRBC.

A moving or stationary rigid 3D RBC of axisymmetric cross-section in the
Direction Splitting solver.

@author A. Goyal - Pacific project 2022 */

class DS_3DRBC: public DS_ImmersedBoundary
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_3DRBC();

      /** @brief Destructor */
      ~DS_3DRBC();
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
      
      /** @brief writes one point of RBC mesh to .vtu file **/
      void write_mesh_to_vtk_file( size_t IB_number, double const& time,
                                   size_t const& cyclenum );

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
      @param copy copied DS_3DRBC object */
      DS_3DRBC( DS_3DRBC const& copy );
      //@}
};

#endif
