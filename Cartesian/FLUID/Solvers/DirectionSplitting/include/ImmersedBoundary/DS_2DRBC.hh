#ifndef _DS_2DRBC__
#define _DS_2DRBC__

#include <DS_ImmersedBoundary.hh>
#include <FV_Mesh.hh>
#include <FV_DiscreteField.hh>
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

      /** @brief Allocates memory & initializes node properties & attributes */
      void initialize_node_properties();

      /** @brief Sets all nodes - number, coordinates, neighbours */
      void set_all_nodes();
      
      /** @brief Projects shape of the membrane to sphere or biconcave **/
      void project_membrane_shape();
      
      /** @brief Sets all trielements - nodes, centers of mass, normals */
      void set_all_trielements();
      
      /** @brief Allocates memory & initializes edge properties & attributes */
      void initialize_edge_properties();
      
      /** @brief Sets all edges and its nodes and normals */
      void set_all_edges();
      
      /** @brief Compute current (& initial) spring lengths */
      void compute_spring_lengths(bool init);
      
      /** @brief Compute unit normal of each edge */
      void compute_edge_normals();
      
      /** @brief Computes angle between edges */
      void compute_edge_angle(bool init);

      /** @brief Computes norm of a 2D vector or array variable */
      double norm(double const* v);
      
      /** @brief Scalar product of two vectors */
      double scalar( double const* v0, double const* v1 );
      
      /** @brief Cross product/vector product of two vectors */
      double cross_2D( geomVector const v0, geomVector const v1 );
      
      /** @brief Computes node based spring and bending constants **/
      void preprocess_membrane_parameters(string const& case_type,
                                          size_t const& num_subtimesteps_RBC);
      
      /** @brief Applies periodic boundary condition to each membrane **/
      void apply_periodic_boundary_conditions(FV_Mesh const* MESH,
                                              size_t const& dim);

      /** @brief IBM: Eulerian velocity to Lagrangian velocity interpolation */
      void eul_to_lag(FV_DiscreteField const* FF, size_t const& dim, 
                      size_t const& comp);

      /** @brief Computes RBC deformation using spring-dashpot model */
      void rbc_dynamics_solver(size_t const& dim, 
                               double const& dt_fluid,
                               string const& case_type);
      
      /** @brief Computes spring force */
      void compute_spring_force( size_t const& dim, 
                                 double const& spring_constant );
      
      /** @brief Computes spring force */
      void compute_linear_spring_force( size_t const& dim, 
                                 double const& spring_constant );
      
      /** @brief Computes bending force */
      void compute_bending_resistance( size_t const& dim, 
                                       double const& bending_spring_constant,
                                       double const& bending_viscous_constant, 
                                       double const& dt );

      /** @brief Computes viscous force */
      void compute_viscous_drag_force( size_t const& dim, 
                                       double const& viscous_drag_constant );

      /** @brief Lagrangian to Eulerian force spreading **/
      void lag_to_eul(FV_DiscreteField* FF, FV_DiscreteField* FF_tag,
                      size_t const& dim, size_t const& comp);
      
      /** @brief Computes the perimeter */
      double perimeter();
      
      /** @brief Writes one point of RBC mesh to .vtu file **/
      void write_mesh_to_vtk_file( size_t IB_number, double const& time,
                                   size_t const& cyclenum );
      
      /** @brief Writes the rbc.pvd file for each RBC membrane */
      void write_rbc_dot_pvd_file();
      
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
