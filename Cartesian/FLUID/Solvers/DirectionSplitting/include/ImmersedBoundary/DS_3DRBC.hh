#ifndef _DS_3DRBC__
#define _DS_3DRBC__

#include <DS_ImmersedBoundary.hh>
#include <FV_Mesh.hh>
#include <FV_DiscreteField.hh>
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

      /** @brief Computes norm of a 3D vector or array variable */
      double norm(double const* v);
      
      /** @brief Scalar product of two vectors */
      double scalar( double const* v0, double const* v1 );
      
      /** @brief Cross product/vector product of two vectors */
      void cross_3D( double const* v0, double const* v1, double* res ); 
      
      /** @brief Computes node based spring and bending constants **/
      void preprocess_membrane_parameters(string const& case_type,
                                          size_t const& num_subtimesteps_RBC);
      
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
      
      /** @brief Computes the statistics of each RBC membrane */
      void compute_stats(string const& directory, string const& filename, 
                         size_t const& dim, double const& time, 
                         size_t const& cyclenum);
      
      /** @brief Computes the centroid of each immersed boundary */
      void compute_centroid(size_t const& dim);
      
      /** @brief Computes axial diameter of RBC: x_max - x_min */
      double compute_axial_diameter();
      
      /** @brief Computes transverse diameter of RBC: max(sqrt(y^2 + z^2)) */
      double compute_transverse_diameter();
      
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
      @param copy copied DS_3DRBC object */
      DS_3DRBC( DS_3DRBC const& copy );
      //@}
};

#endif
