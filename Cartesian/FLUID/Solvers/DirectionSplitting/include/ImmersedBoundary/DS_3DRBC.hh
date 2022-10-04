#ifndef _DS_3DRBC__
#define _DS_3DRBC__

#include <DS_ImmersedBoundary.hh>
#include <FV_Mesh.hh>
#include <FV_DiscreteField.hh>
#include <string>
#include <fstream>
#include <sstream>
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
      /** @brief Allocates memory & initializes all variables of 'Struct Node' */
      void initialize_node_properties(string const& mesh_filename,
                                      size_t const& dim);

      /** @brief Sets all nodes - number, coordinates, neighbours */
      void set_all_nodes(istream& fileIN, size_t const& dim);
      
      /** @brief Computes the radius of RBC from coordinates ASSUMING centroid
      of membrane to be (0, 0, 0) */
      void set_radius(size_t const& num_nodes, size_t const& dim);
      
      /** @brief Scales the node coordinates and RBC radius */
      void scale_all_node_coordinates(size_t const& num_nodes, 
                                      size_t const& dim);
      
      /** @brief Computes the area, normals and centre of mass of each triangle */
      void compute_triangle_area_normals_centre_of_mass(bool init,
                                                        size_t const& dim);
      
      /** @brief Sets the neighbors of each node from triangle nodes */
      void set_all_node_neighbors();
      
      /** @brief Projects shape of the membrane to sphere or biconcave **/
      void project_membrane_shape();
      
      /** @brief Allocates memory & initializes all variables of 'Struct TriElements' */
      void initialize_triangle_properties(size_t const& num_triangles,
                                          size_t const& dim);
      
      /** @brief Sets all trielements - nodes, centers of mass, normals */
      void set_all_trielements(istream& fileIN, size_t const& dim,
                               bool const& MatlabNumb);
      
      /** @brief Allocates memory & initializes edge properties & attributes */
      void initialize_edge_properties(size_t const& dim);
      
      /** @brief Returns a pointer to the edge if it already exists in the list,
      otherwise returns NULL */
      Edge* does_edge_exist( Node const* n2_, Node const* n3_ );

      /** @brief Returns whether 2 sets of coordinates are identical up to eps */
      bool compare( geomVector const& v0, geomVector const& v1, double const& eps );
      
      /** @brief Sets all edges and its nodes and normals */
      void set_all_edges(size_t const& dim);
      
      /** @brief Compute current (& initial) spring lengths */
      void compute_spring_lengths(bool init, size_t const& dim);
      
      /** @brief Compute unit normal of each edge */
      void compute_edge_normals();

      /** @brief Computes angle between edges */
      void compute_edge_angle(bool init);
      
      /** @brief Computes the normal to each triangle = 0.5 * area of triangle */
      void compute_twice_area_vector
      (vector< pair<Node const*,Node const*> > const& ppnodes, geomVector& res );
      
      /** @brief Computes norm of a 3D vector or array variable */
      double norm(geomVector& v);
      
      /** @brief Scalar product of two vectors */
      double scalar( geomVector const& v0, geomVector const& v1 );
      
      /** @brief Cross product/vector product of two vectors */
      void cross_3D( geomVector const& v0, geomVector const& v1, geomVector& res );
      
      /** @brief Computes node based spring and bending constants **/
      void preprocess_membrane_parameters(string const& case_type,
                                          double const& mu,
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
                                 double const& spring_constant,
                                 string const& force_type = "Fedosov" );
      
      /** @brief Computes spring force */
      void compute_linear_spring_force( size_t const& dim, 
                                 double const& spring_constant );
      
      /** @brief Computes bending force */
      void compute_bending_resistance( size_t const& dim, 
                                       double const& bending_spring_constant,
                                       double const& bending_viscous_constant, 
                                       double const& dt,
                                       string const& force_type = "Fedosov" );

      /** @brief Computes viscous force */
      void compute_viscous_drag_force( size_t const& dim, 
                                       double const& viscous_drag_constant );

      /** @brief Lagrangian to Eulerian force spreading **/
      void lag_to_eul(FV_DiscreteField* FF, FV_DiscreteField* FF_tag,
                      size_t const& dim, size_t const& comp);
      
      /** @brief Computes total area & total volume of RBC membrane */
      void compute_total_surface_area_total_volume( bool init, 
                                                    size_t const& dim );

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
      
      /** @brief Computes Taylor Deformation Parameter & Orientation angle */
      void compute_tdp_orientation_angle();
      
      /** @brief Computes the average tangential velocity */
      double compute_avg_tangential_velocity();
      
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
