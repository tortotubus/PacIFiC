#ifndef _DS_ALLIMMERSEDBOUNDARY__
#define _DS_ALLIMMERSEDBOUNDARY__

#include <geomVector.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <doubleArray2D.hh>
#include <boolArray2D.hh>
#include <MAC_Communicator.hh>
#include <MAC_DoubleVector.hh>
#include <FV_Mesh.hh>
#include <FV_TimeIterator.hh>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using std::vector;
using std::istream ;
using std::ostream ;
using std::istringstream ;
using std::string;
class DS_ImmersedBoundary;
class FS_AllImmersedBoundary;
class FV_DiscreteField;
class FV_Mesh;
class FV_TimeIterator;

/** @brief The class DS_AllImmersedBoundary.

The array of all rigid bodies in the Direction Splitting solver.


@author A. Goyal - Pacific project 2022 */

class DS_AllImmersedBoundary
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_AllImmersedBoundary(size_t const& space_dimension
                           , string const& IB_file
                           , size_t const& N_IB
                           , string const& case_type
                           , FV_DiscreteField const* arb_UF
                           , FV_DiscreteField* arb_EulF
                           , FV_DiscreteField* arb_EulF_tag
                           , size_t const& n_RBC_timesteps
                           , string const& dirac_type
                           , size_t const& periodic_dir);

      /** @brief Destructor */
      ~DS_AllImmersedBoundary();
      //@}


   //-- Get methods

      /**@name Get methods */
      //@{
      /** @brief Returns a pointer to the ith DS rigid body */
      DS_ImmersedBoundary* get_ptr_immersed_body( size_t i );

      /** @brief Returns number of immersed boundaries */
      size_t get_number_of_immersed_boundaries() const;

      //@}


   //-- Set methods

      /**@name Set methods */
      //@{

      //@}


   //-- Methods

      /**@name Methods */
      //@{

        /** @brief Function which calls RBC and IBM functions along with
        periodic boundary conditions and parallelisation temporary variables */
        void do_one_inner_iteration( FV_TimeIterator const* t_it );
        
      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{

      //@}

   private: //----------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      size_t m_space_dimension; /**< Space dimension */
      size_t m_nIB; /**< number of immersed boundaries */
      string m_IB_file; /** input file name containing RBC location & shape */
      string m_IB_case_type; /* = Breyannis2000 or other cases */
      size_t m_subtimesteps_RBC; /* number of subtimesteps for RBC iterations */
      string m_dirac_type; /* type of Dirac delta - Balogh, Archer, Roma */
      size_t m_periodic_dir; /* periodic direction - 0 for x, 1 for y & 2 for z*/
      vector<DS_ImmersedBoundary*> m_allDSimmersedboundary; /** pointer of objects
      of DS_ImmersedBoundary class */

      // Pointers to the constant fields and primary grid
      FV_DiscreteField const* UF ;
      FV_DiscreteField* Eul_F ;
      FV_DiscreteField* F_Eul_tag;
      FV_Mesh const* MESH ;

      //@}


    //-- Constructors & Destructor

       /**@name Constructors & Destructor */
       //@{
       /** @brief Copy constructor
       @param copy copied DS_AllImmersedBoundary object */
       DS_AllImmersedBoundary( DS_AllImmersedBoundary const& copy );
       //@}

     //-- Methods

        /**@name Methods */
        //@{
        /** @brief Read the CSV file with RBC parameters */
        void read_shape_and_membrane_parameters();

        /** @brief Read the CSV file with RBS parameters */
        void initialize_variables();

        /** @brief Reads number of lines in csv input data file */
        size_t get_num_lines_in_IB_file();

        /** @brief Generates the mesh for 2D and 3D immersed body */
        void generate_immersed_body_mesh();
        
        /** @brief Projects the shape of 2D and 3D immersed body to
        circle/sphere or biconcave disk/discoid shape */
        void project_shape_of_immersed_body();
        
        /** @brief Positions the immersed body from centers of mass 
        given in csv file fed as input */
        void position_immersed_body();
        
        /** @brief Rotates the immersed body based on xroll, ypitch
        and zyaw angles from csv file fed as input */
        void rotate_immersed_body();
        
        /** @brief Writes the immersed body mesh to .vtu file */
        void write_immersed_body_mesh_to_vtk_file();
        
        /** @brief Computes node based spring, bending constants */
        void preprocess_immersed_body_parameters(string const& case_type, 
                                            size_t const& num_subtimesteps_RBC);
        
        /** @brief Sets the IBM parameters for all immersed bodies */
        void set_IBM_parameters(string const& dirac_type
                              , size_t const& periodic_dir);
        
        /** @brief IBM:Eulerian velocity to Lagrangian velocity interpolation */
        void eul_to_lag_velocity_interpolate();
        //@}
};

#endif
