#ifndef _DS_ALLIMMERSEDBOUNDARY__
#define _DS_ALLIMMERSEDBOUNDARY__

#include <geomVector.hh>
#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <doubleArray2D.hh>
#include <boolArray2D.hh>
#include <MAC_Communicator.hh>
#include <MAC_DoubleVector.hh>
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
                           , size_t const& N_IB);

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
      vector<DS_ImmersedBoundary*> m_allDSimmersedboundary; /** pointer of objects
      of DS_ImmersedBoundary class */

      // Pointers to the constant fields and primary grid
      FV_DiscreteField const* UF ;
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
        /** @brief Read the CSV file with RBS parameters */
        void read_shape_parameters();

        /** @brief Read the CSV file with RBS parameters */
        void initialize_variables();

        /** @brief Reads number of lines in csv input data file */
        size_t get_num_lines_in_IB_file();

        /** @brief Generates the mesh for 2D and 3D immersed body */
        void generate_immersed_body_mesh();
        //@}
};

#endif
