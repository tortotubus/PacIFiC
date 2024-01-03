#ifndef _DS_ALLIMMERSEDBOUNDARIES__
#define _DS_ALLIMMERSEDBOUNDARIES__

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
class FS_AllRigidBodies;
class FV_DiscreteField;
class FV_Mesh;
class FV_TimeIterator;

/** @brief The class DS_AllImmersedBoundaries.

The array of all rigid bodies in the Direction Splitting solver.


@author A. Goyal - Pacific project 2022 */

class DS_AllImmersedBoundaries
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_AllImmersedBoundaries();

      /** @brief Constructor with arguments
      @param dimens number of space dimensions
      @param in input stream where features of rigid bodies are read
      @param b_IB_as_fixed_obstacles treat all IB as fixed obstacles
      @param arb_UF Pointer to flow field UF
      @param arb_scs scale of cell on the rigid body surface as
      compared with the cell of computational grid
      @param arb_macCOMM communicator for MPI communications */
      DS_AllImmersedBoundaries(size_t &dimens
                              , istream &in
                              , bool const &b_IB_as_fixed_obstacles
                              , FV_DiscreteField const *arb_UF
                              , double const& arb_scs 
                              , MAC_Communicator const *arb_macCOMM);

      /** @brief Destructor */
      ~DS_AllImmersedBoundaries();
      //@}


   //-- Get methods

      /**@name Get methods */
      //@{
      /** @brief Returns a pointer to the ith DS rigid body */
      DS_ImmersedBoundary* get_ptr_rigid_body( size_t i );

      /** @brief Returns number of immersed boundaries */
      size_t get_number_immersed_boundaries() const;

      //@}


   //-- Set methods

      /**@name Set methods */
      //@{

      //@}


   //-- Methods

      /**@name Methods */
      //@{

      /** @brief Intialize the surface variables for all rigid bodies */
      void initialize_surface_variables_for_all_IB();

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
      vector<DS_ImmersedBoundary*> m_allDSimmersedboundaries; /**< the vector of all
    	Direction Splitting rigid bodies */

      // Pointers to the constant fields and primary grid
      FV_DiscreteField const* UF ;
      FV_Mesh const* MESH ;

      double surface_cell_scale; /**< a variable to store the scale of surface
      cell on the RB as compared with computational grid cell size */
      MAC_Communicator const *m_macCOMM; /**< Variable for communication
      between processors */
      FS_AllRigidBodies* m_FSallrigidbodies; /**< the pointer to the
    	FS_AllRigidBodies object that contains the vector of all
    	corresponding geometric rigid bodies */

      //@}

      //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Copy constructor
      @param copy copied DS_AllImmersedBoundaries object */
      DS_AllImmersedBoundaries(DS_AllImmersedBoundaries const &copy);
      //@}
};

#endif
