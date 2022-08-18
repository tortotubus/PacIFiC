#ifndef _DS_IMMERSEDBOUNDARY_BUILDERFACTORY__
#define _DS_IMMERSEDBOUNDARY_BUILDERFACTORY__

#include <iostream>
using std::istream;
class DS_ImmersedBoundary;

/** @brief The Class DS_ImmersedBoundary_BuilderFactory.

Immersed Boundary builder factory, creates an instance of immersed body using an
unsigned integer (numebr of corners for a polyhedron or a code for other shapes)
as the main parameter.

@author A. Goyal - Pacific project 2022 */

class DS_ImmersedBoundary_BuilderFactory
{
   public: //-----------------------------------------------------------------

   //-- Static methods

      /** @name Static methods */
      //@{
      /** @brief Creates a Direction Splitting immersed boundary`
      @param pgrb pointer to the corresponding geometric rigid body */
      static DS_ImmersedBoundary* create(size_t const& m_space_dimension);


   protected: //--------------------------------------------------------------


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_ImmersedBoundary_BuilderFactory();

      /** @brief Destructor */
      ~DS_ImmersedBoundary_BuilderFactory();
      //@}
};

#endif
