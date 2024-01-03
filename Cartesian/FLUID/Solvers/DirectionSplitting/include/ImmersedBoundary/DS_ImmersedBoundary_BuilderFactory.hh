#ifndef _DS_IMMERSEDBOUNDARY_BUILDERFACTORY__
#define _DS_IMMERSEDBOUNDARY_BUILDERFACTORY__

#include <iostream>
using std::istream;
class DS_ImmersedBoundary;
class FS_RigidBody;

/** @brief The Class DS_ImmersedBoundary_BuilderFactory.

Immersed Boundary builder factory, creates an instance of rigid body using a unsigned
integer (numebr of corners for a polyhedron or a code for other shapes) as
the main parameter.

@author A. Goyal - Pacific project 2024 */

class DS_ImmersedBoundary_BuilderFactory
{
   public: //-----------------------------------------------------------------

   //-- Static methods

      /** @name Static methods */
      //@{
      /** @brief Creates a Direction Splitting immersed boundary`
      @param pgrb pointer to the corresponding geometric rigid body */
      static DS_ImmersedBoundary* create(FS_RigidBody* pgrb);


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
