#ifndef _DLMFD_RIGIDBODY_BUILDERFACTORY__
#define _DLMFD_RIGIDBODY_BUILDERFACTORY__

#include <FV_DiscreteField.hh>
#include <iostream>
using namespace std;
class FS_RigidBody;
class DLMFD_RigidBody;
class FV_Mesh;

/** @brief The Class DLMFD_RigidBody_BuilderFactory.

Rigid body builder factory, creates an instance of rigid body using a unsigned
integer (number of corners for a polyhedron or a code for other shapes) as
the main parameter.

@author A. Wachs - Pacific project 2025 */

class DLMFD_RigidBody_BuilderFactory
{
  public: //-----------------------------------------------------------------
          //-- Static methods
    /** @name Static methods */
    //@{
    /** @brief Creates a Fictitious Domain rigid body
    @param ptr_geom_rb pointer to the corresponding geometric rigid body */
    static DLMFD_RigidBody *create(FS_RigidBody *ptr_geom_rb,
                                   const bool &are_particles_fixed,
                                   FV_DiscreteField *pField_,
                                   double const critical_distance_);

  protected: //--------------------------------------------------------------
  private:   //----------------------------------------------------------------
             //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{
    /** @brief Default constructor */
    DLMFD_RigidBody_BuilderFactory();

    /** @brief Destructor */
    ~DLMFD_RigidBody_BuilderFactory();
    //@}
};

#endif
