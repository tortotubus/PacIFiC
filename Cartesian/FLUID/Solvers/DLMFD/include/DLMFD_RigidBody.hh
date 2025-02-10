#ifndef _DLMFD_RIGIDBODY__
#define _DLMFD_RIGIDBODY__

#include <FS_RigidBody.hh>

/** @brief The class DLMFD_RigidBody.

A moving or stationary rigid body in the Fictitious Domain solver.

@author M. Houlette - Pacific project 2025 */

class DLMFD_RigidBody
{
public: //-----------------------------------------------------------------
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Default constructor */
    DLMFD_RigidBody();

    /** @brief Constructor with arguments
    @param pgrb Pointer to the geometric rigid body class */
    DLMFD_RigidBody(FS_RigidBody *pgrb);

    /** @brief Destructor */
    ~DLMFD_RigidBody();

    //@}

protected: //--------------------------------------------------------------
    //-- Attributes

    FS_RigidBody *ptr_FSrigidbody;

private: //----------------------------------------------------------------
};

#endif