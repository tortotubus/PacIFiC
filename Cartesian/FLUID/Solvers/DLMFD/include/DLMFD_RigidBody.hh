#ifndef _DLMFD_RIGIDBODY__
#define _DLMFD_RIGIDBODY__

#include <FS_RigidBody.hh>

/** @brief The class DLMFD_RigidBody.

A moving or stationary rigid body in the Fictitious Domain solver.

@author M. Houlette - Pacific project 2025 */

class DLMFD_RigidBody
{
public: //-----------------------------------------------------------------
    DLMFD_RigidBody();

    DLMFD_RigidBody(FS_RigidBody *pgrb);

    ~DLMFD_RigidBody();

protected: //--------------------------------------------------------------
    FS_RigidBody *ptr_FSrigidbody;

private: //----------------------------------------------------------------
};

#endif