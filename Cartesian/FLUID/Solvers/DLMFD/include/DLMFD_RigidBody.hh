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

        ~DLMFD_RigidBody();

        geomVector get_rigid_body_velocity(geomVector const& point);

    
    protected: //--------------------------------------------------------------
   
    private: //----------------------------------------------------------------

        FS_RigidBody* ptr_FSrigidbody;
};

#endif