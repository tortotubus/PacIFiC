#ifndef _DLMFD_SPHERE__
#define _DLMFD_SPHERE__

#include <DLMFD_RigidBody.hh>
#include <string>
using std::string;

/** @brief The class DLMFD_Sphere.

A moving or stationary rigid sphere in the Fictitious Domain solver.

@author A. Wachs - Pacific project 2025 */

class DLMFD_Sphere : public DLMFD_RigidBody
{
public: //------------------------------------------------------------------
    DLMFD_Sphere();

    DLMFD_Sphere(FS_RigidBody *pgrb);

    ~DLMFD_Sphere();

    void update();

    geomVector get_rigid_body_velocity(geomVector const &point) const;

protected: //----------------------------------------------------------------
private:   //----------------------------------------------------------------
};

#endif
