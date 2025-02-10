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
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Default constructor */
    DLMFD_Sphere();

    /** @brief Constructor with arguments
    @param pgrb Pointer to the geometric rigid body class */
    DLMFD_Sphere(FS_RigidBody *pgrb);

    /** @brief Destructor */
    ~DLMFD_Sphere();

    //@}

    //-- Get methods
    /** @name Get methods */
    //@{

    /** @brief Get rigid body velocity
    @param point Point at which the velocity is computed */
    geomVector get_rigid_body_velocity(geomVector const &point) const;

    //@}

    //-- Updating methods
    /** @name Updating methods */
    //@{

    /** @brief Updating method */
    void update();

    //@}

protected: //----------------------------------------------------------------
private:   //----------------------------------------------------------------
};

#endif
