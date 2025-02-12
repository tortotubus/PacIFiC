#ifndef _DLMFD_SPHERE__
#define _DLMFD_SPHERE__

#include <DLMFD_RigidBody.hh>
#include <string>
using namespace std;

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

    //-- Set methods
    /** @name Set methods */
    //@{

    /** @brief Set DLMFD boundary and interior points
    @param critical_distance Critical distance */
    void set_all_points(double critical_distance);

    //@}

    //-- Get methods
    /** @name Get methods */
    //@{

    /** @brief Returns pointer to the rigid body gravity center */
    geomVector const *get_ptr_to_gravity_centre() const;

    /** @brief Returns circumscribed radius */
    double get_circumscribed_radius() const;

    /** @brief Get rigid body velocity
    @param point Point at which the velocity is computed */
    geomVector get_rigid_body_velocity(geomVector const &point) const;

    //@}

    //-- Updating methods
    /** @name Updating methods */
    //@{

    /** @brief Updating method */
    void update();

    /** @brief Updating method
    @param pos updated position
    @param vel updated translation velocity */
    void update_RB_position_and_velocity(geomVector const &pos,
                                         geomVector const &vel,
                                         geomVector const &ang_vel,
                                         vector<geomVector> const &periodic_directions, double const &time_step);

    //@}

protected: //----------------------------------------------------------------
private:   //----------------------------------------------------------------
};

#endif
