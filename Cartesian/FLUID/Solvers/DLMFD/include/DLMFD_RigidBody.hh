#ifndef _DLMFD_RIGIDBODY__
#define _DLMFD_RIGIDBODY__

#include <FS_RigidBody.hh>
#include <geomVector.hh>
#include <DLMFD_BoundaryMultiplierPoint.hh>
#include <DLMFD_InteriorMultiplierPoint.hh>
using namespace std;

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

    //-- Set methods
    /** @name Set methods */
    //@{

    /** @brief Set DLMFD boundary and interior points
    @param critical_distance Critical distance */
    virtual void set_all_points(double critical_distance) = 0;

    //@}

    //-- Get methods
    /** @name Get methods */
    //@{

    /** @brief Returns pointer to the rigid body gravity center */
    virtual geomVector const *get_ptr_to_gravity_centre() const = 0;

    /** @brief Returns circumscribed radius */
    virtual double get_circumscribed_radius() const = 0;

    //@}

    //-- Add methods
    /** @name Add methods */
    //@{

    /** @brief Add a boundary point to the list of boundary points
    @param point Point
    @param critical_distance Critical distance */
    void add_boundary_point(const geomVector &point, double critical_distance);

    //@}

    //-- Update methods
    /** @name Update methods */
    //@{

    /** @brief Update the RB position and velocity
    @param pos updated position
    @param vel updated translation velocity */
    virtual void update_RB_position_and_velocity(geomVector const &pos,
                                                 geomVector const &vel,
                                                 geomVector const &ang_vel,
                                                 vector<geomVector> const &periodic_directions,
                                                 double const &time_step) = 0;

    //@}

protected: //--------------------------------------------------------------
    //-- Attributes

    FS_RigidBody *ptr_FSrigidbody; /* Pointer to geometric Rigid Body */

    // Geometric attributes
    list<DLMFD_BoundaryMultiplierPoint *> boundarypoints; /* List of pointers to boundary points */
    list<DLMFD_InteriorMultiplierPoint *> interiorpoints; /* List of pointers to interior points */

    geomVector gravity_center; /**< center of gravity coordinates */
    double radius;             /**< external radius */

private: //----------------------------------------------------------------
};

#endif