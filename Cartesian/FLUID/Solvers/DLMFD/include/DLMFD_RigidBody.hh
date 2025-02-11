#ifndef _DLMFD_RIGIDBODY__
#define _DLMFD_RIGIDBODY__

#include <FS_RigidBody.hh>
#include <DLMFD_BoundaryMultiplierPoint.hh>
#include <DLMFD_InteriorMultiplierPoint.hh>

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

    /** @brief Set DLMFD boundary and interior points */
    virtual void set_all_points() = 0;

protected: //--------------------------------------------------------------
    //-- Attributes

    FS_RigidBody *ptr_FSrigidbody; /* Pointer to geometric Rigid Body */

    list<DLMFD_BoundaryMultiplierPoint *> boundarypoints; /* List of pointers to boundary points */
    list<DLMFD_InteriorMultiplierPoint *> interiorpoints; /* List of pointers to interior points */

private: //----------------------------------------------------------------
};

#endif