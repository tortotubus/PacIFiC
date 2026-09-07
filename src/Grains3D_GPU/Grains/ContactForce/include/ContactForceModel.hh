#ifndef _CONTACTFORCEMODEL_HH_
#define _CONTACTFORCEMODEL_HH_

#include "ContactInfo.hh"
#include "ContactTable.hh"
#include "Torce.hh"
#include "Vector3.hh"
#include <ostream>
#include <string>

// ContactForceModel types
enum ContactForceModelType
{
    HOOKE,
    HOOKEMEMORY,
    CF_VIRTUAL  // Virtual dispatch fallback (kernel template parameter only)
};

// =================================================================================================
/** @brief The class ContactForceModel.

    Defines the contact forces between two colliding components and computes
    these contact forces.

    @author A.YAZDANI - 2024 - Construction */
// =================================================================================================
template <typename T>
class ContactForceModel
{
protected:
    /**@name Contructors */
    //@{
    /** @brief Default constructor (forbidden except in derived classes) */
    __HOSTDEVICE__
    ContactForceModel();

    /** @brief Copy constructor
        @param cf ContactForceModel object to be copied */
    __HOSTDEVICE__
    ContactForceModel(ContactForceModel<T> const& cf);
    //@}

public:
    /**@name Contructors */
    //@{
    /** @brief Destructor */
    __HOSTDEVICE__
    virtual ~ContactForceModel();
    //@}

    /** @name Get methods */
    //@{
    /** @brief Returns the ContactForceModel type */
    __HOSTDEVICE__
    virtual ContactForceModelType getContactForceModelType() const = 0;
    //@}

    /** @name Methods */
    //@{
    /** @brief Returns a torce based on the contact information
        @param contactInfos geometric contact features
        @param relVelocityAtContact relative velocity at the contact point
        @param relAngVelocity relative angular velocity
        @param vA position of the first component
        @param vB position of the second component
        @param contactHistory pointer to contact history
        @param torceA computed force and torque for the first component
        @param torceB computed force and torque for the second component */
    __HOSTDEVICE__
    virtual void computeForces(const ContactInfo<T>& contactInfos,
                               const Vector3<T>&     relVelocityAtContact,
                               const Vector3<T>&     relAngVelocity,
                               const Vector3<T>&     vA,
                               const Vector3<T>&     vB,
                               ContactHistory<T>*    contactHistory,
                               Torce<T>&             torceA,
                               Torce<T>&             torceB) const
        = 0;

    /** @brief Computes and prints contact force parameter estimates for a binary collision.
        Analytically solves the spring-dashpot ODE to compute the contact time Tc and maximum
        penetration depth delta_max for a head-on collision at relative velocity v0. Also reports
        the maximum elastic force and its ratio to particle weight, and warns if delta_max exceeds
        the crust budget (which would cause GJK to fail or produce incorrect results).
        @param massA         mass of the first component [kg] (use 1e20 for obstacles)
        @param massB         mass of the second component [kg] (use 1e20 for obstacles)
        @param crustA        crust thickness of the first component [m]
        @param crustB        crust thickness of the second component [m]
        @param v0            pre-collisional relative velocity magnitude [m/s]
        @param dt            simulation time step [s]
        @param labelA        display label for the first component
        @param labelB        display label for the second component
        @param out           output stream (e.g. std::cout) */
    __HOST__
    virtual void computeEstimates(T             massA,
                                  T             massB,
                                  T             crustA,
                                  T             crustB,
                                  T             v0,
                                  T             dt,
                                  std::string   labelA,
                                  std::string   labelB,
                                  std::ostream& out) const
        = 0;
    //@}
};

#endif
