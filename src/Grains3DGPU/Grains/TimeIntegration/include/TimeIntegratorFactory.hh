#ifndef _TIMEINTEGRATORFACTORY_HH_
#define _TIMEINTEGRATORFACTORY_HH_

#include "GrainsMemBuffer.hh"
#include "TimeIntegrator.hh"

// =================================================================================================
/** @brief The class TimeIntegratorFactory.

    Creates the numerical scheme for the time integration of the Newton's law
    and the kinematic equations.

    @author A.WACHS - Institut Francais du Petrole - 2011 - Creation
    @author A.WACHS - 2021 - Major cleaning & refactoring
    @author A.YAZDANI - 2024 - Major cleaning for porting to GPU */
// =================================================================================================
template <typename T>
class TimeIntegratorFactory
{
private:
    /**@name Contructors & Destructor */
    //@{
    /** @brief Default constructor (forbidden) */
    TimeIntegratorFactory();

    /** @brief Destructor (forbidden) */
    ~TimeIntegratorFactory();
    //@}

public:
    /**@name Methods */
    //@{
    /** @brief Creates and returns the time integration scheme.
        @param root XML node
        @param dt time step
        @param TI time integrator buffer */
    __HOST__
    static void create(grains::xml::OptionalElement root, T dt,
                       GrainsMemBuffer<TimeIntegrator<T>*, MemType::HOST>& TI);

    //@}
};

#endif
