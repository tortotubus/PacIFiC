#include "TimeIntegratorFactory.hh"
#include "FirstOrderExplicit.hh"
#include "SecondOrderLeapFrog.hh"

/* ============================================================================================== */
/* High-Level Methods                                                                             */
/* ============================================================================================== */
// Creates and stores a TimeIntegrator object in the host memory.
template <typename T>
__HOST__ void TimeIntegratorFactory<T>::create(
    grains::xml::OptionalElement root, T dt,
    GrainsMemBuffer<TimeIntegrator<T>*, MemType::HOST>& TI)
{
    TI.initialize(1);  // Initialize memory for one time integrator

    std::string type
        = grains::xml::requireElement(root).requireAttribute<std::string>("Type");
    if(type == "FirstOrderExplicit")
        TI[0] = new FirstOrderExplicit<T>(dt);
    else if(type == "SecondOrderLeapFrog")
        TI[0] = new SecondOrderLeapFrog<T>(dt);
    else
        GAbort("Unknown time integration! Aborting Grains!");
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class TimeIntegratorFactory<float>;
template class TimeIntegratorFactory<double>;
