#ifndef GRAINS_FORCE_COMPUTATION_CPU_HH
#define GRAINS_FORCE_COMPUTATION_CPU_HH

template <typename T>
struct ContactForceView;
template <typename T>
struct ExternalForceView;
template <typename T>
struct CompositeForceView;

namespace grains::cpu
{
    // Host arrays; sequential accumulation completes before these functions return.
    void computeContactForces(ContactForceView<float> view);
    void computeContactForces(ContactForceView<double> view);
    void addExternalForces(ExternalForceView<float> view);
    void addExternalForces(ExternalForceView<double> view);
    void assembleCompositeForces(CompositeForceView<float> view);
    void assembleCompositeForces(CompositeForceView<double> view);
}
#endif
