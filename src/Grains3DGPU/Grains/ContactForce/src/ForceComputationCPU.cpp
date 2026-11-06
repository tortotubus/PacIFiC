#include "ForceComputationCPU.hh"
#include "ForceComputation.hh"
#include "ContactTableHostOps.hh"

namespace
{
    template <typename T>
    void compute(ContactForceView<T> view)
    {
        for(uint i = 0; i < view.count; ++i)
            view.template compute<grains::cpu::ContactTableOps>(i);
    }
    template <typename T>
    void add(ExternalForceView<T> view)
    {
        for(uint i = view.begin; i < view.end; ++i)
            view.add(i);
    }
    template <typename T>
    void assemble(CompositeForceView<T> view)
    {
        for(uint i = view.begin; i < view.end; ++i)
        {
            uint       master;
            Vector3<T> force, torque;
            if(view.contribution(i, master, force, torque))
            {
                view.torces[master].addForce(force);
                view.torces[master].addTorque(torque);
                view.torces[i].reset();
            }
        }
    }
}
namespace grains::cpu
{
    void computeContactForces(ContactForceView<float> view)
    {
        compute(view);
    }
    void computeContactForces(ContactForceView<double> view)
    {
        compute(view);
    }
    void addExternalForces(ExternalForceView<float> view)
    {
        add(view);
    }
    void addExternalForces(ExternalForceView<double> view)
    {
        add(view);
    }
    void assembleCompositeForces(CompositeForceView<float> view)
    {
        assemble(view);
    }
    void assembleCompositeForces(CompositeForceView<double> view)
    {
        assemble(view);
    }
}
