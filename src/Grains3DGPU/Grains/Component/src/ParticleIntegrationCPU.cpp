#include "ParticleIntegrationCPU.hh"
#include "ParticleIntegration.hh"

namespace
{
    template <typename T>
    void move(ParticleIntegrationView<T> view)
    {
        for(uint i = view.begin; i < view.end; ++i)
            view.move(i);
    }

    template <typename T>
    void advance(ParticleVelocityView<T> view)
    {
        for(uint i = view.begin; i < view.end; ++i)
            view.advanceVelocity(i);
    }

    template <typename T>
    void update(SubBodyTransformView<T> view)
    {
        for(uint i = view.begin; i < view.end; ++i)
            view.update(i, T(1e-12));
    }
}

namespace grains::cpu
{
    void moveParticles(ParticleIntegrationView<float> view)
    {
        move(view);
    }
    void moveParticles(ParticleIntegrationView<double> view)
    {
        move(view);
    }
    void advanceVelocity(ParticleVelocityView<float> view)
    {
        advance(view);
    }
    void advanceVelocity(ParticleVelocityView<double> view)
    {
        advance(view);
    }
    void updateSubBodyPositions(SubBodyTransformView<float> view)
    {
        update(view);
    }
    void updateSubBodyPositions(SubBodyTransformView<double> view)
    {
        update(view);
    }
}
