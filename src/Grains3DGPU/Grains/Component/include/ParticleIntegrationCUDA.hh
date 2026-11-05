#ifndef GRAINS_PARTICLE_INTEGRATION_CUDA_HH
#define GRAINS_PARTICLE_INTEGRATION_CUDA_HH

// This interface needs neither CUDA headers nor the full geometry definitions.
template <typename T>
struct ParticleIntegrationView;
template <typename T>
struct ParticleVelocityView;
template <typename T>
struct SubBodyTransformView;

namespace grains::cuda
{
    void moveParticles(ParticleIntegrationView<float> view);
    void moveParticles(ParticleIntegrationView<double> view);
    void advanceVelocity(ParticleVelocityView<float> view);
    void advanceVelocity(ParticleVelocityView<double> view);
    void updateSubBodyPositions(SubBodyTransformView<float> view);
    void updateSubBodyPositions(SubBodyTransformView<double> view);
}

#endif
