#ifndef GRAINS_COLLISION_COMPUTATION_CUDA_HH
#define GRAINS_COLLISION_COMPUTATION_CUDA_HH

// Ordinary declarations: no CUDA or geometry headers are needed by the interface.
struct NarrowPhaseOptions;
template <typename T>
struct RelativeTransformView;
template <typename T>
struct RelativeCollisionView;
template <typename T>
struct GlobalCollisionView;
template <typename T>
struct ContactTransformView;

namespace grains::cuda
{
    // Device arrays on the default stream. Preserve the previous stage completion
    // semantics: these operations synchronize before returning; empty work is a no-op.
    void computeRelativeTransforms(RelativeTransformView<float> view);
    void computeRelativeTransforms(RelativeTransformView<double> view);
    void detectRelativeCollisions(RelativeCollisionView<float> view, NarrowPhaseOptions options);
    void detectRelativeCollisions(RelativeCollisionView<double> view, NarrowPhaseOptions options);
    void detectGlobalCollisions(GlobalCollisionView<float> view, NarrowPhaseOptions options);
    void detectGlobalCollisions(GlobalCollisionView<double> view, NarrowPhaseOptions options);
    void transformContacts(ContactTransformView<float> view);
    void transformContacts(ContactTransformView<double> view);
}
#endif
