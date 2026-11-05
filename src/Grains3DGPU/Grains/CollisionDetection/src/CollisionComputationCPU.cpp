#include "CollisionComputationCPU.hh"
#include "CollisionComputation.hh"
#include "CollisionDispatch.hh"

namespace
{
    template <typename T>
    void relative(RelativeTransformView<T> view)
    {
        for(uint i = 0; i < view.count; ++i)
            view.compute(i);
    }
    template <typename View>
    void detect(View view, NarrowPhaseOptions options)
    {
        grains::detail::dispatchCollision(
            options,
            [&](auto variant, auto accelerated, auto bounds) {
                for(uint i = 0; i < view.data.count; ++i)
                    view.template detect<decltype(variant)::value,
                                         decltype(accelerated)::value,
                                         decltype(bounds)::value>(i, options.prebuiltShapes);
            });
    }
    template <typename T>
    void transform(ContactTransformView<T> view)
    {
        for(uint i = 0; i < view.count; ++i)
            view.transform(i);
    }
}
namespace grains::cpu
{
    void computeRelativeTransforms(RelativeTransformView<float> view)
    {
        relative(view);
    }
    void computeRelativeTransforms(RelativeTransformView<double> view)
    {
        relative(view);
    }
    void detectRelativeCollisions(RelativeCollisionView<float> view, NarrowPhaseOptions options)
    {
        detect(view, options);
    }
    void detectRelativeCollisions(RelativeCollisionView<double> view, NarrowPhaseOptions options)
    {
        detect(view, options);
    }
    void detectGlobalCollisions(GlobalCollisionView<float> view, NarrowPhaseOptions options)
    {
        detect(view, options);
    }
    void detectGlobalCollisions(GlobalCollisionView<double> view, NarrowPhaseOptions options)
    {
        detect(view, options);
    }
    void transformContacts(ContactTransformView<float> view)
    {
        transform(view);
    }
    void transformContacts(ContactTransformView<double> view)
    {
        transform(view);
    }
}
