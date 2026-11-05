#ifndef GRAINS_COLLISION_DISPATCH_HH
#define GRAINS_COLLISION_DISPATCH_HH

#include "CollisionComputation.hh"
#include <type_traits>

namespace grains::detail
{
    template <typename Function>
    void dispatchGJK(NarrowPhaseOptions options, Function&& function)
    {
        using Johnson      = std::integral_constant<GJKType, GJKType::JOHNSON>;
        using SignedVolume = std::integral_constant<GJKType, GJKType::SIGNEDVOLUME>;
        if(options.algorithm == NarrowPhaseType::GJK && !options.acceleration)
            function(Johnson{}, std::false_type{});
        else if(options.algorithm == NarrowPhaseType::GJK)
            function(Johnson{}, std::true_type{});
        else if(!options.acceleration)
            function(SignedVolume{}, std::false_type{});
        else
            function(SignedVolume{}, std::true_type{});
    }
    template <typename Function>
    void dispatchCollision(NarrowPhaseOptions options, Function&& function)
    {
        auto withBounds = [&](auto bounds) {
            dispatchGJK(options, [&](auto variant, auto accelerated) {
                function(variant, accelerated, bounds);
            });
        };
        using Off      = std::integral_constant<BoundingVolumeType, BoundingVolumeType::OFF>;
        using Box      = std::integral_constant<BoundingVolumeType, BoundingVolumeType::OBB>;
        using Cylinder = std::integral_constant<BoundingVolumeType, BoundingVolumeType::OBC>;
        if(options.boundingVolume == BoundingVolumeType::OBB)
            withBounds(Box{});
        else if(options.boundingVolume == BoundingVolumeType::OBC)
            withBounds(Cylinder{});
        else
            withBounds(Off{});
    }
}
#endif
