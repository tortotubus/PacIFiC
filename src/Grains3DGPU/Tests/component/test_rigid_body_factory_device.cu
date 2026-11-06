#include <cuda_runtime.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <type_traits>

#include "Box.hh"
#include "Cone.hh"
#include "Cylinder.hh"
#include "Rectangle.hh"
#include "RigidBodyFactoryCUDA.cuh"
#include "Sphere.hh"
#include "Superquadric.hh"

namespace
{
template <typename T>
struct DeviceRigidBodySnapshot
{
    int  type;
    T    crustThickness;
    T    mass;
    uint material;
    uint exists;
};

template <typename T>
__GLOBAL__ void captureRigidBodySnapshots(RigidBody<T>* const*       bodies,
                                          DeviceRigidBodySnapshot<T>* snapshots,
                                          uint                        count)
{
    const uint idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx >= count)
        return;

    RigidBody<T>* body = bodies[idx];
    if(body == nullptr)
    {
        snapshots[idx].exists = 0;
        return;
    }

    const auto properties          = body->getPropertiesSnapshot();
    snapshots[idx].type            = static_cast<int>(properties.convex->getConvexType());
    snapshots[idx].crustThickness  = properties.crustThickness;
    snapshots[idx].mass            = properties.mass;
    snapshots[idx].material        = properties.material;
    snapshots[idx].exists          = 1;
}

template <typename T>
class HostRigidBodyOwner
{
public:
    static constexpr uint bodyCount = 6;

    HostRigidBodyOwner(T crustThickness, T density, uint material)
        : bodies(bodyCount)
    {
        bodies[0] = new RigidBody<T>(new Sphere<T>(T(1.25)),
                                     crustThickness,
                                     density,
                                     material);
        bodies[1]
            = new RigidBody<T>(new Box<T>(T(1.0), T(1.5), T(2.0)),
                               crustThickness,
                               density,
                               material);
        bodies[2]
            = new RigidBody<T>(new Cylinder<T>(T(0.75), T(2.5)),
                               crustThickness,
                               density,
                               material);
        bodies[3]
            = new RigidBody<T>(new Cone<T>(T(0.8), T(2.25)),
                               crustThickness,
                               density,
                               material);
        bodies[4]
            = new RigidBody<T>(new Rectangle<T>(T(1.75), T(2.5)),
                               crustThickness,
                               density,
                               material);
        bodies[5]
            = new RigidBody<T>(new Superquadric<T>(T(1.0),
                                                    T(1.25),
                                                    T(1.5),
                                                    T(3.0),
                                                    T(4.0)),
                               crustThickness,
                               density,
                               material);
    }

    ~HostRigidBodyOwner()
    {
        bodies.freePointedObjects();
    }

    GrainsMemBuffer<RigidBody<T>*, MemType::HOST> bodies;
};

template <typename T>
class DeviceRigidBodyOwner
{
public:
    ~DeviceRigidBodyOwner()
    {
        if(ownsObjects)
            grains::cuda::freeRigidBodiesOnDevice(bodies);
    }

    void copyFrom(GrainsMemBuffer<RigidBody<T>*, MemType::HOST>& hostBodies)
    {
        grains::cuda::copyRigidBodiesToDevice(hostBodies, bodies);
        ownsObjects = true;
    }

    void releaseObjects()
    {
        grains::cuda::freeRigidBodiesOnDevice(bodies);
        ownsObjects = false;
    }

    GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE> bodies;

private:
    bool ownsObjects = false;
};

template <typename T>
class RigidBodyFactoryDeviceTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        int               deviceCount = 0;
        const cudaError_t error       = cudaGetDeviceCount(&deviceCount);
        if(error != cudaSuccess || deviceCount == 0)
            GTEST_SKIP() << "No CUDA device available";
    }
};

using RealTypes = ::testing::Types<float, double>;
TYPED_TEST_SUITE(RigidBodyFactoryDeviceTest, RealTypes);

TYPED_TEST(RigidBodyFactoryDeviceTest, CopiesAllSupportedShapesAndDeletesThem)
{
    using T = TypeParam;

    constexpr uint count          = HostRigidBodyOwner<T>::bodyCount;
    const T        crustThickness = T(0.025);
    const T        density        = T(2.75);
    const uint     material       = 7;

    HostRigidBodyOwner<T>   hostBodies(crustThickness, density, material);
    DeviceRigidBodyOwner<T> deviceBodies;
    deviceBodies.copyFrom(hostBodies.bodies);

    GrainsMemBuffer<DeviceRigidBodySnapshot<T>, MemType::DEVICE> deviceSnapshots(count);
    GrainsMemBuffer<DeviceRigidBodySnapshot<T>, MemType::HOST>   hostSnapshots;

    captureRigidBodySnapshots<<<1, count>>>(deviceBodies.bodies.getData(),
                                            deviceSnapshots.getData(),
                                            count);
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
    hostSnapshots.copyFrom(deviceSnapshots);
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);

    constexpr std::array<ConvexType, count> expectedTypes{
        SPHERE, BOX, CYLINDER, CONE, RECTANGLE, SUPERQUADRIC};
    const T relativeTolerance = std::is_same_v<T, float> ? T(1e-5) : T(1e-12);

    for(uint i = 0; i < count; ++i)
    {
        const auto& actual   = hostSnapshots[i];
        const auto  expected = hostBodies.bodies[i]->getPropertiesSnapshot();
        const T     massTolerance
            = relativeTolerance * std::max(T(1), std::abs(expected.mass));

        EXPECT_EQ(actual.exists, 1u) << "body " << i;
        EXPECT_EQ(actual.type, static_cast<int>(expectedTypes[i])) << "body " << i;
        EXPECT_EQ(actual.crustThickness, expected.crustThickness) << "body " << i;
        EXPECT_NEAR(actual.mass, expected.mass, massTolerance) << "body " << i;
        EXPECT_EQ(actual.material, expected.material) << "body " << i;
    }

    deviceBodies.releaseObjects();

    GrainsMemBuffer<RigidBody<T>*, MemType::HOST> releasedPointers;
    releasedPointers.copyFrom(deviceBodies.bodies);
    ASSERT_EQ(cudaDeviceSynchronize(), cudaSuccess);
    for(uint i = 0; i < count; ++i)
        EXPECT_EQ(releasedPointers[i], nullptr) << "body " << i;
}
}  // namespace
