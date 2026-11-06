#include "RigidBodyFactoryCUDA.cuh"

#include "Box.hh"
#include "Cone.hh"
#include "Convex.hh"
#include "Cylinder.hh"
#include "GrainsCUDAError.hh"
#include "Rectangle.hh"
#include "Sphere.hh"
#include "Superquadric.hh"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace
{
    template <typename T, typename... Arguments>
    __GLOBAL__ void createRigidBodiesBatchKernel(RigidBody<T>** bodies,
                                                 const uint*    indices,
                                                 T              crustThickness,
                                                 T              density,
                                                 uint           material,
                                                 ConvexType     convexType,
                                                 uint           count,
                                                 Arguments... args)
    {
        const uint idx = blockIdx.x * blockDim.x + threadIdx.x;
        if(idx >= count)
            return;

        Convex<T>* convex = nullptr;

        if constexpr(sizeof...(args) == 1)
        {
            if(convexType == SPHERE)
                convex = new Sphere<T>(args...);
        }
        else if constexpr(sizeof...(args) == 2)
        {
            if(convexType == CYLINDER)
                convex = new Cylinder<T>(args...);
            else if(convexType == CONE)
                convex = new Cone<T>(args...);
            else if(convexType == RECTANGLE)
                convex = new Rectangle<T>(args...);
        }
        else if constexpr(sizeof...(args) == 3)
        {
            if(convexType == BOX)
                convex = new Box<T>(args...);
        }
        else if constexpr(sizeof...(args) == 5)
        {
            if(convexType == SUPERQUADRIC)
                convex = new Superquadric<T>(args...);
        }

        if(convex != nullptr)
        {
            const uint targetIdx = indices[idx];
            bodies[targetIdx]    = new RigidBody<T>(convex, crustThickness, density, material);
        }
    }

    template <typename T>
    __GLOBAL__ void deleteRigidBodiesKernel(RigidBody<T>** bodies, uint count)
    {
        const uint idx = blockIdx.x * blockDim.x + threadIdx.x;
        if(idx < count)
        {
            delete bodies[idx];
            bodies[idx] = nullptr;
        }
    }
}  // namespace

namespace grains::cuda
{
    template <typename T>
    void copyRigidBodiesToDevice(GrainsMemBuffer<RigidBody<T>*, MemType::HOST>&   hostBodies,
                                 GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>& deviceBodies)
    {
        deviceBodies.initialize(hostBodies.getSize());

        const uint numBodies = static_cast<uint>(hostBodies.getSize());
        if(numBodies == 0)
            return;

        struct RigidBodyProperties
        {
            ConvexType type;
            T          density;
            T          crustThickness;
            uint       material;
            T          params[5];
            uint       numParams;

            bool operator==(const RigidBodyProperties& other) const
            {
                if(type != other.type || numParams != other.numParams)
                    return false;
                for(uint i = 0; i < numParams; ++i)
                    if(std::abs(params[i] - other.params[i]) > T(1e-9))
                        return false;
                return true;
            }
        };

        std::vector<RigidBodyProperties> uniqueTypes;
        std::vector<std::vector<uint>>   indicesPerType;

        for(uint i = 0; i < numBodies; ++i)
        {
            RigidBodyProperties properties;
            Convex<T>*          convex = hostBodies[i]->getConvex();
            properties.type            = convex->getConvexType();

            properties.numParams = 0;
            if(properties.type == SPHERE)
            {
                auto* sphere         = dynamic_cast<Sphere<T>*>(convex);
                properties.params[0] = sphere->getRadius();
                properties.numParams = 1;
            }
            else if(properties.type == BOX)
            {
                auto*      box       = dynamic_cast<Box<T>*>(convex);
                Vector3<T> extent    = box->getExtent();
                properties.params[0] = extent[X];
                properties.params[1] = extent[Y];
                properties.params[2] = extent[Z];
                properties.numParams = 3;
            }
            else if(properties.type == CYLINDER)
            {
                auto* cylinder       = dynamic_cast<Cylinder<T>*>(convex);
                properties.params[0] = cylinder->getRadius();
                properties.params[1] = cylinder->getHeight();
                properties.numParams = 2;
            }
            else if(properties.type == CONE)
            {
                auto* cone           = dynamic_cast<Cone<T>*>(convex);
                properties.params[0] = cone->getRadius();
                properties.params[1] = cone->getHeight();
                properties.numParams = 2;
            }
            else if(properties.type == RECTANGLE)
            {
                auto*      rectangle = dynamic_cast<Rectangle<T>*>(convex);
                Vector3<T> extent    = rectangle->getExtent();
                properties.params[0] = T(2) * extent[X];
                properties.params[1] = T(2) * extent[Y];
                properties.numParams = 2;
            }
            else if(properties.type == SUPERQUADRIC)
            {
                auto*      superquadric = dynamic_cast<Superquadric<T>*>(convex);
                Vector3<T> extent       = superquadric->getExtent();
                Vector3<T> exponent     = superquadric->getExponent();
                properties.params[0]    = extent[X];
                properties.params[1]    = extent[Y];
                properties.params[2]    = extent[Z];
                properties.params[3]    = exponent[X];
                properties.params[4]    = exponent[Y];
                properties.numParams    = 5;
            }

            int typeIndex = -1;
            for(uint j = 0; j < uniqueTypes.size(); ++j)
            {
                if(uniqueTypes[j] == properties)
                {
                    typeIndex = static_cast<int>(j);
                    break;
                }
            }

            const auto snapshot       = hostBodies[i]->getPropertiesSnapshot();
            properties.density        = snapshot.mass / snapshot.convex->computeVolume();
            properties.crustThickness = snapshot.crustThickness;
            properties.material       = snapshot.material;

            if(typeIndex >= 0)
            {
                indicesPerType[typeIndex].push_back(i);
            }
            else
            {
                uniqueTypes.push_back(properties);
                indicesPerType.push_back({i});
            }
        }

        const size_t requiredHeap = numBodies * 1024;
        size_t       currentHeap;
        cudaDeviceGetLimit(&currentHeap, cudaLimitMallocHeapSize);
        if(requiredHeap > currentHeap)
            cudaDeviceSetLimit(cudaLimitMallocHeapSize, requiredHeap);

        for(uint typeIndex = 0; typeIndex < uniqueTypes.size(); ++typeIndex)
        {
            const RigidBodyProperties& properties = uniqueTypes[typeIndex];
            const std::vector<uint>&   indices    = indicesPerType[typeIndex];
            const uint                 count      = static_cast<uint>(indices.size());

            uint* deviceIndices;
            cudaErrCheck(cudaMalloc(&deviceIndices, count * sizeof(uint)));
            cudaErrCheck(cudaMemcpy(deviceIndices,
                                    indices.data(),
                                    count * sizeof(uint),
                                    cudaMemcpyHostToDevice));

            constexpr uint maxBatchSize = 2048;
            const uint     numBatches   = (count + maxBatchSize - 1) / maxBatchSize;

            for(uint batch = 0; batch < numBatches; ++batch)
            {
                const uint batchStart = batch * maxBatchSize;
                const uint batchSize  = std::min(maxBatchSize, count - batchStart);
                const uint numThreads = 256;
                const uint numBlocks  = (batchSize + numThreads - 1) / numThreads;

                if(properties.type == SPHERE)
                {
                    createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(
                        deviceBodies.getData(),
                        deviceIndices + batchStart,
                        properties.crustThickness,
                        properties.density,
                        properties.material,
                        SPHERE,
                        batchSize,
                        properties.params[0]);
                }
                else if(properties.type == BOX)
                {
                    createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(
                        deviceBodies.getData(),
                        deviceIndices + batchStart,
                        properties.crustThickness,
                        properties.density,
                        properties.material,
                        BOX,
                        batchSize,
                        properties.params[0],
                        properties.params[1],
                        properties.params[2]);
                }
                else if(properties.type == CYLINDER)
                {
                    createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(
                        deviceBodies.getData(),
                        deviceIndices + batchStart,
                        properties.crustThickness,
                        properties.density,
                        properties.material,
                        CYLINDER,
                        batchSize,
                        properties.params[0],
                        properties.params[1]);
                }
                else if(properties.type == CONE)
                {
                    createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(
                        deviceBodies.getData(),
                        deviceIndices + batchStart,
                        properties.crustThickness,
                        properties.density,
                        properties.material,
                        CONE,
                        batchSize,
                        properties.params[0],
                        properties.params[1]);
                }
                else if(properties.type == RECTANGLE)
                {
                    createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(
                        deviceBodies.getData(),
                        deviceIndices + batchStart,
                        properties.crustThickness,
                        properties.density,
                        properties.material,
                        RECTANGLE,
                        batchSize,
                        properties.params[0],
                        properties.params[1]);
                }
                else if(properties.type == SUPERQUADRIC)
                {
                    createRigidBodiesBatchKernel<<<numBlocks, numThreads>>>(
                        deviceBodies.getData(),
                        deviceIndices + batchStart,
                        properties.crustThickness,
                        properties.density,
                        properties.material,
                        SUPERQUADRIC,
                        batchSize,
                        properties.params[0],
                        properties.params[1],
                        properties.params[2],
                        properties.params[3],
                        properties.params[4]);
                }
                else
                {
                    cudaFree(deviceIndices);
                    GAbort("Convex type is not implemented for GPU! Aborting Grains!");
                }

                cudaDeviceSynchronize();
                const cudaError_t error = cudaGetLastError();
                if(error != cudaSuccess)
                {
                    cudaFree(deviceIndices);
                    GAbort(("CUDA error in RigidBody batch creation: "
                            + std::string(cudaGetErrorString(error)))
                               .c_str());
                }
            }

            cudaFree(deviceIndices);
        }

        cudaErrCheck(cudaDeviceSynchronize());
    }

    template <typename T>
    void freeRigidBodiesOnDevice(GrainsMemBuffer<RigidBody<T>*, MemType::DEVICE>& deviceBodies)
    {
        const uint count = static_cast<uint>(deviceBodies.getSize());
        if(count == 0)
            return;

        constexpr uint numThreads = 256;
        const uint     numBlocks  = (count + numThreads - 1) / numThreads;
        deleteRigidBodiesKernel<<<numBlocks, numThreads>>>(deviceBodies.getData(), count);
        cudaErrCheck(cudaDeviceSynchronize());
    }

    template void copyRigidBodiesToDevice(GrainsMemBuffer<RigidBody<float>*, MemType::HOST>&,
                                          GrainsMemBuffer<RigidBody<float>*, MemType::DEVICE>&);
    template void copyRigidBodiesToDevice(GrainsMemBuffer<RigidBody<double>*, MemType::HOST>&,
                                          GrainsMemBuffer<RigidBody<double>*, MemType::DEVICE>&);
    template void freeRigidBodiesOnDevice(GrainsMemBuffer<RigidBody<float>*, MemType::DEVICE>&);
    template void freeRigidBodiesOnDevice(GrainsMemBuffer<RigidBody<double>*, MemType::DEVICE>&);
}  // namespace grains::cuda
