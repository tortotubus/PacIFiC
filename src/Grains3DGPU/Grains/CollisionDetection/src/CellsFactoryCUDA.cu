#include "CellsFactoryCUDA.cuh"

#include "Cells.hh"
#include "GrainsCUDAUtils.hh"
#include "Vector3.hh"

namespace
{
    template <typename T, CellOrdering OrderingScheme>
    __GLOBAL__ void constructCellKernel(Cells<T, OrderingScheme>** cells,
                                        uint                       index,
                                        T                          minX,
                                        T                          minY,
                                        T                          minZ,
                                        T                          maxX,
                                        T                          maxY,
                                        T                          maxZ,
                                        T                          cellSize,
                                        uint*                      numCells)
    {
        if(blockIdx.x != 0 || threadIdx.x != 0)
            return;

        cells[index] = new Cells<T, OrderingScheme>(Vector3<T>(minX, minY, minZ),
                                                    Vector3<T>(maxX, maxY, maxZ),
                                                    cellSize);
        *numCells    = cells[index]->getNumCells();
    }
}  // namespace

namespace grains::cuda
{
    template <typename T, CellOrdering OrderingScheme>
    unsigned int constructCellOnDevice(Cells<T, OrderingScheme>** cells,
                                       unsigned int               index,
                                       T                          minX,
                                       T                          minY,
                                       T                          minZ,
                                       T                          maxX,
                                       T                          maxY,
                                       T                          maxZ,
                                       T                          cellSize)
    {
        uint* deviceNumCells;
        cudaErrCheck(cudaMalloc(&deviceNumCells, sizeof(uint)));

        constructCellKernel<<<1, 1>>>(cells,
                                      index,
                                      minX,
                                      minY,
                                      minZ,
                                      maxX,
                                      maxY,
                                      maxZ,
                                      cellSize,
                                      deviceNumCells);
        cudaErrCheck(cudaGetLastError());

        uint numCells = 0;
        cudaErrCheck(cudaMemcpy(&numCells, deviceNumCells, sizeof(uint), cudaMemcpyDeviceToHost));
        cudaErrCheck(cudaFree(deviceNumCells));
        return numCells;
    }

    template unsigned int
        constructCellOnDevice<float, CellOrdering::LINEAR>(Cells<float, CellOrdering::LINEAR>**,
                                                           unsigned int,
                                                           float,
                                                           float,
                                                           float,
                                                           float,
                                                           float,
                                                           float,
                                                           float);
    template unsigned int
        constructCellOnDevice<double, CellOrdering::LINEAR>(Cells<double, CellOrdering::LINEAR>**,
                                                            unsigned int,
                                                            double,
                                                            double,
                                                            double,
                                                            double,
                                                            double,
                                                            double,
                                                            double);
    template unsigned int
        constructCellOnDevice<float, CellOrdering::MORTON>(Cells<float, CellOrdering::MORTON>**,
                                                           unsigned int,
                                                           float,
                                                           float,
                                                           float,
                                                           float,
                                                           float,
                                                           float,
                                                           float);
    template unsigned int
        constructCellOnDevice<double, CellOrdering::MORTON>(Cells<double, CellOrdering::MORTON>**,
                                                            unsigned int,
                                                            double,
                                                            double,
                                                            double,
                                                            double,
                                                            double,
                                                            double,
                                                            double);
}  // namespace grains::cuda
