#ifndef CELLS_FACTORY_CUDA_CUH
#define CELLS_FACTORY_CUDA_CUH

#include "CellsFwd.hh"

namespace grains::cuda
{
    /** Construct one Cells object in device-heap memory and return its total cell count. */
    template <typename T, CellOrdering OrderingScheme>
    unsigned int constructCellOnDevice(Cells<T, OrderingScheme>** cells,
                                       unsigned int               index,
                                       T                          minX,
                                       T                          minY,
                                       T                          minZ,
                                       T                          maxX,
                                       T                          maxY,
                                       T                          maxZ,
                                       T                          cellSize);
}

#endif
