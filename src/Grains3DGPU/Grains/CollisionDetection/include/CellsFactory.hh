#ifndef _CELLSFACTORY_HH_
#define _CELLSFACTORY_HH_

#include "Cells.hh"
#include "CellsFactoryCUDA.cuh"
#include "GrainsMemBuffer.hh"

// =================================================================================================
/** @brief The class CellsFactory.

    Creates cells objects for the simulation. Host construction remains ordinary C++; device-heap
    construction is delegated to the CUDA backend.

    @author A.YAZDANI - 2025 - Construction */
// =================================================================================================
template <typename T, CellOrdering OrderingScheme = CellOrdering::LINEAR>
class CellsFactory
{
private:
    /**@name Constructors & Destructor */
    //@{
    /** @brief Default constructor (forbidden) */
    CellsFactory() = default;

    /** @brief Destructor (forbidden) */
    ~CellsFactory() = default;
    //@}

public:
    /**@name Methods */
    //@{
    /** @brief Creates one Cells object in the requested memory space.
        @param minCorner Minimum corner of the domain
        @param maxCorner Maximum corner of the domain
        @param cellSize Size of each cell
        @param cells Memory buffer for storing the Cells pointer */
    template <MemType M>
    static uint create(const Vector3<T>&                              minCorner,
                       const Vector3<T>&                              maxCorner,
                       T                                              cellSize,
                       GrainsMemBuffer<Cells<T, OrderingScheme>*, M>& cells)
    {
        static_assert(M == MemType::HOST || M == MemType::DEVICE,
                      "CellsFactory only supports HOST or DEVICE memory");

        cells.reserve(1);
        if constexpr(M == MemType::HOST)
            return create_host(minCorner, maxCorner, cellSize, cells);
        else
            return create_device(minCorner, maxCorner, cellSize, cells);
    }

    /** @brief Creates one Cells object on the host. */
    static uint create_host(const Vector3<T>&                                          minCorner,
                            const Vector3<T>&                                          maxCorner,
                            T                                                          cellSize,
                            GrainsMemBuffer<Cells<T, OrderingScheme>*, MemType::HOST>& cells)
    {
        GAssert(cellSize > 0, "Cell size must be positive! Aborting Grains!");
        cells.initialize(1);
        cells[0] = new Cells<T, OrderingScheme>(minCorner, maxCorner, cellSize);
        return cells[0]->getNumCells();
    }

    /** @brief Creates one Cells object in device-heap memory. */
    static uint create_device(const Vector3<T>& minCorner,
                              const Vector3<T>& maxCorner,
                              T                 cellSize,
                              GrainsMemBuffer<Cells<T, OrderingScheme>*, MemType::DEVICE>& cells)
    {
        GrainsMemBuffer<Cells<T, OrderingScheme>*, MemType::HOST> hostCells(1);
        const uint numCells = create_host(minCorner, maxCorner, cellSize, hostCells);
        copyHostToDevice(hostCells, cells);
        delete hostCells[0];
        hostCells[0] = nullptr;
        return numCells;
    }

    /** @brief Reconstructs host Cells objects in device-heap memory. */
    static void
        copyHostToDevice(GrainsMemBuffer<Cells<T, OrderingScheme>*, MemType::HOST>&   hostCells,
                         GrainsMemBuffer<Cells<T, OrderingScheme>*, MemType::DEVICE>& deviceCells)
    {
        deviceCells.initialize(hostCells.getSize());
        for(uint i = 0; i < hostCells.getSize(); ++i)
        {
            if(hostCells[i] == nullptr)
                continue;

            const Vector3<T>& minCorner = hostCells[i]->getMinCorner();
            const Vector3<T>& maxCorner = hostCells[i]->getMaxCorner();
            const T           cellSize  = hostCells[i]->getCellSize();
            GAssert(cellSize > 0, "Cell size must be positive! Aborting Grains!");

            grains::cuda::constructCellOnDevice<T, OrderingScheme>(deviceCells.getData(),
                                                                   i,
                                                                   minCorner[X],
                                                                   minCorner[Y],
                                                                   minCorner[Z],
                                                                   maxCorner[X],
                                                                   maxCorner[Y],
                                                                   maxCorner[Z],
                                                                   cellSize);
        }
    }
    //@}
};

#endif
