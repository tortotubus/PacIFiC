target_sources(Grains3DGPU_obj PRIVATE
  src/Cells.cpp
  src/CollisionComputationCPU.cpp
  src/CollisionDetectionModule.cpp
  src/CollisionShapeDataCPU.cpp
  src/GJK.cpp
  src/NeighborList_Kernels.cpp
)

if(PACIFIC_USE_CUDA)
  target_sources(Grains3DGPU_obj PRIVATE
    src/CellsCUDA.cu
    src/CellsFactoryCUDA.cu
    src/CollisionCompactionCUDA.cu
    src/CollisionComputationCUDA.cu
    src/CollisionDetectionModuleCUDA.cu
    src/CollisionDetectionModule_Kernels.cu
    src/CollisionFilteringCUDA.cu
    src/CollisionShapeDataCUDA.cu
    src/GJKCUDA.cu
    src/LinkedCellCUDA.cu
    src/LinkedCellFactoryCUDA.cu
    src/NeighborListLinkedCellCUDA.cu
    src/NeighborList_Kernels.cu
  )
endif()
