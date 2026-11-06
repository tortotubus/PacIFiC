target_sources(Grains3DGPU_obj PRIVATE
  src/ContactForceModel.cpp
  src/ContactForceModelFactory.cpp
  src/ContactTable.cpp
  src/ForceComputationCPU.cpp
  src/ForceModule.cpp
  src/HookeContactForceModel.cpp
  src/HookeMemoryContactForceModel.cpp
)

if(PACIFIC_USE_CUDA)
  target_sources(Grains3DGPU_obj PRIVATE
    src/ContactForceModelCUDA.cu
    src/ContactForceModelFactoryCUDA.cu
    src/ContactTableCUDA.cu
    src/ForceCompactionCUDA.cu
    src/ForceComputationCUDA.cu
    src/ForceModuleCUDA.cu
    src/HookeContactForceModelCUDA.cu
    src/HookeMemoryContactForceModelCUDA.cu
  )
endif()
