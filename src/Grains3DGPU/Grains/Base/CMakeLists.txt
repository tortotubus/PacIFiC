target_sources(Grains3DGPU_obj PRIVATE
  src/Grains.cpp
  src/GrainsCPU.cpp
  src/GrainsFactory.cpp
  src/GrainsParameters.cpp
)

if(PACIFIC_USE_CUDA)
  target_sources(Grains3DGPU_obj PRIVATE
    src/GrainsFactoryCUDA.cu
    src/GrainsGPU.cu
    src/GrainsMemBufferCUDA.cu
    src/GrainsPostProcessCUDA.cu
    src/GrainsTestDev.cu
    src/GrainsTimerCUDA.cu
  )
else()
  target_sources(Grains3DGPU_obj PRIVATE
    src/GrainsFactoryNoCUDA.cpp
    src/GrainsMemBufferNoCUDA.cpp
    src/GrainsTimerNoCUDA.cpp
  )
endif()
