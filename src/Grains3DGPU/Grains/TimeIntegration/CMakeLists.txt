target_sources(Grains3DGPU_obj PRIVATE
  src/FirstOrderExplicit.cpp
  src/SecondOrderLeapFrog.cpp
  src/TimeIntegrator.cpp
  src/TimeIntegratorFactory.cpp
)

if(PACIFIC_USE_CUDA)
  target_sources(Grains3DGPU_obj PRIVATE
    src/FirstOrderExplicitCUDA.cu
    src/SecondOrderLeapFrogCUDA.cu
    src/TimeIntegratorCUDA.cu
    src/TimeIntegratorFactoryCUDA.cu
  )
endif()
