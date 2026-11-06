target_sources(Grains3DGPU_obj PRIVATE
  src/ComponentManager.cpp
  src/Insertion.cpp
  src/InsertionWindow.cpp
  src/Kinematics.cpp
  src/ParticleIntegrationCPU.cpp
  src/ParticleSortingCPU.cpp
  src/RigidBodyFactory.cpp
  src/Torce.cpp
)

if(PACIFIC_USE_CUDA)
  target_sources(Grains3DGPU_obj PRIVATE
    src/ComponentManagerCUDA.cu
    src/KinematicsCUDA.cu
    src/ParticleIntegrationCUDA.cu
    src/ParticleSortCUDA.cu
    src/ParticleSortingCUDA.cu
    src/RigidBodyFactoryCUDA.cu
  )
endif()
