#
# Subdirectories
#

add_subdirectory(Dtd)
add_subdirectory(Tools)

#
#
#

# ---- Generate Version.hh ----
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/Utility/include/Version.hh.in"
  "${CMAKE_CURRENT_BINARY_DIR}/Utility/include/Version.hh"
  @ONLY
)


#
# Collect Grains3DGPU sources
#

set(GRAINS3D_GPU_INCLUDE_DIRS "")

# Base
list(APPEND GRAINS3D_GPU_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Base/include>")

# CollisionDetection
list(APPEND GRAINS3D_GPU_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CollisionDetection/include>")

# Component
list(APPEND GRAINS3D_GPU_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Component/include>")

# ContactForce
list(APPEND GRAINS3D_GPU_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/ContactForce/include>")

# Geometry
list(APPEND GRAINS3D_GPU_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Geometry/include>")

# PostProcessing
list(APPEND GRAINS3D_GPU_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/PostProcessing/include>")

# TimeIntegration
list(APPEND GRAINS3D_GPU_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/TimeIntegration/include>")

# Utility
list(APPEND GRAINS3D_GPU_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Utility/include>")

# Backend-neutral XML handles. Native libxml2 headers are private to XML/src.
list(APPEND GRAINS3D_GPU_INCLUDE_DIRS "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/XML/include>")



#
# Build Grains3DGPU libraries
#

set(CMAKE_CXX_VISIBILITY_PRESET default)
set(CMAKE_VISIBILITY_INLINES_HIDDEN  OFF)

add_library(Grains3DGPU_obj OBJECT)
add_subdirectory(Base)
add_subdirectory(CollisionDetection)
add_subdirectory(Component)
add_subdirectory(ContactForce)
add_subdirectory(Geometry)
add_subdirectory(PostProcessing)
add_subdirectory(TimeIntegration)
add_subdirectory(Utility)
add_subdirectory(XML)
add_library(PacIFiC::Grains3DGPU_obj ALIAS Grains3DGPU_obj)
add_library(Grains3DGPU SHARED $<TARGET_OBJECTS:Grains3DGPU_obj>)
add_library(PacIFiC::Grains3DGPU ALIAS Grains3DGPU)

# Shared usage requirements for both the reusable objects and the installed library.
foreach(grains_target Grains3DGPU_obj Grains3DGPU)
  target_include_directories(${grains_target} PUBLIC
    ${GRAINS3D_GPU_INCLUDE_DIRS}
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
  target_link_libraries(${grains_target} PUBLIC ZLIB::ZLIB)
  target_link_libraries(${grains_target} PRIVATE LibXml2::LibXml2)
  if(PACIFIC_USE_MPI)
    target_link_libraries(${grains_target} PUBLIC MPI::MPI_C)
  endif()
  if(PACIFIC_USE_OPENMP)
    target_link_libraries(${grains_target} PUBLIC OpenMP::OpenMP_CXX)
  endif()
  if(PACIFIC_USE_CUDA)
    set_target_properties(${grains_target} PROPERTIES
      CUDA_STANDARD ${GRAINS3D_GPU_CUDA_STANDARD}
      CUDA_STANDARD_REQUIRED ${GRAINS3D_GPU_CUDA_STANDARD_REQUIRED}
      CUDA_ARCHITECTURES ${GRAINS3D_GPU_CUDA_ARCHITECTURES}
      CUDA_SEPARABLE_COMPILATION ${GRAINS3D_GPU_CUDA_SEPARABLE_COMPILATION}
      CUDA_RESOLVE_DEVICE_SYMBOLS ${GRAINS3D_GPU_CUDA_RESOLVE_DEVICE_SYMBOLS}
    )
    target_link_libraries(${grains_target} PUBLIC CUDA::cudart CUDA::curand)
    foreach(cuda_opt IN LISTS GRAINS3D_GPU_CUDA_COMPILE_OPTIONS_LIST)
      target_compile_options(${grains_target} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:${cuda_opt}>)
    endforeach()
  endif()
endforeach()

set_target_properties(Grains3DGPU_obj PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CXX_VISIBILITY_PRESET default
)
set_target_properties(Grains3DGPU PROPERTIES
  OUTPUT_NAME Grains3DGPU
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  CXX_VISIBILITY_PRESET default
)
