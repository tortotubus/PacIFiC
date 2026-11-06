target_sources(Grains3DGPU_obj PRIVATE
  src/BoundingBox.cpp
  src/Box.cpp
  src/Cone.cpp
  src/ContactInfo.cpp
  src/Convex.cpp
  src/ConvexFactory.cpp
  src/Cylinder.cpp
  src/Matrix3.cpp
  src/Quaternion.cpp
  src/Rectangle.cpp
  src/RigidBody.cpp
  src/Sphere.cpp
  src/Superquadric.cpp
  src/Transform3.cpp
  src/Vector3.cpp
)

if(PACIFIC_USE_CUDA)
  target_sources(Grains3DGPU_obj PRIVATE
    src/BoundingBoxCUDA.cu
    src/BoxCUDA.cu
    src/ConeCUDA.cu
    src/ContactInfoCUDA.cu
    src/ConvexCUDA.cu
    src/CylinderCUDA.cu
    src/Matrix3CUDA.cu
    src/QuaternionCUDA.cu
    src/RectangleCUDA.cu
    src/RigidBodyCUDA.cu
    src/SphereCUDA.cu
    src/SuperquadricCUDA.cu
    src/Transform3CUDA.cu
  )
endif()
