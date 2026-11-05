
#
# Add Grains3D
#

add_subdirectory(Grains3D)


#
# Add Grains3DGPU
#

if (PACIFIC_BUILD_GRAINS3DGPU)
    add_subdirectory(Grains3DGPU)
endif()

#
# Add Cartesian 
#

add_subdirectory(Cartesian)

#
# Add Octree
#

add_subdirectory(Octree)