add_executable(Grains_rotationmatrix
  src/rotationmatrix.cpp
)

set_target_properties(Grains_rotationmatrix PROPERTIES OUTPUT_NAME rotationmatrix)

target_link_libraries(Grains_rotationmatrix
  PRIVATE
    Grains
)

