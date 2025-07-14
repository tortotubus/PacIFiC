add_executable(rotationmatrix
  src/rotationmatrix.cpp
)

target_link_libraries(rotationmatrix
  PRIVATE
    grains_headers
    grains
)