add_executable(Grains_main src/main.cpp)
target_link_libraries(Grains_main PRIVATE Grains)
set_target_properties(Grains_main PROPERTIES OUTPUT_NAME grains)