#ifndef GRAINS_COLLISION_SHAPE_DATA_CPU_HH
#define GRAINS_COLLISION_SHAPE_DATA_CPU_HH

template <typename T>
struct ShapeTableView;

namespace grains::cpu
{
    // All arrays, including representatives, are host-accessible.
    void buildShapeTables(ShapeTableView<float> view, const unsigned int* representatives);
    void buildShapeTables(ShapeTableView<double> view, const unsigned int* representatives);
}
#endif
