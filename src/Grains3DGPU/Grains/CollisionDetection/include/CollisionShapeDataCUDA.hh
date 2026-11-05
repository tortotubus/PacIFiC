#ifndef GRAINS_COLLISION_SHAPE_DATA_CUDA_HH
#define GRAINS_COLLISION_SHAPE_DATA_CUDA_HH

template <typename T>
struct ShapeTableView;

namespace grains::cuda
{
    // View arrays are device pointers; representatives is a HOST array of view.count
    // entries, uploaded privately. Both fill kernels complete before returning.
    // Empty tables do not access arrays or launch kernels.
    void buildShapeTables(ShapeTableView<float> view, const unsigned int* representatives);
    void buildShapeTables(ShapeTableView<double> view, const unsigned int* representatives);
}
#endif
