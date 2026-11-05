#include "CollisionShapeDataCPU.hh"
#include "CollisionShapeData.hh"

namespace grains::cpu
{
    namespace
    {
        template <typename T>
        void build(ShapeTableView<T> view, const uint* representatives)
        {
            for(uint i = 0; i < view.count; ++i)
            {
                view.fillShape(i, representatives);
                view.fillBounds(i, representatives);
            }
        }
    }
    void buildShapeTables(ShapeTableView<float> view, const uint* representatives)
    {
        build(view, representatives);
    }
    void buildShapeTables(ShapeTableView<double> view, const uint* representatives)
    {
        build(view, representatives);
    }
}
