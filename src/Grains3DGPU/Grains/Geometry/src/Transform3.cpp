#include "Transform3Impl.hh"

// -------------------------------------------------------------------------------------------------
// Constructor using an XML node
template <typename T>
__HOST__ Transform3<T>::Transform3(grains::xml::OptionalElement root)
{
    // Initialization of the 1D array containing all the transformation
    // coefficients ordered as Mxx, Mxy, Mxz, Myx, Myy, Myz, Mzx, Mzy, Mzz,
    // Ox, Oy, Oz
    T t[12] = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};

    // Origin
    const auto& element = grains::xml::requireElement(root);
    auto pos = element.getChild("Centre");
    if(pos)
    {
        T x = T(pos->requireAttribute<double>("X"));
        T y = T(pos->requireAttribute<double>("Y"));
        T z = T(pos->requireAttribute<double>("Z"));
        setOrigin(Vector3<T>(x, y, z));
    }

    // Rotation matrix
    // We only allow rotation matrices, stretching/scaling is prohibited
    auto angPos = element.getChild("AngularPosition");
    const auto& angularPosition = grains::xml::requireElement(angPos);
    string type = angularPosition.requireAttribute<std::string>("Type");
    if(type == "Matrix")
    {
        Matrix3<T>         mat;
        std::string        values = angularPosition.text();
        std::istringstream inValues(values.c_str());
        inValues >> mat;
        setBasis(mat);
        GAssert(isRotation(mat),
                "Input matrix is not a valid rotation matrix in "
                "Quaternion::setQuaternion!");
    }
    else if(type == "Angles")
    {
        // read in degree
        T aX = T(angularPosition.requireAttribute<double>("aX"));
        T aY = T(angularPosition.requireAttribute<double>("aY"));
        T aZ = T(angularPosition.requireAttribute<double>("aZ"));
        // change to Radian
        setBasis(RADS_PER_DEG<T> * aX, RADS_PER_DEG<T> * aY, RADS_PER_DEG<T> * aZ);
    }
    else
    {
        // set basis to zero rotation
        setBasis(T(0), T(0), T(0));
    }
}
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
// Output operator
template <typename T>
__HOST__ std::ostream& operator<<(std::ostream& fileOut, const Transform3<T>& t)
{
    // Orientation first, followed by the position
    fileOut << t.getBasis() << " " << t.getOrigin() << std::endl;
    return (fileOut);
}
// -------------------------------------------------------------------------------------------------
// Input operator
template <typename T>
__HOST__ std::istream& operator>>(std::istream& fileIn, Transform3<T>& t)
{
    Vector3<T> v;
    Matrix3<T> m;

    fileIn >> v;
    fileIn >> m;

    t.setOrigin(v);
    t.setBasis(m);

    return (fileIn);
}
// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class Transform3<float>;
template class Transform3<double>;

#define X(T)                                                                               \
    template std::ostream& operator<< <T>(std::ostream & fileOut, const Transform3<T>& t); \
                                                                                           \
    template std::istream& operator>> <T>(std::istream & fileIn, Transform3<T> & t);
X(float)
X(double)
#undef X
