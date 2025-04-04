#include <DLMFD_XYZBoxPlan.hh>
#include <geomVector.hh>

/* Constructor with arguments
-----------------------------*/
DLMFD_XYZBoxPlan::DLMFD_XYZBoxPlan(const string &type_, const size_t &coordinate_direction_,
                       const double &value_) : DLMFD_GeomBoundary(type_)
{
    MAC_LABEL("DLMFD_XYZBoxPlan::DLMFD_XYZBoxPlan");

    coordinate_direction = coordinate_direction_;
    value = value_;
}

/* Destructor
-------------*/
DLMFD_XYZBoxPlan::~DLMFD_XYZBoxPlan()
{
    MAC_LABEL("DLMFD_XYZBoxPlan::~DLMFD_XYZBoxPlan");
}

/* Search the intersection normal
---------------------------------*/
pair<bool, geomVector> DLMFD_XYZBoxPlan::intersection_normal(
    const geomVector &line_point) const
{
    MAC_LABEL("DLMFD_XYZBoxPlan::intersection_normal");

    pair<bool, geomVector> intersection(true, line_point);

    switch (coordinate_direction)
    {
    case 0:
        intersection.second(0) = value;
        break;
    case 1:
        intersection.second(1) = value;
        break;
    case 2:
        intersection.second(2) = value;
        break;
    }

    return intersection;
}

/* Operator <<
--------------*/
ostream &operator<<(ostream &f, const DLMFD_XYZBoxPlan &G)
{
    MAC_LABEL("DLMFD_XYZBoxPlan::operator <<");

    f << (DLMFD_GeomBoundary &)G << endl;
    f << "Normal coordinate direction = " << G.coordinate_direction << endl;
    if (G.coordinate_direction == 0)
        f << "X = " << G.value;
    else if (G.coordinate_direction == 1)
        f << "Y = " << G.value;
    else
        f << "Z = " << G.value;

    return f;
}

/* Display
----------*/
void DLMFD_XYZBoxPlan::display(ostream &f) const
{
    MAC_LABEL("DLMFD_XYZBoxPlan::display");

    f << "Type = " << geomtype << endl;
    f << "Normal coordinate direction = " << coordinate_direction << endl;
    if (coordinate_direction == 0)
        f << "X = " << value;
    else if (coordinate_direction == 1)
        f << "Y = " << value;
    else
        f << "Z = " << value;
}

/* Translate the geometric boundary
-----------------------------------*/
void DLMFD_XYZBoxPlan::translate(const geomVector &translation_vector,
                           const size_t &translation_direction)
{
    MAC_LABEL("DLMFD_XYZBoxPlan::translate");

    if (coordinate_direction == translation_direction)
        value += translation_vector(translation_direction);
}
