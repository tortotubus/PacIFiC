#include <DLMFD_ParticlePoint.hh>

//---------------------------------------------------------------------------
DLMFD_ParticlePoint::DLMFD_ParticlePoint(const geomVector &point,
                                         const geomVector &gravity_center)
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::DLMFD_ParticlePoint");

    component_number = 0;
    pointCoordinates = point;
    GCPointVector = point - gravity_center;
    valid = true;
}




//---------------------------------------------------------------------------
DLMFD_ParticlePoint::~DLMFD_ParticlePoint()
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::DLMFD_ParticlePoint");
}




//---------------------------------------------------------------------------
geomVector DLMFD_ParticlePoint::get_coordinates() const
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::get_coordinates");

    return pointCoordinates;
}




//---------------------------------------------------------------------------
geomVector const *DLMFD_ParticlePoint::get_ptr_coordinates() const
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::get_ptr_coordinates");

    return &pointCoordinates;
}




//---------------------------------------------------------------------------
double DLMFD_ParticlePoint::get_oneCoordinate(const size_t &dir) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("ParticlePoint:: get_oneCoordinate");

    return pointCoordinates(dir);
}




//---------------------------------------------------------------------------
geomVector DLMFD_ParticlePoint::get_GCPointVector() const
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::get_GCPointVector");

    return GCPointVector;
}




//---------------------------------------------------------------------------
double
DLMFD_ParticlePoint::get_oneCoordinate_GCPointVector(const size_t &dir) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("ParticlePoint:: get_oneCoordinate_GCPointVector");

    return GCPointVector(dir);
}




//---------------------------------------------------------------------------
size_t DLMFD_ParticlePoint::get_compNumber() const
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::get_compNumber");

    return component_number;
}




//---------------------------------------------------------------------------
bool DLMFD_ParticlePoint::isValid() const
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::isValid");

    return valid;
}




//---------------------------------------------------------------------------
void DLMFD_ParticlePoint::set(const size_t &comp, const geomVector &point,
                              const geomVector &gravity_center)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint:: set");

    component_number = comp;
    pointCoordinates = point;
    GCPointVector = pointCoordinates - gravity_center;
    valid = true;
}




//---------------------------------------------------------------------------
void DLMFD_ParticlePoint::set_validity(bool const &valid_)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::set_validity");

    valid = valid_;
}
