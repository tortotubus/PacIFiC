#include <DLMFD_ParticlePoint.hh>

//---------------------------------------------------------------------------
DLMFD_ParticlePoint::DLMFD_ParticlePoint(const geomVector &point, const geomVector &gravity_center)
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::DLMFD_ParticlePoint");

    component_number = 0;
    pointCoordinates = point;
    GCPointVector = point - gravity_center;
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
void DLMFD_ParticlePoint::set(const size_t &comp, const geomVector &point, const geomVector &gravity_center)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint:: set");

    component_number = comp;
    pointCoordinates = point;
    GCPointVector = pointCoordinates - gravity_center;
}

//---------------------------------------------------------------------------
double DLMFD_ParticlePoint::get_oneCoordinate(const size_t &dir) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("ParticlePoint:: get_oneCoordinate");

    return pointCoordinates(dir);
}
