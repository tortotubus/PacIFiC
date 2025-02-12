#include <DLMFD_ParticlePoint.hh>

//---------------------------------------------------------------------------
DLMFD_ParticlePoint::DLMFD_ParticlePoint(const geomVector &point, const geomVector &gravity_center)
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::DLMFD_ParticlePoint");

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
