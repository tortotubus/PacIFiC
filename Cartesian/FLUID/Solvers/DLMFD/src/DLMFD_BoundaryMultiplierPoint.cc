#include <DLMFD_BoundaryMultiplierPoint.hh>
using namespace std;

//---------------------------------------------------------------------------
DLMFD_BoundaryMultiplierPoint::DLMFD_BoundaryMultiplierPoint(
    const geomVector &point, const geomVector &gravity_center)
    : DLMFD_ParticlePoint(point, gravity_center)
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::DLMFD_ParticlePoint");
}




//---------------------------------------------------------------------------
DLMFD_BoundaryMultiplierPoint::~DLMFD_BoundaryMultiplierPoint()
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::DLMFD_ParticlePoint");
}
