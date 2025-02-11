#include <DLMFD_InteriorMultiplierPoint.hh>
using namespace std;

//---------------------------------------------------------------------------
DLMFD_InteriorMultiplierPoint::DLMFD_InteriorMultiplierPoint(const geomVector &position, const geomVector &gravity_center)
    : DLMFD_ParticlePoint(position, gravity_center)
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::DLMFD_ParticlePoint");
}

//---------------------------------------------------------------------------
DLMFD_InteriorMultiplierPoint::~DLMFD_InteriorMultiplierPoint()
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ParticlePoint::DLMFD_ParticlePoint");
}