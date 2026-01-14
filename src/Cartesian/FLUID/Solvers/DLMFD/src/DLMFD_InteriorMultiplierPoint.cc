#include <DLMFD_InteriorMultiplierPoint.hh>
using namespace std;

//---------------------------------------------------------------------------
DLMFD_InteriorMultiplierPoint::DLMFD_InteriorMultiplierPoint(
    const size_t &comp, const geomVector &position, size_t i, size_t j,
    size_t k, const geomVector &gravity_center)
    : DLMFD_ParticlePoint(position, gravity_center)
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_InteriorMultiplierPoint:: DLMFD_InteriorMultiplierPoint");

    localNodeTriplet = new FV_TRIPLET;
    localNodeTriplet->i = i;
    localNodeTriplet->j = j;
    localNodeTriplet->k = k;
}




//---------------------------------------------------------------------------
DLMFD_InteriorMultiplierPoint::~DLMFD_InteriorMultiplierPoint()
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_InteriorMultiplierPoint:: ~DLMFD_InteriorMultiplierPoint");
}




//---------------------------------------------------------------------------
void DLMFD_InteriorMultiplierPoint::set(const size_t &comp,
                                        const geomVector &point, size_t i,
                                        size_t j, size_t k,
                                        const geomVector &gravity_center)
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_InteriorMultiplierPoint:: set");

    DLMFD_ParticlePoint::set(comp, point, gravity_center);
    localNodeTriplet->i = i;
    localNodeTriplet->j = j;
    localNodeTriplet->k = k;
}




//---------------------------------------------------------------------------
FV_TRIPLET *DLMFD_InteriorMultiplierPoint::get_localNodeTriplet()
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_InteriorMultiplierPoint:: get_localNodeTriplet");

    return localNodeTriplet;
}
