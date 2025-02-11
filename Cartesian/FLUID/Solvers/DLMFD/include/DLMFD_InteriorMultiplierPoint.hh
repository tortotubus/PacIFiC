#ifndef DLMFD_INTERIORMULTIPLIERPOINT_HH
#define DLMFD_INTERIORMULTIPLIERPOINT_HH

#include <geomVector.hh>
#include <DLMFD_ParticlePoint.hh>
#include <FV_DiscreteField.hh>
using namespace std;

/** @brief The Class DLMFD_InteriorMultiplierPoint.

Use for the definition of an interior multiplier point of a Rigid Body.

@author M. Houlette - Pacific project 2025 */

class DLMFD_InteriorMultiplierPoint : public DLMFD_ParticlePoint
{
public: //----------------------------------------------------------------
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Constructor with arguments
    @param position point position
    @param gravity_center gravity center of the particle the point belongs to */
    DLMFD_InteriorMultiplierPoint(const geomVector &position, const geomVector &gravity_center);

    /** @brief Destructor */
    ~DLMFD_InteriorMultiplierPoint();

    //@}

private:                          //----------------------------------------------------------------
protected:                        //--------------------------------------------------------------
    FV_TRIPLET *localNodeTriplet; /**< local node MAC triplet */
};

#endif