#ifndef DLMFD_BOUNDARYMULTIPLIERPOINT_HH
#define DLMFD_BOUNDARYMULTIPLIERPOINT_HH

#include <DLMFD_ParticlePoint.hh>
#include <geomVector.hh>
using namespace std;

/** @brief The Class DLMFD_BoundaryMultiplierPoint.

Use for the definition of a boundary multiplier point of a Rigid Body.

@author M. Houlette - Particulate flow project 2025 */

class DLMFD_BoundaryMultiplierPoint : public DLMFD_ParticlePoint
{
  public: //----------------------------------------------------------------
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Constructor with arguments
    @param position point position
    @param gravity_center gravity center of the particle the point belongs to */
    DLMFD_BoundaryMultiplierPoint(const geomVector &point,
                                  const geomVector &gravity_center);

    /** @brief Destructor */
    ~DLMFD_BoundaryMultiplierPoint();

    //@}

  private:   //----------------------------------------------------------------
  protected: //----------------------------------------------------------------
};

#endif