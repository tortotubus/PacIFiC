#ifndef DLMFD_INTERIORMULTIPLIERPOINT_HH
#define DLMFD_INTERIORMULTIPLIERPOINT_HH

#include <DLMFD_ParticlePoint.hh>
#include <FV_DiscreteField.hh>
#include <geomVector.hh>
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
    DLMFD_InteriorMultiplierPoint(const size_t &comp,
                                  const geomVector &position, size_t i,
                                  size_t j, size_t k,
                                  const geomVector &gravity_center);

    /** @brief Destructor */
    ~DLMFD_InteriorMultiplierPoint();

    //@}

    //-- Set methods
    /** @name Set methods */
    //@{

    /** @brief Set
    @param comp Component number
    @param point Point to set
    @param i x index of FV triplet
    @param j y index of FV triplet
    @param k z index of FV triplet
    @param gravity_center Gravity center */
    void set(const size_t &comp, const geomVector &point, size_t i, size_t j,
             size_t k, const geomVector &gravity_center);

    //@}

    //-- Get methods
    /** @name Get methods */
    //@{

    /** @brief Get local Node Triplet */
    FV_TRIPLET *get_localNodeTriplet();

    //@}

  private:   //----------------------------------------------------------------
  protected: //--------------------------------------------------------------
    FV_TRIPLET *localNodeTriplet; /**< local node MAC triplet */
};

#endif