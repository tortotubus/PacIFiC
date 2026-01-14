#ifndef DLMFD_PARTICLEPOINT_HH
#define DLMFD_PARTICLEPOINT_HH

#include <geomVector.hh>
using namespace std;

/** @brief The Class DLMFD_ParticlePoint.

Use for the definition of a multiplier point of a Rigid Body.

@author M. Houlette - Particulate flow project 2025 */

class DLMFD_ParticlePoint
{
  protected: //----------------------------------------------------------------
    geomVector pointCoordinates; /**< coordinates of the point */
    geomVector GCPointVector;    /**< vector from gravity center to the point */
    size_t component_number;     /**< the component number for MAC FVM */
    bool valid;                  /**< validity of the point */

  public: //----------------------------------------------------------------
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Constructor with arguments
    @param point Point position
    @param gravity_center Gravity center of the particle the point belongs to */
    DLMFD_ParticlePoint(const geomVector &point,
                        const geomVector &gravity_center);

    /** @brief Destructor */
    ~DLMFD_ParticlePoint();

    //@}

    //-- Set methods
    /** @name Set methods */
    //@{

    /** @brief Set
    @param comp Component number
    @param point Point to set
    @param gravity_center Gravity center */
    void set(const size_t &comp, const geomVector &point,
             const geomVector &gravity_center);
    /** @brief Set validity */
    void set_validity(bool const &valid_);
    //@}

    //-- Get methods
    /** @name Get methods */
    //@{

    /** @brief Get the coordinates of the point */
    geomVector get_coordinates() const;

    /** @brief Get the pointer to coordinates of the point */
    geomVector const *get_ptr_coordinates() const;

    /** @brief Get the coordinates of the point
    @param dir Direction */
    double get_oneCoordinate(const size_t &dir) const;

    /** @brief Get vector GCPoint */
    geomVector get_GCPointVector() const;

    /** @brief Get the coordinates of the GCPoint vector
    @param dir Direction */
    double get_oneCoordinate_GCPointVector(const size_t &dir) const;

    /** @brief Get the component of the point */
    size_t get_compNumber() const;

    /** @brief Is the point valid ? */
    bool isValid() const;

    //@}

  private: //----------------------------------------------------------------
};

#endif