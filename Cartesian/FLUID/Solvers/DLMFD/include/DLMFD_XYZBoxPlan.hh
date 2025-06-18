#ifndef DLMFD_XYZBOXPLAN_HH
#define DLMFD_XYZBOXPLAN_HH

#include <DLMFD_GeomBoundary.hh>
#include <string>
using namespace std;

/** @brief The Class DLMFD_XYZBoxPlan.

Special type of geometric box: plan of a closed box

@author A. Wachs - Particulate flow project 2007-2009 */

class DLMFD_XYZBoxPlan : public DLMFD_GeomBoundary
{
  protected:
    size_t coordinate_direction; /**< coordinate direction (0,1 or 2) */
    double value;                /**< coordinate value */

  private:
    /** @name Constructors & Destructor */
    //@{
    /** @brief Constructor without argument */
    DLMFD_XYZBoxPlan() {}

    /** @brief Copy constructor */
    DLMFD_XYZBoxPlan(const DLMFD_XYZBoxPlan &M) : DLMFD_GeomBoundary(M) {}
    //@}

  public:
    /** @name Constructors & Destructor */
    //@{
    /** @brief Constructor with arguments
    @param type_ type of geometric boundary
    @param coordinate_direction_ coordinate direction
    @param value_ a value to define the geometric boundary */
    DLMFD_XYZBoxPlan(const string &type_, const size_t &coordinate_direction_,
                     const double &value_);

    /** @brief Destructor */
    ~DLMFD_XYZBoxPlan();
    //@}

    /** @brief Search the intersection normal
    @param line_point point on the normal */
    pair<bool, geomVector>
    intersection_normal(const geomVector &line_point) const;

    /** @brief Operator << */
    friend ostream &operator<<(ostream &f, const DLMFD_XYZBoxPlan &G);

    /** @brief Display */
    void display(ostream &f) const;

    /** @brief Translate the geometric boundary
    @param translation_vector translation vector
    @param translation_direction translation direction */
    void translate(const geomVector &translation_vector,
                   const size_t &translation_direction);
};

#endif
