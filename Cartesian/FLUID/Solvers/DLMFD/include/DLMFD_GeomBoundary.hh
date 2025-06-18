#ifndef DLMFD_GEOMBOUNDARY_HH
#define DLMFD_GEOMBOUNDARY_HH

#include <MAC.hh>
#include <MAC_assertions.hh>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
using namespace std;

class geomVector;

/** @brief The Class DLMFD_GeomBoundary.

Use for the definition of a geometric boundary for the determination of DLMFD
points.

@author A. Wachs - Particulate flow project 2007-2009 */

class DLMFD_GeomBoundary
{
  protected:
    string geomtype; /**< Type of geometric boundary */

  protected:
    /** @name Constructors & Destructor */
    //@{
    /** @brief Constructor without argument */
    DLMFD_GeomBoundary() {}

    /** @brief Copy constructor */
    DLMFD_GeomBoundary(const DLMFD_GeomBoundary &M);
    //@}

  public:
    /** @name Constructors & Destructor */
    //@{
    /** @brief Constructor with arguments
    @param geomtype_ type of geometric boundary */
    DLMFD_GeomBoundary(const string &geomtype_);

    /** @brief Destructor */
    virtual ~DLMFD_GeomBoundary();
    //@}

    /** @brief Operator << */
    friend ostream &operator<<(ostream &f, const DLMFD_GeomBoundary &G);

    /** @brief Display */
    virtual void display(ostream &f) const;

    /** @brief Search the intersection normal
    @param line_point point on the normal */
    virtual pair<bool, geomVector>
    intersection_normal(const geomVector &line_point) const = 0;

    /** @brief Translate geometric boundaries
    @param translation_vector translation vector
    @param translation_direction translation direction */
    virtual void translate(const geomVector &translation_vector,
                           const size_t &translation_direction) = 0;
};

#endif
