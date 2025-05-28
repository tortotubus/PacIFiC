#ifndef _DLMFD_AllGeomBoundaries_HH
#define _DLMFD_AllGeomBoundaries_HH

#include <DLMFD_GeomBoundary.hh>
#include <list>
#include <string>
using namespace std;

class FV_DiscreteField;

/** @brief The Class DLMFD_AllGeomBoundaries.

For the definition of a list of geometric boundaries.

@author A. Wachs - Particulate flow project 2007-2009 */

class DLMFD_AllGeomBoundaries
{
  protected:
    list<DLMFD_GeomBoundary *>
        all_geomboundaries; /**< list of geometric boundaries */

  private:
    /** @name Constructors & Destructor */
    //@{
    /** @brief Constructor without argument */
    DLMFD_AllGeomBoundaries() {}

    /** @brief Copy constructor */
    DLMFD_AllGeomBoundaries(const DLMFD_AllGeomBoundaries &M) {}
    //@}

  public:
    /** @name Constructors & Destructor */
    //@{
    /**
      @brief Constructor with arguments
      @param pField pointer to the field
    */
    DLMFD_AllGeomBoundaries(FV_DiscreteField *pField);

    /** @brief Destructor */
    ~DLMFD_AllGeomBoundaries();
    //@}

    /** @name Get methods */
    //@{

    /** @brief Return a pointer to the list of geometric boundaries */
    const list<DLMFD_GeomBoundary *> *get_ptr_list_geomboundaries() const;
    
    //@}

    /** @brief Operator << */
    friend ostream &operator<<(ostream &f, const DLMFD_AllGeomBoundaries &G);

    /** @brief Display */
    void display(ostream &f) const;

    /** @brief Translate geometric boundaries
    @param translation_vector translation vector
    @param translation_direction translation direction */
    void translate(const geomVector &translation_vector,
                   const size_t &translation_direction);
};

#endif
