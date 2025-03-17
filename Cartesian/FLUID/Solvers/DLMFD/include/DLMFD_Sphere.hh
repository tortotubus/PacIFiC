#ifndef _DLMFD_SPHERE__
#define _DLMFD_SPHERE__

#include <DLMFD_RigidBody.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <string>
using namespace std;

/** @brief The class DLMFD_Sphere.

A moving or stationary rigid sphere in the Fictitious Domain solver.

@author A. Wachs - Pacific project 2025 */

class DLMFD_Sphere : public DLMFD_RigidBody
{
public: //------------------------------------------------------------------
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Default constructor */
    DLMFD_Sphere();

    /** @brief Constructor with arguments
    @param pgrb Pointer to the geometric rigid body class */
    DLMFD_Sphere(FS_RigidBody *pgrb,
                 const bool &are_particles_fixed,
                 FV_DiscreteField *pField_,
                 double const critical_distance_);

    /** @brief Destructor */
    ~DLMFD_Sphere();

    //@}

    //-- Set methods
    /** @name Set methods */
    //@{

    /** @brief Set DLMFD boundary and interior points
    @param critical_distance Critical distance
    @param pField Pointer to constrained field */
    void set_all_points(FV_DiscreteField *pField, double critical_distance);

    /** @brief Set DLMFD boundary points
    @param critical_distance Critical distance
    @param pField Pointer to constrained field */
    void set_boundary_points_list(FV_DiscreteField *pField, double critical_distance);

    /** @brief Set DLMFD interior points
    @param critical_distance Critical distance */
    void set_interior_points_list(FV_DiscreteField *pField, double critical_distance);

    //@}

    //-- Geometric methods
    /** @name Geometric methods */
    //@{

    /** @brief isIn method
    @param x x-component
    @param y y-component
    @param z z-component */
    bool isIn(const geomVector &point) const;

    /** @brief Translate the component from current gravity center to new
    gravity center
    @param newg new gravity center */
    void translateGeometricFeatures(geomVector const &newg);

    //@}

    //-- Add methods
    /** @name Add methods */
    //@{

    void allocate_default_listOfPointsAndVectors_Sphere(const double &critical_distance, FV_DiscreteField *pField);

    //@}

    //-- Add methods
    /** @name Add methods */
    //@{

    /** @brief Remove critical interior points
    @param critical_distance Critical_distance */
    void erase_critical_interior_points_PerProc(double critical_distance);

    //@}

protected: //----------------------------------------------------------------
private:   //----------------------------------------------------------------
};

#endif
