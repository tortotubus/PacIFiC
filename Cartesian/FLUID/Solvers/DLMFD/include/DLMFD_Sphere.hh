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
    DLMFD_Sphere(FS_RigidBody *pgrb);

    /** @brief Destructor */
    ~DLMFD_Sphere();

    //@}

    //-- Set methods
    /** @name Set methods */
    //@{

    /** @brief Set DLMFD boundary and interior points
    @param critical_distance Critical distance
    @param pField Pointer to constrained field */
    void set_all_points(FV_DiscreteField const *pField, double critical_distance);

    /** @brief Set DLMFD boundary points
    @param critical_distance Critical distance
    @param pField Pointer to constrained field */
    void set_boundary_points_list(FV_DiscreteField const *pField, double critical_distance);

    /** @brief Set DLMFD interior points
    @param critical_distance Critical distance */
    void set_interior_points_list(FV_DiscreteField const *pField, double critical_distance);

    //@}

    //-- Get methods
    /** @name Get methods */
    //@{

    /** @brief Returns pointer to the rigid body gravity center */
    geomVector const *get_ptr_to_gravity_centre() const;

    /** @brief Returns circumscribed radius */
    double get_circumscribed_radius() const;

    //@}

    //-- Updating methods
    /** @name Updating methods */
    //@{

    /** @brief Update the RB */
    void update();

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

    /** @brief Allocate default boundary points
    @param nbBPdef Number of default boundary points */
    void allocate_default_boundary_points_sphere(size_t const &nbBPdef);

    /** @brief Allocate default halozone boundary points
    @param nbBPHZdef Number of default halozone boundary points */
    void allocate_default_halozone_boundary_points_sphere(size_t const &nbBPHZdef);

    /** @brief Allocate default interior points
    @param nbIPdef Number of default interior points */
    void allocate_default_interior_points_sphere(size_t const &nbIPdef);

    /** @brief Allocate default halozone interior points
    @param nbIPHZdef Number of default halozone interior points */
    void allocate_default_halozone_interior_points_sphere(size_t const &nbIPHZdef);

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
