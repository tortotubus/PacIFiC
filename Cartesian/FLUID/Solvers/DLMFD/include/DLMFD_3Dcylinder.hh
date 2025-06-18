#ifndef _DLMFD_3DCYLINDER__
#define _DLMFD_3DCYLINDER__

#include <DLMFD_RigidBody.hh>
#include <FS_3Dcylinder.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <string>
using namespace std;

/** @brief The class DLMFD_3Dcylinder.

A moving or stationary rigid 3D cylinder in the Fictitious Domain solver.

@author A. Wachs - Pacific project 2025 */

class DLMFD_3Dcylinder : public DLMFD_RigidBody
{
  public: //------------------------------------------------------------------
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Default constructor */
    DLMFD_3Dcylinder();

    /** @brief Constructor with arguments
    @param pgrb Pointer to the geometric rigid body class
    @param are_particles_fixed True if all particles are fixed
    @param pField_ Constrained field
    @param critical_distance_ Critical distance */
    DLMFD_3Dcylinder(FS_RigidBody *pgrb, const bool &are_particles_fixed,
                     FV_DiscreteField *pField_,
                     double const critical_distance_);

    /** @brief Destructor */
    ~DLMFD_3Dcylinder();

    //@}

    //-- Set methods
    /** @name Set methods */
    //@{

    /** @brief Set the structure of geometrical parameters of the cylinder */
    void set_ptr_FS_3Dcylinder_Additional_Param();

    /** @brief Allocate and set DLMFD boundary and interior points
    @param critical_distance Critical distance
    @param pField Pointer to constrained field */
    void set_all_MAC(FV_DiscreteField *pField, double critical_distance);

    /** @brief Set DLMFD boundary and interior points
    @param critical_distance Critical distance
    @param pField Pointer to constrained field */
    void set_all_points(FV_DiscreteField *pField, double critical_distance);

    /** @brief Set DLMFD boundary points
    @param critical_distance Critical distance
    @param pField Pointer to constrained field */
    void set_boundary_points_list(FV_DiscreteField *pField,
                                  double critical_distance);

    /** @brief Set DLMFD interior points
    @param critical_distance Critical distance */
    void set_interior_points_list(FV_DiscreteField *pField,
                                  double critical_distance);

    //@}

    //-- Update methods
    /** @name Update methods */
    //@{

    /** @brief Update the attributes in case of moving rigid bodies */
    void update();

    //@}

    //-- Get methods
    /** @name Get methods */
    //@{
    /** @brief Get the structure of geometrical parameters of the cylinder */
    FS_3Dcylinder_Additional_Param const *
    get_ptr_FS_3Dcylinder_Additional_Param();

    //@}

    //-- Geometric methods
    /** @name Geometric methods */
    //@{

    /** @brief isIn method
    @param point Point */
    bool isIn(const geomVector &point) const;

    /** @brief Translate the component from current gravity center to new
    gravity center
    @param newg new gravity center */
    void translateGeometricFeatures(geomVector const &newg);

    /** @brief Orientation of the cylinder */
    geomVector particle_orientation_vector() const;

    //@}

    //-- Add methods
    /** @name Add methods */
    //@{

    /** @brief Allocate default list of points and vectors for the cylinder
    @param critical_distance Critical distance
    @param pField Constrained field */
    void allocate_default_listOfPointsAndVectors_3Dcylinder(
        const double &critical_distance, FV_DiscreteField *pField);

    //@}

  protected: //----------------------------------------------------------------
             //-- Attributes
             /** @name Parameters */
    //@{

    FS_3Dcylinder_Additional_Param const
        *pagp; /** Pointer to FS 3Dcylinder parameters */

    geomVector BottomCenter;   /**< Center of the bottom disk */
    geomVector TopCenter;      /**< Center of the top disk */
    geomVector BottomToTopVec; /**< Axial vector from the bottom to the top disk
                                  center */
    geomVector RadialRefVec;   /**< Radial reference vector */
    double cylinder_radius;    /**< Cylinder radius */
    double cylinder_height;    /**< Cylinder height */

    //@}

  private: //-------------------------------------------------------------------
};

#endif
