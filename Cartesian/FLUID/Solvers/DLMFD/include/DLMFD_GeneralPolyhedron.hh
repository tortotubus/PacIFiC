#ifndef _DLMFD_GENERALPOLYHEDRON__
#define _DLMFD_GENERALPOLYHEDRON__

#include <DLMFD_RigidBody.hh>
#include <FS_GeneralPolyhedron.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <string>
using namespace std;

/** @brief The class DLMFD_GeneralPolyhedron.

A moving or stationary rigid 3D box of axisymmetric cross-section in the
Fictitious Domain solver.

@author A. Wachs - Pacific project 2025 */

class DLMFD_GeneralPolyhedron : public DLMFD_RigidBody
{
  public: //------------------------------------------------------------------
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Default constructor */
    DLMFD_GeneralPolyhedron();

    /** @brief Constructor with arguments
    @param pgrb Pointer to the geometric rigid body class */
    DLMFD_GeneralPolyhedron(FS_RigidBody *pgrb, const bool &are_particles_fixed,
                            FV_DiscreteField *pField_,
                            double const critical_distance_);

    /** @brief Destructor */
    ~DLMFD_GeneralPolyhedron();

    //@}

    //-- Set methods
    /** @name Set methods */
    //@{

    void set_ptr_FS_GeneralPolyhedron_Additional_Param();

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

    /** @brief Set boundary points on faces of an isometric cube */
    void setBndPts_IsometricCube(
        vector<geomVector> const &corners_, double spacing,
        FV_Mesh const *primary_grid,
        list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp,
        list<DLMFD_BoundaryMultiplierPoint *>::iterator &bphz);

    /** @brief Set boundary points on faces of an isometric tetrahedron */
    void setBndPts_IsometricTetrahedron(
        vector<geomVector> const &corners_, double spacing,
        FV_Mesh const *primary_grid,
        list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp,
        list<DLMFD_BoundaryMultiplierPoint *>::iterator &bphz);

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

    FS_GeneralPolyhedron_Additional_Param const *
    get_ptr_FS_GeneralPolyhedron_Additional_Param();

    //@}

    //-- Geometric methods
    /** @name Geometric methods */
    //@{

    /** @brief Calculate determinant 4 X 4 for checking a point in tetrahedron
     */
    double calcPointDeterm4by4(const geomVector &pointOne,
                               const geomVector &pointTwo,
                               const geomVector &pointThree,
                               const geomVector &pointFour) const;

    /** @brief Check whether a point is inside a tetrahedron */
    bool checkPointInTetrahedron(const geomVector &pointOne,
                                 const geomVector &pointTwo,
                                 const geomVector &pointThree,
                                 const geomVector &pointFour,
                                 const geomVector &pointToCheck) const;

    /** @brief isIn method
    @param point Point */
    bool isIn(const geomVector &point) const;

    /** @brief Translate the component from current gravity center to new
    gravity center
    @param newg new gravity center */
    void translateGeometricFeatures(geomVector const &newg);

    //@}

    //-- Add methods
    /** @name Add methods */
    //@{

    /** @brief Allocate default list of points and vectors for the polyhedron
    @param critical_distance Critical distance
    @param pField Constrained field */
    void allocate_default_listOfPointsAndVectors_GeneralPolyhedron(
        const double &critical_distance, FV_DiscreteField *pField);

    //@}

    //-- Add methods
    /** @name Add methods */
    //@{

    /** @brief Remove critical interior points
    @param critical_distance Critical_distance */
    void erase_critical_interior_points_PerProc(double critical_distance);

    //@}

  protected: //----------------------------------------------------------------
             //-- Attributes
             /** @name Parameters */
    //@{

    FS_GeneralPolyhedron_Additional_Param const
        *pagp; /** Pointer to FS GeneralPolyhedron parameters */

    vector<geomVector> corners; /**< corner coordinates of the polyhedron */
    vector<vector<size_t>> facesVec; /**< polygonal faces numbering */
    geomVector *g2; /**< slightly randomly translated gravity center */

    //@}

  private: //-------------------------------------------------------------------
};

#endif
