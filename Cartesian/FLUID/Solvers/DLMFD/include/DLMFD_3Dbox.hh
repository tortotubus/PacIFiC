#ifndef _DLMFD_3DBOX__
#define _DLMFD_3DBOX__

#include <DLMFD_RigidBody.hh>
#include <FS_3Dbox.hh>
#include <FV_DiscreteField.hh>
#include <FV_Mesh.hh>
#include <string>
using namespace std;

/** @brief The class DLMFD_3Dbox.

A moving or stationary rigid 3D box of axisymmetric cross-section in the
Fictitious Domain solver.

@author A. Wachs - Pacific project 2025 */

class DLMFD_3Dbox : public DLMFD_RigidBody
{
  public: //------------------------------------------------------------------
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Default constructor */
    DLMFD_3Dbox();

    /** @brief Constructor with arguments
    @param pgrb Pointer to the geometric rigid body class
    @param are_particles_fixed True if all particles are fixed
    @param pField_ Constrained field
    @param critical_distance_ Critical distance */

    DLMFD_3Dbox(FS_RigidBody *pgrb, const bool &are_particles_fixed,
                FV_DiscreteField *pField_, double const critical_distance_);

    /** @brief Destructor */
    ~DLMFD_3Dbox();

    //@}

    //-- Set methods
    /** @name Set methods */
    //@{

    /** @brief Set the structure of geometrical parameters of the cube */
    void set_ptr_FS_3Dbox_Additional_Param();

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

    /** @brief Get the structure of geometrical parameters of the cube */
    FS_3Dbox_Additional_Param const *get_ptr_FS_3Dbox_Additional_Param();

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

    /** @brief Proximity query based on bounding box */
    bool proximityQuery(DLMFD_RigidBody const *second_component,
                        const double &distance) const;

    /** @brief Translate the component from current gravity center to new
    gravity center
    @param newg new gravity center */
    void translateGeometricFeatures(geomVector const &newg);

    //@}

    //-- Add methods
    /** @name Add methods */
    //@{

    /** @brief Allocate default list of points and vectors for the box
    @param critical_distance Critical distance
    @param pField Constrained field */
    void allocate_default_listOfPointsAndVectors_3Dbox(
        const double &critical_distance, FV_DiscreteField *pField);

    //@}

  protected: //----------------------------------------------------------------
             //-- Attributes
             /** @name Parameters */
    //@{

    FS_3Dbox_Additional_Param const *pagp; /** Pointer to FS 3Dbox parameters */

    vector<geomVector> corners; /**< corner coordinates of the polyhedron */
    vector<vector<size_t>> facesVec; /**< polygonal faces numbering */
    geomVector coor_min; /**< minimal coordinates of the bounding box */
    geomVector coor_max; /**< maximal coordinates of the bounding box */
    geomVector *g2;      /**< slightly randomly translated gravity center */

    vector<double> box_min, box_max;
    vector<size_t> Npoints; /** Surface points on each direction
                                of 3D box.*/
                            //@}

  private: //-------------------------------------------------------------------
};

#endif
