#ifndef _DLMFD_RIGIDBODY__
#define _DLMFD_RIGIDBODY__

#include <FS_RigidBody.hh>
#include <geomVector.hh>
#include <DLMFD_BoundaryMultiplierPoint.hh>
#include <DLMFD_InteriorMultiplierPoint.hh>
using namespace std;

/** @brief The class DLMFD_RigidBody.

A moving or stationary rigid body in the Fictitious Domain solver.

@author M. Houlette - Pacific project 2025 */

class DLMFD_RigidBody
{
public: //-----------------------------------------------------------------
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Default constructor */
    DLMFD_RigidBody();

    /** @brief Constructor with arguments
    @param pgrb Pointer to the geometric rigid body class */
    DLMFD_RigidBody(FS_RigidBody *pgrb);

    /** @brief Destructor */
    ~DLMFD_RigidBody();

    //@}

    //-- Set methods
    /** @name Set methods */
    //@{

    /** @brief Set DLMFD boundary and interior points
    @param critical_distance Critical distance
    @param pField Constrained field */
    virtual void set_all_points(FV_DiscreteField const *pField, double critical_distance) = 0;

    /** @brief Set DLMFD boundary point in the list
    @param TO Write */
    void setBndPoint(double const &x, double const &y,
                     double const &z, FV_Mesh const *primary_grid,
                     list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp);

    //@}

    //-- Get methods
    /** @name Get methods */
    //@{

    /** @brief Returns pointer to the rigid body gravity center */
    virtual geomVector const *get_ptr_to_gravity_centre() const = 0;

    /** @brief Returns circumscribed radius */
    virtual double get_circumscribed_radius() const = 0;

    /** @brief Returns the number of output points
    @param withIntPts True if the interior points are taken into account */
    size_t get_npts_output(bool const &withIntPts);

    //@}

    //-- Add methods
    /** @name Add methods */
    //@{

    /** @brief Extend the list of boundary points
    @param np Number of points */
    void extent_bp_list(size_t const &np);

    /** @brief Extend the list of interior points
    @param np Number of points */
    void extent_ip_list(size_t const &np);

    //@}

    //-- Update methods
    /** @name Update methods */
    //@{

    /** @brief Update the RB position and velocity
    @param pos updated position
    @param vel updated translation velocity */
    virtual void update_RB_position_and_velocity(geomVector const &pos,
                                                 geomVector const &vel,
                                                 geomVector const &ang_vel,
                                                 vector<geomVector> const &periodic_directions,
                                                 double const &time_step) = 0;

    //@}

    //-- Geometric methods
    /** @name Geometric methods */
    //@{

    /** @brief isIn method
    @param x x-component
    @param y y-component
    @param z z-component */
    virtual bool isIn(double const &x, double const &y, double const &z) const = 0;

    /** @brief Print the particule points coordinates */
    void print_partPointsCoordinates(ofstream &f,
                                     geomVector const *translated_distance_vector,
                                     const string &text2write_before,
                                     const string &text2write_after,
                                     bool const &withIntPts) const;

    //@}

    //-- Clear methods
    /** @name Clear methods */
    //@{

    /** @brief Clear the lists of points */
    void clear_listOfPointsAndVectors();

    //@}

    //-- Boolean methods
    /** @name Boolean methods */
    //@{

    /** @brief True if interior_points is empty */
    bool is_interior_points_empty();

    //@}

protected: //--------------------------------------------------------------
    //-- Attributes

    FS_RigidBody *ptr_FSrigidbody; /* Pointer to geometric Rigid Body */

    // Geometric attributes
    list<DLMFD_BoundaryMultiplierPoint *> boundary_points; /* List of pointers to boundary points */
    size_t nBP;                                            /**< number of boundary points on proc */
    list<DLMFD_InteriorMultiplierPoint *> interior_points; /* List of pointers to interior points */
    size_t nIP;                                            /**< number of interior points on proc */

    size_t_array2D *index_min; /**< lower bound index in the mesh related to the
 solid component circumscribing box */
    size_t_array2D *index_max; /**< upper bound index in the mesh related to the
     solid component circumscribing box */

    geomVector gravity_center; /**< Center of gravity coordinates */
    double radius;             /**< External radius */

private: //----------------------------------------------------------------
};

#endif