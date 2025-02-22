#ifndef _DLMFD_RIGIDBODY__
#define _DLMFD_RIGIDBODY__

#include <FS_RigidBody.hh>
#include <geomVector.hh>
#include <DLMFD_BoundaryMultiplierPoint.hh>
#include <DLMFD_InteriorMultiplierPoint.hh>
#include <size_t_array3D.hh>
#include <doubleVector.hh>
using namespace std;

/**
  @brief The Structure ULBD_RHSInfos.
  @details Information to compute vectors at both the particle and matrix
    system levels in the Uzawa solving algorithm.
  @author A. Wachs - Particulate flow project 2007-2009
*/
struct ULBD_RHSInfos
{
    size_t compIdx;                       /**< component number */
    list<size_t> NodeNoU;                 /**< Velocity node number */
    list<FV_TRIPLET> MacTripletU;         /**< Velocity MAC triplets */
    list<size_t> positionU;               /**< Velocity position in the matrix system */
    list<double> omega_delta;             /**< scalar product omega*delta at the point */
    const DLMFD_ParticlePoint *ptr_point; /**< pointer to the DLM/FD point */
};

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
    void set_boundary_point(const geomVector &point, list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp);

    /** @brief Set DLMFD interior point in the list
    @param TO Write */
    void set_interior_point(const size_t &comp,
                            const geomVector &point,
                            size_t i, size_t j, size_t k,
                            list<DLMFD_InteriorMultiplierPoint *>::iterator &ip);

    /** @brief Set points infos
    @param pField Constrained field */
    void set_points_infos(FV_DiscreteField const *pField);

    /** @brief Set coupling factor
    @param rho_f Fluid density
    @param explicit_treatment Is explicit treated */
    void set_coupling_factor(double const &rho_f, bool const &explicit_treatment);

    /** @brief Set mass, density and inertia of RB */
    void set_mass_and_density_and_inertia();

    /** @brief Set translational velocity of RB */
    void set_translational_velocity();

    /** @brief Set angular velocity of RB */
    void set_angular_velocity();

    /** @brief Set t_tran of RB */
    void set_Tu(const geomVector &ttran);

    /** @brief Set t_rot of RB */
    void set_Trot(const geomVector &trot);

    //@}

    //-- Get methods
    /** @name Get methods */
    //@{

    /** @brief Returns the component type of the RB */
    string get_component_type() const;

    /** @brief Returns pointer to the rigid body gravity center */
    virtual geomVector const *get_ptr_to_gravity_centre() const = 0;

    /** @brief Returns circumscribed radius */
    virtual double get_circumscribed_radius() const = 0;

    /** @brief Returns the number of output points
    @param withIntPts True if the interior points are taken into account */
    size_t get_npts_output(bool const &withIntPts);

    /** @brief Get the rigid body velocity at the given point
    @param point Point */
    geomVector get_rigid_body_velocity(geomVector const &point) const;

    /** @brief Get the rigid body angular velocity*/
    geomVector get_rigid_body_angular_velocity() const;

    /** @brief Returns the translational velocity of the RB */
    geomVector get_rigid_body_translational_velocity() const;

    /** @brief Returns a tuple of mass and density of RB */
    tuple<double, double> get_mass_and_density() const;

    /** @brief Returns a tuple of mass and density of RB */
    vector<vector<double>> get_inertia() const;

    /** @brief Get q_tran */
    geomVector const get_Qu() const;

    /** @brief Get q_rot_3D */
    geomVector const get_Qrot() const;

    /** @brief Get t_tran */
    geomVector const get_Tu() const;

    /** @brief Get t_rot_3D */
    geomVector const get_Trot() const;

    //@}

    //-- Add methods
    /** @name Add methods */
    //@{

    /** @brief Extend the list of boundary points
    @param np Number of points */
    void extend_bp_list(size_t const &np);

    /** @brief Extend the list of interior points
    @param np Number of points */
    void extend_ip_list(size_t const &np);

    //@}

    //-- Update methods
    /** @name Update methods */
    //@{

    /** @brief Update the RB */
    virtual void update() = 0;

    /** @brief Update the RB position and velocity
    @param pos updated position
    @param vel updated translation velocity */
    void update_RB_position_and_velocity(geomVector const &pos,
                                         geomVector const &vel,
                                         geomVector const &ang_vel,
                                         vector<geomVector> const &periodic_directions,
                                         double const &time_step);

    //@}

    //-- Geometric methods
    /** @name Geometric methods */
    //@{

    /** @brief isIn method
    @param point Point */
    virtual bool isIn(const geomVector &point) const = 0;

    /** @brief Print the particule points coordinates */
    void print_partPointsCoordinates(ofstream &f,
                                     geomVector const *translated_distance_vector,
                                     const string &text2write_before,
                                     const string &text2write_after,
                                     bool const &withIntPts) const;

    /** @brief Translate the component from current gravity center to new
    gravity center
    @param newg new gravity center */
    virtual void translateGeometricFeatures(geomVector const &newg) = 0;

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

    // -- DLMFD computing methods
    /** @name DLMFD computing methods */
    //@{

    /** @brief Compute the q_U quantities */
    void compute_Qu(bool init);

    /** @brief Compute the q_w quantities */
    void compute_Qrot(bool init);

    /** @brief Add qtran to q_tran
    @param qtran To add to q_tran */
    void add_to_Qu(const geomVector &qtran);

    /** @brief Add qrot to q_rot_3D
    @param qrot To add to q_rot_3D */
    void add_to_Qrot(const geomVector &qrot);

    /** @brief Correct q_tran and q_rot vectors for initial Uzawa iteration i.e.
    set q_tran = (1-rho_f/rho_s)*M*U(n-1) + (1-rho_f/rho_s)*M*g_split - q_tran */
    void correctQvectorsAndInitUzawa_Velocity(const double &rho_f,
                                              const double &timestep,
                                              const geomVector &gravity_vector_split,
                                              const geomVector &gravity_vector);
    /** @brief Invert the inertia_3D matrix */
    vector<vector<double>> calcInvers3by3Matrix(const vector<vector<double>> &oldMatrix);

    /** @brief (1-rho_f/rho_s)*M*t_tran = q_tran and
   (1-rho_f/rho_s)*(I*t_rot+t_rot x I*t_rot) = q_rot
    @param rho_f Fluid density
    @param timestep Time step */
    void calculateParticleVelocities(double const &rho_f, double const &timestep);

    /** @brief Set U = t_tran*alpha+U*beta and omega = t_rot*alpha+omega*beta */
    void updateSolutionVectors_Velocity(const double &alpha, const double &beta);

    /** @brief Set the translational t_tran and rotational t_rot vectors, in
   case of constant velocity over the particle i.e. ttran and trot are
   initialized with the velocity to be imposed  */
    void setTVectors_constant();

    //@}

    // -- DLMFD other methods
    /** @name DLMFD other methods */
    //@{

    /** @brief Fill the DLMFD points infos */
    void fill_DLMFD_pointInfos(FV_DiscreteField const *pField,
                               struct ULBD_RHSInfos &OnePointInfos,
                               bool const &SecondOrderBP);

    double compute_weight_generic(double const &x,
                                  size_t const &i, size_t const &order);

    double Q2weighting(size_t const &i, double const &x,
                       double const &y, double const &z);

    //@}

protected: //--------------------------------------------------------------
    //-- Attributes

    FS_RigidBody *ptr_FSrigidbody; /* Pointer to geometric Rigid Body */

    // Geometric attributes
    size_t dim;
    size_t ndof;            /**< number of degrees of freedom */
    size_t ntotal_fieldunk; /**< total number of field unknowns for this solid component on the proc */

    list<DLMFD_BoundaryMultiplierPoint *> boundary_points; /* List of pointers to boundary points */
    size_t nBP;                                            /**< number of boundary points on proc */
    list<DLMFD_InteriorMultiplierPoint *> interior_points; /* List of pointers to interior points */
    size_t nIP;                                            /**< number of interior points on proc */

    string component_type; /**< type of component : particle, obstacle, ... */

    list<struct ULBD_RHSInfos> points_infos; /**< list of points infos to
        compute vectors at both the particle and matrix system levels
        in the Uzawa solving algorithm */

    doubleVector VEC_lambda = doubleVector(0, 0.); /**< DLM vector */
    doubleVector VEC_w = doubleVector(0, 0.);      /**< DLM work vector for the Uzawa/DLMFD algorithm */

    geomVector t_tran;   /**< translational velocity unknown */
    geomVector q_tran;   /**< translational velocity right hand side */
    geomVector t_rot_3D; /**< rotational velocity unknown */
    geomVector q_rot_3D; /**< rotational velocity right hand side */

    size_t_array2D *index_min; /**< lower bound index in the mesh related to the
 solid component circumscribing box */
    size_t_array2D *index_max; /**< upper bound index in the mesh related to the
     solid component circumscribing box */

    geomVector translational_velocity; /**< translational velocity */
    geomVector angular_velocity_3D;    /**< angular velocity */

    double fluidsolid_coupling_factor; /**< fluid/solid coupling factor */
    double rho_s;                      /**< Rigid body density */
    double mass;                       /**< mass */
    double volume;                     /**< volume */
    geomVector gravity_center;         /**< Center of gravity coordinates */
    double radius;                     /**< External radius */

    vector<vector<double>> inertia_3D;         /**< inertia tensor */
    vector<vector<double>> inversedInertia_3D; /**< inverse of the inertia tensor */

private: //----------------------------------------------------------------
    //-- Private attributes

    static size_t_array3D *Q2numb; /**< Q2 cubic interpolation (i,j,k)->basis */
};

#endif