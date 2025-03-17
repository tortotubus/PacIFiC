#ifndef _DLMFD_RIGIDBODY__
#define _DLMFD_RIGIDBODY__

#include <FS_RigidBody.hh>
#include <geomVector.hh>
#include <DLMFD_BoundaryMultiplierPoint.hh>
#include <DLMFD_InteriorMultiplierPoint.hh>
#include <DLMFD_ProjectionNavierStokesSystem.hh>
#include <size_t_array3D.hh>
#include <doubleVector.hh>
using namespace std;

/** @brief The Structure ULBD_RHSInfos
    @details Information to compute vectors at both the particle and matrix
    system levels in the Uzawa solving algorithm. */
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
    DLMFD_RigidBody(FS_RigidBody *pgrb,
                    const bool &are_particles_fixed,
                    FV_DiscreteField *pField_);

    /** @brief Destructor */
    ~DLMFD_RigidBody();

    //@}

    //-- Set methods
    /** @name Set methods */
    //@{

    /** @brief Set component type */
    void set_component_type();

    void set_ptr_constrained_field(FV_DiscreteField *pField__);

    void initialize_listOfDLMFDPoints();

    /** @brief Set DLMFD boundary and interior points
    @param critical_distance Critical distance
    @param pField Constrained field */
    virtual void set_all_points(FV_DiscreteField *pField, double critical_distance) = 0;

    /** @brief Set DLMFD boundary point in the list
    @param TO Write */
    void set_boundary_point(const geomVector &point, list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp);

    /** @brief Set DLMFD halozone boundary point in the list
    param TO Write */
    void set_halozone_boundary_point(const geomVector &point, list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp);

    /** @brief Set DLMFD interior point in the list
    @param TO Write */
    void set_interior_point(const size_t &comp,
                            const geomVector &point,
                            size_t i, size_t j, size_t k,
                            list<DLMFD_InteriorMultiplierPoint *>::iterator &ip);

    /** @brief Set DLMFD halozone interior point in the list
    @param TO Write */
    void set_halozone_interior_point(const size_t &comp,
                                     const geomVector &point,
                                     size_t i, size_t j, size_t k,
                                     list<DLMFD_InteriorMultiplierPoint *>::iterator &ip);

    void setBndPoint(geomVector const &point,
                     FV_Mesh const *primary_grid,
                     list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp,
                     list<DLMFD_BoundaryMultiplierPoint *>::iterator &bphz);

    /** @brief Set boundary points located on a segment in the boundary_points list */
    void setPtsOnEdge(const geomVector &firstPoint,
                      const geomVector &secondPoint,
                      const double &critical_distance, const bool &addBeginPnt,
                      FV_Mesh const *primary_grid,
                      list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp,
                      list<DLMFD_BoundaryMultiplierPoint *>::iterator &bphz);

    /** @brief Remove critical interior points
    @param critical_distance Critical_distance */
    virtual void erase_critical_interior_points_PerProc(double critical_distance);

    /** @brief Remove critical boundary points too close from an other particle
    @param second_component Other particle
    @param critical_distance Critical_distance */
    void erase_critical_boundary_points_ptp_PerProc(DLMFD_RigidBody *second_component, const double &critical_distance);

    /** @brief Remove critical boundary points too close from an other particle
    @param second_component Other particle
    @param critical_distance Critical_distance */
    void erase_critical_boundary_points_ptp_PerProc_oneGC(DLMFD_RigidBody *second_component, const double &critical_distance);

    /** @brief Set points infos
    @param pField Constrained field */
    void set_points_infos(FV_DiscreteField *pField);

    void check_allocation_DLMFD_Cvectors();

    void fill_DLMFD_Cvectors();

    /** @brief Set coupling factor
    @param rho_f Fluid density
    @param explicit_treatment Is explicit treated */
    void set_coupling_factor(double const &rho_f, bool const &explicit_treatment);

    /** @brief Set mass, density and inertia of RB */
    void set_mass_and_density_and_inertia();

    /** @brief Set volume of RB */
    void set_volume();

    /** @brief Set translational velocity of RB from FS */
    void set_translational_velocity();

    /** @brief Set translational velocity of RB
    @param vtran Velocity to set */
    void set_translational_velocity(const geomVector &vtran);

    /** @brief Set angular velocity of RB from FS */
    void set_angular_velocity_3D();

    /** @brief Set angular velocity of RB
    @param vtran Angular velocity to set */
    void set_angular_velocity_3D(const geomVector &vrot);

    /** @brief Set q_tran of RB */
    void set_Qu(const geomVector &qtran);

    /** @brief Set q_rot_3D of RB */
    void set_Qrot(const geomVector &qrot);

    /** @brief Set t_tran of RB */
    void set_Tu(const geomVector &ttran);

    /** @brief Set t_rot of RB */
    void set_Trot(const geomVector &trot);

    void set_ttran_ncomp(const size_t &ncomp_);

    //@}

    //-- Get methods
    /** @name Get methods */
    //@{

    int get_number_periodicClones() const;

    vector<geomVector> *get_periodicClones_DirectionsVector() const;

    /** @brief Returns the component type of the RB */
    string get_component_type() const;

    GEOMETRICSHAPE get_GeometricObjectType() const;

    /** @brief Returns pointer to the rigid body gravity center */
    geomVector const *get_ptr_to_gravity_centre() const;

    /** @brief Returns circumscribed radius */
    double get_circumscribed_radius() const;

    /** @brief Returns the number of output points
    @param withIntPts True if the interior points are taken into account */
    size_t get_npts_output(bool const &withIntPts);

    /** @brief Get the rigid body velocity at the given point
    @param point Point */
    geomVector get_rigid_body_velocity(geomVector const &point) const;

    /** @brief Get the rigid body angular velocity from FS */
    geomVector get_rigid_body_angular_velocity() const;

    /** @brief Get the rigid body angular velocity object */
    geomVector get_angular_velocity_3D() const;

    /** @brief Get the rigid body translational velocity from FS */
    geomVector get_rigid_body_translational_velocity() const;

    /** @brief Get the rigid body translational velocity object */
    geomVector get_translational_velocity() const;

    /** @brief Returns the mass of RB from FS */
    double get_rigid_body_mass() const;

    /** @brief Returns the density of RB from FS */
    double get_rigid_body_density() const;

    /** @brief Returns the density object of RB */
    double get_density() const;

    /** @brief Returns the volume object of RB */
    double get_volume() const;

    /** @brief Returns the inertia tensor of RB */
    vector<vector<double>> get_rigid_body_inertia() const;

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

    /** @brief Extend the list of halozone boundary points
    @param np Number of points */
    void extend_bphz_list(size_t const &np);

    /** @brief Extend the list of interior points
    @param np Number of points */
    void extend_ip_list(size_t const &np);

    /** @brief Extend the list of halozone interior points
    @param np Number of points */
    void extend_iphz_list(size_t const &np);

    void allocate_default_listOfPointsAndVectors(size_t const &nbIPdef,
                                                 size_t const &nbBPdef,
                                                 size_t const &nbIPHZdef,
                                                 size_t const &nbBPHZdef);

    //@}

    //-- Update methods
    /** @name Update methods */
    //@{

    /** @brief Update the RB */
    virtual void update();

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

    /** @brief Proximity query */
    virtual bool proximityQuery(DLMFD_RigidBody const *second_component, const double &distance) const;

    void DLMFDPoints_in_ContactRegion(list<DLMFD_BoundaryMultiplierPoint *> &BP_contactRegion,
                                      list<DLMFD_BoundaryMultiplierPoint *> &BPHZ_contactRegion,
                                      geomVector const &refPoint,
                                      double const &distance_contactRegion);

    /** @brief True if there is at least one valid point on proc */
    bool hasDLMFDPointsOnProc() const;

    /** @brief True if there is at least one halozone point on proc */
    bool hasDLMFDPoints_inHalozone_OnProc() const;

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

    void nullify_Uzawa_vectors();

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

    void compute_fluid_rhs(DLMFD_ProjectionNavierStokesSystem *GLOBAL_EQ, bool init);

    /** @brief Compute residuals x=<alpha,tu-(tU+tomega^GM)>_P */
    void compute_x_residuals_Velocity(FV_DiscreteField *pField);

    /** @brief Set r = - x and w = r */
    void compute_r_and_w_FirstUzawaIteration();

    /** @brief Return the VEC_r.VEC_r scalar product */
    double compute_r_dot_r() const;

    /** @brief Return the VEC_r.VEC_r scalar product */
    double compute_w_dot_x() const;

    /** @brief Update lambda and r i.e. do: lambda -= alpha * w and r -= alpha * x */
    void update_lambda_and_r(const double &alpha);

    /** @brief Update w i.e. do: w = r + beta * w */
    void update_w(const double &beta);

    /** @brief Compute the explicit momentum equations right hand side
    <lambda,v>_P */
    void compute_fluid_DLMFD_explicit(DLMFD_ProjectionNavierStokesSystem *GLOBAL_EQ,
                                      bool bulk,
                                      FV_DiscreteField *pField);

    //@}

    // -- DLMFD other methods
    /** @name DLMFD other methods */
    //@{

    /** @brief Fill the DLMFD points infos */
    void fill_DLMFD_pointInfos(FV_DiscreteField *pField,
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
    FV_DiscreteField *pField_;

    bool is_particle_fixed; /**< is the solid component a particle treated as a fixed obstacle */

    // Geometric attributes
    size_t dim;
    size_t ndof;            /**< number of degrees of freedom */
    size_t ntotal_fieldunk; /**< total number of field unknowns for this solid component on the proc */

    list<DLMFD_BoundaryMultiplierPoint *> boundary_points;          /* List of pointers to boundary points */
    size_t nBP;                                                     /**< number of boundary points on proc */
    list<DLMFD_InteriorMultiplierPoint *> interior_points;          /* List of pointers to interior points */
    size_t nIP;                                                     /**< number of interior points on proc */
    list<DLMFD_BoundaryMultiplierPoint *> halozone_boundary_points; /**< list of boundary points in the halo zone of this process */
    size_t nBPHZ;                                                   /**< number of boundary points in halozone */
    list<DLMFD_InteriorMultiplierPoint *> halozone_interior_points; /**< list of interior points in the halo zone of this process */
    size_t nIPHZ;                                                   /**< number of interior points in halozone */

    string component_type; /**< type of component : particle, obstacle, ... */

    list<struct ULBD_RHSInfos> points_infos; /**< list of points infos to
        compute vectors at both the particle and matrix system levels
        in the Uzawa solving algorithm */

    doubleVector VEC_r;      /**< residual vector for the Uzawa/DLMFD algorithm */
    doubleVector VEC_x;      /**< work vector for the Uzawa/DLMFD algorithm */
    doubleVector VEC_lambda; /**< DLM vector */
    doubleVector VEC_w;      /**< DLM work vector for the Uzawa/DLMFD */

    geomVector t_tran;   /**< translational velocity unknown */
    geomVector q_tran;   /**< translational velocity right hand side */
    geomVector t_rot_3D; /**< rotational velocity unknown */
    geomVector q_rot_3D; /**< rotational velocity right hand side */

    size_t_array2D *index_min; /**< lower bound index in the mesh related to the
    solid component circumscribing box */
    size_t_array2D *index_max; /**< upper bound index in the mesh related to the
     solid component circumscribing box */

    size_t *NDOF_comp;       /**< vector of component number of constrained DOF */
    double *NDOF_leverage;   /**< vector of leverage of constrained velocity
      DOF */
    size_t *NDOF_globalpos;  /**< vector of global position in the linear system
     of constrained DOF */
    size_t *NDOF_nfieldUNK;  /**< vector of number of field unknowns associated
         to a constrained DOF */
    double *NDOF_deltaOmega; /**< vector of delta*omega weights for each field
        unknown associated to a constrained DOF */
    size_t *NDOF_FVTriplet;  /**< vector of MAC triplet (i,j,k) for each field
     unknown associated to a constrained DOF */

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

    vector<geomVector> *periodic_directions; /**< periodic clones whose position is gravity_center + (*periodic_directions)[i] */

private: //----------------------------------------------------------------
    //-- Private attributes

    static size_t_array3D *Q2numb; /**< Q2 cubic interpolation (i,j,k)->basis */
};

#endif