#ifndef _DLMFD_ALLRIGIDBODIES__
#define _DLMFD_ALLRIGIDBODIES__

#include <DLMFD_AllGeomBoundaries.hh>
#include <DLMFD_ProjectionNavierStokesSystem.hh>
#include <DLMFD_RigidBody.hh>
#include <FS_AllRigidBodies.hh>
#include <FV_DiscreteField.hh>
#include <MAC_Communicator.hh>
#include <doubleArray2D.hh>
#include <size_t_array2D.hh>
#include <size_t_vector.hh>
#include <DLMFD_System.hh>
#include <FS_AllRigidBodies.hh>
#include <FV_DiscreteField.hh>
#include <MAC_Communicator.hh>
#include <doubleArray2D.hh>
#include <size_t_array2D.hh>
#include <size_t_vector.hh>
using namespace std;

/** @brief The class DLMFD_AllRigidBodies.

The array of all rigid bodies in the Fictitious Domain solver.

@author M. Houlette - Pacific project 2025 */

class DLMFD_AllRigidBodies
{
  public: //-----------------------------------------------------------------
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Constructor with arguments
    @param dim Number of space dimensions
    @param time Time
    @param solidFluid_transferStream Input stream where features of rigid bodies
    are read
    @param are_particles_fixed_ True if all the particles are fixed
    @param UU Pointer to flow field UF
    @param PP Pointer to flow field PF
    @param critical_distance Critical distance */
    DLMFD_AllRigidBodies(size_t &dim, double const &time,
                         istringstream &solidFluid_transferStream,
                         bool const &are_particles_fixed_, FV_DiscreteField *UU,
                         FV_DiscreteField *PP, double const critical_distance);

    /** @brief Destructor */
    ~DLMFD_AllRigidBodies();

    //@}

    // -- Set methods
    /** @name Set methods */
    //@{

    void set_components_type();

    void set_b_output_hydro_forceTorque(bool const &is_output);

    /** @brief Allocate and set the multiplier points
    @param critical_distance Critical distance */
    void set_all_MAC(double critical_distance);

    /** @brief Set the multiplier points
    @param critical_distance Critical distance */
    void set_all_points(double critical_distance);

    /** @brief Erase critical DLMFD points */
    void eraseCriticalDLMFDPoints(const double &time, double critical_distance);

    /** @brief Set the list of IDs of on proc */
    void set_listIdOnProc();

    /** @brief Set constrained field
    @param pField_ Constrained field to set */
    void set_ptr_constrained_field(FV_DiscreteField *pField_);

    void set_ptr_constrained_field_in_all_particles();

    void set_geometric_boundaries();

    /** @brief Set points infos for all rigid bodies */
    void set_points_infos();

    /** @brief Allocate DLMFD C vectors */
    void check_allocation_DLMFD_Cvectors();

    /** @brief Fill DLMFD C vectors */
    void fill_DLMFD_Cvectors();

    /** @brief Set MPI data */
    void set_MPI_data();

    /** @brief Set coupling factor of each rigid body
    @param rho_f Fluid density
    @param explicit_treatment If true, we solve by treating the added mass term
    (rho_f/rho_s) * (dU/dt) as explicit in the Newton's law */
    void set_coupling_factor(double const &rho_f,
                             bool const &explicit_treatment);

    /** @brief Set mass, density and inertia of all RBs */
    void set_mass_and_density_and_volume_and_inertia(
        istringstream &solidFluid_transferStream);

    /** @brief Set t_tran of RB i */
    void set_Tu(const size_t i, const geomVector &ttran);

    /** @brief Set t_rot of RB i */
    void set_Trot(const size_t i, const geomVector &trot);

    /** @brief Set the number of moving RBs */
    void const set_npart();

    /** @brief Set translational velocity of RB i */
    void set_translational_velocity(const size_t i, const geomVector &vtran);

    /** @brief Set angular velocity of RB i */
    void set_angular_velocity_3D(const size_t i, const geomVector &vrot);

    /** @brief Set output frequency */
    void set_output_frequency(size_t const output_frequency_);

    void set_ttran_ncomp(const size_t &ncomp_);

    //@}

    // -- Get methods
    /** @name Get methods */
    //@{

    /** @brief Get the total number of rigid bodies */
    size_t get_number_rigid_bodies() const;

    /** @brief Get shared rigid bodies on procs */
    list<int> const *get_SharedOnProc() const;

    /** @brief Get the pointer to the processes shared solid components */
    vector<size_t_vector> const *get_v_AllSharedOnProcs() const;

    /** @brief Get q_tran of rigid body i
    @param i Rigid body index */
    geomVector const get_Qu(const size_t i) const;

    /** @brief Get q_rot_3D of rigid body i
    @param i Rigid body index */
    geomVector const get_Qrot(const size_t i) const;

    /** @brief Get t_tran of rigid body i
    @param i Rigid body index */
    geomVector const get_Tu(const size_t i) const;

    /** @brief Get t_rot_3D of rigid body i
    @param i Rigid body index */
    geomVector const get_Trot(const size_t i) const;

    size_t const get_npart();

    list<int> const *get_onProc() const;

    /** @brief Get translational velocity of rigid body i
    @param i Rigid body index */
    geomVector get_translational_velocity(const size_t i) const;

    /** @brief Get angular velocity of rigid body i
    @param i Rigid body index */
    geomVector get_angular_velocity_3D(const size_t i) const;

    double const get_volume(const size_t i) const;

    //@}

    // -- Update methods
    /** @name Update methods */
    //@{

    /** @brief Update method
    @param critical_distance Critical distance */
    void update(istringstream &solidFluid_transferStream);

    /** @brief Store the RBs velocities to transfer to Grains3D
    @param vecVel Storage vector */
    void particles_velocities_output(vector<vector<double>> &vecVel) const;

    //@}

    // -- Output methods
    /** @name Output methods */
    //@{

    /** @brief Output DLMFD points in Paraview
    @param filename File name */
    void
    output_DLMFDPoints_PARAVIEW(const string &filename,
                                geomVector const *translated_distance_vector,
                                const bool &withIntPts) const;

    bool is_hydro_forceTorque_postprocessed() const;

    /** @brief Output force and torque
    @param nothing File name */
    void particles_hydrodynamic_force_output(
        const string &path_name, const bool &b_restart, const double &time,
        const double &timestep, const double &rho_f,
        vector<vector<double>> const *Iw_Idw);

    /** @brief MPI routine to sum the contributions of the hydrodynamic force of
    all procs
    @param b_restart True if the simulation is reloaded */
    void sum_DLM_hydrodynamic_force_output(const bool &b_restart);

    /** @brief Read the particles outputs from the explorer
    @param exp Explorer */
    void read_particles_outputs(MAC_ModuleExplorer const *exp);

    /** @brief Write in the output files
    @param path_name Path where to write
    @param b_restart True if the simulation is reloaded
    @param time Time */
    void particles_features_output(const string &path_name,
                                   const bool &b_restart,
                                   const double &time) const;

    //@}

    // -- DLMFD solving methods
    /** @name DLMFD solving methods */
    //@{

    /** @brief Nullify all the vectors used in Uzawa algorithm */
    void nullify_all_Uzawa_vectors();

    /** @brief Compute the q_U quantities */
    void compute_all_Qu(bool init);

    /** @brief Compute the q_w quantities */
    void compute_all_Qrot(bool init);

    /** @brief Add qtran to q_tran of rigid body i */
    void add_to_Qu(const size_t i, const geomVector &qtran);

    /** @brief Add qrot to q_rot_3D of rigid body i */
    void add_to_Qrot(const size_t i, const geomVector &qrot);

    /** @brief Complete the initialization of Uzawa solving algorithm */
    void
    completeFirstUzawaIteration_Velocity(const double &rho_f,
                                         const double &timestep,
                                         const geomVector &gravity_vector_split,
                                         const geomVector &gravity_vector);

    /** @brief Correct the particle velocities after broadcast of t_tran and
   t_rot i.e. do on each particle: U = t_tran and omega = t_rot */
    void update_ParticlesVelocities_afterBcast_of_T();

    /** @brief Solve one Uzawa iteration for the particles system
    @param rho_f Fluid density
    @param timestep Time step */
    void solve_Particles_OneUzawaIter_Velocity(double const &rho_f,
                                               double const &timestep);

    /** @brief Compute the momentum equations right hand side <w,v>_P
    as a list of position and entries in the right hand side vector of the
    global matrix system
    @param GLOBAL_EQ Linear system object
    @param init True if very first Uzawa iteration */
    void compute_fluid_rhs(DLMFD_System *GLOBAL_EQ,
                           bool init);

    /** @brief Compute residuals x=<alpha,tu-(tU+tomega^GM)>_P */
    void compute_x_residuals_Velocity();

    /** @brief Set r = - x and w = r */
    void compute_r_and_w_FirstUzawaIteration();

    /** @brief Return the VEC_r.VEC_r scalar product on all particles */
    double compute_r_dot_r() const;

    /** @brief Return the VEC_w.VEC_x scalar product on all particles */
    double compute_w_dot_x() const;

    /** @brief Update lambda and r i.e. do on each particle: lambda -= alpha * w
     * and r -= alpha * x */
    void update_lambda_and_r(const double &alpha);

    /** @brief Correct the particle velocities i.e. do on each particle: U +=
     * t_tran*alpha and omega += t_rot*alpha  */
    void update_ParticlesVelocities(const double &alpha);

    /** @brief Update w i.e. do on each particle: w = r + beta * w */
    void update_w(const double &beta);

    /** @brief Compute the momentum equations right hand side <w,v>_P */
    void
    compute_fluid_DLMFD_explicit(DLMFD_System *GLOBAL_EQ,
                                 bool bulk);

    //@}

    // -- Allocation methods
    /** @name Allocation methods */
    //@{

    /** @brief Allocate initial array for n-1 translational and angular velocity
     */
    void allocate_translational_angular_velocity_array();

    //@}

    // -- Geometric methods
    /** @name Geometric methods */
    //@{

    /** @brief Compute : min[ abs(gravity_center_i(direction) - coordinate) ],
    where i is a rigid body
    @param coordinate Coordinate
    @param direction Direction */
    double compute_minimum_distance_to_bottom(const double &coordinate,
                                              const size_t &direction) const;
    /** @brief Translate the geometric boundaries
    @param translation_vector Translation vector
    @param translation_direction Translation direction */
    void translate_geometricBoundaries(const geomVector &translation_vector,
                                       const size_t &translation_direction);

    //@}

  protected: //----------------------------------------------------------------
    bool b_output_hydro_forceTorque;
    bool are_particles_fixed;

  private: //----------------------------------------------------------------
    //-- Attributes
    /** @name Parameters */
    //@{

    // MPI data
    MAC_Communicator const *pelCOMM;
    size_t size_procs;
    size_t rank;
    size_t master;

    size_t dim = 3;          //*< Dimension */
    size_t RBs_number;       /**< Number of rigid bodies */
    size_t particles_number; /**< Number of particles */
    FS_AllRigidBodies
        *ptr_FSallrigidbodies; /**< Pointer to the geometric rigid bodies */
    vector<DLMFD_RigidBody *>
        vec_ptr_DLMFDallrigidbodies; /**<  Vector of DLMFD rigid
                                        bodies pointers*/
    list<int> entireOnProc; /**< list of ids of rigid bodies entirely located
                             on this process */

    list<int> SharedOnProc; /**< list of ids of rigid bodies partially
                               located on this process */
    list<int>
        onProc; /**< list of ids of rigid bodies located on this process */
    vector<size_t_vector>
        *v_AllSharedOnProcs;       /**< vector of size nprocs stored on master
                                    process, element i is       the list of ids of
                                    rigid bodies partially located on process
                                    i */
    list<int> *l_AllSharedOnProcs; /**< list, stored on master process, of ids
                                  of rigid bodies shared on all processes,
                                  except those which are shared by master */

    string
        output_x_velocity; /**< X-velocity output type: yes, no, mean or all */
    string
        output_y_velocity; /**< Y-velocity output type: yes, no, mean or all */
    string
        output_z_velocity; /**< Z-velocity output type: yes, no, mean or all */
    string
        output_x_position; /**< X-position output type: yes, no, mean or all */
    string
        output_y_position; /**< Y-position output type: yes, no, mean or all */
    string
        output_z_position; /**< Z-position output type: yes, no, mean or all */
    string output_x_omega; /**< X-Omega output type: yes, no, mean or all */
    string output_y_omega; /**< Y-Omega output type: yes, no, mean or all */
    string output_z_omega; /**< Z-Omega output type: yes, no, mean or all */
    string output_orientation; /**< Orientation output type: yes, no, mean or
                                  all */

    doubleArray2D *translational_velocity_nm1; /**< array containing the
    translational velocity of all particules at previous time */
    doubleArray2D *angular_velocity_3D_nm1;    /**< array containing the
           angular velocity of all particules at previous time */

    // Pointers to the constant fields and primary grid
    FV_DiscreteField *pField; /**< Pointer to constrained field*/
    FV_Mesh const *MESH;

    DLMFD_AllGeomBoundaries
        *GeoBoundaries; /**< geometric boundaries of the domain */

    // Output
    size_t output_frequency;

    //@}
};

#endif