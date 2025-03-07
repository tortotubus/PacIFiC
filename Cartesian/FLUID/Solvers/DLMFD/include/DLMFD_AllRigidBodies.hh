#ifndef _DLMFD_ALLRIGIDBODIES__
#define _DLMFD_ALLRIGIDBODIES__

#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <doubleArray2D.hh>
#include <FS_AllRigidBodies.hh>
#include <DLMFD_RigidBody.hh>
#include <FV_DiscreteField.hh>
#include <MAC_Communicator.hh>
#include <DLMFD_ProjectionNavierStokesSystem.hh>
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
    @param pelCOMM_ Communicator
    @param solidFluid_transferStream Input stream where features of rigid bodies are read
    @param are_particles_fixed_ True if the all the particles are
    obstacles
    @param UU Pointer to flow field UF
    @param PP Pointer to flow field PF */
    DLMFD_AllRigidBodies(size_t &dim,
                         MAC_Communicator const *pelCOMM_,
                         istringstream &solidFluid_transferStream,
                         bool const &are_particles_fixed_,
                         FV_DiscreteField *UU,
                         FV_DiscreteField *PP);

    /** @brief Destructor */
    ~DLMFD_AllRigidBodies();

    //@}

    // -- Set methods
    /** @name Set methods */
    //@{

    void set_b_output_hydro_forceTorque(bool const &is_output);

    /** @brief Return the total number of rigid bodies
    @param critical_distance Critical distance */
    void set_all_points(double critical_distance);

    /** @brief Return the total number of rigid bodies */
    void eraseCriticalDLMFDPoints(const double &time, double critical_distance);

    /** @brief Set the list of IDs of on proc */
    void set_listIdOnProc();

    /** @brief Set constrained field
    @param pField Constrained field */
    void set_ptr_constrained_field(FV_DiscreteField *pField_);

    /** @brief Set points infos for all rigid bodies */
    void set_points_infos();

    void check_allocation_DLMFD_Cvectors();

    void fill_DLMFD_Cvectors();

    /** @brief Set MPI data
    @param pelCOMM_ Communicator */
    void set_MPI_data(MAC_Communicator const *pelCOMM_);

    /** @brief Set coupling factor of each rigid body
    @param rho_f Fluid density
    @param explicit_treatment Is explicit treated */
    void set_coupling_factor(double const &rho_f, bool const &explicit_treatment);

    /** @brief Set mass, density and inertia of all RBs */
    void set_mass_and_density_and_volume_and_inertia();

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

    //@}

    // -- Get methods
    /** @name Get methods */
    //@{

    /** @brief Return the total number of rigid bodies */
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
    /** @name Get methods */
    //@{

    /** @brief Update method
    @param critical_distance Critical distance */
    void update(double const &time, double critical_distance, istringstream &solidFluid_transferStream);

    //@}

    // -- Output methods
    /** @name Output methods */
    //@{

    /** @brief Output DLMFD points in Paraview
    @param filename File name */
    void output_DLMFDPoints_PARAVIEW(const string &filename,
                                     geomVector const *translated_distance_vector,
                                     const bool &withIntPts,
                                     size_t rank) const;

    /** @brief Output force and torque
    @param nothing File name */
    void particles_hydrodynamic_force_output(const string &path_name,
                                             const bool &b_restart,
                                             const double &time,
                                             const double &timestep,
                                             const double &rho_f,
                                             vector<vector<double>> const *Iw_Idw);

    /** @brief TODO
    @param nothing File name */
    void sum_DLM_hydrodynamic_force_output(const bool &b_restart);

    //@}

    // -- DLMFD solving methods
    /** @name DLMFD solving methods */
    //@{

    void nullify_all_Uzawa_vectors();

    /** @brief Allocate Uzawa vectors */
    void allocate_initialize_Uzawa_vectors();

    /** @brief Compute the q_U quantities */
    void compute_all_Qu(bool init);

    /** @brief Compute the q_w quantities */
    void compute_all_Qrot(bool init);

    /** @brief Add qtran to q_tran of rigid body i */
    void add_to_Qu(const size_t i, const geomVector &qtran);

    /** @brief Add qrot to q_rot_3D of rigid body i */
    void add_to_Qrot(const size_t i, const geomVector &qrot);

    /** @brief Complete the initialization of Uzawa solving algorithm */
    void completeFirstUzawaIteration_Velocity(const double &rho_f,
                                              const double &timestep,
                                              const geomVector &gravity_vector_split,
                                              const geomVector &gravity_vector);

    /** @brief Correct the particle velocities after broadcast of t_tran and
   t_rot i.e. do on each particle: U = t_tran and omega = t_rot */
    void update_ParticlesVelocities_afterBcast_of_T();

    /** @brief Solve one Uzawa iteration for the particles system
    @param rho_f Fluid density
    @param timestep Time step */
    void solve_Particles_OneUzawaIter_Velocity(double const &rho_f, double const &timestep);

    /** @brief Compute the momentum equations right hand side <w,v>_P
   as a list of position and entries in the right hand side vector of the
   global matrix system */
    void compute_fluid_rhs(DLMFD_ProjectionNavierStokesSystem *GLOBAL_EQ, bool init);

    /** @brief Compute residuals x=<alpha,tu-(tU+tomega^GM)>_P */
    void compute_x_residuals_Velocity();

    /** @brief Set r = - x and w = r */
    void compute_r_and_w_FirstUzawaIteration();

    /** @brief Return the VEC_r.VEC_r scalar product on all particles */
    double compute_r_dot_r() const;

    /** @brief Return the VEC_w.VEC_x scalar product on all particles */
    double compute_w_dot_x() const;

    /** @brief Update lambda and r i.e. do on each particle: lambda -= alpha * w and r -= alpha * x */
    void update_lambda_and_r(const double &alpha);

    /** @brief Correct the particle velocities i.e. do on each particle: U += t_tran*alpha and omega += t_rot*alpha  */
    void update_ParticlesVelocities(const double &alpha);

    /** @brief Update w i.e. do on each particle: w = r + beta * w */
    void update_w(const double &beta);

    /** @brief Compute the momentum equations right hand side <w,v>_P */
    void compute_fluid_DLMFD_explicit(DLMFD_ProjectionNavierStokesSystem *GLOBAL_EQ, bool bulk);

    //@}

    // -- Side methods
    /** @name Side methods */
    //@{

    /** @brief Allocate initial array for n-1 translational and angular velocity */
    void allocate_translational_angular_velocity_array();

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

    size_t dim = 3;                                        //*< Dimension */
    size_t RBs_number;                                     /**< Number of rigid bodies */
    size_t particles_number;                               /**< Number of particles */
    FS_AllRigidBodies *ptr_FSallrigidbodies;               /**< Pointer to the geometric rigid bodies */
    vector<DLMFD_RigidBody *> vec_ptr_DLMFDallrigidbodies; /**<  Pointer to the vector of DLMFD rigid bodies */
    list<int> entireOnProc;                                /**< list of ids of solid component entirely located
                                                            on this process */

    list<int> SharedOnProc;                    /**< list of ids of solid component partially located on this process */
    list<int> onProc;                          /**< list of ids of solid component located on this process */
    vector<size_t_vector> *v_AllSharedOnProcs; /**< vector of size nprocs stored on master process, element i is
                                                the list of ids of solid component partially located on process i */
    list<int> *l_AllSharedOnProcs;             /**< list, stored on master process, of ids
                                              of solid component shared on all processes, except those which are
                                              shared by master */

    doubleArray2D *translational_velocity_nm1; /**< array containing the
    translational velocity of all particules at previous time */
    doubleArray2D *angular_velocity_3D_nm1;    /**< array containing the
           angular velocity of all particules at previous time */

    // Pointers to the constant fields and primary grid
    FV_DiscreteField *pField; /**< Pointer to constrained field*/
    FV_Mesh const *MESH;

    // Output
    size_t output_frequency;

    //@}
};

#endif