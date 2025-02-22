#ifndef _DLMFD_ALLRIGIDBODIES__
#define _DLMFD_ALLRIGIDBODIES__

#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <FS_AllRigidBodies.hh>
#include <DLMFD_RigidBody.hh>
#include <FV_DiscreteField.hh>
#include <MAC_Communicator.hh>
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
    @param are_particles_fixed True if the all the particles are
    obstacles
    @param UF Pointer to flow field UF
    @param PF Pointer to flow field PF */
    DLMFD_AllRigidBodies(size_t &dim,
                         MAC_Communicator const *pelCOMM_,
                         istringstream &solidFluid_transferStream,
                         bool const &are_particles_fixed,
                         FV_DiscreteField const *UF,
                         FV_DiscreteField const *PF);

    /** @brief Destructor */
    ~DLMFD_AllRigidBodies();
    //@}

    // -- Set methods
    /** @name Set methods */
    //@{

    /** @brief Return the total number of rigid bodies
    @param critical_distance Critical distance */
    void set_all_points(double critical_distance);

    /** @brief Set the list of IDs of on proc */
    void set_listIdOnProc();

    /** @brief Set constrained field
    @param pField Constrained field */
    void set_ptr_constrained_field(FV_DiscreteField const *pField_);

    /** @brief Set points infos for all rigid bodies
    @param pField Constrained field */
    void set_points_infos(FV_DiscreteField const *pField_);

    /** @brief Set MPI data
    @param pelCOMM_ Communicator */
    void set_MPI_data(MAC_Communicator const *pelCOMM_);

    /** @brief Set coupling factor of each rigid body
    @param rho_f Fluid density
    @param explicit_treatment Is explicit treated */
    void set_coupling_factor(double const &rho_f, bool const &explicit_treatment);

    /** @brief Set mass, density and inertia of all RBs */
    void set_mass_and_density_and_inertia();

    /** @brief Set t_tran of RB i */
    void set_Tu(const size_t i, const geomVector &ttran);

    /** @brief Set t_rot of RB i */
    void set_Trot(const size_t i, const geomVector &trot);

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

    //@}

    // -- Update methods
    /** @name Get methods */
    //@{

    /** @brief Update method
    @param critical_distance Critical distance */
    void update(double critical_distance, istringstream &solidFluid_transferStream);

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

    //@}

    // -- DLMFD solving methods
    /** @name DLMFD solving methods */
    //@{

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

    //@}

protected: //----------------------------------------------------------------
private:   //----------------------------------------------------------------
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
    FS_AllRigidBodies *ptr_FSallrigidbodies;               /**< Pointer to the geometric rigid bodies */
    vector<DLMFD_RigidBody *> vec_ptr_DLMFDallrigidbodies; /**<  Pointer to the vector of DLMFD rigid bodies */
    list<int> entireOnProc;                                /**< list of ids of solid component entirely located
                                                            on this process */

    list<int> SharedOnProc;                          /**< list of ids of solid component partially located on this process */
    list<int> onProc;                                /**< list of ids of solid component located on this process */
    vector<size_t_vector> const *v_AllSharedOnProcs; /**< vector of size nprocs stored on master process, element i is
                                                the list of ids of solid component partially located on process i */
    list<int> *l_AllSharedOnProcs;                   /**< list, stored on master process, of ids
                                                    of solid component shared on all processes, except those which are
                                                    shared by master */

    //@}

    // Pointers to the constant fields and primary grid
    FV_DiscreteField const *pField; /**< Pointer to constrained field*/
    FV_Mesh const *MESH;
};

#endif