#ifndef DLMFD_FictitiousDomain_HH
#define DLMFD_FictitiousDomain_HH

#include <DLMFD_AllRigidBodies.hh>
#include <DLMFD_ProjectionNavierStokesSystem.hh>
#include <FS_SolidPlugIn.hh>
#include <FV_DomainAndFields.hh>
#include <FV_TimeIterator.hh>
#include <MAC_Communicator.hh>
#include <MAC_DoubleVector.hh>
#include <MAC_ModuleExplorer.hh>
#include <PAC_computingtime.hh>
#include <PAC_solvercomputingtime.hh>
using namespace std;

/** @brief The Class DLMFD_FictitiousDomain.

Solver for the coupling with particles using a Distributed Lagrange
Multiplier/Fictitious Domain method.

@author A. Wachs & M. Houlette - Pacific project 2024-2025 */

struct NavierStokes2FluidSolid
{
    // Output
    string solid_resDir;
    size_t output_frequency;

    // Parameters
    double rho_f;
    geomVector gravity_vector;
    geomVector split_gravity_vector;

    // Linear resolution
    DLMFD_ProjectionNavierStokesSystem *GLOBAL_EQ;
    size_t velocitylevelDiscrField;

    // Booleans
    bool b_restart;

    // Fields
    FV_DiscreteField *UU;
    FV_DiscreteField *PP;
};

class DLMFD_FictitiousDomain : public MAC_Object,
                               public PAC_ComputingTime,
                               public PAC_SolverComputingTime
{
  public: //-----------------------------------------------------------------
    //-- Public class attributes
    static bool b_SecondOrderInterpol;
    static doubleVector *dbnull;
    static double BoundaryPointsSpacing_coef;

    //-- Substeps of the step by step progression

    /** @brief Create and initialize an instance of DLMFD_FictitiousDomain
    @param a_owner The MAC-based object
    @param dom Domain
    @param exp To read the data file */
    static DLMFD_FictitiousDomain *create(MAC_Object *a_owner,
                                          FV_DomainAndFields const *dom,
                                          MAC_ModuleExplorer const *exp,
                                          NavierStokes2FluidSolid transfert);

    /** @name Substeps of the step by step progression */
    //@{

    /** @brief Tasks performed before the main loop
    @param t_it Time iterator */
    void do_before_time_stepping(FV_TimeIterator const *t_it);

    /** @brief Tasks performed at the main loop
    @param t_it Time iterator
    @param sub_prob_number Sub problem number */
    void do_one_inner_iteration(FV_TimeIterator const *t_it,
                                size_t &sub_prob_number);

    /** @brief Tasks performed just after the main loop
    @param t_it Time iterator */
    void do_after_inner_iterations_stage(FV_TimeIterator const *t_it);

    /** @brief Tasks performed for additional savings
    @param t_it Time iterator */
    void do_additional_savings(int const &cycleNumber,
                               FV_TimeIterator const *t_it,
                               const double &translated_distance,
                               const size_t &translation_direction);

    //@}

    //-- Set methods
    /** @name Set methods */
    //@{

    /** @brief Setting the critical distance attribute
    @param critical_distance_ Critical distance to set */
    void set_critical_distance(double critical_distance_);

    /** @brief Set the Paraview translated distance vector
    @param translated_distance translated distance magnitude
    @param translation_direction translation direction */
    void set_Paraview_translated_distance_vector(
        const double &translated_distance, const size_t &translation_direction);

    /** @brief Sets the Paraview post-processing translation vector in case of
    projection-translation */
    void setParaviewPostProcessingTranslationVector(const double &tvx,
                                                    const double &tvy,
                                                    const double &tvz);

    //@}

    //-- Get methods
    /** @name Get methods */
    //@{

    /** @brief Get DLMFD explicit boolean */
    bool const get_explicit_DLMFD() const;

    //@}

    //-- DLMFD solver methods
    /** @name DLMFD solver methods */
    //@{

    /** @brief Purely granular Newton's law
    @param t_it Time iterator
    @param sub_prob_number Sub problem number */
    void update_rigid_bodies(FV_TimeIterator const *t_it,
                             size_t &sub_prob_number);

    /** @brief Uzawa saddle point problem that solves the rigid body constraint
    @param t_it Time iterator
    @param sub_prob_number Sub problem number */
    void run_DLMFD_UzawaSolver(FV_TimeIterator const *t_it,
                               size_t &sub_prob_number);

    /** @brief Initialization of the DLMFD problem to solve
    @param t_it Time iterator */
    void DLMFD_construction(FV_TimeIterator const *t_it);

    /** @brief Solve the DLMFD problem
    @param t_it Time iterator */
    void DLMFD_solving(FV_TimeIterator const *t_it);

    //@}

    // -- Geometric methods
    /** @name Geometric methods */
    //@{

    double Compute_distance_to_bottom(const double &coordinate,
                                      const size_t &direction) const;

    void translate_all(const geomVector &translation_vector,
                       const size_t &translation_direction);

    //@}

  protected: //----------------------------------------------------------------
    //-- Class attributes
    size_t levelDiscrField;

  private: //----------------------------------------------------------------
    /** @name Constructors & Destructor */
    //@{

    /** @brief Constructor with arguments
    @param a_owner The MAC-based object
    @param dom Domain
    @param exp To read the data file */
    DLMFD_FictitiousDomain(MAC_Object *a_owner, FV_DomainAndFields const *dom,
                           MAC_ModuleExplorer const *exp,
                           NavierStokes2FluidSolid transfert);

    /** @brief Finalize construction */
    void finalize_construction(MAC_ModuleExplorer const *exp);

    /** @brief Constructor with arguments */
    ~DLMFD_FictitiousDomain();

    //@}

    //-- Parallel operations for DLM/FD problem
    /** @name Parallel operations for DLM/FD problem */
    //@{

    /** @brief Compute sum of q_tran=-<lambda,V>_P, where lambda
    is the Lagrange multipliers field and V the test function for the particle
    translational velocity, and q_rot = -<lambda,xi^GM>_P, where xi is the
    test function for the particle rotational velocity on a particle, of each
    process (MPI Reduction operation) on master process
    @param t_it time iterator */
    void calculate_ReductionFor_qtranAndqrot(FV_TimeIterator const *t_it);

    /** @brief Complete the initialization of DLM/FD-Uzawa solving algorithm
    i.e. do  on each particle of master process:
    <ul>
    <li> Set q_tran = (1-rho_f/rho_s)*M*U(n-1) - q_tran
    and q_rot = (1-rho_f/rho_s)*I*omega(n-1) - q_rot.
    <li> Solve (1-rho_f/rho_s)*M*t_tran = q_tran
    and (1-rho_f/rho_s)*I*t_rot = q_rot.
    <li> Set U = t_tran and omega = t_rot.
    </ul>
    then for each particle, broadcast t_tran and t_rot from master process
    on each process,
    and finally do for each particle of each process: U = t_tran and
    omega = t_rot. At the end of the method, each instance of a particle on
    every process has the same (and right) t_tran and t_rot.
    @param t_it time iterator */
    void velocity_broadcast_andUpdate_First(FV_TimeIterator const *t_it);

    /** @brief Solve one DLM/FD-Uzawa iteration for the particles system on
    master process i.e. do on
    each particle: Solve (1-rho_f/rho_s)*M*t_tran = q_tran and
    (1-rho_f/rho_s)*I*t_rot = q_rot,
    then for each particle, broadcast t_tran and t_rot from master process
    on each process. At the end of the method, each instance of a particle on
    every process has the same (and right) t_tran and t_rot.
    @param t_it time iterator */
    void velocity_broadcast_andUpdateInOneIt(FV_TimeIterator const *t_it);

    /** @brief Broadcast t_tran and t_rot of all particles shared by processes
    from master to processes that own these particles
    @param t_it time iterator */
    void
    Broadcast_tVectors_sharedParticles_MasterToAll(FV_TimeIterator const *t_it);

    /** @brief Set converged velocity of all particles on master process
    @param t_it time iterator */
    void Set_Velocity_AllParticles_Master(FV_TimeIterator const *t_it);

    /** @brief Compute the momentum equations right hand side
    coef*<lambda,v>_P and set it in vector q of the matrix system
    @param t_it time iterator
    @param init true if initialization step of Uzawa algorithm */
    void compute_fluid_LBD_rhs(FV_TimeIterator const *t_it, bool init);

    //@}

    //-- Output methods
    /** @name Output methods */
    //@{

    /** @brief Writing PVTU
    @param filename File name */
    void write_PVTU_multiplier_file(string const &filename) const;

    //@}

    //-- Attributes

    FV_DiscreteField *UU;
    FV_DiscreteField *PP;

    DLMFD_ProjectionNavierStokesSystem *GLOBAL_EQ;

    // Physical Parameters
    size_t dim;
    double rho_f;
    geomVector gravity_vector;
    geomVector split_gravity_vector;

    // Numerical parameters
    size_t sub_prob_number;
    double critical_distance;
    string coupling_scheme;

    // Grains3D variables
    string solidSolverType;
    FS_SolidPlugIn *solidSolver;
    string solidSolver_insertionFile;
    string solidSolver_simulationFile;
    istringstream *solidFluid_transferStream;
    DLMFD_AllRigidBodies *allrigidbodies;

    // Uzawa algorithm
    double Uzawa_DLMFD_precision;
    int Uzawa_DLMFD_maxiter;

    // Booleans
    bool b_restart;
    bool b_explicit_added_mass;
    bool are_particles_fixed;
    bool b_solidSolver_parallel;
    bool b_correct_particle_acceleration;
    bool b_ExplicitDLMFD;

    // MPI data
    MAC_Communicator const *macCOMM;
    size_t size_proc;
    size_t my_rank;
    size_t is_master;
    string *transferString;

    // Output
    string SolidSolverResultsDirectory;
    ostringstream Paraview_saveMultipliers_pvd;
    geomVector Paraview_translated_distance_vector;
    bool b_particles_verbose;

    // Transfer of velocity between solid solver and fluid solver
    vector<vector<double>> velocitiesVecGrains;
    vector<vector<double>> const *Iw_Idw;
};

#endif
