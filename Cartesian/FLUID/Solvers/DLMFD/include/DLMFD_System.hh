#ifndef _DLMFD_SYSTEM_BUILDERFACTORY__
#define _DLMFD_SYSTEM_BUILDERFACTORY__

#include <iostream>
#include <FV_DiscreteField.hh>
#include <MAC_Object.hh>
#include <MAC_ModuleExplorer.hh>
using namespace std;

/** @brief The Class DLMFD_System.

Linear system resolution class, that is either equal to DLMFD_ProjectionNavierStokesSystem
or DLMFD_DirectionSplittingSystem, depending on the flow solver that is actually used.

@author M. Houlette - Pacific project 2025 */

class DLMFD_System
{
public: //-----------------------------------------------------------------
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Default constructor */
    DLMFD_System();

    /** @brief Copy constructor */
    DLMFD_System(DLMFD_System const &other);

    /** @brief Constructor with arguments
    @param pgrb Pointer to the geometric rigid body class */
    DLMFD_System(MAC_Object *a_owner,
                 MAC_ModuleExplorer const *exp,
                 FV_DiscreteField *mac_UF,
                 FV_DiscreteField *mac_PF,
                 size_t const &NS_Viscous_TimeAccuracy_,
                 size_t const &NS_Advection_TimeAccuracy_,
                 bool const &b_pressure_rescaling_,
                 bool const &b_ExplicitPressureGradient_,
                 bool const &b_HighOrderPressureCorrection_,
                 bool is_stressCal_);
    //@});

    /** @brief Destructor */
    ~DLMFD_System();

    //@}

    //-- Access

    /** @name Access */
    //@{
    /** @brief Return the velocity solution vector */
    virtual LA_SeqVector const *get_solution_velocity(void) const;

    /** @brief Return the pressure solution vector */
    virtual LA_SeqVector const *get_solution_pressure(void) const;

    /** @brief Return the pressure Dirichlet BC vector */
    virtual LA_Vector *get_pressure_DirichletBC_vector(void);

    /** @brief Return the periodic pressure drop vector */
    virtual LA_Vector *get_unitary_periodic_pressure_drop_vector(void);
    //@}

    /** @name Translation projection from peligriff */
    //@{

    virtual void set_velocity_unknown(size_t i_row, double xx);

    virtual void set_pressure_unknown(size_t i_row, double xx);

    virtual void synchronize_velocity_unknown_vector();

    virtual void synchronize_pressure_unknown_vector();

    virtual void nullify_DLMFD_Nm1_rhs();

    virtual void set_rhs_DLMFD_Nm1(size_t i_row, double xx);

    virtual LA_SeqVector const *get_rhs_DLMFD_Nm1();

    virtual void synchronize_rhs_DLMFD_Nm1_vector();

    //@}

    virtual void synchronize_rhs_periodic_pressure_vector();

    /** @name Set */
    //@{

    virtual void set_periodic_pressure_rhs_item(size_t i_row, double xx);

    //@}

    //-- Basic operations on matrices & vectors

    /** @name Basic operations on matrices & vectors */
    //@{
    /** @brief Initialize the velocity unknown vector with field values */
    virtual void initialize_velocity(void);

    /** @brief Initialize the pressure unknown vector with field values */
    virtual void initialize_pressure(void);

    /** @brief Store velocity vector at previous time step */
    virtual void at_each_time_step(void);

    /** @brief Compute velocity change from one time step to the next one */
    virtual double compute_velocity_change(void);

    /** @brief Compute velocity divergence norm */
    virtual double compute_velocity_divergence_norm(void);

    /** @brief Nullify velocity advection (inertia) right hand side */
    virtual void nullify_velocity_advection_rhs(void);

    /** @brief Assemble the velocity advection rhs term coef*u*grad(u)
    @param AdvectionScheme advection scheme ("Upwind" or "TVD")
    @param advecting_level level of advecting velocity
    @param coef prefactor
    @param advected_level level of advected velocity */
    virtual void assemble_velocity_advection(string const &AdvectionScheme,
                                             size_t advecting_level, double const &coef, size_t advected_level);

    /** @brief Finalize constant matrices */
    virtual void finalize_constant_matrices(void);

    /** @brief Compute second order velocity advection-diffusion rhs
    @param restart if the simulation is a restart
    @param iteration_number iteration number
    @param b_with_advection add advection term
    @param dpdl periodic pressure gradient source term (is 0 if no periodic
    pressure drop is imposed) */
    virtual void compute_velocityAdvectionDiffusion_rhs(bool const &b_restart,
                                                        size_t const &iteration_number,
                                                        bool const &b_with_advection,
                                                        double const &dpdl);

    /** @brief Assemble velocity viscous matrix and rhs
    @param coef_lap laplacian coefficient */
    virtual void assemble_velocity_viscous_matrix_rhs(double const &coef_lap);

    /** @brief Assemble velocity unsteady matrix
    @param coef_lap mass coefficient */
    virtual void assemble_velocity_unsteady_matrix(double const &coef);

    /** @brief Assemble velocity divergence matrix and rhs
    @param coef coefficient (either -1 or 1 ) */
    virtual void assemble_pdivv_matrix_rhs(double const &coef);

    /** @brief Assemble pressure laplacian matrix and rhs
    @param coef_lap laplacian coefficient */
    virtual void assemble_pressure_laplacian_matrix_rhs(double const &coef_lap);

    /** @brief Correct the pressure laplacian operator in case of all non
    Dirichlet BCs i.e. add coefdiag to one diagonal coefficient to ensure that
    pressure is 0 at this node */
    virtual void pressure_laplacian_correction(void);

    /** @brief Store u.grad(u)^(n-2) at first sub-iteration when CFL condition
    is not met and 2nd order scheme degenerates to 1st order */
    virtual void store_ugradu_Nm2(size_t const &n_advection_subtimesteps);
    //@}

    //-- Solver

    /** @name Solvers */
    //@{
    /** @brief Velocity Diffusion (viscous) solver with velocity advection
    in the rhs */
    virtual bool VelocityDiffusion_solver(void);

    /** @brief Velocity-pressure correction solver: (i) solve a pressure
    Poisson problem and (ii) update velocity and pressure
    @param density fluid density
    @param viscosity fluid viscosity
    @param timestep time step magnitude */
    virtual double VelocityPressure_correction_solver(
        double const &density, double const &viscosity,
        double const &timestep);

    /** @brief Velocity Advection solver */
    virtual bool VelocityAdvection_solver(void);
    //@}

    //-- Input / Output methods

    /** @name Input / Output methods */
    //@{
    //@}

    //-- Persistence

    // Add the vectors to be stored for restart
    virtual void add_storable_objects(MAC_ListIdentity *list) const;

    // --------------------------------------------------------------------------
    // --------------------------------------------------------------------------
    // --------------------------------------------------------------------------
    // --------------------------- DLMFD FRAMEWORK ------------------------------

    /** @brief Initialize the velocity right hand side of momentum equations
    in the DLM/FD problem i.e. compute (ro/dt)*U(n-1) before the DLM/FD
    solution */
    virtual void updateFluid_DLMFD_rhs();

    /** @brief Nullify q vector */
    virtual void nullify_QUvector();

    /** @brief Add the contribution coef*<lambda,v> to the q
    vector i.e. the DLM right hand side of momentum equations in the DLM/FD
    problem
    @param transferVec the entry
    @param index the position in the vector
    @param coef parameter */
    virtual void assemble_inQUvector(double transferVal, size_t index, double coef);

    /** @brief Solve the fluid system at the matrix level */
    virtual void solve_FluidVel_DLMFD_Init(const double &time);

    virtual LA_SeqVector const *get_solution_U() const;

    virtual double const get_DLMFD_convergence_criterion() const;

    virtual int const get_DLMFD_maxiter() const;

    /** @brief Initialize Qu vector as Bt*w, with w the pressure descent direction */
    virtual void initialize_QUvector_with_divv_rhs();

    /** @brief solve A.tu = quf = <w,v> */
    virtual void solve_FluidVel_DLMFD_Iter(const double &time);

    /** @brief Return the velocity solution vector t */
    virtual LA_SeqVector const *get_tVector_U() const;

    /** @brief Update u+=alpha.t */
    virtual void update_FluidVel_OneUzawaIter(const double &alpha);

    virtual void store_DLMFD_rhs();

    virtual void re_initialize_explicit_DLMFD(bool const &restart, string const &rootfilename_dlm);

    /** @brief Reload other data than FV fields for restart  */
    virtual void do_additional_savings(string const &rootfilename_dlm);

private:   //-----------------------------------------------------------------
protected: //-----------------------------------------------------------------
};

#endif
