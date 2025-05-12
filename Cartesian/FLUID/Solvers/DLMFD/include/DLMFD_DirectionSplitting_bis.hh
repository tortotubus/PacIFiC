#ifndef DLMFD_DirectionSplitting_bis_HH
#define DLMFD_DirectionSplitting_bis_HH

#include <mpi.h>
#include <FV_OneStepIteration.hh>
#include <geomVector.hh>
#include <PAC_computingtime.hh>
#include <boolVector.hh>
#include <PAC_solvercomputingtime.hh>
#include <MAC_DoubleVector.hh>
#include <DS_PID.hh>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <DLMFD_FictitiousDomain.hh>
using namespace std;

class MAC_Communicator;
class FV_DiscreteField;
class LA_Vector;
class LA_SeqVector;
class DLMFD_DirectionSplittingSystem_bis;
class LA_SeqMatrix;

/** @brief The Class DLMFD_DirectionSplitting_bis.

Server for the resolution of the unsteady momentum equation by a second order
time integrator and a Finite Volume MAC scheme on rectangular grids.

@author A. Goyal - Pacific project 2025 */

/** @brief MPIVar include all vectors required while message passing */
struct MPIVarNS
{
    size_t *size;
    double ***send;
    double ***receive;
};

class DLMFD_DirectionSplitting_bis : public FV_OneStepIteration, public PAC_ComputingTime, public PAC_SolverComputingTime
{
public: //-----------------------------------------------------------------
        //-- Substeps of the step by step progression
    /** @name Substeps of the step by step progression */
    //@{
    /** @brief Tasks performed at initialization of the algorithm, before
    starting the time stepping loop
    @param t_it time iterator */
    virtual void do_before_time_stepping(FV_TimeIterator const *t_it,
                                         std::string const &basename);

    /** @brief Perform one time step
    @param t_it time iterator */
    virtual void do_one_inner_iteration(FV_TimeIterator const *t_it);

    /** @brief Tasks performed at initialization of each time step
    @param t_it time iterator */
    virtual void do_before_inner_iterations_stage(
        FV_TimeIterator const *t_it);

    /** @brief Tasks performed after of each time step
    @param t_it time iterator */
    virtual void do_after_inner_iterations_stage(
        FV_TimeIterator const *t_it);

    /** @brief Tasks performed at the end of the time stepping loop */
    virtual void do_after_time_stepping(void);

    /** @brief Save additional data than fields
    @param t_it time iterator
    @param cycleNumber cycle number */
    virtual void do_additional_savings(FV_TimeIterator const *t_it,
                                       int const &cycleNumber);
    //@}

    //-- Projection-Translation methods

    /** @name Projection-Translation methods */
    //@{
    /** @brief Projection of the field on the translated position of the grid
     */
    void fields_projection();

    //@}

protected: //--------------------------------------------------------------
private:   //----------------------------------------------------------------
    static DLMFD_DirectionSplitting_bis const *PROTOTYPE;

    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{
    /** @brief Destructor */
    virtual ~DLMFD_DirectionSplitting_bis();

    /** @brief Copy constructor */
    DLMFD_DirectionSplitting_bis(DLMFD_DirectionSplitting_bis const &other);

    /** @brief Constructor with arguments
    @param a_owner the MAC-based object
    @param exp to read the data file */
    DLMFD_DirectionSplitting_bis(MAC_Object *a_owner,
                             FV_DomainAndFields const *dom,
                             MAC_ModuleExplorer const *exp);

    /** @brief Constructor without argument */
    DLMFD_DirectionSplitting_bis();

    /** @brief Create a clone
    @param a_owner the MAC-based object
    @param dom mesh and fields
    @param prms set of parameters
    @param exp to read the data file */
    virtual DLMFD_DirectionSplitting_bis *create_replica(MAC_Object *a_owner,
                                                     FV_DomainAndFields const *dom,
                                                     MAC_ModuleExplorer *exp) const;

    //@}

    //-- Basic discrete system building

    /** @name Basic discrete system building */
    //@{
    /** @brief Call the function to assemble 1D matrices
    for both velocity and pressure field*/
    void assemble_1D_matrices(FV_DiscreteField const *FF, FV_TimeIterator const *t_it);

    /** @brief Assemble 1D matrices for both velocity
    and pressure field in all directions */
    void assemble_field_matrix(FV_DiscreteField const *FF, FV_TimeIterator const *t_it);

    /** @brief Assemble 1D schur matrices for both
    velocity and pressure field in all directions */
    void assemble_field_schur_matrix(FV_DiscreteField const *FF);

    /** @brief Assemble advection term for Upwind spacial scheme */
    double assemble_advection_Upwind(size_t const &advecting_level, double const &coef, size_t const &advected_level, size_t const &i, size_t const &j, size_t const &k, size_t const &component);

    /** @brief Assemble advection term for Centered spacial scheme */
    double assemble_advection_Centered(size_t const &advecting_level, double const &coef, size_t const &advected_level, size_t const &i, size_t const &j, size_t const &k, size_t const &component);

    /** @brief Assemble advection term for Centered spacial scheme with FD corrections*/
    double assemble_advection_Centered_FD(double const &coef, size_t const &i, size_t const &j, size_t const &k, size_t const &comp, size_t const &level);

    // /** @brief Assemble advection term for Centered spacial scheme with CutCell corrections */
    // double assemble_advection_Centered_CutCell(FV_TimeIterator const *t_it, double const &coef, size_t const &i, size_t const &j, size_t const &k, size_t const &comp, size_t const &level);

    /** @brief Assemble advection term for TVD spacial scheme */
    double assemble_advection_TVD(size_t const &advecting_level, double const &coef, size_t const &advected_level, size_t const &i, size_t const &j, size_t const &k, size_t const &component) const;

    /** @brief Assemble rhs for velocity in any direction */
    double velocity_local_rhs(size_t const &j, size_t const &k, double const &gamma, FV_TimeIterator const *t_it, size_t const &comp, size_t const &dir);

    /** @brief Compute first derivative */
    std::tuple<double, double>
    compute_first_derivative(size_t const &comp, size_t const &i, size_t const &j, size_t const &k, size_t const &dir, size_t const &level);

    /** @brief Compute diffusive term of velocity field from previous timestep */
    double compute_un_component(size_t const &comp, size_t const &i, size_t const &j, size_t const &k, size_t const &dir, size_t const &level);
    /** @brief Compute diffusive term of velocity field from previous timestep */
    double compute_un_component_FD(size_t const &comp, size_t const &i, size_t const &j, size_t const &k, size_t const &dir, size_t const &level);
    // /** @brief Compute diffusive term of velocity field from previous timestep */
    // double compute_un_component_FV(size_t const &comp, size_t const &i, size_t const &j, size_t const &k, size_t const &dir, size_t const &level);
    /** @brief Compute diffusive term of pressure field from previous timestep */
    double compute_p_component(size_t const &comp, size_t const &i, size_t const &j, size_t const &k);
    /** @brief Compute advective term based on either Upwind or TVD spacial scheme */
    double compute_adv_component(FV_TimeIterator const *t_it, size_t const &comp, size_t const &i, size_t const &j, size_t const &k);
    /** @brief Assemble rhs term after calling compute_**_component */
    void assemble_DS_un_at_rhs(FV_TimeIterator const *t_it, double const &gamma);

    void assemble_velocity_diffusion_terms();

    void assemble_velocity_advection_terms(FV_TimeIterator const *t_it);

    void calculate_row_indexes(FV_DiscreteField const *FF);

    void compute_velocity_divergence(FV_DiscreteField const *FF);

    double divergence_of_U_noCorrection(size_t const &i, size_t const &j, size_t const &k, size_t const &comp, size_t const &level);

    /** @brief Call functions to assemble rhs for pressure or velocity fields in any direction */
    double assemble_local_rhs(size_t const &j, size_t const &k, double const &gamma, FV_TimeIterator const *t_it, size_t const &comp, size_t const &dir, size_t const &field);

    /** @brief Assemble divergence for penelty step */
    double calculate_velocity_divergence_FD(size_t const &i, size_t const &j, size_t const &k, size_t const &level);

    // /** @brief Assemble divergence for penelty step */
    // double calculate_velocity_divergence_cutCell(size_t const &i, size_t const &j, size_t const &k, size_t const &level);

    double assemble_velocity_gradients(class doubleVector &grad, size_t const &i, size_t const &j, size_t const &k, size_t const &level);

    // void detect_fresh_cells_and_neighbours();

    void calculate_divergence_weighting(FV_TimeIterator const *t_it);

    double pressure_local_rhs(size_t const &j, size_t const &k, FV_TimeIterator const *t_it, size_t const &dir);

    // /** @brief Initialize the velocity on the velocity nodes in MAC grid*/
    // void initialize_grid_nodes_on_rigidbody(vector<size_t> const &list);

    void ugradu_initialization();

    // void correct_pressure_1st_layer_solid(size_t const &level);

    // void correct_pressure_2nd_layer_solid(size_t const &level);

    void correct_mean_pressure(size_t const &level);

    void predicted_pressure_drop(FV_TimeIterator const *t_it);

    double get_current_mean_flow_speed();

    /** @brief Solve interface unknowns for
    both fields in any particular direction */
    void solve_interface_unknowns(FV_DiscreteField *FF, double const &gamma, FV_TimeIterator const *t_it, size_t const &comp, size_t const &dir, size_t const &level);
    /** @brief Unpack the interface variable sent
    by master processor to slave processor */
    void unpack_ue(size_t const &comp, double *received_data, size_t const &dir, size_t const &p, size_t const &field);
    /** @brief Unpack the data sent by "data_packing" and compute
    the interface unknown; and pack ue for sending to slave processor */
    void unpack_compute_ue_pack(size_t const &comp, size_t const &dir, size_t const &p, size_t const &field);

    /** @brief Pack Aei*(Aii)-1*fi and fe for sending to master processor */
    void data_packing(FV_DiscreteField const *FF, size_t const &p, double const &fe, size_t const &comp, size_t const &dir);
    /** @brief Compute Aei*(Aii)-1*fi required to compute interface unknown */
    void compute_Aei_ui(struct TDMatrix *arr, struct LocalVector *VEC, size_t const &comp, size_t const &dir, size_t const &r_index);

    /** @brief Solve interface unknown for all cases */
    void DS_interface_unknown_solver(LA_SeqVector *rhs, size_t const &comp, size_t const &dir, size_t const &field, size_t const &r_index);

    /** @brief Solve i in j and k; e.g. solve x in y ank z */
    void Solve_i_in_jk(FV_DiscreteField *FF, FV_TimeIterator const *t_it, size_t const &dir_i, size_t const &dir_j, size_t const &dir_k, double const &gamma, size_t const &level);

    /** Pressure predictor */
    void NS_first_step(FV_TimeIterator const *t_it);

    /** Velocity update */
    void NS_velocity_update(FV_TimeIterator const *t_it);

    /** Pressure update (penalty step) */
    void NS_pressure_update(FV_TimeIterator const *t_it);

    /** Pressure correction */
    void NS_final_step(FV_TimeIterator const *t_it);

    void write_output_field(FV_DiscreteField const *FF);

    void output_L2norm_divergence();

    void output_L2norm_pressure(size_t const &level);

    void output_L2norm_velocity(size_t const &level);

    doubleVector compute_outOfDomain_Pressure(size_t const &level);

    /** @brief Compute velocity change from one time step to the
    next one with the direction splitting solution method */
    double compute_DS_velocity_change(void);
    //@}

    //-- Direction splitting communicators

    /** @name Direction splitting communicators */
    /** @brief Create the sub-communicators */
    void create_DS_subcommunicators(void);
    void processor_splitting(int const &color, int const &key, size_t const &dir);

    void allocate_mpi_variables(FV_DiscreteField const *FF);
    void deallocate_mpi_variables();

    /** @brief Free the sub-communicators */
    void free_DS_subcommunicators(void);

    //@}

    //-- Projection-Translation methods

    /** @name Projection-Translation methods */
    //@{

    /** @brief Build the translation vector */
    void set_translation_vector();

    /** @brief Build the field projection-translation interpolations */
    void build_links_translation();
    //@}

private: //----------------------------------------------------------------
         //-- Class attributes
    //-- Attributes

    FV_DiscreteField *UF;
    FV_DiscreteField *PF;

    DLMFD_DirectionSplittingSystem_bis *GLOBAL_EQ;
    DS_PID *controller;

    size_t sub_prob_number;
    double imposed_CFL;

    size_t nb_procs;
    size_t my_rank;
    size_t is_master;
    size_t dim;
    // 0th element for P and 1st element for U
    size_t nb_comps[2];

    MAC_Communicator const *macCOMM;
    MPI_Comm DS_Comm_i[3];

    int rank_in_i[3];
    int nb_ranks_comm_i[3];

    // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
    struct MPIVarNS first_pass[2][3];
    // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
    struct MPIVarNS second_pass[2][3];
    // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
    struct MPIVarNS data_for_S[2][3];

    double mu;
    double kai;
    string AdvectionScheme;
    string StencilCorrection;
    bool is_CConlyDivergence;
    double FluxRedistThres;
    size_t AdvectionTimeAccuracy;
    double rho;
    bool b_restart;
    bool is_solids;

    bool is_stressCal;
    string ViscousStressOrder;
    string PressureStressOrder;
    double surface_cell_scale;
    size_t stressCalFreq;
    bool is_par_motion;

    // Grid motion
    bool b_projection_translation;
    int outOfDomain_boundaryID;
    FV_Mesh const *primary_grid;
    double critical_distance_translation;
    geomVector MVQ_translation_vector;
    size_t translation_direction;
    double bottom_coordinate;
    double translated_distance;

    double Qold;

    boolVector const *P_periodic_comp;
    boolVector const *U_periodic_comp;
    MAC_DoubleVector *gravity_vector;
    geomVector split_gravity_vector;
    geomVector gravity_vector_geom;
    bool is_periodic[2][3];

    bool exceed, turn;

    // DLMFD solver
    DLMFD_FictitiousDomain *dlmfd_solver;
    bool b_ExplicitDLMFD;

    // Post-processing
    string resultsDirectory;
};

#endif
