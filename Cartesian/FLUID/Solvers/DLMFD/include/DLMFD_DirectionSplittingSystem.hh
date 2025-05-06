#ifndef REG_HEAT_EQUATION_SYSTEM_HH
#define REG_HEAT_EQUATION_SYSTEM_HH
#include <MAC_Object.hh>
#include <DLMFD_System.hh>
#include <utility>
#include <boolVector.hh>
using namespace std;

class MAC_ModuleExplorer;
class MAC_Communicator;
class MAC_Timer;
class size_t_vector;
class intVector;
class doubleVector;
class LA_Matrix;
class LA_SeqMatrix;
class LA_Vector;
class LA_SeqVector;
class LA_Scatter;
class LA_Solver;
class LA_CRSmatrix;
class FV_SystemNumbering;
class FV_DiscreteField;
class FV_TimeIterator;

/** @brief TDMatrix include all elements of block matrices (ii,ie,ei,ee) */
struct TDMatrix
{
    LA_SeqVector ***ii_main;
    LA_SeqVector ***ii_super;
    LA_SeqVector ***ii_sub;
    LA_SeqMatrix ***ie;
    LA_SeqMatrix ***ei;
    LA_SeqMatrix ***ee;
};
/** @brief Product matrix is composed of products elements of block matrices (ii,ie,ei,ee) */
struct ProdMatrix
{
    LA_SeqMatrix **ei_ii_ie;
    LA_SeqVector **ii_ie;
    LA_SeqVector **result;
};
/** @brief LocalVector to be used in storing the local values and interface values of DOF */
struct LocalVector
{
    LA_SeqVector **local_T;
    LA_SeqVector **local_solution_T;
    LA_SeqVector **T;
    LA_SeqVector **interface_T;
};
/** @brief The Class DLMFD_DirectionSplittingSystem.
Matrix systems for the resolution of the heat equation.
@author A. Wachs - Pacific project 2017 */
class DLMFD_DirectionSplittingSystem : public MAC_Object, public DLMFD_System
{
private: //----------------------------------------------------------
         //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{
    /** @brief Constructor without argument */
    DLMFD_DirectionSplittingSystem(void);
    /** @brief Destructor */
    ~DLMFD_DirectionSplittingSystem(void);

    /** @brief Copy constructor */
    DLMFD_DirectionSplittingSystem(DLMFD_DirectionSplittingSystem const &other);

    /** @brief Operator ==
    @param other the right hand side */
    DLMFD_DirectionSplittingSystem &operator=(DLMFD_DirectionSplittingSystem const &other);

    /** @brief Constructor with arguments
    @param a_owner the MAC-based object
    @param exp to read the data file
    @param mac_UF FV velocity field */
    DLMFD_DirectionSplittingSystem(MAC_Object *a_owner,
                                   MAC_ModuleExplorer const *exp,
                                   FV_DiscreteField *mac_UF,
                                   FV_DiscreteField *mac_PF);
    //@}

public: //-----------------------------------------------------------
        //-- Instance delivery and initialization
    /** @name Instance delivery and initialization */
    //@{
    /** @brief Create and initialize an instance of DLMFD_DirectionSplittingSystem
    @param a_owner the MAC-based object
    @param exp to read the data file
    @param mac_UF FV velocity field */
    static DLMFD_DirectionSplittingSystem *create(MAC_Object *a_owner,
                                                  MAC_ModuleExplorer const *exp,
                                                  FV_DiscreteField *mac_UF,
                                                  FV_DiscreteField *mac_PF);
    //@}

    //-- Access

    /** @name Access */
    //@{
    /** @brief Return the DS velocity solution vector DS_UF */
    LA_SeqVector const *get_solution_DS_velocity(void) const;

    /** @brief Return the DS pressure solution vector DS_PF */
    LA_SeqVector const *get_solution_DS_pressure(void) const;

    /** @brief Return the matrix system of spacial discretization */
    TDMatrix *get_A(size_t const &field);
    /** @brief Return the Schur complement of spacial discretization */
    TDMatrix *get_Schur(size_t const &field);
    /** @brief Return the Schur complement of Schur complement in case of periodic domain */
    TDMatrix *get_DoubleSchur(size_t const &field);
    /** @brief Return the product matrix of Schur complement */
    ProdMatrix *get_SchurP(size_t const &field);
    /** @brief Return RHS for the Schur complement */
    LocalVector *get_Schur_VEC(size_t const &field);
    /** @brief Return the product matrix of spacial discretization */
    ProdMatrix *get_Ap(size_t const &field);
    /** @brief Return RHS for the matrix system of spacial discretization */
    LocalVector *get_VEC(size_t const &field);

    //-- Basic operations on matrices & vectors

    /** @name Basic operations on matrices & vectors */
    //@{
    /** @brief Initialize the velocity unknown vector with field values */
    void initialize_DS_velocity(void);
    /** @brief Initialize the pressure unknown vector with field values */
    void initialize_DS_pressure(void);

    /** @brief Store velocity vector at previous time step */
    void at_each_time_step(void);
    /** @brief Compute velocity change from one time step to the
    next one with the direction splitting solution method */
    double compute_DS_velocity_change(void);

    //-- Solver

    /** @name Solvers */
    //@{
    /** @brief Solve the DS splitting problem in x by performing the
    matrix-vector product A_x^-1.Vx and transfer in the distributed vector */
    void DS_NavierStokes_solver(FV_DiscreteField *FF, size_t const &j, size_t const &k, size_t const &min_i, size_t const &comp, size_t const &dir, size_t const &field, size_t const &r_index);

    //@}

    /** @brief Synchronize the DS solution vector for velocity*/
    void synchronize_DS_solution_vec(void);
    /** @brief Synchronize the DS solution vector for pressure*/
    void synchronize_DS_solution_vec_P(void);

    //-- Output methods

    /** @name Output methods */
    //@{
    //@}

    /** @brief Calls interior function for different conditions to compute the product of Aei*inv(Aii)*Aie */
    void compute_product_matrix(struct TDMatrix *arr, struct ProdMatrix *prr, size_t const &comp, size_t const &dir, size_t const &field, size_t const &r_index);
    /** @brief Compute the product of Aei*inv(Aii)*Aie in any direction for any field*/
    void compute_product_matrix_interior(struct TDMatrix *arr, struct ProdMatrix *prr, size_t const &comp, size_t const &column, size_t const &dir, size_t const &r_index);

    //-- Utilities

    /** @name Utilities */
    //@{
    /** @brief Solve Linear system mat_A*x = rhs with only three vectors of mat_A(x,y,z) using thomas algorithm  */
    void mod_thomas_algorithm(TDMatrix *arr, LA_SeqVector *rhs, size_t const &comp, size_t const &dir, size_t const &r_index);
    /** @brief Compute the modified super diagonal for thomas algorithm  */
    void pre_thomas_treatment(size_t const &comp, size_t const &dir, struct TDMatrix *arr, size_t const &r_index);
    //@}

    // --------------------------------------------------------------------------
    // --------------------------------------------------------------------------
    // --------------------------------------------------------------------------
    // --------------------------- DLMFD FRAMEWORK ------------------------------

    /** @brief Finalize constant matrices */
    void finalize_constant_matrices();

    /** @brief Assemble velocity unsteady matrix for DLMFD workflow
    @param coef_lap mass coefficient */
    void assemble_velocity_unsteady_matrix(double const &coef);

    /** @brief Initialize the velocity right hand side of momentum equations
    in the DLM/FD problem i.e. compute (ro/dt)*U(n-1) before the DLM/FD
    solution */
    void updateFluid_DLMFD_rhs();

    /** @brief Nullify q vector */
    void nullify_QUvector();

    /** @brief Nullify Explicit DLMFD vector */
    void nullify_Explicit_DLMFD_Cvector();

    /** @brief Add the contribution coef*<lambda,v> to the q
    vector i.e. the DLM right hand side of momentum equations in the DLM/FD
    problem
    @param transferVec the entry
    @param index the position in the vector
    @param coef parameter */
    void assemble_inQUvector(double transferVal, size_t index, double coef);

    /** @brief Assemble the explicit DLMFD C vector
    @param transferVec the entry
    @param index the position in the vector
    @param coef parameter */
    void assemble_inExplicit_DLMFD_Cvector(double transferVal, size_t index, double coef);

    /** @brief Solve the fluid system at the matrix level */
    void solve_FluidVel_DLMFD_Init(const double &time);

    LA_SeqVector const *get_solution_U() const;

    /** @brief solve A.tu = quf = <w,v> */
    void solve_FluidVel_DLMFD_Iter(const double &time);

    /** @brief Return the velocity solution vector t */
    LA_SeqVector const *get_tVector_U() const;

    /** @brief Update u+=alpha.t */
    void update_FluidVel_OneUzawaIter(const double &alpha);

    void store_DLMFD_rhs();

    void re_initialize_explicit_DLMFD(bool const &restart);

    double get_explicit_DLMFD_at_index(size_t index) const;

protected: //--------------------------------------------------------
private:   //----------------------------------------------------------
    /** @name Initialize matrices & vectors */
    //@{
    /** @brief Create matrices & vectors (without allocating memory)
    @param exp to read the data file */
    void build_system(MAC_ModuleExplorer const *exp);

    /** @brief Allocate memory for matrices & vectors */
    void re_initialize(void);
    //@}

    //-- Attributes
    FV_DiscreteField *UF;
    FV_DiscreteField *PF;

    // Local vectors
    LA_SeqVector *UF_DS_LOC;
    LA_SeqVector *PF_DS_LOC;
    LA_SeqVector *T_LOC;

    // Global velocity solution vectors
    LA_Vector *VEC_DS_UF;
    LA_Vector *VEC_DS_UF_previoustime;
    LA_Vector *VEC_DS_UF_timechange;

    // Global pressure solution vectors
    LA_Vector *VEC_DS_PF;

    // Work vectors
    LA_Vector *VEC_q;
    double *vector_rhs_VelocityDLMFD_Nm1;
    LA_Vector *VEC_t;
    LA_Vector *VEC_r;
    LA_Vector *VEC_w;

    // Solvers
    LA_Solver *SOLVER_A_VelocityUnsteady;

    // Matrices & rhs
    LA_Matrix *MAT_D_velocityUnsteadyPlusDiffusion;
    LA_Matrix *MAT_A_VelocityUnsteady;

    LA_Vector *VEC_rhs_D_velocityDiffusionPlusBodyTerm;
    LA_Vector *VEC_rhs_A_Velocity;

    // Unknowns numbering
    FV_SystemNumbering *UF_NUM;
    FV_SystemNumbering *PF_NUM;

    // Direction splitting matrices
    LA_SeqMatrix *MAT_velocityUnsteadyPlusDiffusion_1D;

    // Spacitial discretization matrices
    struct TDMatrix A[2][3];      // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
    struct ProdMatrix Ap[2][3];   // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
    struct LocalVector VEC[2][3]; // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions

    // Schur complement matrices
    struct TDMatrix Schur[2][3]; // [0,1] are for pressure and velocity;[0,1,2] are for x, y and z directions
    struct ProdMatrix SchurP[2][3];
    struct LocalVector Schur_VEC[2][3];

    // Schur complement of Schur complement
    struct TDMatrix DoubleSchur[2][3];

    size_t dim;
    MAC_Communicator const *pelCOMM;
    size_t nb_comps[2];

    /** Processor positions in x,y,z */
    size_t proc_pos_in_i[3];
    /** Number of Processors in x,y,z */
    size_t nb_procs_in_i[3];

    bool is_periodic[2][3];
    boolVector const *U_periodic_comp;
    boolVector const *P_periodic_comp;

    //-- DLMFD objects

    // Explicit treatment
    LA_Vector *VEC_rhs_VelocityDLMFD_Nm1;
    bool b_NS_ExplicitDLMFD;
};

#endif
