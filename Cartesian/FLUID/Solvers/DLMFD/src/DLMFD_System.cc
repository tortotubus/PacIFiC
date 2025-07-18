#include <DLMFD_System.hh>
#include <MAC_Error.hh>

//----------------------------------------------------------------------
DLMFD_System::DLMFD_System(MAC_Object *a_owner, MAC_ModuleExplorer const *exp,
                           FV_DiscreteField *mac_UF, FV_DiscreteField *mac_PF,
                           size_t const &NS_Viscous_TimeAccuracy_,
                           size_t const &NS_Advection_TimeAccuracy_,
                           bool const &b_pressure_rescaling_,
                           bool const &b_ExplicitPressureGradient_,
                           bool const &b_HighOrderPressureCorrection_,
                           bool is_stressCal_)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: DLMFD_System");

    // DLMFD Uzawa convergence criterions
    Uzawa_DLMFD_precision = 1.e-5;
    if (exp->has_entry("Uzawa_DLMFD_precision"))
        Uzawa_DLMFD_precision = exp->double_data("Uzawa_DLMFD_precision");

    Uzawa_DLMFD_maxiter = 100;
    if (exp->has_entry("Uzawa_DLMFD_maxiter"))
        Uzawa_DLMFD_maxiter = exp->int_data("Uzawa_DLMFD_maxiter");
}




//----------------------------------------------------------------------
DLMFD_System::DLMFD_System()
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: DLMFD_System");
}




//----------------------------------------------------------------------
DLMFD_System::~DLMFD_System()
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: DLMFD_System");
}




//----------------------------------------------------------------------
void DLMFD_System::at_each_time_step(void)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: at_each_time_step");

    string error_message = "DLMFD_System::at_each_time_step ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
double DLMFD_System::compute_velocity_change(void)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: compute_velocity_change");

    string error_message = "DLMFD_System::compute_velocity_change ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
double DLMFD_System::compute_velocity_divergence_norm(void)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: compute_velocity_divergence_norm");

    string error_message = "DLMFD_System::compute_velocity_divergence_norm ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::nullify_velocity_advection_rhs(void)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: nullify_velocity_advection_rhs");

    string error_message = "DLMFD_System::nullify_velocity_advection_rhs ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
LA_SeqVector const *DLMFD_System::get_solution_velocity(void) const
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: get_solution_U");

    string error_message = "DLMFD_System::get_solution_velocity ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
LA_SeqVector const *DLMFD_System::get_solution_pressure(void) const
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: get_solution_pressure");

    string error_message = "DLMFD_System::get_solution_pressure ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::initialize_velocity(void)
//----------------------------------------------------------------------
{
    MAC_LABEL("MAC_NavierStokesSystem:: initialize_velocity");

    string error_message = "DLMFD_System::initialize_velocity ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::initialize_pressure(void)
//----------------------------------------------------------------------
{
    MAC_LABEL("MAC_NavierStokesSystem:: initialize_pressure");

    string error_message = "DLMFD_System::initialize_pressure ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::finalize_constant_matrices(void)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: finalize_constant_matrices");

    string error_message = "DLMFD_System::finalize_constant_matrices ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::compute_velocityAdvectionDiffusion_rhs(
    bool const &b_restart, size_t const &iteration_number,
    bool const &b_with_advection, double const &dpdl)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: "
              "compute_velocityAdvectionDiffusion_rhs");

    string error_message =
        "DLMFD_System::compute_velocityAdvectionDiffusion_rhs ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::assemble_velocity_viscous_matrix_rhs(double const &coef_lap)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: assemble_velocity_viscous_matrix_rhs");

    string error_message =
        "DLMFD_System::assemble_velocity_viscous_matrix_rhs ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::assemble_velocity_unsteady_matrix(double const &coef)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: assemble_velocity_unsteady_matrix");

    string error_message = "DLMFD_System::assemble_velocity_unsteady_matrix ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::assemble_pdivv_matrix_rhs(double const &coef)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: assemble_pdivv_matrix_rhs");

    string error_message = "DLMFD_System::assemble_pdivv_matrix_rhs ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::assemble_velocity_advection(string const &AdvectionScheme,
                                               size_t advecting_level,
                                               double const &coef,
                                               size_t advected_level)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: assemble_velocity_advection");

    string error_message = "DLMFD_System::assemble_velocity_advection ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::assemble_pressure_laplacian_matrix_rhs(
    double const &coef_lap)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: assemble_pressure_laplacian_matrix_rhs");

    string error_message =
        "DLMFD_System::assemble_pressure_laplacian_matrix_rhs ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::pressure_laplacian_correction(void)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: pressure_laplacian_correction");

    string error_message = "DLMFD_System::pressure_laplacian_correction ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::assemble_inExplicit_DLMFD_Cvector(double transferVal,
                                                     size_t index, double coef)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::assemble_inExplicit_DLMFDvector");
}




//----------------------------------------------------------------------
bool DLMFD_System::VelocityDiffusion_solver(void)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: VelocityDiffusion_solver");

    string error_message = "DLMFD_System::VelocityDiffusion_solver ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
double DLMFD_System::VelocityPressure_correction_solver(double const &density,
                                                        double const &viscosity,
                                                        double const &timestep)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: VelocityPressure_correction_solver");

    string error_message = "DLMFD_System::VelocityPressure_correction_solver ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
LA_Vector *DLMFD_System::get_pressure_DirichletBC_vector(void)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: get_pressure_DirichletBC_vector");

    string error_message = "DLMFD_System::get_pressure_DirichletBC_vector ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
LA_Vector *DLMFD_System::get_unitary_periodic_pressure_drop_vector(void)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: "
              "get_unitary_periodic_pressure_drop_vector");

    string error_message =
        "DLMFD_System::get_unitary_periodic_pressure_drop_vector ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::set_velocity_unknown(size_t i_row, double xx)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::set_velocity_unknown");

    string error_message = "DLMFD_System::set_velocity_unknown ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::set_pressure_unknown(size_t i_row, double xx)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::set_pressure_unknown");

    string error_message = "DLMFD_System::set_pressure_unknown ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::synchronize_velocity_unknown_vector()
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::synchronize_velocity_unknown_vector");

    string error_message = "DLMFD_System::synchronize_velocity_unknown_vector ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::synchronize_pressure_unknown_vector()
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::synchronize_pressure_unknown_vector");

    string error_message = "DLMFD_System::synchronize_pressure_unknown_vector ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::nullify_DLMFD_Nm1_rhs()
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::nullify_DLMFD_Nm1_rhs");

    string error_message = "DLMFD_System::nullify_DLMFD_Nm1_rhs ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::set_rhs_DLMFD_Nm1(size_t i_row, double xx)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::set_rhs_DLMFD_Nm1");

    string error_message = "DLMFD_System::set_rhs_DLMFD_Nm1 ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
LA_SeqVector const *DLMFD_System::get_rhs_DLMFD_Nm1()
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::get_rhs_DLMFD_Nm1");

    string error_message = "DLMFD_System::get_rhs_DLMFD_Nm1 ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::synchronize_rhs_DLMFD_Nm1_vector()
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::synchronize_rhs_DLMFD_Nm1_vector");

    string error_message = "DLMFD_System::synchronize_rhs_DLMFD_Nm1_vector ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::synchronize_rhs_periodic_pressure_vector()
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::synchronize_rhs_periodic_pressure_vector");

    string error_message =
        "DLMFD_System::synchronize_rhs_periodic_pressure_vector ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::set_periodic_pressure_rhs_item(size_t i_row, double xx)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::set_periodic_pressure_rhs_item");

    string error_message = "DLMFD_System::set_periodic_pressure_rhs_item ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::store_ugradu_Nm2(size_t const &n_advection_subtimesteps)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: store_ugradu_Nm2");

    string error_message = "DLMFD_System::store_ugradu_Nm2 ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
bool DLMFD_System::VelocityAdvection_solver(void)
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System:: VelocityAdvection_solver");

    string error_message = "DLMFD_System::VelocityAdvection_solver ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}




//----------------------------------------------------------------------
void DLMFD_System::initialize_QUvector_with_divv_rhs()
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::initialize_QUvector_with_divv_rhs");

    string error_message = "DLMFD_System::initialize_QUvector_with_divv_rhs ";
    error_message += "should not be called !! Check implementation";
    MAC_Error::object()->raise_plain(error_message);
}

// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
// --------------------------- DLMFD FRAMEWORK ------------------------------

//----------------------------------------------------------------------
double const DLMFD_System::get_DLMFD_convergence_criterion() const
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::get_DLMFD_convergence_criterion");

    return Uzawa_DLMFD_precision;
}




//----------------------------------------------------------------------
int const DLMFD_System::get_DLMFD_maxiter() const
//----------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_System::get_DLMFD_convergence_criterion");

    return Uzawa_DLMFD_maxiter;
}
