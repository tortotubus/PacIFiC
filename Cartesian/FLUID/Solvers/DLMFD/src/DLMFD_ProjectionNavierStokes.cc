#include <DLMFD_AllRigidBodies.hh>
#include <DLMFD_ProjectionNavierStokes.hh>
#include <FS_SolidPlugIn.hh>
#include <FS_SolidPlugIn_BuilderFactory.hh>
#include <FV.hh>
#include <FV_DiscreteField.hh>
#include <FV_DiscreteField_Centered.hh>
#include <FV_DomainAndFields.hh>
#include <FV_DomainBuilder.hh>
#include <FV_Mesh.hh>
#include <FV_SystemNumbering.hh>
#include <FV_TimeIterator.hh>
#include <GrainsExec.hh>
#include <LA_Vector.hh>
#include <MAC.hh>
#include <MAC_Application.hh>
#include <MAC_Communicator.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <MAC_ListIdentity.hh>
#include <MAC_ModuleExplorer.hh>
#include <MAC_ModuleIterator.hh>
#include <MAC_Root.hh>
#include <MAC_Vector.hh>
#include <PAC_Misc.hh>
#include <cstdlib>
#include <intVector.hh>
#include <math.h>

DLMFD_ProjectionNavierStokes const *DLMFD_ProjectionNavierStokes::PROTOTYPE =
    new DLMFD_ProjectionNavierStokes();
bool DLMFD_ProjectionNavierStokes::b_pressure_drop_each_time = false;

//---------------------------------------------------------------------------
DLMFD_ProjectionNavierStokes::DLMFD_ProjectionNavierStokes(void)
    //--------------------------------------------------------------------------
    : FV_OneStepIteration("DLMFD_ProjectionNavierStokes"),
      PAC_ComputingTime("Solver")
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: DLMFD_ProjectionNavierStokes");
}




//---------------------------------------------------------------------------
DLMFD_ProjectionNavierStokes *
DLMFD_ProjectionNavierStokes::create_replica(MAC_Object *a_owner,
                                             FV_DomainAndFields const *dom,
                                             MAC_ModuleExplorer *exp) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: create_replica");
    MAC_CHECK(create_replica_PRE(a_owner, dom, exp));

    DLMFD_ProjectionNavierStokes *result =
        new DLMFD_ProjectionNavierStokes(a_owner, dom, exp);

    MAC_CHECK(create_replica_POST(result, a_owner, dom, exp));
    return (result);
}




//---------------------------------------------------------------------------
DLMFD_ProjectionNavierStokes::DLMFD_ProjectionNavierStokes(
    MAC_Object *a_owner, FV_DomainAndFields const *dom,
    MAC_ModuleExplorer const *exp)
    //---------------------------------------------------------------------------
    : FV_OneStepIteration(a_owner, dom, exp), PAC_ComputingTime("Solver"),
      UU(dom->discrete_field("velocity")), PP(dom->discrete_field("pressure")),
      GLOBAL_EQ(0), density(exp->double_data("Density")),
      viscosity(exp->double_data("Viscosity")), imposed_CFL(0.5),
      restartFileName_0("savingB"), restartFileName_1("savingA"),
      AdvectionScheme("TVD"), resultsDirectory("Res"), ViscousTimeAccuracy(1),
      AdvectionTimeAccuracy(1), b_ExplicitPressureGradient(false),
      b_HighOrderPressureCorrection(false), b_restart(false),
      b_pressure_rescaling(false), b_ExplicitDLMFD(false),
      explicitDLMFD_restartFilename_Prefix("DLM_"),
      b_projection_translation(dom->primary_grid()->is_translation_active()),
      primary_grid(dom->primary_grid()), critical_distance_translation(0.),
      translation_direction(0), bottom_coordinate(0.), translated_distance(0.),
      compute_flow_rate_on("none")
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: DLMFD_ProjectionNavierStokes");
    MAC_ASSERT(PP->discretization_type() == "centered");
    MAC_ASSERT(UU->discretization_type() == "staggered");

    // Call of MAC_Communicator routine to set the rank of each proces and
    // the number of processes during execution of DLMFD_ProjectionNavierStokes
    macCOMM = MAC_Exec::communicator();
    my_rank = macCOMM->rank();
    nb_ranks = macCOMM->nb_ranks();
    is_master = 0;

    // Timing routines
    if (my_rank == is_master)
    {
        CT_set_start();
        SCT_insert_app("Objects_Creation");
        SCT_set_start("Objects_Creation");
    }

    // Is the run a follow up of a previous job
    b_restart = MAC_Application::is_follow();

    // Results directory
    if (exp->has_entry("Results_directory"))
        resultsDirectory = exp->string_data("Results_directory");

    // Clear results directory in case of a new run
    if (!b_restart)
        PAC_Misc::clearAllFiles(resultsDirectory, "Savings", my_rank);

    // Get space dimension
    dim = UU->primary_grid()->nb_space_dimensions();
    if (dim == 1)
    {
        string error_message = "Space dimension should either 2 or 3";
        MAC_Error::object()->raise_bad_data_value(exp, "nb_space_dimensions",
                                                  error_message);
    }

    // Pressure rescaling
    b_pressure_rescaling = PP->all_BCs_nonDirichlet(0);

    // Imposed CFL
    if (exp->has_entry("Imposed_CFL"))
        imposed_CFL = exp->double_data("Imposed_CFL");

    // Advection scheme
    if (exp->has_entry("AdvectionScheme"))
        AdvectionScheme = exp->string_data("AdvectionScheme");
    if (AdvectionScheme != "Upwind" && AdvectionScheme != "TVD")
    {
        string error_message = "   - Upwind\n   - TVD";
        MAC_Error::object()->raise_bad_data_value(exp, "AdvectionScheme",
                                                  error_message);
    }
    if (AdvectionScheme == "TVD" &&
        UU->primary_grid()->get_security_bandwidth() < 2)
    {
        string error_message = "   >= 2 with TVD scheme";
        MAC_Error::object()->raise_bad_data_value(exp, "security_bandwidth",
                                                  error_message);
    }

    // Viscous term time accuracy
    if (exp->has_entry("ViscousTimeAccuracy"))
        ViscousTimeAccuracy = exp->int_data("ViscousTimeAccuracy");
    if (ViscousTimeAccuracy != 1 && ViscousTimeAccuracy != 2)
    {
        string error_message = "   - 1\n   - 2\n   ";
        MAC_Error::object()->raise_bad_data_value(exp, "ViscousTimeAccuracy",
                                                  error_message);
    }

    // Advection term time accuracy
    if (exp->has_entry("AdvectionTimeAccuracy"))
        AdvectionTimeAccuracy = exp->int_data("AdvectionTimeAccuracy");
    if (AdvectionTimeAccuracy != 1 && AdvectionTimeAccuracy != 2)
    {
        string error_message = "   - 1\n   - 2\n   ";
        MAC_Error::object()->raise_bad_data_value(exp, "AdvectionTimeAccuracy",
                                                  error_message);
    }

    // Use explicit pressure gradient in advection-diffusion equation
    if (exp->has_entry("ExplicitPressureGradient"))
        b_ExplicitPressureGradient = exp->bool_data("ExplicitPressureGradient");

    // High order pressure correction
    if (exp->has_entry("HighOrderPressureCorrection"))
        b_HighOrderPressureCorrection =
            exp->bool_data("HighOrderPressureCorrection");

    // Pressure drop handled each time step
    if (exp->has_entry("Pressure_drop_handled_each_time_step"))
        b_pressure_drop_each_time =
            exp->bool_data("Pressure_drop_handled_each_time_step");

    // Set the split gravity vector
    split_gravity_vector.resize(dim);
    split_gravity_vector.setVecZero();
    if (exp->has_entry("Split_gravity_vector"))
    {
        doubleVector *grav = new doubleVector(dim);
        *grav = exp->doubleVector_data("Split_gravity_vector");
        for (size_t idx = 0; idx < dim; ++idx)
            split_gravity_vector(idx) = (*grav)(idx);
        delete grav;
    }

    // Set the gravity vector
    gravity_vector.resize(dim);
    gravity_vector.setVecZero();
    if (exp->has_entry("Gravity_vector"))
    {
        doubleVector *grav = new doubleVector(dim);
        *grav = exp->doubleVector_data("Gravity_vector");
        for (size_t idx = 0; idx < dim; ++idx)
            gravity_vector(idx) = (*grav)(idx);
        delete grav;
    }

    // Set the output frequency
    size_t output_frequency = 1;
    if (exp->has_entry("Output_frequency"))
        output_frequency = size_t(exp->int_data("Output_frequency"));

    // Critical distance
    if (b_projection_translation)
    {
        if (exp->has_entry("Critical_Distance_Translation"))
            critical_distance_translation =
                exp->double_data("Critical_Distance_Translation");
        else
        {
            string error_message = " Projection-Translation is active but ";
            error_message += "Critical_Distance_Translation is NOT defined.";
            MAC_Error::object()->raise_bad_data_value(
                exp, "Projection_Translation", error_message);
        }
    }

    // NS parameters
    if (my_rank == is_master)
    {
        MAC::out() << endl << "*** Navier & Stokes algorithm" << endl << endl;
        MAC::out() << "   Name = Projection" << endl;
        MAC::out() << "   Viscous term time accuracy = " << ViscousTimeAccuracy
                   << endl;
        MAC::out() << "   Advection term time accuracy = "
                   << AdvectionTimeAccuracy << endl;
        MAC::out() << "   Imposed CFL = " << imposed_CFL << endl;
        MAC::out() << "   Explicit Pressure Gradient = "
                   << (b_ExplicitPressureGradient ? "true" : "false") << endl;
        MAC::out() << "   High Order Pressure Correction = "
                   << (b_HighOrderPressureCorrection ? "true" : "false")
                   << endl;
        MAC::out() << endl << endl;
    }

    // Flow rate
    if (exp->has_entry("Flow_rate"))
    {
        compute_flow_rate_on = exp->string_data("Flow_rate");
        if (compute_flow_rate_on != "none")
            PAC_Misc::is_main_boundary(compute_flow_rate_on, dim, "Flow_rate",
                                       exp);

        if (exp->has_entry("Flow_rate_Frequency"))
            flow_rate_frequency = exp->int_data("Flow_rate_Frequency");
        else
            flow_rate_frequency = 1;
    }

    // Build the matrix system
    MAC_ModuleExplorer *se =
        exp->create_subexplorer(0, "DLMFD_ProjectionNavierStokesSystem");
    GLOBAL_EQ = DLMFD_ProjectionNavierStokesSystem::create(
        this, se, UU, PP, ViscousTimeAccuracy, AdvectionTimeAccuracy,
        b_pressure_rescaling, b_ExplicitPressureGradient,
        b_HighOrderPressureCorrection);
    se->destroy();

    // Timing routines
    if (my_rank == is_master)
    {
        SCT_insert_app("Matrix_Assembly&Initialization");
        SCT_insert_app("VelocityPressure_CorrectionStep");
        SCT_insert_app("AdvectionDiffusion_PredictionStep");
        SCT_insert_app("Advection_PredictionStep");
        SCT_insert_app("Diffusion_PredictionStep");
        SCT_insert_app("FluidSolid_CorrectionStep");
        SCT_get_elapsed_time("Objects_Creation");
    }

    // Create Fictitious Domain solver
    struct NavierStokes2FluidSolid transfert;
    transfert.solid_resDir = resultsDirectory;
    transfert.rho_f = density;
    transfert.gravity_vector = gravity_vector;
    transfert.split_gravity_vector = split_gravity_vector;
    transfert.GLOBAL_EQ = GLOBAL_EQ;
    transfert.velocitylevelDiscrField = 0;
  transfert.nb_levels = 1;
    transfert.b_restart = b_restart;
    transfert.output_frequency = output_frequency;
    transfert.UU = UU;
    transfert.PP = PP;

    dlmfd_solver = DLMFD_FictitiousDomain::create(a_owner, dom, exp, transfert);

  // Use explicit DLMFD in N&S advection-diffusion equation
  b_ExplicitDLMFD = dlmfd_solver->get_explicit_DLMFD();
  if (b_ExplicitDLMFD)
    GLOBAL_EQ->re_initialize_explicit_DLMFD(b_restart);
}




//---------------------------------------------------------------------------
DLMFD_ProjectionNavierStokes::~DLMFD_ProjectionNavierStokes(void)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: ~DLMFD_ProjectionNavierStokes");
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::do_one_inner_iteration(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: do_one_inner_iteration");
    MAC_CHECK_PRE(do_one_inner_iteration_PRE(t_it));

    start_total_timer("DLMFD_ProjectionNavierStokes:: do_one_inner_iteration");
    start_solving_timer();
    macCOMM->barrier();

    // Solve the constrained problem momentum + mass conservation
    // by a fractional step projection algorithm
    NavierStokes_Projection(t_it);

    // Use the Fictitious Domain method to solve the coupling
    // between the fluid and the solid
    if (my_rank == is_master)
        SCT_set_start("FluidSolid_CorrectionStep");

    dlmfd_solver->do_one_inner_iteration(t_it, sub_prob_number);

    if (my_rank == is_master)
        SCT_get_elapsed_time("FluidSolid_CorrectionStep");

    stop_solving_timer();
    stop_total_timer();
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::do_before_time_stepping(
    FV_TimeIterator const *t_it, std::string const &basename)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: do_before_time_stepping");

    start_total_timer("DLMFD_ProjectionNavierStokes:: do_before_time_stepping");

    FV_OneStepIteration::do_before_time_stepping(t_it, basename);

    if (!b_pressure_drop_each_time)
        // Assemble constant matrix of periodic pressure drop
        if (UU->primary_grid()->is_periodic_pressure_drop())
            assemble_periodic_pressure_rhs();

    // Projection-Translation
    if (b_projection_translation)
    {
        set_translation_vector();

        if (MVQ_translation_vector(translation_direction) < 0.)
            bottom_coordinate =
                (*primary_grid
                      ->get_global_main_coordinates())[translation_direction](
                    0);
        else
            bottom_coordinate =
                (*primary_grid
                      ->get_global_main_coordinates())[translation_direction](
                    (*primary_grid->get_global_max_index())(
                        translation_direction));

        geomVector vt(dim);
        vt(translation_direction) = primary_grid->get_translation_distance();
        dlmfd_solver->setParaviewPostProcessingTranslationVector(-vt(0), -vt(1),
                                                                 -vt(2));

        build_links_translation();
    }

    if (my_rank == is_master)
        SCT_set_start("Matrix_Assembly&Initialization");

    // Assemble matrices
    if (my_rank == is_master)
        MAC::out() << "         Assemble matrices & vectors" << endl;

    // Velocity unsteady matrix
    if (my_rank == is_master)
        MAC::out() << "            Velocity unsteady matrix" << endl;
    GLOBAL_EQ->assemble_velocity_unsteady_matrix(density / t_it->time_step());

    // Velocity viscous matrix and rhs
    // Note: we assemble here the total viscous matrix, in case of 2nd order
    // Crank-Nicholson scheme for the viscous term, half viscosity for the
    // velocity operator is taken care of at the matrix level in
    // DLMFD_ProjectionNavierStokesSystem:: finalize_constant_matrices
    if (my_rank == is_master)
        MAC::out() << "            Velocity viscous matrix & rhs" << endl;
    GLOBAL_EQ->assemble_velocity_viscous_matrix_rhs(-viscosity);

    // Velocity divergence matrix and rhs
    if (my_rank == is_master)
        MAC::out() << "            Velocity divergence matrix & rhs" << endl;
    GLOBAL_EQ->assemble_pdivv_matrix_rhs(-1.);

    // Pressure laplacian matrix and rhs
    if (my_rank == is_master)
        MAC::out() << "            Pressure laplacian matrix & rhs" << endl;
    GLOBAL_EQ->assemble_pressure_laplacian_matrix_rhs(-1.);
    GLOBAL_EQ->pressure_laplacian_correction();

    // Pressure Dirichlet BC
    assemble_pressure_DirichletBC_in_momentumEquation(
        GLOBAL_EQ->get_pressure_DirichletBC_vector());

    if (b_pressure_drop_each_time)
        // Unitary periodic pressure gradient rhs
        if (UU->primary_grid()->is_periodic_flow())
            assemble_unitary_periodic_pressure_gradient_rhs(
                GLOBAL_EQ->get_unitary_periodic_pressure_drop_vector());

    // Synchronize and finalize matrices
    GLOBAL_EQ->finalize_constant_matrices();

    // Initialize velocity & pressure
    GLOBAL_EQ->initialize_velocity();
    GLOBAL_EQ->initialize_pressure();

    // Do additional reload
    do_additional_reload(basename);

    if (my_rank == is_master)
        SCT_get_elapsed_time("Matrix_Assembly&Initialization");

    // Do DLMFD solid/fluid set-up
    dlmfd_solver->do_before_time_stepping(t_it);

    // Flow rate
    if (compute_flow_rate_on != "none")
    {
        if (my_rank == is_master)
        {
            MAC::out() << endl << "   +++ Flow rate ++++++" << endl << endl;
            MAC::out() << "      Output every " << flow_rate_frequency
                       << " time steps" << endl;
            MAC::out() << "      Boundary name = " << compute_flow_rate_on
                       << endl;
            MAC::out() << endl << endl;
        }
        compute_flow_rate(t_it, b_restart);
    }

    stop_total_timer();
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::do_after_time_stepping(void)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: do_after_time_stepping");

    // Elapsed time by sub-problems
    if (my_rank == is_master)
    {
        double cputime = CT_get_elapsed_time();
        MAC::out() << endl << "Full problem" << endl;
        write_elapsed_time_smhd(MAC::out(), cputime, "Computing time");
        SCT_get_summary(MAC::out(), cputime);
    }
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::do_before_inner_iterations_stage(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL(
        "DLMFD_ProjectionNavierStokes:: do_before_inner_iterations_stage");

    start_total_timer(
        "DLMFD_ProjectionNavierStokes:: do_before_inner_iterations_stage");

    FV_OneStepIteration::do_before_inner_iterations_stage(t_it);

    // Perform matrix level operations before each time step
    GLOBAL_EQ->at_each_time_step();

    if (b_pressure_drop_each_time)
        // Update pressure drop in case of periodic imposed flow rate
        if (UU->primary_grid()->is_periodic_flow_rate())
            update_pressure_drop_imposed_flow_rate(t_it);

    stop_total_timer();
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::do_after_inner_iterations_stage(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: do_after_inner_iterations_stage");

    start_total_timer(
        "DLMFD_ProjectionNavierStokes:: do_after_inner_iterations_stage");

    FV_OneStepIteration::do_after_inner_iterations_stage(t_it);

    // Compute velocity change over the time step
    double velocity_time_change =
        GLOBAL_EQ->compute_velocity_change() / t_it->time_step();
    if (my_rank == is_master)
        MAC::out() << "         Velocity change = "
                   << MAC::doubleToString(ios::scientific, 14,
                                          velocity_time_change)
                   << endl;

    // Compute Norm(div(u))
    double normdivu = GLOBAL_EQ->compute_velocity_divergence_norm();
    if (my_rank == is_master)
        MAC::out() << "         L2 Norm of distribvec(div(u)) = "
                   << MAC::doubleToString(ios::scientific, 14, normdivu)
                   << endl;
    compute_and_print_divu_norm();

    // Flow rate
    compute_flow_rate(t_it, false);

    // DLMFD computation
    dlmfd_solver->do_after_inner_iterations_stage(t_it);

    // Projection translation
    if (b_projection_translation)
    {
        double distance_to_bottom = dlmfd_solver->Compute_distance_to_bottom(
            bottom_coordinate, translation_direction);

        if (my_rank == is_master)
            MAC::out() << "         Distance to bottom = "
                       << MAC::doubleToString(ios::scientific, 5,
                                              distance_to_bottom)
                       << endl;

        if (distance_to_bottom < critical_distance_translation)
        {

            if (my_rank == is_master)
                MAC::out() << "         -> -> -> -> -> -> -> -> -> -> -> ->"
                           << endl
                           << "         !!!     Domain Translation      !!!"
                           << endl
                           << "         -> -> -> -> -> -> -> -> -> -> -> ->"
                           << endl;

            translated_distance +=
                MVQ_translation_vector(translation_direction);
            if (my_rank == is_master)
                MAC::out() << "         Translated distance = "
                           << translated_distance << endl;

            fields_projection();

            dlmfd_solver->translate_all(MVQ_translation_vector,
                                        translation_direction);

            if (MVQ_translation_vector(translation_direction) < 0.)
                bottom_coordinate =
                    (*primary_grid->get_global_main_coordinates())
                        [translation_direction](0);
            else
                bottom_coordinate =
                    (*primary_grid->get_global_main_coordinates())
                        [translation_direction](
                            (*primary_grid->get_global_max_index())(
                                translation_direction));
        }
    }

    stop_total_timer();
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::do_additional_savings(
    FV_TimeIterator const *t_it, int const &cycleNumber)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: do_additional_savings");

    start_total_timer("DLMFD_ProjectionNavierStokes:: do_additional_savings");

    // DLMFD additional savings
    dlmfd_solver->do_additional_savings(cycleNumber, t_it, translated_distance,
                                        translation_direction);

    // Elapsed time by sub-problems
    if (my_rank == is_master)
    {
        double cputime = CT_get_elapsed_time();
        MAC::out() << endl << "Full problem" << endl;
        write_elapsed_time_smhd(MAC::out(), cputime, "Computation time");
        SCT_get_summary(MAC::out(), cputime);
    }

    stop_total_timer();
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::do_additional_save_for_restart(
    FV_TimeIterator const *t_it, size_t const &restartCycleNumber,
    std::string const &basename)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: do_additional_save_for_restart");

    start_total_timer(
        "DLMFD_ProjectionNavierStokes:: do_additional_save_for_restart");

    stop_total_timer();
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::set_translation_vector()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: set_translation_vector");

    MVQ_translation_vector.resize(primary_grid->nb_space_dimensions());
    translation_direction = primary_grid->get_translation_direction();
    MVQ_translation_vector(translation_direction) =
        primary_grid->get_translation_magnitude();
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::fields_projection()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: fields_projection");

    FV_Mesh *pmesh = const_cast<FV_Mesh *>(primary_grid);
    pmesh->translation();

    UU->translation_projection(0, 1, !b_ExplicitDLMFD);
    synchronize_velocity_field(0);

    PP->translation_projection(0, 1);
    synchronize_pressure_field(0);

    if (b_ExplicitDLMFD)
    {
        // Copy back DLMFD vector on velocity field
        UU->update_free_DOFs_value(0, GLOBAL_EQ->get_rhs_DLMFD_Nm1());

        // Projection-translation at the level of the field
        UU->translation_projection(0, 1, b_ExplicitDLMFD);

        // Synchronize the DLMFD forcing term at previous time
        synchronize_DLMFD_Nm1_rhs(0);

        // Re-assign velocity field from vector
        UU->update_free_DOFs_value(0, GLOBAL_EQ->get_solution_U());
    }
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::synchronize_velocity_field(size_t level)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: synchronize_velocity_field");

    size_t nb_comps = UU->nb_components();

    // Transfer field values owned by the proc only to the global vector
    for (size_t comp = 0; comp < nb_comps; ++comp)
    {
        // Get local min and max indices
        size_t_vector min_unknown_index(dim, 0);
        for (size_t l = 0; l < dim; ++l)
            min_unknown_index(l) =
                UU->get_min_index_unknown_handled_by_proc(comp, l);
        size_t_vector max_unknown_index(dim, 0);
        for (size_t l = 0; l < dim; ++l)
            max_unknown_index(l) =
                UU->get_max_index_unknown_handled_by_proc(comp, l);

        // Set value in global vector
        size_t k = 0;
        for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i)
            for (size_t j = min_unknown_index(1); j <= max_unknown_index(1);
                 ++j)
                if (dim == 2)
                    GLOBAL_EQ->set_velocity_unknown(
                        UU->DOF_global_number(i, j, k, comp),
                        UU->DOF_value(i, j, k, comp, level));
                else
                    for (k = min_unknown_index(2); k <= max_unknown_index(2);
                         ++k)
                        GLOBAL_EQ->set_velocity_unknown(
                            UU->DOF_global_number(i, j, k, comp),
                            UU->DOF_value(i, j, k, comp, level));
    }

    // Transfer back values from the global vector to the field
    GLOBAL_EQ->synchronize_velocity_unknown_vector();
    UU->update_free_DOFs_value(level, GLOBAL_EQ->get_solution_U());
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::synchronize_pressure_field(size_t level)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: synchronize_pressure_field");

    // Transfer field values owned by the proc only to the global vector
    // Get local min and max indices
    size_t_vector min_unknown_index(dim, 0);
    for (size_t l = 0; l < dim; ++l)
        min_unknown_index(l) = PP->get_min_index_unknown_handled_by_proc(0, l);
    size_t_vector max_unknown_index(dim, 0);
    for (size_t l = 0; l < dim; ++l)
        max_unknown_index(l) = PP->get_max_index_unknown_handled_by_proc(0, l);

    // Set value in global vector
    size_t k = 0;
    for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i)
        for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j)
            if (dim == 2)
                GLOBAL_EQ->set_pressure_unknown(
                    PP->DOF_global_number(i, j, k, 0),
                    PP->DOF_value(i, j, k, 0, level));
            else
                for (k = min_unknown_index(2); k <= max_unknown_index(2); ++k)
                    GLOBAL_EQ->set_pressure_unknown(
                        PP->DOF_global_number(i, j, k, 0),
                        PP->DOF_value(i, j, k, 0, level));

    // Transfer back values from the global vector to the field
    GLOBAL_EQ->synchronize_pressure_unknown_vector();
    PP->update_free_DOFs_value(level, GLOBAL_EQ->get_solution_pressure());
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::synchronize_DLMFD_Nm1_rhs(size_t level)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: synchronize_DLMFD_Nm1_rhs");

    size_t nb_comps = UU->nb_components();

    // Nullify the DLMFD forcing term at previous time at the matrix level
    GLOBAL_EQ->nullify_DLMFD_Nm1_rhs();

    // Transfer field values owned by the proc only to the global vector
    for (size_t comp = 0; comp < nb_comps; ++comp)
    {
        // Get local min and max indices
        size_t_vector min_unknown_index(dim, 0);
        for (size_t l = 0; l < dim; ++l)
            min_unknown_index(l) =
                UU->get_min_index_unknown_handled_by_proc(comp, l);
        size_t_vector max_unknown_index(dim, 0);
        for (size_t l = 0; l < dim; ++l)
            max_unknown_index(l) =
                UU->get_max_index_unknown_handled_by_proc(comp, l);

        // Set value in global vector
        size_t k = 0;
        for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i)
            for (size_t j = min_unknown_index(1); j <= max_unknown_index(1);
                 ++j)
                if (dim == 2)
                    GLOBAL_EQ->set_rhs_DLMFD_Nm1(
                        UU->DOF_global_number(i, j, k, comp),
                        UU->DOF_value(i, j, k, comp, level));
                else
                    for (k = min_unknown_index(2); k <= max_unknown_index(2);
                         ++k)
                        GLOBAL_EQ->set_rhs_DLMFD_Nm1(
                            UU->DOF_global_number(i, j, k, comp),
                            UU->DOF_value(i, j, k, comp, level));
    }

    // Transfer back values from the global vector to the field
    GLOBAL_EQ->synchronize_rhs_DLMFD_Nm1_vector();
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::NavierStokes_Projection(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: NavierStokes_Projection");

    // Solve advection-diffusion velocity prediction step
    NavierStokes_AdvectionDiffusion_PredictionStep(t_it);

    // Project velocity on a divergence free space by
    // solving a pressure Poisson problem, and update velocity and pressure
    NavierStokes_VelocityPressure_CorrectionStep(t_it);

    // Copy back velocity solution in field
    UU->update_free_DOFs_value(0, GLOBAL_EQ->get_solution_velocity());

    if (!b_pressure_drop_each_time)
        // Correct velocity in case of periodic imposed flow rate
        if (UU->primary_grid()->is_periodic_flow_rate())
            periodic_flow_rate_update(t_it);
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::
    NavierStokes_AdvectionDiffusion_PredictionStep(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: "
              "NavierStokes_AdvectionDiffusion_PredictionStep");

    sub_prob_number = 1;

    // Compute CFL
    double computed_CFL = UU->compute_CFL(t_it, 0);
    size_t n_advection_subtimesteps = unsigned(computed_CFL / imposed_CFL) + 1;
    double dpdl = 0.;

    // Calculate periodic pressure gradient from periodic pressure drop
    // Is 0 if the flow is not periodic
    if (UU->primary_grid()->is_periodic_flow())
    {
        size_t comp = UU->primary_grid()->get_periodic_flow_direction();
        dpdl = -UU->primary_grid()->get_periodic_pressure_drop() /
               (UU->primary_grid()->get_main_domain_max_coordinate(comp) -
                UU->primary_grid()->get_main_domain_min_coordinate(comp));
    }

    if (n_advection_subtimesteps == 1)
    {
        if (my_rank == is_master)
        {
            MAC::out()
                << "------------------------------------------------------"
                << endl;
            MAC::out() << "Sub-problem " << sub_prob_number
                       << " : Advection-diffusion prediction problem" << endl;
            MAC::out()
                << "------------------------------------------------------"
                << endl;
            MAC::out() << "CFL : imposed = " << imposed_CFL
                       << " computed = " << computed_CFL
                       << "  Nb of sub time steps = "
                       << n_advection_subtimesteps << endl;
            SCT_set_start("AdvectionDiffusion_PredictionStep");
        }

        // Compute velocity advection rhs
        GLOBAL_EQ->assemble_velocity_advection(AdvectionScheme, 0, -density, 0);

        // Compute the advection-diffusion rhs
        // The method adds the periodic pressure drop in case of periodic flow
        GLOBAL_EQ->compute_velocityAdvectionDiffusion_rhs(
            b_restart, t_it->iteration_number(), true, dpdl);

        // If the explicit pressure gradient is added to the rhs, this is done
        // by DLMFD_ProjectionNavierStokesSystem::VelocityDiffusion_solver at
        // the matrix level. However here we want to add the gradient of the
        // actual pressure, while at the matrix level at that stage, the
        // pressure vector contains the pressure correction, not the pressure,
        // hence we re-initialize the pressure vector with the pressure
        if (b_ExplicitPressureGradient)
            GLOBAL_EQ->initialize_pressure();

        // Solve the diffusion system with advection in rhs and the explicit
        // pressure gradient in case b_ExplicitPressureGradient is true
        GLOBAL_EQ->VelocityDiffusion_solver();
        ++sub_prob_number;
        if (my_rank == is_master)
            SCT_get_elapsed_time("AdvectionDiffusion_PredictionStep");
    }
    else
    {
        // Sub-problem: Advection problem
        if (my_rank == is_master)
        {
            MAC::out()
                << "------------------------------------------------------"
                << endl;
            MAC::out() << "Sub-problem " << sub_prob_number
                       << " : Advection prediction problem" << endl;
            MAC::out()
                << "------------------------------------------------------"
                << endl;
            MAC::out() << "CFL : imposed = " << imposed_CFL
                       << " computed = " << computed_CFL
                       << "  Nb of sub time steps = "
                       << n_advection_subtimesteps << endl;
            if (AdvectionTimeAccuracy == 2)
                MAC::out() << "Degenerates as order 1 scheme, affects the whole"
                              " splitting algorithm"
                           << endl;
            SCT_set_start("Advection_PredictionStep");
        }

        for (size_t i = 0; i < n_advection_subtimesteps; ++i)
        {
            // Copy back velocity solution in field for i != 0 as for i = 0, the
            // field already has the right values of the velocity
            if (i)
                UU->update_free_DOFs_value(0,
                                           GLOBAL_EQ->get_solution_velocity());

            // Compute velocity advection rhs
            GLOBAL_EQ->assemble_velocity_advection(
                AdvectionScheme, 0, -density / double(n_advection_subtimesteps),
                0);

            // Solve advection problem
            GLOBAL_EQ->VelocityAdvection_solver();

            // Store velocity advection rhs at t^n-1 in case of 2nd order scheme
            // Then if at next time the scheme is again 2nd (no sub-time
            // stepping) the explicit Adams-Bashforth term is properly computed
            if (AdvectionTimeAccuracy == 2 && i == 0)
                GLOBAL_EQ->store_ugradu_Nm2(n_advection_subtimesteps);
        }
        ++sub_prob_number;
        if (my_rank == is_master)
            SCT_get_elapsed_time("Advection_PredictionStep");

        // Sub-problem: Diffusion (viscous) problem
        // Compute the diffusion rhs (note the flase boolean as the last
        // parameter) The method adds the periodic pressure drop in case of
        // periodic flow
        GLOBAL_EQ->compute_velocityAdvectionDiffusion_rhs(
            b_restart, t_it->iteration_number(), false, dpdl);

        // If the explicit pressure gradient is added to the rhs, this is done
        // by DLMFD_ProjectionNavierStokesSystem::VelocityDiffusion_solver at
        // the matrix level. However here we want to add the gradient of the
        // actual pressure, while at the matrix level at that stage, the
        // pressure vector contains the pressure correction, not the pressure,
        // hence we re-initialize the pressure vector with the pressure
        if (b_ExplicitPressureGradient)
            GLOBAL_EQ->initialize_pressure();

        // Solve the diffusion (viscous) problem
        if (my_rank == is_master)
        {
            MAC::out()
                << "------------------------------------------------------"
                << endl;
            MAC::out() << "Sub-problem " << sub_prob_number
                       << " : Viscous problem" << endl;
            MAC::out()
                << "------------------------------------------------------"
                << endl;
            SCT_set_start("Diffusion_PredictionStep");
        }
        GLOBAL_EQ->VelocityDiffusion_solver();
        ++sub_prob_number;
        if (my_rank == is_master)
            SCT_get_elapsed_time("Diffusion_PredictionStep");
    }
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::NavierStokes_VelocityPressure_CorrectionStep(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: "
              "NavierStokes_VelocityPressure_CorrectionStep");

    if (my_rank == is_master)
    {
        MAC::out() << "------------------------------------------------------"
                   << endl;
        MAC::out() << "Sub-problem " << sub_prob_number
                   << " : Velocity-pressure correction problem" << endl;
        MAC::out() << "------------------------------------------------------"
                   << endl;
        SCT_set_start("VelocityPressure_CorrectionStep");
    }

    // Solve the pressure Poisson problem, correct pressure and update velocity
    GLOBAL_EQ->VelocityPressure_correction_solver(density, viscosity,
                                                  t_it->time_step());
    ++sub_prob_number;

    // Copy back pressure solution in field
    if (b_ExplicitPressureGradient)
        PP->add_to_free_DOFs_value(0, GLOBAL_EQ->get_solution_pressure());
    else
        PP->update_free_DOFs_value(0, GLOBAL_EQ->get_solution_pressure());

    if (my_rank == is_master)
        SCT_get_elapsed_time("VelocityPressure_CorrectionStep");
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::do_additional_reload(string const &basename)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: do_additional_reload");
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::add_storable_objects(MAC_ListIdentity *list)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: add_storable_objects");

    GLOBAL_EQ->add_storable_objects(list);
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::
    assemble_pressure_DirichletBC_in_momentumEquation(LA_Vector *VEC_rhs)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: "
              "assemble_pressure_BC_in_momentumEquation");

    if (my_rank == is_master)
        MAC::out() << "            Pressure BC momentum equation" << endl;

    // Parameters
    size_t_vector min_unknown_index(dim, 0);
    size_t_vector max_unknown_index(dim, 0);
    size_t nb_comps = UU->nb_components(), center_pos_in_matrix = 0, comp;
    double dxC, dyC, dzC, ai;

    FV_SHIFT_TRIPLET shift;

    for (comp = 0; comp < nb_comps; ++comp)
    {
        // Get local min and max indices
        for (size_t l = 0; l < dim; ++l)
            min_unknown_index(l) =
                UU->get_min_index_unknown_handled_by_proc(comp, l);
        for (size_t l = 0; l < dim; ++l)
            max_unknown_index(l) =
                UU->get_max_index_unknown_handled_by_proc(comp, l);

        shift = UU->shift_staggeredToStaggered(comp);

        // Perform assembling
        for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i)
        {
            dxC = UU->get_cell_size(i, comp, 0);
            for (size_t j = min_unknown_index(1); j <= max_unknown_index(1);
                 ++j)
            {
                dyC = UU->get_cell_size(j, comp, 1);

                if (dim == 2)
                {
                    size_t k = 0;
                    center_pos_in_matrix = UU->DOF_global_number(i, j, k, comp);

                    switch (comp)
                    {
                    // The First Component (u)
                    case 0:
                        // Right (X)
                        if (PP->DOF_has_imposed_Dirichlet_value(i + shift.i, j,
                                                                k, 0))
                        {
                            ai = -dyC;
                            VEC_rhs->add_to_item(
                                center_pos_in_matrix,
                                ai * PP->DOF_value(i + shift.i, j, k, 0, 0));
                        }
                        // Left (X)
                        if (PP->DOF_has_imposed_Dirichlet_value(i + shift.i - 1,
                                                                j, k, 0))
                        {
                            ai = dyC;
                            VEC_rhs->add_to_item(
                                center_pos_in_matrix,
                                ai *
                                    PP->DOF_value(i + shift.i - 1, j, k, 0, 0));
                        }
                        break;

                    // The second Component (v)
                    case 1:
                        // Top (Y)
                        if (PP->DOF_has_imposed_Dirichlet_value(i, j + shift.j,
                                                                k, 0))
                        {
                            ai = -dxC;
                            VEC_rhs->add_to_item(
                                center_pos_in_matrix,
                                ai * PP->DOF_value(i, j + shift.j, k, 0, 0));
                        }
                        // Bottom (Y)
                        if (PP->DOF_has_imposed_Dirichlet_value(
                                i, j + shift.j - 1, k, 0))
                        {
                            ai = dxC;
                            VEC_rhs->add_to_item(
                                center_pos_in_matrix,
                                ai *
                                    PP->DOF_value(i, j + shift.j - 1, k, 0, 0));
                        }
                        break;
                    }
                }
                else
                {
                    for (size_t k = min_unknown_index(2);
                         k <= max_unknown_index(2); ++k)
                    {
                        dzC = UU->get_cell_size(k, comp, 2);
                        center_pos_in_matrix =
                            UU->DOF_global_number(i, j, k, comp);

                        switch (comp)
                        {
                        // The First Component (u)
                        case 0:
                            // Right (X)
                            if (PP->DOF_has_imposed_Dirichlet_value(i + shift.i,
                                                                    j, k, 0))
                            {
                                ai = -dyC * dzC;
                                VEC_rhs->add_to_item(
                                    center_pos_in_matrix,
                                    ai *
                                        PP->DOF_value(i + shift.i, j, k, 0, 0));
                            }
                            // Left (X)
                            if (PP->DOF_has_imposed_Dirichlet_value(
                                    i + shift.i - 1, j, k, 0))
                            {
                                ai = dyC * dzC;
                                VEC_rhs->add_to_item(
                                    center_pos_in_matrix,
                                    ai * PP->DOF_value(i + shift.i - 1, j, k, 0,
                                                       0));
                            }
                            break;

                        // The Second Component (v)
                        case 1:
                            // Top (Y)
                            if (PP->DOF_has_imposed_Dirichlet_value(
                                    i, j + shift.j, k, 0))
                            {
                                ai = -dxC * dzC;
                                VEC_rhs->add_to_item(
                                    center_pos_in_matrix,
                                    ai *
                                        PP->DOF_value(i, j + shift.j, k, 0, 0));
                            }
                            // Bottom (Y)
                            if (PP->DOF_has_imposed_Dirichlet_value(
                                    i, j + shift.j - 1, k, 0))
                            {
                                ai = dxC * dzC;
                                VEC_rhs->add_to_item(
                                    center_pos_in_matrix,
                                    ai * PP->DOF_value(i, j + shift.j - 1, k, 0,
                                                       0));
                            }
                            break;

                        // The Third Component (w)
                        case 2:
                            // Front (Z)
                            if (PP->DOF_has_imposed_Dirichlet_value(
                                    i, j, k + shift.k, 0))
                            {
                                ai = -dxC * dyC;
                                VEC_rhs->add_to_item(
                                    center_pos_in_matrix,
                                    ai *
                                        PP->DOF_value(i, j, k + shift.k, 0, 0));
                            }
                            // Behind (Z)
                            if (PP->DOF_has_imposed_Dirichlet_value(
                                    i, j, k + shift.k - 1, 0))
                            {
                                ai = dxC * dyC;
                                VEC_rhs->add_to_item(
                                    center_pos_in_matrix,
                                    ai * PP->DOF_value(i, j, k + shift.k - 1, 0,
                                                       0));
                            }
                            break;
                        }
                    }
                }
            }
        }
    }

    VEC_rhs->synchronize();
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::
    assemble_unitary_periodic_pressure_gradient_rhs(LA_Vector *VEC_rhs)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: "
              "assemble_unitary_periodic_pressure_gradient_rhs");

    cout << VEC_rhs << endl;

    // Parameters
    size_t comp = UU->primary_grid()->get_periodic_flow_direction();
    double dxC, dyC, dzC;
    size_t center_pos_in_matrix = 0;

    // Get local min and max indices
    size_t_vector min_unknown_index(dim, 0);
    for (size_t l = 0; l < dim; ++l)
        min_unknown_index(l) =
            UU->get_min_index_unknown_handled_by_proc(comp, l);
    size_t_vector max_unknown_index(dim, 0);
    for (size_t l = 0; l < dim; ++l)
        max_unknown_index(l) =
            UU->get_max_index_unknown_handled_by_proc(comp, l);

    // Perform assembling
    for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i)
    {
        dxC = UU->get_cell_size(i, comp, 0);
        for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j)
        {
            dyC = UU->get_cell_size(j, comp, 1);

            if (dim == 2)
            {
                size_t k = 0;
                center_pos_in_matrix = UU->DOF_global_number(i, j, k, comp);
                VEC_rhs->set_item(center_pos_in_matrix, dxC * dyC);
            }
            else
            {
                for (size_t k = min_unknown_index(2); k <= max_unknown_index(2);
                     ++k)
                {
                    dzC = UU->get_cell_size(k, comp, 2);
                    center_pos_in_matrix = UU->DOF_global_number(i, j, k, comp);
                    VEC_rhs->set_item(center_pos_in_matrix, dxC * dyC * dzC);
                }
            }
        }
    }

    VEC_rhs->synchronize();
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::assemble_periodic_pressure_rhs()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: assemble_periodic_pressure_rhs");

    // Parameters
    size_t comp = UU->primary_grid()->get_periodic_flow_direction();
    double dpdl = -UU->primary_grid()->get_periodic_pressure_drop() /
                  (UU->primary_grid()->get_main_domain_max_coordinate(comp) -
                   UU->primary_grid()->get_main_domain_min_coordinate(comp));
    double dxC, dyC, dzC;
    size_t center_pos_in_matrix = 0;

    // Get local min and max indices
    size_t_vector min_unknown_index(dim, 0);
    for (size_t l = 0; l < dim; ++l)
        min_unknown_index(l) =
            UU->get_min_index_unknown_handled_by_proc(comp, l);
    size_t_vector max_unknown_index(dim, 0);
    for (size_t l = 0; l < dim; ++l)
        max_unknown_index(l) =
            UU->get_max_index_unknown_handled_by_proc(comp, l);

    // Perform assembling
    for (size_t i = min_unknown_index(0); i <= max_unknown_index(0); ++i)
    {
        dxC = UU->get_cell_size(i, comp, 0);
        for (size_t j = min_unknown_index(1); j <= max_unknown_index(1); ++j)
        {
            dyC = UU->get_cell_size(j, comp, 1);

            if (dim == 2)
            {
                size_t k = 0;
                center_pos_in_matrix = UU->DOF_global_number(i, j, k, comp);
                GLOBAL_EQ->set_periodic_pressure_rhs_item(center_pos_in_matrix,
                                                          dpdl * dxC * dyC);
            }
            else
            {
                for (size_t k = min_unknown_index(2); k <= max_unknown_index(2);
                     ++k)
                {
                    dzC = UU->get_cell_size(k, comp, 2);
                    center_pos_in_matrix = UU->DOF_global_number(i, j, k, comp);
                    GLOBAL_EQ->set_periodic_pressure_rhs_item(
                        center_pos_in_matrix, dpdl * dxC * dyC * dzC);
                }
            }
        }
    }

    GLOBAL_EQ->synchronize_rhs_periodic_pressure_vector();
}

// ------ NOTE TO THE DEVELOPER ------
// The following draft method has been added by Anthony in order to implement a 
// prediction of the periodic pressure gradient in the case of imposed 
// flow rate. It is kept here as a history, and is not used at all in the 
// Projection/DLMFD workflow.


// //---------------------------------------------------------------------------
// void DLMFD_ProjectionNavierStokes::update_pressure_drop_imposed_flow_rate(
//     FV_TimeIterator const *t_it)
// //---------------------------------------------------------------------------
// {
//   MAC_LABEL(
//       "DLMFD_ProjectionNavierStokes::
//       update_pressure_drop_imposed_flow_rate");

//   FV_Mesh const *primary_mesh = UU->primary_grid();

//   // Get the periodic flow direction
//   size_t periodic_flow_direction =
//   primary_mesh->get_periodic_flow_direction();

//   // Get the boundary name
//   string boundary_name = periodic_flow_direction == 0 ? "right" :
//   periodic_flow_direction == 1 ? "top"
//                                                                                                : "front";
//   // Compute the flow rate
//   double flow_rate = PAC_Misc::compute_flow_rate(UU, 0, boundary_name);

//   // Get the imposed flow rate
//   double imposed_flow_rate = primary_mesh->get_periodic_flow_rate();

//   // Correction
//   double kp = 1.e1;
//   //   double ki = 1.e-1, ks = 1.e0;
//   double old_ppd = primary_mesh->get_periodic_pressure_drop();
//   //    static double Ik = 0;
//   //    static double flow_rate_previoustime = 0 ;
//   //    static double ppd_previoustime = 0 ;

//   //    double delta_ppd = kp * ( imposed_flow_rate - flow_rate ) + Ik ;
//   //    double new_ppd = old_ppd - delta_ppd ;
//   //    Ik = Ik + 0.1 * t_it->time_step() * ( ki * ( imposed_flow_rate -
//   flow_rate )
//   //    	+ ks * ( delta_ppd - kp * ( imposed_flow_rate - flow_rate ) -Ik
//   ) );

//   double delta_ppd = kp * (imposed_flow_rate - flow_rate);

//   double new_ppd = old_ppd - delta_ppd;

//   //    if ( t_it->iteration_number() == 1 )
//   //      new_ppd = - density * imposed_flow_rate * 4. / ( t_it->time_step()
//   * 1. );

//   const_cast<FV_Mesh *>(primary_mesh)->set_periodic_pressure_drop(new_ppd);
//   //    MAC::out() << "Flow rate = " << flow_rate << endl;
//   //    MAC::out() << "Ppd = " << new_ppd << endl;

//   string fr_filename = resultsDirectory + "/flowrate.dat";
//   // MAC::out() << fr_filename << endl;
//   ofstream FILEOUT(fr_filename.c_str(), ios::app);
//   FILEOUT << t_it->time() << " " << new_ppd << " " << flow_rate << endl;
//   FILEOUT.close();
// }

//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::update_pressure_drop_imposed_flow_rate(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: "
              "update_pressure_drop_imposed_flow_rate");

    FV_Mesh const *primary_mesh = UU->primary_grid();

    // Get the periodic flow direction
    size_t periodic_flow_direction =
        primary_mesh->get_periodic_flow_direction();

    // Area of boundary face perpendicular to periodic flow direction
    double face_area =
        primary_mesh->get_main_domain_boundary_perp_to_direction_measure(
            periodic_flow_direction);

    // Length of domain in the periodic direction
    double periodic_length = UU->primary_grid()->get_main_domain_max_coordinate(
                                 periodic_flow_direction) -
                             UU->primary_grid()->get_main_domain_min_coordinate(
                                 periodic_flow_direction);

    // Get the boundary name
    string boundary_name = periodic_flow_direction == 0   ? "right"
                           : periodic_flow_direction == 1 ? "top"
                                                          : "front";
    // Compute the flow rate
    double flow_rate = PAC_Misc::compute_flow_rate(UU, 0, boundary_name);

    // Get the imposed flow rate
    double imposed_periodic_flow_rate = primary_mesh->get_periodic_flow_rate();

    // Velocity correction and nominal flow rate
    double per_vel = -t_it->time_step() / (density * periodic_length);
    double nominal_flow_rate = per_vel * face_area;

    double delta_q = imposed_periodic_flow_rate - flow_rate;

    double ppd = delta_q / nominal_flow_rate;

    if (my_rank == is_master)
        MAC::out() << "         Periodic pressure drop = "
                   << MAC::doubleToString(ios::scientific, 14, ppd) << endl;

    const_cast<FV_Mesh *>(primary_mesh)->set_periodic_pressure_drop(ppd);
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::periodic_flow_rate_update(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: periodic_flow_rate_update");

    FV_Mesh const *primary_mesh = UU->primary_grid();
    size_t periodic_flow_direction =
        primary_mesh->get_periodic_flow_direction();

    // Area of boundary face perpendicular to periodic flow direction
    double face_area =
        primary_mesh->get_main_domain_boundary_perp_to_direction_measure(
            periodic_flow_direction);

    // Length of domain in the periodic direction
    double periodic_length = UU->primary_grid()->get_main_domain_max_coordinate(
                                 periodic_flow_direction) -
                             UU->primary_grid()->get_main_domain_min_coordinate(
                                 periodic_flow_direction);

    // Imposed flow rate
    double imposed_periodic_flow_rate = primary_mesh->get_periodic_flow_rate();

    // Velocity correction and nominal flow rate
    double per_vel = -t_it->time_step() / (density * periodic_length);
    double nominal_flow_rate = per_vel * face_area;

    // Compute flow rate difference and updated pressure drop
    string boundary_name = periodic_flow_direction == 0   ? "right"
                           : periodic_flow_direction == 1 ? "top"
                                                          : "front";
    double delta_q = imposed_periodic_flow_rate -
                     PAC_Misc::compute_flow_rate(UU, 0, boundary_name);

    double ppd = delta_q / nominal_flow_rate;
    if (my_rank == is_master)
        MAC::out() << "Periodic pressure drop = "
                   << MAC::doubleToString(ios::scientific, 14, ppd) << endl;
    const_cast<FV_Mesh *>(primary_mesh)->set_periodic_pressure_drop(ppd);

    // Correct velocity field
    UU->add_to_DOFs_value(periodic_flow_direction, 0, ppd * per_vel);
    GLOBAL_EQ->initialize_velocity();
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::compute_and_print_divu_norm(void)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: compute_and_print_divu_norm");

    size_t i, j, k;
    size_t_vector min_unknown_index(dim, 0);
    size_t_vector max_unknown_index(dim, 0);
    double dux, duy, dx, dy, dz, duz = 0.;
    double div_velocity = 0.;
    double cell_div = 0., max_divu = 0., divu = 0.;

    FV_SHIFT_TRIPLET shift = PP->shift_staggeredToCentered();
    for (size_t l = 0; l < dim; ++l)
        min_unknown_index(l) = PP->get_min_index_unknown_handled_by_proc(0, l);
    for (size_t l = 0; l < dim; ++l)
        max_unknown_index(l) = PP->get_max_index_unknown_handled_by_proc(0, l);

    for (i = min_unknown_index(0); i <= max_unknown_index(0); ++i)
    {
        dx = PP->get_cell_size(i, 0, 0);
        for (j = min_unknown_index(1); j <= max_unknown_index(1); ++j)
        {
            dy = PP->get_cell_size(j, 0, 1);
            if (dim == 2)
            {
                k = 0;

                // Divergence of u (x component)
                dux = UU->DOF_value(shift.i + i, j, k, 0, 0) -
                      UU->DOF_value(shift.i + i - 1, j, k, 0, 0);

                // Divergence of u (y component)
                duy = UU->DOF_value(i, shift.j + j, k, 1, 0) -
                      UU->DOF_value(i, shift.j + j - 1, k, 1, 0);

                cell_div = dux * dy + duy * dx;
                max_divu = MAC::max(MAC::abs(cell_div) / (dx * dy), max_divu);
                div_velocity += cell_div * cell_div / (dx * dy);
            }
            else
            {
                for (k = min_unknown_index(2); k <= max_unknown_index(2); ++k)
                {
                    dz = PP->get_cell_size(k, 0, 2);

                    // Divergence of u (x component)
                    dux = UU->DOF_value(shift.i + i, j, k, 0, 0) -
                          UU->DOF_value(shift.i + i - 1, j, k, 0, 0);

                    // Divergence of u (y component)
                    duy = UU->DOF_value(i, shift.j + j, k, 1, 0) -
                          UU->DOF_value(i, shift.j + j - 1, k, 1, 0);

                    // Divergence of u(z component)
                    duz = UU->DOF_value(i, j, shift.k + k, 2, 0) -
                          UU->DOF_value(i, j, shift.k + k - 1, 2, 0);

                    cell_div = dux * dy * dz + duy * dx * dz + duz * dx * dy;
                    max_divu =
                        MAC::max(MAC::abs(cell_div) / (dx * dy * dz), max_divu);
                    div_velocity += cell_div * cell_div / (dx * dy * dz);
                }
            }
        }
    }

    FV_Mesh const *primary_mesh = UU->primary_grid();
    double domain_measure =
        dim == 2
            ? primary_mesh->get_main_domain_boundary_perp_to_direction_measure(
                  0) *
                  primary_mesh
                      ->get_main_domain_boundary_perp_to_direction_measure(1)
            : primary_mesh->get_main_domain_boundary_perp_to_direction_measure(
                  0) *
                  (primary_mesh->get_main_domain_max_coordinate(2) -
                   primary_mesh->get_main_domain_min_coordinate(2));

    div_velocity = macCOMM->sum(div_velocity);
    div_velocity = MAC::sqrt(div_velocity);
    max_divu = macCOMM->max(max_divu);
    if (my_rank == is_master)
        MAC::out() << "         Norm of div(u) on grid: L2 = " << div_velocity
                   << "  Linf = " << max_divu << endl;
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::compute_flow_rate(
    FV_TimeIterator const *t_it, bool check_restart)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: compute_flow_rate");

    static size_t flow_rate_counter = 0;

    if (compute_flow_rate_on != "none")
    {
        if (!flow_rate_counter)
        {
            flow_rate =
                PAC_Misc::compute_flow_rate(UU, 0, compute_flow_rate_on);
            if (my_rank == is_master)
            {
                MAC::out() << "         Flow rate through "
                           << compute_flow_rate_on << " boundary = "
                           << MAC::doubleToString(ios::scientific, 6, flow_rate)
                           << endl;
                string filename = resultsDirectory + "/" + "flowrate_" +
                                  compute_flow_rate_on + ".res";
                if (check_restart)
                    GrainsExec::checkTime_outputFile(filename, t_it->time());
                else
                {
                    ofstream flfile(filename.c_str(), ios::app);
                    flfile << MAC::doubleToString(ios::scientific, 6,
                                                  t_it->time())
                           << " "
                           << MAC::doubleToString(ios::scientific, 6, flow_rate)
                           << endl;
                    flfile.close();
                }
            }
        }
        ++flow_rate_counter;
        if (flow_rate_counter == flow_rate_frequency)
            flow_rate_counter = 0;
    }
}




//---------------------------------------------------------------------------
void DLMFD_ProjectionNavierStokes::build_links_translation()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_ProjectionNavierStokes:: build_links_translation");

    UU->create_transproj_interpolation();
    PP->create_transproj_interpolation();
}
