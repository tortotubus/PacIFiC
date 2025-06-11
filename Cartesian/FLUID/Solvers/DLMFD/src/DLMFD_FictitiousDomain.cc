#include <DLMFD_FictitiousDomain.hh>
#include <FS_SolidPlugIn_BuilderFactory.hh>
#include <FV_Mesh.hh>
#include <MAC.hh>
#include <MAC_Error.hh>
#include <MAC_Exec.hh>
#include <math.h>
using namespace std;

bool DLMFD_FictitiousDomain::b_SecondOrderInterpol = false;
double DLMFD_FictitiousDomain::BoundaryPointsSpacing_coef = 1.;
doubleVector *DLMFD_FictitiousDomain::dbnull = new doubleVector(0, 0.);

//---------------------------------------------------------------------------
DLMFD_FictitiousDomain::DLMFD_FictitiousDomain(
    MAC_Object *a_owner, FV_DomainAndFields const *dom,
    MAC_ModuleExplorer const *exp, NavierStokes2FluidSolid transfert)
    : MAC_Object(a_owner), PAC_ComputingTime("Solver"), solidSolver(NULL),
      solidFluid_transferStream(NULL), transferString(NULL)
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: DLMFD_FictitiousDomain");

    // MPI data
    macCOMM = MAC_Exec::communicator();
    my_rank = macCOMM->rank();
    size_proc = macCOMM->nb_ranks();
    is_master = 0;

    // Timing routines
    if (my_rank == is_master)
    {
        CT_set_start();
        SCT_insert_app("Objects_Creation");
        SCT_set_start("Objects_Creation");
    }

    // -- Read the data from the explorer

    // Read the coupling scheme
    coupling_scheme = "Standard";
    if (exp->has_entry("Coupling_scheme"))
        coupling_scheme = exp->string_data("Coupling_scheme");
    if (coupling_scheme != "Standard" &&
        coupling_scheme != "PredictorCorrector")
    {
        string error_message = "   - Standard\n   - PredictorCorrector";
        MAC_Error::object()->raise_bad_data_value(exp, "Coupling_scheme",
                                                  error_message);
    }

    // Use 2nd order interpolation on boundary points
    if (exp->has_entry("DLMFD_BP_2ndOrder"))
        b_SecondOrderInterpol = exp->bool_data("DLMFD_BP_2ndOrder");

    // Whether to treat added mass explicitly or not
    b_explicit_added_mass = false;
    if (exp->has_entry("Explicit_added_mass"))
        b_explicit_added_mass = exp->bool_data("Explicit_added_mass");

    // Level of verbosity
    b_particles_verbose = true;
    if (exp->has_entry("Particles_verbose"))
        b_particles_verbose = exp->bool_data("Particles_verbose");

    // Explicit DLMFD treatment
    b_ExplicitDLMFD = false;
    if (exp->has_entry("ExplicitDLMFD"))
        b_ExplicitDLMFD = exp->bool_data("ExplicitDLMFD");

    // Are particles fixed
    are_particles_fixed = false;
    if (exp->has_entry("Particles_as_FixedObstacles"))
        are_particles_fixed = exp->bool_data("Particles_as_FixedObstacles");

    if (exp->has_entry("BPSpacingCoef"))
        DLMFD_FictitiousDomain::BoundaryPointsSpacing_coef =
            exp->double_data("BPSpacingCoef");

    // Output objects
    SolidSolverResultsDirectory = transfert.solid_resDir;

    // Instantiate discrete fields
    dim = 3;
    UU = transfert.UU;
    PP = transfert.PP;

    // Allocate constrained DOFs array in the constrained field
    UU->allocate_DLMFDconstrainedDOFs();

    // Set physical parameters
    rho_f = transfert.rho_f;
    gravity_vector = transfert.gravity_vector;
    split_gravity_vector = transfert.split_gravity_vector;

    // Create solid readable files
    solidSolverType = "Grains3D";
    b_solidSolver_parallel = false;
    b_correct_particle_acceleration = true;
    solidSolver_insertionFile = "Grains/Init/insert.xml";
    solidSolver_simulationFile = "Grains/Res/simul.xml";
    int error = 0;

    // Restart
    b_restart = transfert.b_restart;

    // Initiate the GRAINS plugins
    solidSolver = FS_SolidPlugIn_BuilderFactory::create(
        solidSolverType, solidSolver_insertionFile, solidSolver_simulationFile,
        rho_f, b_correct_particle_acceleration, b_restart,
        dom->primary_grid()->get_smallest_constant_grid_size(),
        b_solidSolver_parallel, error);

    // Create the Solid/Fluid transfer stream
    solidFluid_transferStream = new istringstream;
    solidSolver->getSolidBodyFeatures(solidFluid_transferStream);

    // Set critical distance
    double grid_size = UU->primary_grid()->get_smallest_constant_grid_size();
    critical_distance = sqrt(double(dim)) * grid_size;

    // Create the Rigid Bodies
    allrigidbodies = new DLMFD_AllRigidBodies(
        dim, 0., *solidFluid_transferStream, are_particles_fixed, UU, PP,
        critical_distance);

    if (exp->has_entry("Output_hydro_forceTorque"))
        allrigidbodies->set_b_output_hydro_forceTorque(
            exp->bool_data("Output_hydro_forceTorque"));

    allrigidbodies->set_output_frequency(transfert.output_frequency);

    allrigidbodies->allocate_translational_angular_velocity_array();

    // Set the coupling factor for each rigid body
    if (are_particles_fixed)
        allrigidbodies->set_coupling_factor(rho_f, b_explicit_added_mass);

    // Create the instances of resolution objects
    GLOBAL_EQ = transfert.GLOBAL_EQ;
    levelDiscrField = transfert.velocitylevelDiscrField;
    nb_levels = transfert.nb_levels;

    // Get the DLMFD convergence criterion and maximum iterations allowed
    Uzawa_DLMFD_precision = GLOBAL_EQ->get_DLMFD_convergence_criterion();
    Uzawa_DLMFD_maxiter = GLOBAL_EQ->get_DLMFD_maxiter();

    finalize_construction(exp);

    // Timing routines
    if (my_rank == is_master)
    {
        SCT_insert_app("update_rigid_bodies");
        SCT_insert_app("DLMFD_construction");
        SCT_insert_app("DLMFD_solving");
        SCT_get_elapsed_time("Objects_Creation");
    }

    // Do DLMFD before the projection step
    b_DLMFD_before_projection = false;
    if (exp->has_entry("DLMFD_start"))
        b_DLMFD_before_projection = exp->bool_data("DLMFD_start");
}




//---------------------------------------------------------------------------
DLMFD_FictitiousDomain::~DLMFD_FictitiousDomain(void)
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: ~DLMFD_FictitiousDomain");

    if (solidSolver)
        delete solidSolver;
}




//---------------------------------------------------------------------------
DLMFD_FictitiousDomain *DLMFD_FictitiousDomain::create(
    MAC_Object *a_owner, FV_DomainAndFields const *dom,
    MAC_ModuleExplorer const *exp, NavierStokes2FluidSolid transfert)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: create");

    DLMFD_FictitiousDomain *dlmfd_solver =
        new DLMFD_FictitiousDomain(a_owner, dom, exp, transfert);
    return (dlmfd_solver);
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::do_one_inner_iteration(FV_TimeIterator const *t_it,
                                                    size_t &sub_prob_number)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: do_one_inner_iteration");

    // Newtons's law -- Prediction step
    if (my_rank == is_master)
        SCT_set_start("update_rigid_bodies");

    update_rigid_bodies(t_it, sub_prob_number);

    if (my_rank == is_master)
        SCT_get_elapsed_time("update_rigid_bodies");

    // DLMFD process -- Correction step
    run_DLMFD_UzawaSolver(t_it, sub_prob_number);

    // Sum DLM on master process for hydrodynamic force & torque output
    allrigidbodies->sum_DLM_hydrodynamic_force_output(b_restart);
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::do_after_inner_iterations_stage(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain::do_after_inner_iterations_stage");

    if (!are_particles_fixed)
        if (my_rank == is_master)
        {
            // Update velocity on the solid side
            allrigidbodies->particles_velocities_output(velocitiesVecGrains);
            solidSolver->UpdateParticlesVelocities(velocitiesVecGrains,
                                                   b_explicit_added_mass);
        }

    // Write output files
    if (my_rank == is_master && !are_particles_fixed)
        allrigidbodies->particles_features_output(
            SolidSolverResultsDirectory + "/", b_restart, t_it->time());

    allrigidbodies->particles_hydrodynamic_force_output(
        SolidSolverResultsDirectory + "/", b_restart, t_it->time(),
        t_it->time_step(), rho_f, Iw_Idw);
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::do_before_time_stepping(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: do_before_time_stepping");

    Paraview_saveMultipliers_pvd << "<?xml version=\"1.0\"?>" << endl;
    Paraview_saveMultipliers_pvd
        << "<VTKFile type=\"Collection\" version=\"0.1\""
        << " byte_order=\"LittleEndian\">" << endl;
    Paraview_saveMultipliers_pvd << "<Collection>" << endl;

    // Hydrodynamic force and torque
    allrigidbodies->sum_DLM_hydrodynamic_force_output(b_restart);
    allrigidbodies->particles_hydrodynamic_force_output(
        SolidSolverResultsDirectory + "/", b_restart, t_it->time(),
        t_it->time_step(), rho_f, Iw_Idw);

    // Check that corresponding post-processing writer exists in Grains3D
    solidSolver->checkParaviewPostProcessing(SolidSolverResultsDirectory);
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::do_additional_savings(
    int const &cycleNumber, FV_TimeIterator const *t_it,
    const double &translated_distance, const size_t &translation_direction)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: do_additional_savings");

    set_Paraview_translated_distance_vector(translated_distance,
                                            translation_direction);

    string filename;

    filename = "saveMultipliersT" + MAC::intToString(cycleNumber);

    if (my_rank == is_master)
    {
        Paraview_saveMultipliers_pvd << "<DataSet timestep=\"" << t_it->time()
                                     << "\" " << "group=\"\" part=\"0\" file=\""
                                     << filename << ".pvtu\"/>" << endl;

        ofstream f(
            (SolidSolverResultsDirectory + "/saveMultipliers.pvd").c_str(),
            ios::out);

        string str = Paraview_saveMultipliers_pvd.str();
        f << str;

        f << "</Collection>" << endl;
        f << "</VTKFile>" << endl;
        f.close();

        write_PVTU_multiplier_file(filename);
    }

    allrigidbodies->output_DLMFDPoints_PARAVIEW(
        SolidSolverResultsDirectory + "/" + filename,
        &Paraview_translated_distance_vector, true);
    if (my_rank == is_master)
        solidSolver->saveResults("", t_it->time(), cycleNumber);

    // Elapsed time by sub-problems
    if (my_rank == is_master)
    {
        double cputime = CT_get_elapsed_time();
        MAC::out() << endl << "DLMFD problem" << endl;
        write_elapsed_time_smhd(MAC::out(), cputime, "Computation time");
        SCT_get_summary(MAC::out(), cputime);
    }
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::write_PVTU_multiplier_file(
    string const &filename) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("MY_FictitiousDomain:: write_PVTU_multiplier_file");

    ofstream f((SolidSolverResultsDirectory + "/" + filename + ".pvtu").c_str(),
               ios::out);
    f << "<?xml version=\"1.0\"?>" << endl;
    f << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" "
      << "byte_order=\"LittleEndian\">" << endl;
    f << "<PUnstructuredGrid GhostLevel=\"0\">" << endl;
    f << "<PPoints>" << endl;
    f << "<PDataArray NumberOfComponents=\"3\" type=\"Float32\" "
         "format=\"ascii\">"
      << endl;
    f << "</PDataArray>" << endl;
    f << "</PPoints>" << endl;
    f << "<PCells>" << endl;
    f << "<PDataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\">"
      << endl;
    f << "</PDataArray>" << endl;
    f << "<PDataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\">"
      << endl;
    f << "</PDataArray>" << endl;
    f << "<PDataArray Name=\"types\" type=\"Int32\" format=\"ascii\">" << endl;
    f << "</PDataArray>" << endl;
    f << "</PCells>" << endl;
    for (size_t i = 0; i < size_proc; ++i)
    {
        f << "<Piece Source=\"" << filename << "_" << i << ".vtu\">" << endl;
        f << "</Piece>" << endl;
    }
    f << "</PUnstructuredGrid>" << endl;
    f << "</VTKFile>" << endl;
    f.close();
}




//---------------------------------------------------------------------------
double DLMFD_FictitiousDomain::Compute_distance_to_bottom(
    const double &coordinate, const size_t &direction) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: Compute_distance_to_bottom");

    return (allrigidbodies->compute_minimum_distance_to_bottom(coordinate,
                                                               direction));
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::translate_all(const geomVector &translation_vector,
                                           const size_t &translation_direction)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: translate_all");

    allrigidbodies->translate_geometricBoundaries(translation_vector,
                                                  translation_direction);

    double trans_dist = UU->primary_grid()->get_translation_distance();
    if (my_rank == is_master)
        solidSolver->setParaviewPostProcessingTranslationVector(
            translation_direction == 0 ? -trans_dist : 0.,
            translation_direction == 1 ? -trans_dist : 0.,
            translation_direction == 2 ? -trans_dist : 0.);
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::set_Paraview_translated_distance_vector(
    const double &translated_distance, const size_t &translation_direction)
//---------------------------------------------------------------------------
{
    MAC_LABEL(
        "DLMFD_FictitiousDomain:: set_Paraview_translated_distance_vector");

    Paraview_translated_distance_vector(translation_direction) =
        translated_distance;
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::setParaviewPostProcessingTranslationVector(
    const double &tvx, const double &tvy, const double &tvz)
//---------------------------------------------------------------------------
{
    MAC_LABEL(
        "DLMFD_FictitiousDomain:: setParaviewPostProcessingTranslationVector");

    solidSolver->setParaviewPostProcessingTranslationVector(tvx, tvy, tvz);
}




//---------------------------------------------------------------------------
bool const DLMFD_FictitiousDomain::get_explicit_DLMFD() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain::get_explicit_DLMFD");

    return b_ExplicitDLMFD;
}




//---------------------------------------------------------------------------
bool const DLMFD_FictitiousDomain::get_DLMFD_before_projection() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain::get_DLMFD_before_projection");

    return b_DLMFD_before_projection;
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::finalize_construction(
    MAC_ModuleExplorer const *exp)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain::finalize_construction");

    macCOMM->barrier();

    allrigidbodies->read_particles_outputs(exp);

    // Resize the Paraview translated distance vector
    Paraview_translated_distance_vector.resize(dim);

    // Set Iw_Idw and velocitiesVecGrains
    if (my_rank == is_master)
    {
        if (allrigidbodies->is_hydro_forceTorque_postprocessed())
        {
            vector<double> work(6, 0.);
            Iw_Idw =
                new vector<vector<double>>(allrigidbodies->get_npart(), work);
        }

        velocitiesVecGrains.reserve(allrigidbodies->get_npart());
        for (size_t i = 0; i < allrigidbodies->get_npart(); i++)
            velocitiesVecGrains.push_back(vector<double>(6, 0.));
    }

    macCOMM->barrier();
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::update_rigid_bodies(FV_TimeIterator const *t_it,
                                                 size_t &sub_prob_number)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: update_rigid_bodies");

    if (!are_particles_fixed)
    {
        if (my_rank == is_master)
        {
            MAC::out()
                << "------------------------------------------------------"
                << endl;
            MAC::out() << "Sub-problem " << sub_prob_number
                       << " : Solid problem -- Prediction" << endl;
            MAC::out()
                << "------------------------------------------------------"
                << endl;
        }

        solidFluid_transferStream = new istringstream;
        if (my_rank == is_master)
        {
            // Update the Rigid Bodies (Prediction problem)
            solidSolver->Simulation(t_it->time_step(), true,
                                    coupling_scheme == "PredictorCorrector", 1.,
                                    b_explicit_added_mass);
        }
        solidSolver->getSolidBodyFeatures(solidFluid_transferStream);
        if (my_rank == is_master)
            MAC::out() << "Solid components written in stream by solid solver"
                       << endl;

        ++sub_prob_number;

        macCOMM->barrier();
    }
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::run_DLMFD_UzawaSolver(FV_TimeIterator const *t_it,
                                                   size_t &sub_prob_number)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: run_DLMFD_UzawaSolver");

    if (my_rank == is_master)
    {
        if (are_particles_fixed)
        {
            MAC::out()
                << "------------------------------------------------------"
                << endl;
            MAC::out() << "Sub-problem " << sub_prob_number
                       << " : DLMFD problem" << endl;
            MAC::out()
                << "------------------------------------------------------"
                << endl;
        }
        else
        {
            MAC::out()
                << "------------------------------------------------------"
                << endl;
            MAC::out() << "Sub-problem " << sub_prob_number
                       << " : DLMFD problem -- Correction" << endl;
            MAC::out()
                << "------------------------------------------------------"
                << endl;
        }
    }

    // Initialize the DLMFD correction problem
    if (my_rank == is_master)
        SCT_set_start("DLMFD_construction");

    DLMFD_construction(t_it);

    if (my_rank == is_master)
        SCT_get_elapsed_time("DLMFD_construction");

    // Solve the DLMFD correction problem
    if (my_rank == is_master)
        SCT_set_start("DLMFD_solving");

    DLMFD_solving(t_it);

    if (my_rank == is_master)
        SCT_get_elapsed_time("DLMFD_solving");

    ++sub_prob_number;

    if (b_DLMFD_before_projection)
    {
        if (my_rank == is_master)
            MAC::out() << "Uzawa problem completed" << endl;
    }
    else if (my_rank == is_master)
        MAC::out() << "Uzawa problem completed" << endl << endl;
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::DLMFD_construction(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: DLMFD_construction");

    // Update the rigid bodies features before properly solve the 
    // Uzawa algorithm
    if (!are_particles_fixed)
    {
        // Initialize the constrained DOFs array in the constrained field
        UU->initialize_DLMFDconstrainedDOFs();

        allrigidbodies->set_ptr_constrained_field(UU);
        allrigidbodies->set_ptr_constrained_field_in_all_particles();
        allrigidbodies->set_ttran_ncomp(UU->nb_components());

        // Update the rigid bodies
        allrigidbodies->update(*solidFluid_transferStream);

        // Set valid points
        allrigidbodies->set_all_MAC(critical_distance);
        allrigidbodies->eraseCriticalDLMFDPoints(t_it->time(),
                                                 critical_distance);

        // Set the onProc IDs
        allrigidbodies->set_listIdOnProc();

        // Set the fluid/solid coupling factor
        allrigidbodies->set_coupling_factor(rho_f, b_explicit_added_mass);

        // Set points infos
        allrigidbodies->set_points_infos();

        // Fill DLMFD vectors
        allrigidbodies->check_allocation_DLMFD_Cvectors();
        allrigidbodies->fill_DLMFD_Cvectors();
    }

    // Nullify Uzawa vectors in particles
    allrigidbodies->nullify_all_Uzawa_vectors();

    macCOMM->barrier();
    if (my_rank == is_master)
        MAC::out()
            << "Update of DLMFD_AllRigidBodies completed on all processes"
            << endl;
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::DLMFD_solving(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: DLMFD_solving");

    double alpha = 0., nr = 0., nrkm1 = 0., nwx = 0., beta = 0.;
    int iter = 0;

    //-- The implementation is formerly detailed in the paper
    // PeliGRIFF, a parallel DEM-DLM/FD direct numerical simulation tool
    // for 3D particulate flows

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // IMPORTANT REMARKS:
    // * During the iterative process, the fluid velocity field UU contains the
    // values of tu (VEC_t at the matrix level) whereas its actual values
    // are stored in VEC_U.
    // Once the algorithm has converged, we ultimately copy back VEC_U to the
    // field to recover the proper values.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // ----- INITIALISATION OF DLM/FD-UZAWA SOLVING ALGORITHM -----

    struct timezone tz;
    struct timeval total_start, start, start_iter;
    struct timeval total_end, end, end_iter;
    double elapsed_time = 0.;
    if (my_rank == is_master)
        gettimeofday(&total_start, &tz);

    // Assemble fluid velocity RHS vector i.e. compute fu = A.u(n-1) where A is
    // the velocity unsteady matrix and U(n-1) the velocity at previous
    // sub-problem
    // Pay Attention: VEC_rhs_A_Velocity must be set.
    GLOBAL_EQ->updateFluid_DLMFD_rhs();

    // For all particles, compute qtran = -<lambda,V>_P ( = qu )
    // and qrot = -<lambda,xi^GM>_P at the particles/field level
    allrigidbodies->compute_all_Qu(true);
    allrigidbodies->compute_all_Qrot(true);

    // Calculate sum of qtran and qtrot of each process
    // (MPI Reduction operation) on master process
    calculate_ReductionFor_qtranAndqrot(t_it);

    // Complete the initialization of DLM/FD-Uzawa solving algorithm
    // do, on each particle of master process:
    //    Solve (1-rho_f/rho_s)*M*t_tran = q_tran
    //     and (1-rho_f/rho_s)*I*t_rot = q_rot.
    // then for each particle, broadcast t_tran and t_rot from master process
    // on each process. At the end of the method, each instance of
    // a particle on every process has the same (and right) t_tran and t_rot.
    velocity_broadcast_andUpdate_First(t_it);

    // For all particles and all processes, do: U = t_tran and omega = t_rot
    allrigidbodies->update_ParticlesVelocities_afterBcast_of_T();

    // Nullify Qu vector
    GLOBAL_EQ->nullify_QUvector();

    // Compute the DLM right hand side of the momentum equations
    // quf = -<lambda,v> at the particles/field level
    compute_fluid_LBD_rhs(t_it, true);

    // Solve the fluid system at the matrix level
    // * Compute fuf = (ro/dt)*U(n-1)-<lambda,v> = fu + quf
    // * Solve A.u = fuf
    GLOBAL_EQ->solve_FluidVel_DLMFD_Init(t_it->time());

    // Transfer velocity unknown vector u in the UU (velocity) field
    UU->update_free_DOFs_value(levelDiscrField, GLOBAL_EQ->get_solution_U());

    // Compute the residual vector x = <alpha,u-(U+omega^GM)>_P
    allrigidbodies->compute_x_residuals_Velocity();

    // Set r = - x and w = r
    allrigidbodies->compute_r_and_w_FirstUzawaIteration();

    // Compute nr = r.r
    nr = allrigidbodies->compute_r_dot_r();

    if (my_rank == is_master && b_particles_verbose)
        MAC::out() << "Residuals = "
                   << MAC::doubleToString(ios::scientific, 14, sqrt(nr)) << " "
                   << MAC::doubleToString(ios::scientific, 14, nr)
                   << "  Iterations = 0" << endl;

    while (sqrt(nr) > Uzawa_DLMFD_precision && iter < Uzawa_DLMFD_maxiter)
    {
        iter++;

        // For all particles, compute qtran = -<w,V>_P ( = qu )
        // and qrot = -<w,xi^GM>_P at the particles/field level
        allrigidbodies->compute_all_Qu(false);
        allrigidbodies->compute_all_Qrot(false);

        // Calculate sum of qtran and qtrot of each process
        // (MPI Reduction operation) on master process
        calculate_ReductionFor_qtranAndqrot(t_it);

        // do, on each particle of master process:
        //    Solve (1-rho_f/rho_s)*M*t_tran = q_tran
        //     and (1-rho_f/rho_s)*I*t_rot = q_rot.
        // then for each particle, broadcast t_tran and t_rot from master
        // process on each process. At the end of the method, each instance of
        // a particle on every process has the same (and right) t_tran and
        // t_rot.
        velocity_broadcast_andUpdateInOneIt(t_it);

        // Nullify Qu vector
        GLOBAL_EQ->nullify_QUvector();

        // Add the DLM right hand side of the momentum equations
        // <w,v> at the particles/field level to quf
        compute_fluid_LBD_rhs(t_it, false);

        // Solve the fluid system at the matrix level
        // i.e. solve A.tu = quf = <w,v>
        GLOBAL_EQ->solve_FluidVel_DLMFD_Iter(t_it->time());

        // Transfer velocity unknown vector tu in the UU (velocity) field
        UU->update_free_DOFs_value(levelDiscrField, GLOBAL_EQ->get_tVector_U());

        // Compute the residual vector x = <alpha,tu-(tU+tomega^GM)>_P
        allrigidbodies->compute_x_residuals_Velocity();

        // Compute nrkm1 = r.r^(iter-1)
        nrkm1 = nr;

        // Compute the w.x dot product
        nwx = allrigidbodies->compute_w_dot_x();

        // Compute alpha=nrkm1/nwx
        alpha = nrkm1 / nwx;

        // Update lambda-=alpha.w, r-=alpha.x
        allrigidbodies->update_lambda_and_r(alpha);

        // Update particles and fluid velocities
        allrigidbodies->update_ParticlesVelocities(alpha);

        // Update u+=alpha.t
        GLOBAL_EQ->update_FluidVel_OneUzawaIter(alpha);

        // Compute nr = r.r
        nr = allrigidbodies->compute_r_dot_r();

        // Compute beta=nr/nrkm1
        beta = nr / nrkm1;

        // Update w = r + beta.w
        allrigidbodies->update_w(beta);

        if (my_rank == is_master && b_particles_verbose)
            MAC::out() << "Residuals = "
                       << MAC::doubleToString(ios::scientific, 14, sqrt(nr))
                       << " " << MAC::doubleToString(ios::scientific, 14, nr)
                       << "  Iterations = " << iter << endl;
    }

    if (my_rank == is_master)
        MAC::out() << "Residuals = "
                   << MAC::doubleToString(ios::scientific, 14, sqrt(nr))
                   << "  Nb iterations = " << iter << endl;

    // Copy velocity of all particles on master
    Set_Velocity_AllParticles_Master(t_it);

    // Copy back the fluid velocity values to the field
    UU->update_free_DOFs_value(levelDiscrField, GLOBAL_EQ->get_solution_U());

    // Transfer values from the level of computation to additional levels if
    // needed
    if (!b_DLMFD_before_projection)
        for (size_t i = 1; i < nb_levels; i++)
            UU->copy_DOFs_value(levelDiscrField, i);

    if ((my_rank == is_master) && (b_particles_verbose))
    {
        gettimeofday(&total_end, &tz);
        double total_elapsed_time =
            total_end.tv_sec - total_start.tv_sec +
            double(total_end.tv_usec - total_start.tv_usec) / 1e6;
        MAC::out() << "total time = " << total_elapsed_time << endl;
    }

    // Store DLMFD forcing term
    if (b_ExplicitDLMFD)
    {
        GLOBAL_EQ->nullify_QUvector();
        allrigidbodies->compute_fluid_DLMFD_explicit(GLOBAL_EQ, true);
        GLOBAL_EQ->store_DLMFD_rhs();
    }
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::calculate_ReductionFor_qtranAndqrot(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: calculate_ReductionFor_qtranAndqrot");

    size_t i, k, indexcomp, nsolid;
    geomVector tempqtran(dim), tempqrot(dim);

    if (my_rank != is_master)
    {
        list<int> const *plistShared = allrigidbodies->get_SharedOnProc();
        list<int>::const_iterator il;
        nsolid = plistShared->size();
        doubleVector sendbuf(2 * nsolid * dim);
        i = 0;
        for (il = plistShared->begin(); il != plistShared->end(); il++, ++i)
        {
            tempqtran = allrigidbodies->get_Qu(*il);
            tempqrot = allrigidbodies->get_Qrot(*il);

            for (indexcomp = 0; indexcomp < dim; ++indexcomp)
            {
                sendbuf(2 * i * dim + indexcomp) = tempqtran(indexcomp);
                sendbuf((2 * i + 1) * dim + indexcomp) = tempqrot(indexcomp);
            }
        }
        macCOMM->send(is_master, sendbuf);
    }
    else
    {
        vector<size_t_vector> const *pallshared =
            allrigidbodies->get_v_AllSharedOnProcs();
        for (k = 1; k < size_proc; ++k)
        {
            nsolid = (*pallshared)[k].size();
            doubleVector recvbuf(2 * nsolid * dim);
            macCOMM->receive(k, recvbuf);
            for (i = 0; i < nsolid; ++i)
            {
                for (indexcomp = 0; indexcomp < dim; ++indexcomp)
                {
                    tempqtran(indexcomp) = recvbuf(2 * i * dim + indexcomp);
                    tempqrot(indexcomp) =
                        recvbuf((2 * i + 1) * dim + indexcomp);
                }
                allrigidbodies->add_to_Qu((*pallshared)[k](i), tempqtran);
                allrigidbodies->add_to_Qrot((*pallshared)[k](i), tempqrot);
            }
        }
    }
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::velocity_broadcast_andUpdate_First(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: velocity_broadcast_andUpdateInOneIt");

    // Done by all processes
    allrigidbodies->completeFirstUzawaIteration_Velocity(
        rho_f, t_it->time_step(), split_gravity_vector, gravity_vector);

    // Broadcast tU and trot of shared solid components from master to
    // other processes
    Broadcast_tVectors_sharedParticles_MasterToAll(t_it);
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::velocity_broadcast_andUpdateInOneIt(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: velocity_broadcast_andUpdateInOneIt");

    // Done by all processes
    allrigidbodies->solve_Particles_OneUzawaIter_Velocity(rho_f,
                                                          t_it->time_step());

    // Broadcast tU and trot of shared solid components from master to
    // other processes
    Broadcast_tVectors_sharedParticles_MasterToAll(t_it);
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::Broadcast_tVectors_sharedParticles_MasterToAll(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain:: "
              "Broadcast_tVectors_sharedParticles_MasterToAll");

    geomVector ttran(dim);
    geomVector trot(dim);
    size_t i, k, nsolid, indexcomp;

    if (my_rank == is_master)
    {
        vector<size_t_vector> const *pallshared =
            allrigidbodies->get_v_AllSharedOnProcs();
        for (k = 1; k < size_proc; ++k)
        {
            nsolid = (*pallshared)[k].size();
            doubleVector sendbuf(2 * nsolid * dim);
            for (i = 0; i < nsolid; ++i)
            {
                ttran = allrigidbodies->get_Tu((*pallshared)[k](i));
                trot = allrigidbodies->get_Trot((*pallshared)[k](i));

                for (indexcomp = 0; indexcomp < dim; ++indexcomp)
                {
                    sendbuf(2 * i * dim + indexcomp) = ttran(indexcomp);
                    sendbuf((2 * i + 1) * dim + indexcomp) = trot(indexcomp);
                }
            }
            macCOMM->send(k, sendbuf);
        }
    }
    else
    {
        list<int> const *pshared = allrigidbodies->get_SharedOnProc();
        list<int>::const_iterator il;

        nsolid = pshared->size();

        doubleVector recvbuf(2 * nsolid * dim);
        macCOMM->receive(is_master, recvbuf);
        i = 0;
        for (il = pshared->begin(); il != pshared->end(); il++, ++i)
        {
            for (indexcomp = 0; indexcomp < dim; ++indexcomp)
            {
                ttran(indexcomp) = recvbuf(2 * i * dim + indexcomp);
                trot(indexcomp) = recvbuf((2 * i + 1) * dim + indexcomp);
            }
            allrigidbodies->set_Tu(*il, ttran);
            allrigidbodies->set_Trot(*il, trot);
        }
    }
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::compute_fluid_LBD_rhs(FV_TimeIterator const *t_it,
                                                   bool init)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain::compute_fluid_LBD_rhs");

    allrigidbodies->compute_fluid_rhs(GLOBAL_EQ, init);
}




//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::Set_Velocity_AllParticles_Master(
    FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_FictitiousDomain::Set_Velocity_AllParticles_Master");

    size_t i, k, j, indexcomp, npart;
    geomVector vtran(dim), vrot(dim);
    size_t nallpart = allrigidbodies->get_npart();

    // Number of particles sent by each process
    list<int> const *plistOnProc = allrigidbodies->get_onProc();
    list<int>::const_iterator il;
    npart = plistOnProc->size();
    intVector numberOfParticlesSent(size_proc);
    macCOMM->gather(npart, numberOfParticlesSent, is_master);

    if (my_rank != is_master)
    {
        intVector sendbuf_id(npart);
        doubleVector sendbuf_vel(2 * npart * dim);
        i = 0;
        for (il = plistOnProc->begin(); il != plistOnProc->end(); il++, ++i)
        {
            sendbuf_id(i) = *il;
            vtran = allrigidbodies->get_translational_velocity(*il);
            vrot = allrigidbodies->get_angular_velocity_3D(*il);
            for (indexcomp = 0; indexcomp < dim; ++indexcomp)
            {
                sendbuf_vel(2 * i * dim + indexcomp) = vtran(indexcomp);
                sendbuf_vel((2 * i + 1) * dim + indexcomp) = vrot(indexcomp);
            }
        }
        macCOMM->send(is_master, sendbuf_id);
        macCOMM->send(is_master, sendbuf_vel);
    }
    else
    {
        for (k = 1; k < size_proc; ++k)
        {
            npart = numberOfParticlesSent(k);
            intVector recvbuf_id(npart);
            doubleVector recvbuf_vel(2 * npart * dim);
            macCOMM->receive(k, recvbuf_id);
            macCOMM->receive(k, recvbuf_vel);
            for (i = 0; i < npart; ++i)
            {
                for (indexcomp = 0; indexcomp < dim; ++indexcomp)
                {
                    vtran(indexcomp) = recvbuf_vel(2 * i * dim + indexcomp);
                    vrot(indexcomp) =
                        recvbuf_vel((2 * i + 1) * dim + indexcomp);
                }
                allrigidbodies->set_translational_velocity(recvbuf_id(i),
                                                           vtran);
                allrigidbodies->set_angular_velocity_3D(recvbuf_id(i), vrot);
            }
        }
    }
}
