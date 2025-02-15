#include <DLMFD_FictitiousDomain.hh>
#include <FS_SolidPlugIn_BuilderFactory.hh>
#include <FV_Mesh.hh>
#include <MAC_Exec.hh>
#include <MAC.hh>
using namespace std;

//---------------------------------------------------------------------------
DLMFD_FictitiousDomain::DLMFD_FictitiousDomain(MAC_Object *a_owner,
                                               FV_DomainAndFields const *dom,
                                               MAC_ModuleExplorer const *exp,
                                               NavierStokes2FluidSolid transfert) : MAC_Object(a_owner)
//--------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: DLMFD_FictitiousDomain");

   // MPI data
   pelCOMM = MAC_Exec::communicator();
   size_proc = pelCOMM->nb_ranks();

   // Output objects
   SolidSolverResultsDirectory = transfert.solid_resDir;

   // Instantiate discrete fields
   UU = dom->discrete_field("velocity");
   PP = dom->discrete_field("pressure");

   // Create solid readable files
   solidSolverType = "Grains3D";
   b_solidSolver_parallel = false;
   solidSolver_insertionFile = "Grains/Init/insert.xml";
   solidSolver_simulationFile = "Grains/Res/simul.xml";
   int error = 0;

   // Create the Solid/Fluid transfer stream
   solidSolver = FS_SolidPlugIn_BuilderFactory::create(solidSolverType,
                                                       solidSolver_insertionFile, solidSolver_simulationFile, density, false,
                                                       b_restart, dom->primary_grid()->get_smallest_constant_grid_size(),
                                                       b_solidSolver_parallel, error);

   solidFluid_transferStream = NULL;
   solidSolver->getSolidBodyFeatures(solidFluid_transferStream);

   // Create the Rigid Bodies
   allrigidbodies = new DLMFD_AllRigidBodies(dim, *solidFluid_transferStream,
                                             are_particles_fixed, UU, PP);
}

//---------------------------------------------------------------------------
DLMFD_FictitiousDomain::~DLMFD_FictitiousDomain(void)
//--------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: ~DLMFD_FictitiousDomain");
}

//---------------------------------------------------------------------------
DLMFD_FictitiousDomain *DLMFD_FictitiousDomain::create(MAC_Object *a_owner,
                                                       FV_DomainAndFields const *dom,
                                                       MAC_ModuleExplorer const *exp,
                                                       NavierStokes2FluidSolid transfert)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: create");

   DLMFD_FictitiousDomain *dlmfd_solver = new DLMFD_FictitiousDomain(a_owner, dom, exp, transfert);
   return (dlmfd_solver);
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::do_one_inner_iteration(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: do_one_inner_iteration");

   update_rigid_bodies(t_it);
   run_DLMFD_UzawaSolver(t_it);
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::do_additional_savings(int const &cycleNumber,
                                                   FV_TimeIterator const *t_it)
// const double &translated_distance,
// const size_t &translation_direction)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: do_additional_savings");

   string filename;

   filename = "saveMultipliersT" + MAC::intToString(cycleNumber);

   Paraview_saveMultipliers_pvd << "<DataSet timestep=\"" << t_it->time()
                                << "\" " << "group=\"\" part=\"0\" file=\"" << filename << ".pvtu\"/>"
                                << endl;

   ofstream f((SolidSolverResultsDirectory + "/saveMultipliers.pvd").c_str(), ios::out);

   string str = Paraview_saveMultipliers_pvd.str();
   for (string::iterator it = str.begin(); it != str.end(); ++it)
      f << *it;

   f << "</Collection>" << endl;
   f << "</VTKFile>" << endl;
   f.close();

   write_PVTU_multiplier_file(filename);

   allrigidbodies->output_DLMFDPoints_PARAVIEW(SolidSolverResultsDirectory + "/" + filename,
                                               &Paraview_translated_distance_vector,
                                               true, pelCOMM->rank());
   
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::write_PVTU_multiplier_file(string const &filename) const
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
   f << "<PDataArray NumberOfComponents=\"3\" type=\"Float32\" format=\"ascii\">"
     << endl;
   f << "</PDataArray>" << endl;
   f << "</PPoints>" << endl;
   f << "<PCells>" << endl;
   f << "<PDataArray Name=\"connectivity\" type=\"Int32\" format=\"ascii\">"
     << endl;
   f << "</PDataArray>" << endl;
   f << "<PDataArray Name=\"offsets\" type=\"Int32\" format=\"ascii\">" << endl;
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
void DLMFD_FictitiousDomain::update_rigid_bodies(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: update_rigid_bodies");

   sub_prob_number = 3;
   MAC::out() << "-----------------------------------------" << "-------------" << endl;
   MAC::out() << "Sub-problem " << sub_prob_number
              << " : Rigid Bodies updating" << endl;
   MAC::out() << "-----------------------------------------" << "-------------" << endl;

   // Update the Rigid Bodies (Prediction problem)
   solidSolver->Simulation(t_it->time_step());
   MAC::out() << "Solid components written in stream by solid solver" << endl;
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::run_DLMFD_UzawaSolver(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: run_DLMFD_UzawaSolver");

   sub_prob_number = 4;
   MAC::out() << "-----------------------------------------" << "-------------" << endl;
   MAC::out() << "Sub-problem " << sub_prob_number
              << " : DLMFD solving" << endl;
   MAC::out() << "-----------------------------------------" << "-------------" << endl;

   // Initialize the DLMFD correction problem
   DLMFD_construction(t_it);

   // Solve the DLMFD correction problem
   cout << "SOLVING THE DLMFD PROLEM" << endl;
   // allrigidbodies->update_RB_position_and_velocity()
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::DLMFD_construction(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: initialize_DLMFD_problem");

   double critical_distance = 1.0e-3;
   allrigidbodies->set_all_points(critical_distance);
   cout << "INITIALIZING THE DLMFD PROLEM" << endl;
}
