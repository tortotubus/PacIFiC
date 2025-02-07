#include <DLMFD_FictitiousDomain.hh>
#include <FS_SolidPlugIn_BuilderFactory.hh>
#include <FV_Mesh.hh>
using namespace std;

//---------------------------------------------------------------------------
DLMFD_FictitiousDomain::DLMFD_FictitiousDomain(MAC_Object *a_owner,
                                               FV_DomainAndFields const *dom,
                                               MAC_ModuleExplorer const *exp) : MAC_Object(a_owner)
//--------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: DLMFD_FictitiousDomain");

   // Create solid readable files
   solidSolverType = "Grains3D";
   b_solidSolver_parallel = false;
   solidSolver_insertionFile = "Grains/Init/insert.xml";
   solidSolver_simulationFile = "Grains/Res/simul.xml";
   int error = 0;

   // Create rigid bodies object
   solidSolver = FS_SolidPlugIn_BuilderFactory::create(solidSolverType,
                                                       solidSolver_insertionFile, solidSolver_simulationFile, density, false,
                                                       b_restart, dom->primary_grid()->get_smallest_constant_grid_size(),
                                                       b_solidSolver_parallel, error);

   solidFluid_transferStream = NULL;
   solidSolver->getSolidBodyFeatures(solidFluid_transferStream);

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
                                                       MAC_ModuleExplorer const *exp)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: create");
   DLMFD_FictitiousDomain *result = new DLMFD_FictitiousDomain(a_owner, dom, exp);
   return (result);
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::do_one_inner_iteration()
//---------------------------------------------------------------------------
{

   update_rigid_bodies();
   run_DLMFD_UzawaSolver();
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::update_rigid_bodies()
//---------------------------------------------------------------------------
{
   sub_prob_number = 3;
   MAC::out() << "-----------------------------------------" << "-------------" << endl;
   MAC::out() << "Sub-problem " << sub_prob_number
              << " : Rigid Bodies updating" << endl;
   MAC::out() << "-----------------------------------------" << "-------------" << endl;
   allrigidbodies->update(*solidFluid_transferStream);
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::run_DLMFD_UzawaSolver()
//---------------------------------------------------------------------------
{
   sub_prob_number = 4;
   MAC::out() << "-----------------------------------------" << "-------------" << endl;
   MAC::out() << "Sub-problem " << sub_prob_number
              << " : DLMFD processing" << endl;
   MAC::out() << "-----------------------------------------" << "-------------" << endl;
   cout << "Hello World from Uzawa algorithm"
        << endl;
}
