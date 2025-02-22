#include <DLMFD_FictitiousDomain.hh>
#include <FS_SolidPlugIn_BuilderFactory.hh>
#include <FV_Mesh.hh>
#include <MAC_Exec.hh>
#include <MAC.hh>
#include <math.h>
using namespace std;

bool DLMFD_FictitiousDomain::b_SecondOrderInterpol = false;

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
   rank = pelCOMM->rank();
   master = 0;

   // Output objects
   SolidSolverResultsDirectory = transfert.solid_resDir;

   // Instantiate discrete fields
   dim = 3;
   UU = dom->discrete_field("velocity");
   PP = dom->discrete_field("pressure");

   // Set physical parameters
   rho_f = transfert.rho_f;
   gravity_vector = transfert.gravity_vector;
   split_gravity_vector = transfert.split_gravity_vector;

   // Create solid readable files
   solidSolverType = "Grains3D";
   b_solidSolver_parallel = false;
   solidSolver_insertionFile = "Grains/Init/insert.xml";
   solidSolver_simulationFile = "Grains/Res/simul.xml";
   int error = 0;

   // Create the Solid/Fluid transfer stream
   solidSolver = FS_SolidPlugIn_BuilderFactory::create(solidSolverType,
                                                       solidSolver_insertionFile, solidSolver_simulationFile, rho_f, false,
                                                       b_restart, dom->primary_grid()->get_smallest_constant_grid_size(),
                                                       b_solidSolver_parallel, error);

   solidFluid_transferStream = NULL;
   solidSolver->getSolidBodyFeatures(solidFluid_transferStream);

   // Create the Rigid Bodies
   allrigidbodies = new DLMFD_AllRigidBodies(dim, pelCOMM, *solidFluid_transferStream,
                                             are_particles_fixed, UU, PP);

   // Create the instance of DLMFD_ProjectionNavierStokesSystem
   GLOBAL_EQ = transfert.GLOBAL_EQ;

   // Use 2nd order interpolation on boundary points
   if (exp->has_entry("DLMFD_BP_2ndOrder"))
      b_SecondOrderInterpol = exp->bool_data("DLMFD_BP_2ndOrder");

   // Whether to treat added mass explicitly or not
   b_explicit_added_mass = false;
   if (exp->has_entry("Explicit_added_mass"))
      b_explicit_added_mass = exp->bool_data("Explicit_added_mass");

   // Set the coupling factor for each rigid body
   allrigidbodies->set_coupling_factor(transfert.rho_f, b_explicit_added_mass);
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

   // Newtons's law -- Prediction step
   update_rigid_bodies(t_it);

   // DLMFD process -- Correction step
   run_DLMFD_UzawaSolver(t_it);
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::do_before_time_stepping(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: do_before_time_stepping");

   Paraview_saveMultipliers_pvd << "<?xml version=\"1.0\"?>" << endl;
   Paraview_saveMultipliers_pvd << "<VTKFile type=\"Collection\" version=\"0.1\""
                                << " byte_order=\"LittleEndian\">" << endl;
   Paraview_saveMultipliers_pvd << "<Collection>" << endl;
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::do_additional_savings(int const &cycleNumber,
                                                   FV_TimeIterator const *t_it)
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
   f << str;

   f << "</Collection>" << endl;
   f << "</VTKFile>" << endl;
   f.close();

   write_PVTU_multiplier_file(filename);

   allrigidbodies->output_DLMFDPoints_PARAVIEW(SolidSolverResultsDirectory + "/" + filename,
                                               &Paraview_translated_distance_vector,
                                               true, pelCOMM->rank());

   solidSolver->saveResults("", t_it->time(), cycleNumber);
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
void DLMFD_FictitiousDomain::set_critical_distance(double critical_distance_)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: set_critical_distance");

   critical_distance = critical_distance_;
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::update_rigid_bodies(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: update_rigid_bodies");

   sub_prob_number = 3;
   MAC::out() << "-----------------------------------------" << "-------------" << endl;
   MAC::out() << "Sub-problem " << sub_prob_number
              << " : Rigid Bodies updating -- Prediction" << endl;
   MAC::out() << "-----------------------------------------" << "-------------" << endl;

   // Update the Rigid Bodies (Prediction problem)
   solidSolver->Simulation(t_it->time_step());
   solidSolver->getSolidBodyFeatures(solidFluid_transferStream);

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
              << " : DLMFD solving -- Correction" << endl;
   MAC::out() << "-----------------------------------------" << "-------------" << endl;

   // Initialize the DLMFD correction problem
   DLMFD_construction(t_it);

   // Solve the DLMFD correction problem
   DLMFD_solving(t_it);

   MAC::out() << "Uzawa problem completed" << endl;
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::DLMFD_construction(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: DLMFD_construction");

   double grid_size = UU->primary_grid()->get_smallest_constant_grid_size();
   set_critical_distance(sqrt(double(dim)) * grid_size);

   allrigidbodies->update(critical_distance, *solidFluid_transferStream);
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::DLMFD_solving(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: DLMFD_solving");

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

   
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::calculate_ReductionFor_qtranAndqrot(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: calculate_ReductionFor_qtranAndqrot");

   size_t i, k, indexcomp, nsolid;
   geomVector tempqtran(dim), tempqrot(dim);

   if (rank != master)
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
      pelCOMM->send(master, sendbuf);
   }
   else
   {
      vector<size_t_vector> const *pallshared = allrigidbodies->get_v_AllSharedOnProcs();
      for (k = 1; k < size_proc; ++k)
      {
         nsolid = (*pallshared)[k].size();
         doubleVector recvbuf(2 * nsolid * dim);
         pelCOMM->receive(k, recvbuf);
         for (i = 0; i < nsolid; ++i)
         {
            for (indexcomp = 0; indexcomp < dim; ++indexcomp)
            {
               tempqtran(indexcomp) = recvbuf(2 * i * dim + indexcomp);
               tempqrot(indexcomp) = recvbuf((2 * i + 1) * dim + indexcomp);
            }
            allrigidbodies->add_to_Qu((*pallshared)[k](i), tempqtran);
            allrigidbodies->add_to_Qrot((*pallshared)[k](i), tempqrot);
         }
      }
   }
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::velocity_broadcast_andUpdate_First(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: velocity_broadcast_andUpdateInOneIt");

   // Done by all processes
   allrigidbodies->completeFirstUzawaIteration_Velocity(rho_f, t_it->time_step(), split_gravity_vector, gravity_vector);

   // Broadcast tU and trot of shared solid components from master to
   // other processes
   Broadcast_tVectors_sharedParticles_MasterToAll(t_it);
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::velocity_broadcast_andUpdateInOneIt(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: velocity_broadcast_andUpdateInOneIt");

   // Done by all processes
   allrigidbodies->solve_Particles_OneUzawaIter_Velocity(rho_f, t_it->time_step());

   // Broadcast tU and trot of shared solid components from master to
   // other processes
   Broadcast_tVectors_sharedParticles_MasterToAll(t_it);
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::Broadcast_tVectors_sharedParticles_MasterToAll(FV_TimeIterator const *t_it)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: Broadcast_tVectors_sharedParticles_MasterToAll");

   geomVector ttran(dim);
   geomVector trot(dim);
   size_t i, k, nsolid, indexcomp;

   if (rank == master)
   {
      vector<size_t_vector> const *pallshared = allrigidbodies->get_v_AllSharedOnProcs();
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
         pelCOMM->send(k, sendbuf);
      }
   }
   else
   {
      list<int> const *pshared = allrigidbodies->get_SharedOnProc();
      list<int>::const_iterator il;

      nsolid = pshared->size();

      doubleVector recvbuf(2 * nsolid * dim);
      pelCOMM->receive(master, recvbuf);
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
