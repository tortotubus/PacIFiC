#include <DLMFD_AllRigidBodies.hh>
#include <DLMFD_RigidBody_BuilderFactory.hh>
#include <doubleArray2D.hh>
#include <fstream>
using namespace std;

//---------------------------------------------------------------------------
DLMFD_AllRigidBodies::DLMFD_AllRigidBodies(size_t &dim,
                                           MAC_Communicator const *pelCOMM_,
                                           istringstream &in,
                                           bool const &are_particles_fixed,
                                           FV_DiscreteField const *UF,
                                           FV_DiscreteField const *PF)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: DLMFD_AllRigidBodies");

   // Set the rigid bodies
   ptr_FSallrigidbodies = new FS_AllRigidBodies(dim, in, are_particles_fixed);
   RBs_number = ptr_FSallrigidbodies->get_number_rigid_bodies();

   DLMFD_RigidBody *ptr_dlmfd_rb = NULL;

   for (size_t i = 0; i < RBs_number; ++i)
   {
      vec_ptr_DLMFDallrigidbodies.push_back(ptr_dlmfd_rb);
      vec_ptr_DLMFDallrigidbodies[i] = DLMFD_RigidBody_BuilderFactory::create(
          ptr_FSallrigidbodies->get_ptr_rigid_body(i));
   }

   // Set MPI
   set_MPI_data(pelCOMM_);

   // Set the constrained field
   set_ptr_constrained_field(UF);

   // Set the onProc IDs
   set_listIdOnProc();

   // Set the points infos for all rigid bodies
   set_points_infos(pField);

   // Set mass, density and inertia tensor of each rigid body
   set_mass_and_density_and_inertia();
}

//---------------------------------------------------------------------------
DLMFD_AllRigidBodies::~DLMFD_AllRigidBodies()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: ~DLMFD_AllRigidBodies");
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_all_points(double critical_distance)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_all_points");

   for (size_t i = 0; i < RBs_number; i++)
   {
      vec_ptr_DLMFDallrigidbodies[i]->set_all_points(pField, critical_distance);
   }
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_listIdOnProc()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_all_points");

   entireOnProc.clear();
   SharedOnProc.clear();
   onProc.clear();

   for (size_t i = 0; i < RBs_number; i++)
   {
      onProc.push_back(i);
      entireOnProc.push_back(i);
   }

   size_t j = 0, k;
   list<int>::iterator il;

   size_t_vector v_SharedOnProc(SharedOnProc.size());
   for (il = SharedOnProc.begin(); il != SharedOnProc.end(); il++, ++j)
      v_SharedOnProc(j) = *il;

   if (rank != master)
      pelCOMM->send(master, v_SharedOnProc);
   else
   {
      v_AllSharedOnProcs = new vector<size_t_vector>(size_procs, v_SharedOnProc);
      // for (size_t k = 1; k < size_procs; ++k)
      //    pelCOMM->receive(k, (*v_AllSharedOnProcs)[k]);
      l_AllSharedOnProcs = new list<int>;
      for (k = 0; k < size_procs; ++k)
         for (j = 0; j < (*v_AllSharedOnProcs)[k].size(); ++j)
            l_AllSharedOnProcs->push_back((*v_AllSharedOnProcs)[k](j));
      l_AllSharedOnProcs->sort();
      l_AllSharedOnProcs->unique();
   }
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_ptr_constrained_field(FV_DiscreteField const *pField_)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_ptr_constrained_field");

   pField = pField_;
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_points_infos(FV_DiscreteField const *pField_)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_ptr_constrained_field");

   for (size_t i = 0; i < RBs_number; i++)
   {
      vec_ptr_DLMFDallrigidbodies[i]->set_points_infos(pField);
   }
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_MPI_data(MAC_Communicator const *pelCOMM_)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_MPI");

   pelCOMM = pelCOMM_;
   size_procs = pelCOMM->nb_ranks();
   rank = pelCOMM->rank();
   master = 0;
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_coupling_factor(double const &rho_f, bool const &explicit_treatment)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_coupling_factor");

   list<int>::iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->set_coupling_factor(rho_f, explicit_treatment);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_mass_and_density_and_inertia()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_mass_and_density_and_inertia");

   list<int>::iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->set_mass_and_density_and_inertia();
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_Tu(const size_t i, const geomVector &ttran)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_Tu");

   vec_ptr_DLMFDallrigidbodies[i]->set_Tu(ttran);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_Trot(const size_t i, const geomVector &trot)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_Trot");

   vec_ptr_DLMFDallrigidbodies[i]->set_Trot(trot);
}

//---------------------------------------------------------------------------
size_t DLMFD_AllRigidBodies::get_number_rigid_bodies() const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: get_number_rigid_bodies");

   return (RBs_number);
}

//---------------------------------------------------------------------------
list<int> const *DLMFD_AllRigidBodies::get_SharedOnProc() const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: get_SharedOnProc");

   return &SharedOnProc;
}

//---------------------------------------------------------------------------
vector<size_t_vector> const *DLMFD_AllRigidBodies::get_v_AllSharedOnProcs() const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: get_v_AllSharedOnProcs");

   return v_AllSharedOnProcs;
}

//---------------------------------------------------------------------------
geomVector const DLMFD_AllRigidBodies::get_Qu(const size_t i) const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::get_Qu");

   return vec_ptr_DLMFDallrigidbodies[i]->get_Qu();
}

//---------------------------------------------------------------------------
geomVector const DLMFD_AllRigidBodies::get_Qrot(const size_t i) const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::get_Qrot");

   return vec_ptr_DLMFDallrigidbodies[i]->get_Qrot();
}

//---------------------------------------------------------------------------
geomVector const DLMFD_AllRigidBodies::get_Tu(const size_t i) const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::get_Tu");

   return vec_ptr_DLMFDallrigidbodies[i]->get_Tu();
}

//---------------------------------------------------------------------------
geomVector const DLMFD_AllRigidBodies::get_Trot(const size_t i) const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::get_Trot");

   return vec_ptr_DLMFDallrigidbodies[i]->get_Trot();
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::update(double critical_distance, istringstream &solidFluid_transferStream)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: update");

   ptr_FSallrigidbodies->update(solidFluid_transferStream);

   for (size_t i = 0; i < RBs_number; i++)
   {
      vec_ptr_DLMFDallrigidbodies[i]->update();
   }

   set_all_points(critical_distance);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::output_DLMFDPoints_PARAVIEW(const string &filename,
                                                       geomVector const *translated_distance_vector,
                                                       const bool &withIntPts,
                                                       size_t rank) const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: output_DLMFDPoints_PARAVIEW");

   ostringstream ossRK;
   ossRK << rank;

   size_t output_npts = 0;
   list<int>::const_iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      output_npts += vec_ptr_DLMFDallrigidbodies[*il]->get_npts_output(withIntPts);

   ofstream f((filename + "_" + ossRK.str() + ".vtu").c_str(), ios::out);

   f << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" "
     << "byte_order=\"LittleEndian\">" << endl;
   f << "<UnstructuredGrid>" << endl;
   f << "<Piece NumberOfPoints=\"" << output_npts << "\""
     << " NumberOfCells=\"" << output_npts << "\">" << endl;
   f << "<Points>" << endl;
   f << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
     << "Format=\"ascii\">" << endl;

   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->print_partPointsCoordinates(f, translated_distance_vector, "", "", withIntPts);

   f << "</DataArray>" << endl;
   f << "</Points>" << endl;
   f << "<Cells>" << endl;
   f << "<DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">"
     << endl;
   for (size_t i = 0; i < output_npts; ++i)
      f << i << " ";
   f << endl
     << "</DataArray>" << endl;
   f << "<DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">" << endl;
   for (size_t i = 1; i < output_npts + 1; ++i)
      f << i << " ";
   f << endl
     << "</DataArray>" << endl;
   f << "<DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">" << endl;
   for (size_t i = 0; i < output_npts; ++i)
      f << "1 ";
   f << endl
     << "</DataArray>" << endl;
   f << "</Cells>" << endl;
   f << "</Piece>" << endl;
   f << "</UnstructuredGrid>" << endl;
   f << "</VTKFile>" << endl;
   f.close();
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::compute_all_Qu(bool init)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: compute_all_Qu");

   list<int>::iterator il;

   // Compute Qu for particles on this process
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->compute_Qu(init);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::compute_all_Qrot(bool init)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: compute_all_Qrot");

   list<int>::iterator il;

   // Compute Qrot for particles on this process
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->compute_Qrot(init);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::add_to_Qu(const size_t i, const geomVector &qtran)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: add_to_Qu");

   vec_ptr_DLMFDallrigidbodies[i]->add_to_Qu(qtran);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::add_to_Qrot(const size_t i, const geomVector &qrot)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: add_to_Qrot");

   vec_ptr_DLMFDallrigidbodies[i]->add_to_Qrot(qrot);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::completeFirstUzawaIteration_Velocity(const double &rho_f,
                                                                const double &timestep,
                                                                const geomVector &gravity_vector_split,
                                                                const geomVector &gravity_vector)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::completeFirstUzawaIteration_Velocity");

   list<int>::iterator il;
   // Particles entire on this process
   for (il = entireOnProc.begin(); il != entireOnProc.end(); il++)
   {
      if (vec_ptr_DLMFDallrigidbodies[*il]->get_component_type() != "O")
      {
         vec_ptr_DLMFDallrigidbodies[*il]->correctQvectorsAndInitUzawa_Velocity(rho_f, timestep, gravity_vector_split, gravity_vector);
         vec_ptr_DLMFDallrigidbodies[*il]->calculateParticleVelocities(rho_f, timestep);
         vec_ptr_DLMFDallrigidbodies[*il]->updateSolutionVectors_Velocity(1., 0.);
      }
      else
         vec_ptr_DLMFDallrigidbodies[*il]->setTVectors_constant();
   }

   // Particles shared by processes treated on master process
   if (l_AllSharedOnProcs)
   {
      for (il = l_AllSharedOnProcs->begin(); il != l_AllSharedOnProcs->end(); il++)
      {
         if (vec_ptr_DLMFDallrigidbodies[*il]->get_component_type() != "O")
         {
            vec_ptr_DLMFDallrigidbodies[*il]->correctQvectorsAndInitUzawa_Velocity(rho_f, timestep, gravity_vector_split, gravity_vector);
            vec_ptr_DLMFDallrigidbodies[*il]->calculateParticleVelocities(rho_f, timestep);
            vec_ptr_DLMFDallrigidbodies[*il]->updateSolutionVectors_Velocity(1., 0.);
         }
         else
            vec_ptr_DLMFDallrigidbodies[*il]->setTVectors_constant();
      }
   }
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::update_ParticlesVelocities_afterBcast_of_T()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::update_ParticlesVelocities_afterBcast_of_T");

   list<int>::iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->updateSolutionVectors_Velocity(1., 0.);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::solve_Particles_OneUzawaIter_Velocity(double const &rho_f, double const &timestep)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: solve_Particles_OneUzawaIter_Velocity");

   list<int>::iterator il;
   // Particles entire on this process
   for (il = entireOnProc.begin(); il != entireOnProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->calculateParticleVelocities(rho_f, timestep);

   // // Particles shared by processes treated on master process
   // if (l_AllSharedOnProcs)
   //    for (il = l_AllSharedOnProcs->begin(); il != l_AllSharedOnProcs->end(); il++)
   //       solidcomponents[*il]->calculateParticleVelocities(densityFl, timestep);
}
