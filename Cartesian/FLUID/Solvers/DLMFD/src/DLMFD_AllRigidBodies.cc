#include <DLMFD_AllRigidBodies.hh>
#include <DLMFD_RigidBody_BuilderFactory.hh>
#include <doubleArray2D.hh>
#include <fstream>
#include <math.h>
using namespace std;

//---------------------------------------------------------------------------
DLMFD_AllRigidBodies::DLMFD_AllRigidBodies(size_t &dim,
                                           MAC_Communicator const *pelCOMM_,
                                           istringstream &in,
                                           bool const &are_particles_fixed_,
                                           FV_DiscreteField *UU,
                                           FV_DiscreteField *PP) : are_particles_fixed(are_particles_fixed_)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: DLMFD_AllRigidBodies");

   // Set the rigid bodies
   ptr_FSallrigidbodies = new FS_AllRigidBodies(dim, in, are_particles_fixed_);
   RBs_number = ptr_FSallrigidbodies->get_number_rigid_bodies();

   DLMFD_RigidBody *ptr_dlmfd_rb = NULL;

   for (size_t i = 0; i < RBs_number; ++i)
   {
      vec_ptr_DLMFDallrigidbodies.push_back(ptr_dlmfd_rb);
      vec_ptr_DLMFDallrigidbodies[i] = DLMFD_RigidBody_BuilderFactory::create(
          ptr_FSallrigidbodies->get_ptr_rigid_body(i));
   }

   // Set the number of moving RBs
   set_npart();

   // Set MPI
   set_MPI_data(pelCOMM_);

   // Set the constrained field
   set_ptr_constrained_field(UU);

   // Set the onProc IDs
   set_listIdOnProc();

   // Set the points infos for all rigid bodies
   set_points_infos();

   // Set mass, density and inertia tensor of each rigid body
   ptr_FSallrigidbodies->update(in);
   set_mass_and_density_and_volume_and_inertia();

   // Set the force and torque boolean
   b_output_hydro_forceTorque = false;
}

//---------------------------------------------------------------------------
DLMFD_AllRigidBodies::~DLMFD_AllRigidBodies()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: ~DLMFD_AllRigidBodies");
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_b_output_hydro_forceTorque(bool const &is_output)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::set_b_output_hydro_forceTorque");

   b_output_hydro_forceTorque = is_output;
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
void DLMFD_AllRigidBodies::eraseCriticalDLMFDPoints(const double &time, double critical_distance)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_all_points");

   // !!! IMPORTANT REMARK !!!
   // The order in which operations on DLM/FD points are performed
   // is crucial, so do not play with it if it is not purposely

   // TODO: erase critical boundary points

   double coef_cd = 1. / sqrt(double(dim));

   for (size_t i = 0; i < RBs_number; i++)
      vec_ptr_DLMFDallrigidbodies[i]->erase_critical_interior_points_PerProc(coef_cd * critical_distance);
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
void DLMFD_AllRigidBodies::set_ptr_constrained_field(FV_DiscreteField *pField_)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_ptr_constrained_field");

   pField = pField_;
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_points_infos()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_ptr_constrained_field");

   for (size_t i = 0; i < RBs_number; i++)
   {
      vec_ptr_DLMFDallrigidbodies[i]->set_points_infos(pField);
   }
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::check_allocation_DLMFD_Cvectors()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: check_allocation_DLMFD_Cvectors");

   for (size_t i = 0; i < RBs_number; i++)
   {
      vec_ptr_DLMFDallrigidbodies[i]->check_allocation_DLMFD_Cvectors();
   }
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::fill_DLMFD_Cvectors()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: fill_DLMFD_Cvectors");

   for (size_t i = 0; i < RBs_number; i++)
   {
      vec_ptr_DLMFDallrigidbodies[i]->fill_DLMFD_Cvectors();
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
void DLMFD_AllRigidBodies::set_mass_and_density_and_volume_and_inertia()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_mass_and_density_and_volume_and_inertia");

   list<int>::iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
   {
      vec_ptr_DLMFDallrigidbodies[*il]->set_mass_and_density_and_inertia();
      vec_ptr_DLMFDallrigidbodies[*il]->set_volume();
   }
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
void const DLMFD_AllRigidBodies::set_npart()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::set_npart");

   particles_number = 0;

   for (size_t i = 0; i < RBs_number; i++)
      if (vec_ptr_DLMFDallrigidbodies[i]->get_component_type() != "O")
         particles_number++;
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_translational_velocity(const size_t i, const geomVector &vtran)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::set_translational_velocity");

   vec_ptr_DLMFDallrigidbodies[i]->set_translational_velocity(vtran);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_angular_velocity_3D(const size_t i, const geomVector &vrot)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::set_angular_velocity_3D");

   vec_ptr_DLMFDallrigidbodies[i]->set_angular_velocity_3D(vrot);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_output_frequency(size_t output_frequency_)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::set_output_frequency");

   output_frequency = output_frequency_;
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
size_t const DLMFD_AllRigidBodies::get_npart()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::get_npart");

   return particles_number;
}

//---------------------------------------------------------------------------
list<int> const *DLMFD_AllRigidBodies::get_onProc() const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::get_onProc");

   return &onProc;
}

//---------------------------------------------------------------------------
geomVector DLMFD_AllRigidBodies::get_translational_velocity(const size_t i) const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::get_translational_velocity");

   return vec_ptr_DLMFDallrigidbodies[i]->get_translational_velocity();
}

//---------------------------------------------------------------------------
geomVector DLMFD_AllRigidBodies::get_angular_velocity_3D(const size_t i) const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::get_angular_velocity");

   return vec_ptr_DLMFDallrigidbodies[i]->get_angular_velocity_3D();
}

//---------------------------------------------------------------------------
double const DLMFD_AllRigidBodies::get_volume(const size_t i) const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::get_volume");

   return vec_ptr_DLMFDallrigidbodies[i]->get_volume();
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::update(double const &time, double critical_distance, istringstream &solidFluid_transferStream)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: update");

   ptr_FSallrigidbodies->update(solidFluid_transferStream);

   for (size_t i = 0; i < RBs_number; i++)
   {
      vec_ptr_DLMFDallrigidbodies[i]->update();
   }

   // Set valid points
   set_all_points(critical_distance);
   eraseCriticalDLMFDPoints(time, critical_distance);

   // Set points infos
   set_points_infos();

   // Allocate Uzawa vectors
   allocate_initialize_Uzawa_vectors();

   // Fill DLMFD vectors
   check_allocation_DLMFD_Cvectors();
   fill_DLMFD_Cvectors();
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
void DLMFD_AllRigidBodies::particles_hydrodynamic_force_output(const string &path_name,
                                                               const bool &b_restart,
                                                               const double &time,
                                                               const double &timestep,
                                                               const double &rho_f,
                                                               vector<vector<double>> const *Iw_Idw)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::particles_hydrodynamic_force_output");

   static size_t counter = 0;
   static size_t restart_counter = 0;
   ostringstream oss;
   oss << counter;
   restart_counter += !b_restart;

   if (b_output_hydro_forceTorque)
   {
      geomVector vtrans;
      geomVector vrot;
      size_t i, indexcomp;

      // Compute hydro force & torque, write to files and update n-1 velocity
      if (counter == 0 && restart_counter)
      {
         list<int>::const_iterator il;
         double dlmfd_mass = 0., density_ratio = 0., mass = 0.;
         string file_name;
         geomVector omega(3), tempqtran(dim), tempqrot(dim);

         // Sum contributions of each processor to the different particles

         // Compute hydrodynamic force & torque and write to a file
         // Fh = <lambda,V>_P + (rho_f/rho_s)MdU/dt
         // Mh = <lambda,xi^GM>_P + (rho_f/rho_s).( Idom/dt + om x (I.om) )
         // If "extended", compute the rigid-body rotation contributions
         // a) om x (MU) (denoted ocmu below)
         // b) om x (I.om) (denoted ocjo below)
         // and write them in separate files
         // In 2D:
         // a) om x (MU) : x and y components
         // b) om x (I.om) is zero by construction
         // and we write xocmu and yocmu
         // In 3D:
         // a) om x (MU) : x, y, z components
         // b) om x (I.om) is zero by construction for spheres and non-zero
         // otherwise
         // and we write all 6 files
         if (rank == master)
         {
            geomVector hydro_force(dim);
            geomVector hydro_torque(dim);

            // Write hydro force & torque to files
            // -----------------------------------

            file_name = path_name + "x-hydroForce_time.res";
            // Manually remove the content
            ofstream xhf(file_name.c_str(), ios::app);
            if (time < 0.001 * timestep)
            {
               ofstream xhf(file_name.c_str(), ios::out | ios::trunc);
               xhf.close();
               xhf.open(file_name.c_str(), ios::app);
            }
            xhf.setf(ios::scientific, ios::floatfield);
            xhf.precision(6);
            xhf << time;

            file_name = path_name + "y-hydroForce_time.res";
            ofstream yhf(file_name.c_str(), ios::app);
            yhf.setf(ios::scientific, ios::floatfield);
            yhf.precision(6);
            yhf << time;

            file_name = path_name + "z-hydroForce_time.res";
            ofstream zhf(file_name.c_str(), ios::app);
            zhf.setf(ios::scientific, ios::floatfield);
            zhf.precision(6);
            zhf << time;

            file_name = path_name + "x-hydroTorque_time.res";
            ofstream xht(file_name.c_str(), ios::app);
            xht.setf(ios::scientific, ios::floatfield);
            xht.precision(6);
            xht << time;

            file_name = path_name + "y-hydroTorque_time.res";
            ofstream yht(file_name.c_str(), ios::app);
            yht.setf(ios::scientific, ios::floatfield);
            yht.precision(6);
            yht << time;

            file_name = path_name + "z-hydroTorque_time.res";
            ofstream zht(file_name.c_str(), ios::app);
            zht.setf(ios::scientific, ios::floatfield);
            zht.precision(6);
            zht << time;

            // Write hydro force and update velocity array at previous time
            for (i = 0; i < RBs_number; i++)
            {
               // Get the Lagrange multiplier contribution
               tempqtran = get_Qu(i);
               tempqrot = get_Qrot(i);

               // Hydrodynamic force based on translational velocity
               vtrans = get_translational_velocity(i);
               dlmfd_mass = rho_f * get_volume(i);

               for (indexcomp = 0; indexcomp < dim; ++indexcomp)
               {
                  hydro_force(indexcomp) = -tempqtran(indexcomp);

                  if (vec_ptr_DLMFDallrigidbodies[i]->get_component_type() != "O" && !are_particles_fixed)
                  {
                     // hydro_force(indexcomp) += dlmfd_mass * (vtrans(indexcomp) - (*translational_velocity_nm1)(i, indexcomp)) / timestep;
                     (*translational_velocity_nm1)(i, indexcomp) = vtrans(indexcomp);
                  }
               }
               xhf << " " << hydro_force(0);
               yhf << " " << hydro_force(1);
               zhf << " " << hydro_force(2);

               // Hydrodynamic torque based on angular velocity
               vrot = get_angular_velocity_3D(i);
               density_ratio = rho_f / vec_ptr_DLMFDallrigidbodies[i]->get_density();

               for (indexcomp = 0; indexcomp < dim; ++indexcomp)
               {
                  hydro_torque(indexcomp) = -tempqrot(indexcomp);

                  if (vec_ptr_DLMFDallrigidbodies[i]->get_component_type() != "O" && !are_particles_fixed)
                  {
                     hydro_torque(indexcomp) += density_ratio * (*Iw_Idw)[i][dim + indexcomp] / timestep;
                     (*angular_velocity_3D_nm1)(i, indexcomp) = vrot(indexcomp);
                  }
               }
               // Add w x Iw to all torque component
               if (vec_ptr_DLMFDallrigidbodies[i]->get_component_type() != "O" && !are_particles_fixed)
               {
                  hydro_torque(0) += density_ratio * (vrot(1) * (*Iw_Idw)[i][2] - vrot(2) * (*Iw_Idw)[i][1]);
                  hydro_torque(1) += density_ratio * (vrot(2) * (*Iw_Idw)[i][0] - vrot(0) * (*Iw_Idw)[i][2]);
                  hydro_torque(2) += density_ratio * (vrot(0) * (*Iw_Idw)[i][1] - vrot(1) * (*Iw_Idw)[i][0]);
               }

               xht << " " << hydro_torque(0);
               yht << " " << hydro_torque(1);
               zht << " " << hydro_torque(2);
            }
            // Close files
            // Hydro force
            xhf << endl;
            xhf.close();
            yhf << endl;
            yhf.close();
            zhf << endl;
            zhf.close();

            // Hydro torque
            xht << endl;
            xht.close();
            yht << endl;
            yht.close();
            zht << endl;
            zht.close();
         }
      }
      // Only update n-1 velocity
      else
      {
         if (rank == master)
         {
            for (i = 0; i < particles_number; i++)
            {
               vtrans = vec_ptr_DLMFDallrigidbodies[i]->get_translational_velocity();
               vrot = vec_ptr_DLMFDallrigidbodies[i]->get_angular_velocity_3D();
               for (indexcomp = 0; indexcomp < dim; ++indexcomp)
               {
                  (*translational_velocity_nm1)(i, indexcomp) = vtrans(indexcomp);
                  (*angular_velocity_3D_nm1)(i, indexcomp) = vrot(indexcomp);
               }
            }
         }
      }
   }
   ++counter;
   if (counter == output_frequency)
      counter = 0;
   ++restart_counter;
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::sum_DLM_hydrodynamic_force_output(const bool &b_restart)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::sum_DLM_hydrodynamic_force_output");

   static size_t counter = 0;
   static size_t restart_counter = 0;
   restart_counter += !b_restart;
   size_t no = RBs_number;

   if (b_output_hydro_forceTorque)
   {
      size_t i, indexcomp;

      if (counter == 0 && restart_counter)
      {
         list<int>::const_iterator il;
         size_t k, ns, id;
         double intTodouble = 0.1;

         // Sum contributions of each processor to the different particles
         // same as MY_FictitiousDomain::calculate_ReductionFor_qtranAndqrot
         geomVector tempqtran(dim), tempqrot(dim);

         if (rank != master)
         {
            // Compute on each particle on proc
            // q_tran = -<lambda,V>_P
            // q_rot = -<lambda,xi^GM>_P
            compute_all_Qu(true);
            compute_all_Qrot(true);

            ns = onProc.size();
            doubleVector sendbuf(ns * (2 * dim + 1));
            i = 0;
            for (il = onProc.begin(); il != onProc.end(); il++)
            {
               tempqtran = vec_ptr_DLMFDallrigidbodies[*il]->get_Qu();
               tempqrot = vec_ptr_DLMFDallrigidbodies[*il]->get_Qrot();

               sendbuf(i) = double(*il) + intTodouble;
               ++i;
               for (indexcomp = 0; indexcomp < dim; ++indexcomp)
               {
                  sendbuf(i + indexcomp) = tempqtran(indexcomp);
                  sendbuf(i + dim + indexcomp) = tempqrot(indexcomp);
               }
               i += 2 * dim;
            }
            pelCOMM->send(master, sendbuf);
         }
         else
         {
            // Compute q_tran = -<lambda,V>_P and q_rot = -<lambda,xi^GM>_P
            // for all particles ("all" is crucial !!)
            // Do not use compute_all_Qu and compute_all_Qrot that act on particles
            // that are on proc only
            // If a particle is not on master proc, methods
            // SolidComponent::compute_Qu and SolidComponent::compute_Qrot
            // initialize q_tran and q_rot to zero
            for (i = 0; i < no; i++)
            {
               vec_ptr_DLMFDallrigidbodies[i]->compute_Qu(true);
               vec_ptr_DLMFDallrigidbodies[i]->compute_Qrot(true);
            }

            for (k = 1; k < size_procs; ++k)
            {
               doubleVector recvbuf(0); // No need to set the reception buffer size
                                        // Automatically set by the receive method
                                        // see MY_Communicator:: receive
               pelCOMM->receive(k, recvbuf);
               size_t recvsize = recvbuf.size();
               for (i = 0; i < recvsize;)
               {
                  id = size_t(recvbuf(i));
                  ++i;
                  for (indexcomp = 0; indexcomp < dim; ++indexcomp)
                  {
                     tempqtran(indexcomp) = recvbuf(i + indexcomp);
                     tempqrot(indexcomp) = recvbuf(i + dim + indexcomp);
                  }
                  vec_ptr_DLMFDallrigidbodies[id]->add_to_Qu(tempqtran);
                  vec_ptr_DLMFDallrigidbodies[id]->add_to_Qrot(tempqrot);

                  i += 2 * dim;
               }
            }
         }
      }
   }

   ++counter;
   if (counter == output_frequency)
      counter = 0;
   ++restart_counter;
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::nullify_all_Uzawa_vectors()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::nullify_all_Uzawa_vectors");

   list<int>::iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->nullify_Uzawa_vectors();
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::allocate_initialize_Uzawa_vectors()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::allocate_initialize_Uzawa_vectors");

   list<int>::iterator il;

   // Compute Qu for particles on this process
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->allocate_initialize_Uzawa_vectors();
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
   {
      vec_ptr_DLMFDallrigidbodies[*il]->updateSolutionVectors_Velocity(1., 0.);
   }
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

   // Particles shared by processes treated on master process
   if (l_AllSharedOnProcs)
      for (il = l_AllSharedOnProcs->begin(); il != l_AllSharedOnProcs->end(); il++)
         vec_ptr_DLMFDallrigidbodies[*il]->calculateParticleVelocities(rho_f, timestep);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::compute_fluid_rhs(DLMFD_ProjectionNavierStokesSystem *GLOBAL_EQ, bool init)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: compute_fluid_rhs");

   list<int>::const_iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->compute_fluid_rhs(GLOBAL_EQ, init);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::compute_x_residuals_Velocity()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::compute_x_residuals_Velocity");

   list<int>::const_iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->compute_x_residuals_Velocity(pField);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::compute_r_and_w_FirstUzawaIteration()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::compute_x_residuals_Velocity");

   list<int>::const_iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->compute_r_and_w_FirstUzawaIteration();
}

//---------------------------------------------------------------------------
double DLMFD_AllRigidBodies::compute_r_dot_r() const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::compute_r_dot_r");

   double rdotr = 0.;

   list<int>::const_iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      rdotr += vec_ptr_DLMFDallrigidbodies[*il]->compute_r_dot_r();

   double rdotr_collective = pelCOMM->sum(rdotr);

   return rdotr_collective;
}

//---------------------------------------------------------------------------
double DLMFD_AllRigidBodies::compute_w_dot_x() const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::compute_w_dot_x");

   double wdotx = 0.;

   list<int>::const_iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      wdotx += vec_ptr_DLMFDallrigidbodies[*il]->compute_w_dot_x();

   double wdotx_collective = pelCOMM->sum(wdotx);

   return wdotx_collective;
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::update_lambda_and_r(const double &alpha)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::update_lambda_and_r");

   list<int>::iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->update_lambda_and_r(alpha);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::update_ParticlesVelocities(const double &alpha)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::update_ParticlesVelocities");

   list<int>::iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
   {
      if (vec_ptr_DLMFDallrigidbodies[*il]->get_component_type() == "P" || vec_ptr_DLMFDallrigidbodies[*il]->get_component_type() == "PP")
         vec_ptr_DLMFDallrigidbodies[*il]->updateSolutionVectors_Velocity(alpha, 1.);
      else
         vec_ptr_DLMFDallrigidbodies[*il]->updateSolutionVectors_Velocity(0., 1.);
   }
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::update_w(const double &beta)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::update_lambda_and_r");

   list<int>::iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->update_w(beta);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::compute_fluid_DLMFD_explicit(DLMFD_ProjectionNavierStokesSystem *GLOBAL_EQ, bool bulk)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::compute_fluid_DLMFD_explicit");

   list<int>::iterator il;
   for (il = onProc.begin(); il != onProc.end(); il++)
      vec_ptr_DLMFDallrigidbodies[*il]->compute_fluid_DLMFD_explicit(GLOBAL_EQ, bulk, pField);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::allocate_translational_angular_velocity_array()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies::allocate_translational_angular_velocity_array");

   if (rank == master && b_output_hydro_forceTorque)
   {
      translational_velocity_nm1 = new doubleArray2D(particles_number, dim);
      angular_velocity_3D_nm1 = new doubleArray2D(particles_number, dim);
   }
}
