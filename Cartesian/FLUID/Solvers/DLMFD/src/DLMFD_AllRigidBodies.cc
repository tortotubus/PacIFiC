#include <DLMFD_AllRigidBodies.hh>
#include <DLMFD_RigidBody_BuilderFactory.hh>
#include <fstream>
using namespace std;

//---------------------------------------------------------------------------
DLMFD_AllRigidBodies::DLMFD_AllRigidBodies(size_t &dim,
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

   // Set the constrained field
   pField = UF;

   // Set the onProc IDs
   set_listIdOnProc();
}

//---------------------------------------------------------------------------
DLMFD_AllRigidBodies::~DLMFD_AllRigidBodies()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: ~DLMFD_AllRigidBodies");
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::set_all_points(double critical_distance) const
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

   for (size_t i = 0; i < RBs_number; i++)
   {
      onProc.push_back(i);
   }
}

//---------------------------------------------------------------------------
size_t DLMFD_AllRigidBodies::get_number_rigid_bodies() const
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: get_number_rigid_bodies");

   return (RBs_number);
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::update(istringstream &solidFluid_transferStream)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: update");

   ptr_FSallrigidbodies->update(solidFluid_transferStream);
   cout << "Hello World from Rigid body updating"
        << endl;
}

//---------------------------------------------------------------------------
void DLMFD_AllRigidBodies::update_RB_position_and_velocity(geomVector const &pos,
                                                           geomVector const &vel,
                                                           geomVector const &ang_vel,
                                                           vector<geomVector> const &periodic_directions,
                                                           double const &time_step)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: set_all_points");

   for (size_t i = 0; i < RBs_number; i++)
   {
      vec_ptr_DLMFDallrigidbodies[i]->update_RB_position_and_velocity(pos, vel, ang_vel, periodic_directions, time_step);
   }
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
