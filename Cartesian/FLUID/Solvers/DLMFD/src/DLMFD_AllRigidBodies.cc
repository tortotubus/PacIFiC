#include <DLMFD_AllRigidBodies.hh>
#include <DLMFD_RigidBody_BuilderFactory.hh>
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
   ptr_FSallrigidbodies = new FS_AllRigidBodies(dim, in, are_particles_fixed);
   RBs_number = ptr_FSallrigidbodies->get_number_rigid_bodies();

   DLMFD_RigidBody *ptr_dlmfd_rb = NULL;

   for (size_t i = 0; i < RBs_number; ++i)
   {
      vec_ptr_DLMFDallrigidbodies.push_back(ptr_dlmfd_rb);
      vec_ptr_DLMFDallrigidbodies[i] = DLMFD_RigidBody_BuilderFactory::create(
          ptr_FSallrigidbodies->get_ptr_rigid_body(i));
   }
}

//---------------------------------------------------------------------------
DLMFD_AllRigidBodies::~DLMFD_AllRigidBodies()
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: ~DLMFD_AllRigidBodies");
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
