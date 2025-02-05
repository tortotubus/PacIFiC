#include <DLMFD_AllRigidBodies.hh>
#include <FS_AllRigidBodies.hh>

//---------------------------------------------------------------------------
DLMFD_AllRigidBodies:: DLMFD_AllRigidBodies(size_t& dim, istream& in, bool const& are_particles_fixed)
//--------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_AllRigidBodies:: DLMFD_AllRigidBodies");

    ptr_FSallrigidbodies = new FS_AllRigidBodies(dim, in, are_particles_fixed);
    RBs_number = ptr_FSallrigidbodies->get_number_rigid_bodies();
}

//---------------------------------------------------------------------------
DLMFD_AllRigidBodies:: ~DLMFD_AllRigidBodies()
//--------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: ~DLMFD_AllRigidBodies");
}

//---------------------------------------------------------------------------
size_t DLMFD_AllRigidBodies:: get_number_rigid_bodies() const
//--------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_AllRigidBodies:: get_number_rigid_bodies");
   
   return(RBs_number);
}


