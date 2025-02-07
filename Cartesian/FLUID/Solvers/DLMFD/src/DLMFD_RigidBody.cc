#include <DLMFD_RigidBody.hh>
#include <FS_RigidBody.hh>
#include <geomVector.hh>

//---------------------------------------------------------------------------
DLMFD_RigidBody::DLMFD_RigidBody()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: DLMFD_RigidBody");
}

//---------------------------------------------------------------------------
DLMFD_RigidBody::DLMFD_RigidBody(FS_RigidBody *pgrb) : ptr_FSrigidbody(pgrb)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: DLMFD_RigidBody");
}

//---------------------------------------------------------------------------
DLMFD_RigidBody::~DLMFD_RigidBody()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: ~DLMFD_RigidBody");
}
