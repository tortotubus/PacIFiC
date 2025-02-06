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
DLMFD_RigidBody::~DLMFD_RigidBody()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: ~DLMFD_RigidBody");
}

//---------------------------------------------------------------------------
geomVector DLMFD_RigidBody::get_rigid_body_velocity(geomVector const &point)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: get_rigid_body_velocity");
    geomVector virtual_velocity = ptr_FSrigidbody->rigid_body_velocity(point);

    return (virtual_velocity);
}