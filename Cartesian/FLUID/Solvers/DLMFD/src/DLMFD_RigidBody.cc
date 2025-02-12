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

//---------------------------------------------------------------------------
void DLMFD_RigidBody::add_boundary_point(const geomVector &point, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: add_boundary_point");

    DLMFD_BoundaryMultiplierPoint *bmp = new DLMFD_BoundaryMultiplierPoint(point, gravity_center);
    boundarypoints.push_back(bmp);
}