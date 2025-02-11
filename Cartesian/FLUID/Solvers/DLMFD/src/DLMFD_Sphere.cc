#include <DLMFD_Sphere.hh>
#include <DLMFD_RigidBody.hh>
#include <FS_Sphere.hh>
using namespace std;

//---------------------------------------------------------------------------
DLMFD_Sphere::DLMFD_Sphere() : DLMFD_RigidBody()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: DLMFD_Sphere");
}

//---------------------------------------------------------------------------
DLMFD_Sphere::DLMFD_Sphere(FS_RigidBody *pgrb) : DLMFD_RigidBody(pgrb)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: DLMFD_Sphere");
}

//---------------------------------------------------------------------------
DLMFD_Sphere::~DLMFD_Sphere()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: ~DLMFD_Sphere");
}

void DLMFD_Sphere::set_all_points()
{
    MAC_LABEL("DLMFD_Sphere:: set_all_points");

    cout << "Hello World from setting DLMFD points" << endl;
}

//---------------------------------------------------------------------------
geomVector DLMFD_Sphere::get_rigid_body_velocity(geomVector const &point) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: get_rigid_body_velocity");
    return (ptr_FSrigidbody->rigid_body_velocity(point));
}

//---------------------------------------------------------------------------
void DLMFD_Sphere::update()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: update");
}