#include <DLMFD_RigidBody_BuilderFactory.hh>
#include <MAC.hh>
#include <DLMFD_RigidBody.hh>
#include <FS_RigidBody.hh>
#include <DLMFD_Sphere.hh>
using namespace std;

//---------------------------------------------------------------------------
DLMFD_RigidBody *DLMFD_RigidBody_BuilderFactory::create(FS_RigidBody *ptr_geom_rb,
                                                        FV_DiscreteField *pField_,
                                                        double const critical_distance_)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody_BuilderFactory:: create");

    DLMFD_RigidBody *ptr_dlmfd_rb = NULL;

    // Build the Fictitious Domain rigid body
    switch (ptr_geom_rb->get_shape_type())
    {
    case GEOM_SPHERE:
        ptr_dlmfd_rb = new DLMFD_Sphere(ptr_geom_rb, pField_, critical_distance_);
        break;

    default:
        MAC::out() << "Unknown geometric shape in "
                      "DLMFD_RigidBody_BuilderFactory::create"
                   << endl;
    }

    return (ptr_dlmfd_rb);
}

//---------------------------------------------------------------------------
DLMFD_RigidBody_BuilderFactory::DLMFD_RigidBody_BuilderFactory()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody_BuilderFactory:: DLMFD_RigidBody_BuilderFactory");
}

//---------------------------------------------------------------------------
DLMFD_RigidBody_BuilderFactory::~DLMFD_RigidBody_BuilderFactory()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody_BuilderFactory:: ~DLMFD_RigidBody_BuilderFactory");
}
