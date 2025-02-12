#include <DLMFD_Sphere.hh>
#include <DLMFD_RigidBody.hh>
#include <FS_Sphere.hh>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
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

    gravity_center = *DLMFD_Sphere::get_ptr_to_gravity_centre();
    radius = DLMFD_Sphere::get_circumscribed_radius();
}

//---------------------------------------------------------------------------
DLMFD_Sphere::~DLMFD_Sphere()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: ~DLMFD_Sphere");
}

//---------------------------------------------------------------------------
void DLMFD_Sphere::set_all_points(double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: set_all_points");

    //-- Boundary points
    double pi = acos(-1.);
    double quart_perimeter = pi * radius / 2.,
           half_perimeter = pi * radius;
    unsigned nbpts = 2 * int(half_perimeter / critical_distance);
    unsigned nbptsZ = int(quart_perimeter / critical_distance), ipoint = 0;
    double angle = nbptsZ <= 0 ? pi / 2. : pi / (2. * nbptsZ),
           angleLocal = 2 * pi / nbpts,
           gravityCenterX = gravity_center(0),
           gravityCenterY = gravity_center(1),
           gravityCenterZ = gravity_center(2);
    geomVector point(3);
    cout << half_perimeter << endl;

    // Equator
    for (size_t i = 0; i <= nbpts - 1; ++i)
    {
        point(0) = radius * cos(i * angleLocal) + gravityCenterX;
        point(1) = radius * sin(i * angleLocal) + gravityCenterY;
        point(2) = gravityCenterZ;
        add_boundary_point(point, critical_distance);
    }

    // Upper tip
    point(0) = gravityCenterX;
    point(1) = gravityCenterY;
    point(2) = gravityCenterZ + radius;
    add_boundary_point(point, critical_distance);

    // Lower tip
    point(2) = gravityCenterZ - radius;
    add_boundary_point(point, critical_distance);

    // Surface
    double totAngle = angle;
    while (ipoint < nbptsZ - 1)
    {
        vector<double> zcoord(2);
        double local_radius = radius * cos(totAngle);
        zcoord[0] = gravityCenterZ + radius * sin(totAngle);
        zcoord[1] = gravityCenterZ - radius * sin(totAngle);

        half_perimeter = pi * local_radius;
        nbpts = 2 * int(half_perimeter / critical_distance);
        angleLocal = 2 * pi / nbpts;

        for (size_t jindex = 0; jindex < 2; ++jindex)
            for (size_t i = 0; i <= nbpts - 1; ++i)
            {
                point(0) = local_radius * cos(i * angleLocal) + gravityCenterX;
                point(1) = local_radius * sin(i * angleLocal) + gravityCenterY;
                point(2) = zcoord[jindex];
                add_boundary_point(point, critical_distance);
            }

        totAngle += angle;
        ++ipoint;
    }

    //-- Interior points
    cout << "Hello World from setting DLMFD points" << endl;
}

//---------------------------------------------------------------------------
double DLMFD_Sphere::get_circumscribed_radius() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DS_Sphere:: get_circumscribed_radius()");

    return (ptr_FSrigidbody->get_circumscribed_radius());
}

//---------------------------------------------------------------------------
geomVector const *DLMFD_Sphere::get_ptr_to_gravity_centre() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DS_Sphere:: get_ptr_to_gravity_centre( )");

    return (dynamic_cast<FS_RigidBody *>(ptr_FSrigidbody)->get_ptr_to_gravity_centre());
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

//---------------------------------------------------------------------------
void DLMFD_Sphere::update_RB_position_and_velocity(geomVector const &pos,
                                                   geomVector const &vel,
                                                   geomVector const &ang_vel,
                                                   vector<geomVector> const &periodic_directions, double const &time_step)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: update_RB_position_and_velocity");

    return (ptr_FSrigidbody->update_RB_position_and_velocity(pos, vel, ang_vel, periodic_directions, time_step));
}