#include <DLMFD_Sphere.hh>
#include <DLMFD_RigidBody.hh>
#include <DLMFD_InteriorMultiplierPoint.hh>
#include <DLMFD_BoundaryMultiplierPoint.hh>
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

    nIP = 0;
    nBP = 0;
}

//---------------------------------------------------------------------------
DLMFD_Sphere::~DLMFD_Sphere()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: ~DLMFD_Sphere");
}

//---------------------------------------------------------------------------
void DLMFD_Sphere::set_all_points(FV_DiscreteField const *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: set_all_points");

    nIP = 0;
    nBP = 0;

    //-- Boundary points
    set_boundary_points_list(pField, critical_distance);

    //-- Interior points
    set_interior_points_list(pField, critical_distance);
}

//---------------------------------------------------------------------------
void DLMFD_Sphere::set_boundary_points_list(FV_DiscreteField const *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: set_boundary_points_list");

    double pi = acos(-1.);
    double spacing = critical_distance;

    FV_Mesh const *primary_grid = pField->primary_grid();

    size_t nbBPdef = pow(3.809 * radius / spacing, 2.);
    if (boundary_points.empty())
        allocate_default_boundary_points_sphere(nbBPdef);

    list<DLMFD_BoundaryMultiplierPoint *>::iterator bp = boundary_points.begin();

    double spiral_spacing_correction = spacing / sqrt(3.);

    double x, y, z;
    geomVector point;

    // Number of points
    size_t NSpiral = size_t(pow(3.809 * radius / spacing, 2.));

    // Spiral points construction
    double hk, thetak, phik, phikm1 = 0., TwoPi = 2. * pi, Cspiral = 3.6 / sqrt(NSpiral), dphi = 0.;
    for (size_t k = 0; k < NSpiral; ++k)
    {
        hk = -1. + 2. * double(k) / (NSpiral - 1.);
        thetak = acos(hk);
        if (k == 0)
        {
            phik = 0.;
            thetak = pi - 0.5 * spiral_spacing_correction / radius;
        }
        else if (k == NSpiral - 1)
        {
            phik = phikm1 + 1. * dphi;
            if (phik > TwoPi)
                phik -= TwoPi;
            thetak = 0.5 * spiral_spacing_correction / radius;
        }
        else
        {
            dphi = Cspiral / sqrt(1. - hk * hk);
            phik = phikm1 + dphi;
            if (phik > TwoPi)
                phik -= TwoPi;
        }

        phikm1 = phik;

        if (k == 1)
            thetak -= 0.4 * spiral_spacing_correction / radius;

        if (k == NSpiral - 2)
        {
            phik -= 0.1 * dphi;
            thetak += 0.25 * spiral_spacing_correction / radius;
        }

        x = radius * cos(phik) * sin(thetak);
        y = radius * sin(phik) * sin(thetak);
        z = radius * cos(thetak);

        point = geomVector(x, y, z);

        set_boundary_point(gravity_center + point, bp);
    }
}

//---------------------------------------------------------------------------
void DLMFD_Sphere::set_interior_points_list(FV_DiscreteField const *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: set_interior_points_list");

    double pi = acos(-1);

    size_t dim = 3;
    doubleVector coor_min(dim, 0.), coor_max(dim, 0.);
    size_t ncomps = pField->nb_components();
    vector<doubleVector const *> work(dim, NULL);
    vector<vector<doubleVector const *>> field_mesh(ncomps, work);

    double x, y, z;
    geomVector point;

    // Cubic sub-domain occupied by the component
    // Correction by 0.0001*critical_distance avoids rounding errors associated
    // to comparison of doubles
    for (size_t m = 0; m < dim; ++m)
    {
        coor_min(m) = gravity_center(m) - radius;
        coor_max(m) = gravity_center(m) + radius;
        for (size_t comp = 0; comp < ncomps; ++comp)
        {
            field_mesh[comp][m] = pField->get_DOF_coordinates_vector(comp, m);
            (*index_min)(comp, m) = FV_Mesh::min_index(field_mesh[comp][m], coor_min(m) - 0.0001 * critical_distance);
            (*index_max)(comp, m) = FV_Mesh::max_index(field_mesh[comp][m], coor_max(m) + 0.0001 * critical_distance);
        }
    }

    // Interior points
    double mesh_size = critical_distance / sqrt(3.);
    size_t n = size_t(2. * radius / mesh_size) + 1;
    size_t nbIPdef = size_t(3.2 * pi * n * n * n / 6.);

    if (interior_points.empty())
        allocate_default_interior_points_sphere(nbIPdef);

    list<DLMFD_InteriorMultiplierPoint *>::iterator ip = interior_points.begin();

    size_t __nip = interior_points.size();
    for (size_t comp = 0; comp < ncomps; ++comp)
        for (size_t i = (*index_min)(comp, 0); i <= (*index_max)(comp, 0); ++i)
        {
            x = (*field_mesh[comp][0])(i);
            for (size_t j = (*index_min)(comp, 1); j <= (*index_max)(comp, 1); ++j)
            {
                y = (*field_mesh[comp][1])(j);
                for (size_t k = (*index_min)(comp, 2); k <= (*index_max)(comp, 2); ++k)
                {
                    z = (*field_mesh[comp][2])(k);

                    point = geomVector(x, y, z);

                    if (isIn(point))
                    {
                        set_interior_point(comp, point, i, j, k, ip);
                    }
                }
            }
        }
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
void DLMFD_Sphere::update()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: update");

    gravity_center = *DLMFD_Sphere::get_ptr_to_gravity_centre();
    radius = DLMFD_Sphere::get_circumscribed_radius();
}

//---------------------------------------------------------------------------
bool DLMFD_Sphere::isIn(const geomVector &point) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: isIn");

    double x = point(0), y = point(1), z;
    if (point.getVecSize() == 3)
        z = point(2);
    else
        z = 0.;

    return (gravity_center.calcDist(x, y, z) < radius);
}

//---------------------------------------------------------------------------
void DLMFD_Sphere::translateGeometricFeatures(geomVector const &newg)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: translateGeometricFeatures");

    gravity_center = newg;
}

//---------------------------------------------------------------------------
void DLMFD_Sphere::allocate_default_boundary_points_sphere(size_t const &nbBPdef)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: allocate_default_boundary_points");

    DLMFD_BoundaryMultiplierPoint *bmp = NULL;
    for (size_t i = 0; i < nbBPdef; ++i)
    {
        boundary_points.push_back(bmp);

        geomVector virtual_point = geomVector(0., 0., 0.);
        boundary_points.back() = new DLMFD_BoundaryMultiplierPoint(virtual_point, gravity_center);
    }
}

//---------------------------------------------------------------------------
void DLMFD_Sphere::allocate_default_interior_points_sphere(size_t const &nbIPdef)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: allocate_default_interior_points");

    DLMFD_InteriorMultiplierPoint *imp = NULL;
    for (size_t i = 0; i < nbIPdef; ++i)
    {
        interior_points.push_back(imp);

        geomVector virtual_point = geomVector(0., 0., 0.);
        interior_points.back() = new DLMFD_InteriorMultiplierPoint(0, virtual_point, 0, 0, 0, gravity_center);
    }
}