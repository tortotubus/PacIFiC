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
DLMFD_Sphere::DLMFD_Sphere(FS_RigidBody *pgrb,
                           FV_DiscreteField *pField_,
                           double const critical_distance_) : DLMFD_RigidBody(pgrb, pField_)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: DLMFD_Sphere");

    gravity_center = *get_ptr_to_gravity_centre();
    radius = get_circumscribed_radius();

    nIP = 0;
    nBP = 0;
    nIPHZ = 0;
    nBPHZ = 0;

    // Allocate default DLMFD points and DLMFD vectors
    allocate_default_listOfPointsAndVectors_Sphere(critical_distance_, pField_);

    // Set DLMFD points
    set_all_points(pField_, critical_distance_);
}

//---------------------------------------------------------------------------
DLMFD_Sphere::~DLMFD_Sphere()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: ~DLMFD_Sphere");
}

//---------------------------------------------------------------------------
void DLMFD_Sphere::set_all_points(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: set_all_points");

    nIP = 0;
    nBP = 0;
    nIPHZ = 0;
    nBPHZ = 0;

    //-- Boundary points
    set_boundary_points_list(pField, critical_distance);

    //-- Interior points
    set_interior_points_list(pField, critical_distance);
}

//---------------------------------------------------------------------------
void DLMFD_Sphere::set_boundary_points_list(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: set_boundary_points_list");

    double pi = acos(-1.);
    double spacing = critical_distance;

    FV_Mesh const *primary_grid = pField->primary_grid();

    list<DLMFD_BoundaryMultiplierPoint *>::iterator bp = boundary_points.begin();
    list<DLMFD_BoundaryMultiplierPoint *>::iterator bphz = halozone_boundary_points.begin();
    for (size_t i = 0; i < nBP; ++i)
        bp++;
    for (size_t i = 0; i < nBPHZ; ++i)
        bphz++;

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
        geomVector GC_plus_point = gravity_center + point;

        if (primary_grid->is_in_domain_with_halozone(GC_plus_point(0), GC_plus_point(1), GC_plus_point(2)))
        {
            if (primary_grid->is_in_domain_on_current_processor(GC_plus_point(0), GC_plus_point(1), GC_plus_point(2)))
            {
                set_boundary_point(GC_plus_point, bp);
            }
            else
            {
                set_halozone_boundary_point(GC_plus_point, bphz);
            }
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_Sphere::set_interior_points_list(FV_DiscreteField *pField, double critical_distance)
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
    list<DLMFD_InteriorMultiplierPoint *>::iterator ip = interior_points.begin();
    list<DLMFD_InteriorMultiplierPoint *>::iterator iphz = halozone_interior_points.begin();
    for (size_t i = 0; i < nIP; ++i)
        ip++;
    for (size_t i = 0; i < nIPHZ; ++i)
        iphz++;

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
                        // !! Add interior points that are unknowns only !!
                        if (pField->DOF_is_unknown(i, j, k, comp) && !pField->DOF_is_constrained(i, j, k, comp))
                        {
                            pField->set_DOF_constrained(i, j, k, comp);
                            if (pField->DOF_is_unknown_handled_by_proc(i, j, k, comp))
                            {
                                set_interior_point(comp, point, i, j, k, ip);
                            }
                            else
                            {
                                set_halozone_interior_point(comp, point, i, j, k, iphz);
                            }
                        }
                    }
                }
            }
        }
}

//---------------------------------------------------------------------------
void DLMFD_Sphere::allocate_default_listOfPointsAndVectors_Sphere(const double &critical_distance, FV_DiscreteField *pField)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DS_Sphere::allocate_default_listOfPointsAndVectors_Sphere()");

    double pi = acos(-1.);
    double mesh_size = critical_distance / sqrt(3.);
    double spacing = critical_distance;

    size_t security_bandwidth = pField->primary_grid()->get_security_bandwidth();

    size_t n = size_t(2. * radius / mesh_size) + 1;

    size_t nb = size_t(pi * radius / spacing);
    size_t nbIPdef = size_t(3.2 * pi * n * n * n / 6.);
    size_t nbIPHZdef = security_bandwidth * n * n * 3 * 2;

    size_t nbBPdef = size_t(pow(3.809 * radius / spacing, 2.));

    size_t nbBPHZdef = size_t(0.25 * nbBPdef);

    allocate_default_listOfPointsAndVectors(nbIPdef, nbBPdef, nbIPHZdef, nbBPHZdef);
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

    gravity_center = *get_ptr_to_gravity_centre();
    radius = get_circumscribed_radius();
    translational_velocity = get_rigid_body_translational_velocity();
    angular_velocity_3D = get_rigid_body_angular_velocity();
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
void DLMFD_Sphere::erase_critical_interior_points_PerProc(double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: erase_critical_interior_points_PerProc");

    list<DLMFD_BoundaryMultiplierPoint *>::iterator imp, ihp;
    list<DLMFD_InteriorMultiplierPoint *>::iterator ivi;
    double dist = 0., reduced_radius = radius - critical_distance;
    bool erase = false;

    if (nIP)
    {
        ivi = interior_points.begin();
        for (size_t j = 0; j < nIP; ++j, ivi++)
        {
            if ((*ivi)->isValid())
            {
                erase = false;
                dist = gravity_center.calcDist((*ivi)->get_coordinates());
                if (dist > reduced_radius)
                {
                    // Compare to boundary points on this process
                    if (!boundary_points.empty())
                    {
                        imp = boundary_points.begin();
                        for (size_t i = 0; i < nBP && !erase; ++i, imp++)
                            if ((*imp)->isValid())
                                if (((*imp)->get_coordinates()).calcDist((*ivi)->get_coordinates()) < critical_distance)
                                    erase = true;
                    }

                    // Compare to boundary points in halo zone
                    if (!halozone_boundary_points.empty())
                    {
                        ihp = halozone_boundary_points.begin();
                        for (size_t i = 0; i < nBPHZ && !erase; ++i, ihp++)
                            if ((*ihp)->isValid())
                                if (((*ihp)->get_coordinates()).calcDist((*ivi)->get_coordinates()) < critical_distance)
                                    erase = true;
                    }
                }

                if (erase)
                    (*ivi)->set_validity(false);
            }
        }
    }
}
