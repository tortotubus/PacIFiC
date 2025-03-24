#include <DLMFD_3Dcylinder.hh>
#include <DLMFD_RigidBody.hh>
#include <DLMFD_InteriorMultiplierPoint.hh>
#include <DLMFD_BoundaryMultiplierPoint.hh>
#include <DLMFD_FictitiousDomain.hh>
#include <FS_3Dcylinder.hh>
#include <math.h>
#include <iostream>
#include <set>
#include <utility>
#include <sstream>
#include <string>
#include <climits>
using namespace std;
#define THRESHOLD 1.e-7
//---------------------------------------------------------------------------
DLMFD_3Dcylinder::DLMFD_3Dcylinder() : DLMFD_RigidBody()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dcylinder:: DLMFD_3Dcylinder");
}

//---------------------------------------------------------------------------
DLMFD_3Dcylinder::DLMFD_3Dcylinder(FS_RigidBody *pgrb,
                                   const bool &are_particles_fixed,
                                   FV_DiscreteField *pField_,
                                   double const critical_distance_) : DLMFD_RigidBody(pgrb, are_particles_fixed, pField_)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dcylinder:: DLMFD_3Dcylinder");

    // Set the FS 3Dcylinder additional parameters object
    set_ptr_FS_3Dcylinder_Additional_Param();

    gravity_center = *get_ptr_to_gravity_centre();
    radius = get_circumscribed_radius();

    // TODO: maybe no need to define these objects
    BottomCenter = pagp->BottomCenter;
    TopCenter = pagp->TopCenter;
    BottomToTopVec = pagp->BottomToTopVec;
    RadialRefVec = pagp->RadialRefVec;
    cylinder_radius = pagp->cylinder_radius;
    cylinder_height = pagp->cylinder_height;

    nIP = 0;
    nBP = 0;
    nIPHZ = 0;
    nBPHZ = 0;

    // Allocate default DLMFD points and DLMFD vectors
    allocate_default_listOfPointsAndVectors_3Dcylinder(critical_distance_, pField_);

    // Set DLMFD points
    set_all_points(pField_, critical_distance_);
}

//---------------------------------------------------------------------------
DLMFD_3Dcylinder::~DLMFD_3Dcylinder()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dcylinder:: ~DLMFD_3Dcylinder");
}

//---------------------------------------------------------------------------
void DLMFD_3Dcylinder::set_ptr_FS_3Dcylinder_Additional_Param()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dcylinder::set_ptr_FS_3Dcylinder_Additional_Param");

    pagp = get_ptr_FS_3Dcylinder_Additional_Param();
}

//---------------------------------------------------------------------------
void DLMFD_3Dcylinder::set_all_points(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dcylinder:: set_all_points");

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
void DLMFD_3Dcylinder::set_boundary_points_list(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dcylinder::set_boundary_points_list");

    FV_Mesh const *primary_grid = pField->primary_grid();

    list<DLMFD_BoundaryMultiplierPoint *>::iterator bp = boundary_points.begin();
    list<DLMFD_BoundaryMultiplierPoint *>::iterator bphz = halozone_boundary_points.begin();
    for (size_t i = 0; i < nBP; ++i)
        bp++;
    for (size_t i = 0; i < nBPHZ; ++i)
        bphz++;

    double pi = acos(-1.);
    double spacing = critical_distance * DLMFD_FictitiousDomain::BoundaryPointsSpacing_coef;
    size_t npts_radius = size_t(cylinder_radius / spacing) + 1;
    double delta_radius = cylinder_radius / (npts_radius - 1);
    size_t npts_height = size_t(cylinder_height / sqrt(3) * 2 / spacing) + 1;
    double delta_height = cylinder_height / (npts_height - 1);
    double local_radius, local_angle, local_radius_ratio;
    size_t npts_local_radius;
    geomVector unit_axial(BottomToTopVec), newPoint(3);
    unit_axial /= unit_axial.calcNorm();
    geomVector n_cross_rad(unit_axial ^ RadialRefVec);

    // Bottom and Top centers
    setBndPoint(BottomCenter, primary_grid, bp, bphz);
    setBndPoint(TopCenter, primary_grid, bp, bphz);

    // Cylinder height (diamond meshing)
    npts_local_radius = size_t(2. * pi * cylinder_radius / spacing);

    for (size_t i = 0; i < npts_height; ++i)
    {
        double bin;
        // odd or even
        if (i % 2 == 0)
        {
            bin = 0.;
        }
        else
        {
            bin = 1.;
        }
        for (size_t j = 0; j < npts_local_radius; ++j)
        {
            local_angle = 2. * pi * double(j) / double(npts_local_radius) + pi * bin / double(npts_local_radius);
            newPoint = cos(local_angle) * RadialRefVec + sin(local_angle) * n_cross_rad + BottomCenter + i * delta_height * unit_axial;
            setBndPoint(newPoint, primary_grid, bp, bphz);
        }
    }

    // Spiral layout
    double goldenAngle = pi * (3 - sqrt(5));
    size_t nPtsDisk;
    double theoreticalPacking = pi / 2 / sqrt(3);

    // Fixed-point algorithm to obtain the right number of point on the disk
    double nPtsDiski = 4. * cylinder_radius * cylinder_radius / (spacing * spacing) * theoreticalPacking;
    double nPtsDiskOld = 0;

    while (fabs(nPtsDiski - nPtsDiskOld) > 1e-2)
    {
        nPtsDiskOld = nPtsDiski;
        nPtsDiski = 4 * pow(cylinder_radius, 2.0) / pow(spacing, 2.0) * (theoreticalPacking - 1 / (pow(nPtsDiski / theoreticalPacking, 0.5) + 1));
    }
    nPtsDisk = size_t(nPtsDiski);

    for (size_t j = 1; j <= nPtsDisk; ++j)
    {
        local_angle = double(j) * goldenAngle;
        newPoint = cos(local_angle) * RadialRefVec + sin(local_angle) * n_cross_rad;
        newPoint *= sqrt(double(j) / nPtsDisk) * (1. - spacing / 2. / cylinder_radius);

        // Bottom disk
        newPoint += BottomCenter;
        setBndPoint(newPoint, primary_grid, bp, bphz);

        // Top disk
        newPoint += BottomToTopVec;
        setBndPoint(newPoint, primary_grid, bp, bphz);
    }
}

//---------------------------------------------------------------------------
void DLMFD_3Dcylinder::set_interior_points_list(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dcylinder:: set_interior_points_list");

    size_t dim = 3;
    doubleVector coor_min(dim, 0.), coor_max(dim, 0.);
    size_t ncomps = pField->nb_components();
    vector<doubleVector const *> work(dim, NULL);
    vector<vector<doubleVector const *>> field_mesh(ncomps, work);

    double x, y, z;
    geomVector node(3);

    // Cuboid sub-domain occupied by the component
    // Correction by 0.1*critical_distance avoids rounding errors associated
    // to comparison of doubles and truncation errors from Grains3D on corners
    // coordinates
    for (size_t m = 0; m < dim; ++m)
    {
        coor_min(m) = gravity_center(m) - radius;
        coor_max(m) = gravity_center(m) + radius;

        for (size_t comp = 0; comp < ncomps; ++comp)
        {
            field_mesh[comp][m] = pField->get_DOF_coordinates_vector(comp, m);
            (*index_min)(comp, m) = FV_Mesh::min_index(field_mesh[comp][m], coor_min(m) - 0.1 * critical_distance);
            (*index_max)(comp, m) = FV_Mesh::max_index(field_mesh[comp][m], coor_max(m) + 0.1 * critical_distance);
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
            node(0) = x;
            for (size_t j = (*index_min)(comp, 1); j <= (*index_max)(comp, 1); ++j)
            {
                y = (*field_mesh[comp][1])(j);
                node(1) = y;
                for (size_t k = (*index_min)(comp, 2); k <= (*index_max)(comp, 2); ++k)
                {
                    z = (*field_mesh[comp][2])(k);
                    node(2) = z;

                    if (isIn(node))
                    {
                        // !! Add interior points that are unknowns only !!
                        if (pField->DOF_is_unknown(i, j, k, comp) && !pField->DOF_is_constrained(i, j, k, comp))
                        {
                            pField->set_DOF_constrained(i, j, k, comp);
                            if (pField->DOF_is_unknown_handled_by_proc(i, j, k, comp))
                            {
                                set_interior_point(comp, node, i, j, k, ip);
                            }
                            else
                            {
                                set_halozone_interior_point(comp, node, i, j, k, iphz);
                            }
                        }
                    }
                }
            }
        }
}

//---------------------------------------------------------------------------
void DLMFD_3Dcylinder::update()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dcylinder::update");

    gravity_center = *get_ptr_to_gravity_centre();
    radius = get_circumscribed_radius();
    translational_velocity = get_rigid_body_translational_velocity();
    angular_velocity_3D = get_rigid_body_angular_velocity();
    inertia_3D = get_rigid_body_inertia();

    BottomCenter = pagp->BottomCenter;
    TopCenter = pagp->TopCenter;
    BottomToTopVec = pagp->BottomToTopVec;
    RadialRefVec = pagp->RadialRefVec;
    cylinder_radius = pagp->cylinder_radius;
    cylinder_height = pagp->cylinder_height;
}

//---------------------------------------------------------------------------
FS_3Dcylinder_Additional_Param const *DLMFD_3Dcylinder::get_ptr_FS_3Dcylinder_Additional_Param()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dcylinder::get_ptr_FS_3Dcylinder_Additional_Param");

    return (dynamic_cast<FS_3Dcylinder *>(ptr_FSrigidbody)->get_ptr_FS_3Dcylinder_Additional_Param());
}

//---------------------------------------------------------------------------
void DLMFD_3Dcylinder::allocate_default_listOfPointsAndVectors_3Dcylinder(const double &critical_distance, FV_DiscreteField *pField)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DS_3Dcylinder::allocate_default_listOfPointsAndVectors_3Dcylinder()");

    double pi = acos(-1.);
    double spacing = critical_distance;
    double mesh_size = critical_distance / sqrt(3.);
    size_t security_bandwidth = pField->primary_grid()->get_security_bandwidth();
    size_t nSI = size_t(2. * cylinder_radius / mesh_size) + 1;
    size_t nHI = size_t(cylinder_height / mesh_size) + 1;
    size_t nbIPdef = size_t(3.2 * nHI * pi * nSI * nSI / 4.);
    size_t nbIPHZdef = security_bandwidth * 3 * nSI * nSI;

    // Exterior
    size_t npts_radius = size_t(cylinder_radius / critical_distance) + 1;
    double delta_radius = cylinder_radius / (npts_radius - 1);
    size_t npts_height = size_t(cylinder_height / sqrt(3) * 2 / spacing) + 1;
    double local_radius;
    size_t npts_local_radius;
    size_t nbBPdef = 0;

    // Bottom and Top centers
    nbBPdef += 2;

    // Cylinder height
    npts_local_radius = size_t(2. * pi * cylinder_radius / critical_distance);
    nbBPdef += npts_local_radius * npts_height;

    // Spiralling distribution
    double goldenAngle = pi * (3 - sqrt(5));
    size_t nPtsDisk;
    double theoreticalPacking = pi / 2 / sqrt(3);
    double nPtsDiski = 4. * cylinder_radius * cylinder_radius / (spacing * spacing) * theoreticalPacking;
    double nPtsDiskOld = 0;

    // Loop giving the optimal number of points on the disk (converge in 4-5 iterations)
    while (fabs(nPtsDiski - nPtsDiskOld) > 1e-2)
    {
        nPtsDiskOld = nPtsDiski;
        nPtsDiski = 4 * pow(cylinder_radius, 2.0) / pow(spacing, 2.0) * (theoreticalPacking - 1 / (pow(nPtsDiski / theoreticalPacking, 0.5) + 1));
    }
    nbBPdef += 2 * size_t(nPtsDiski);

    size_t nbBPHZdef = size_t(0.5 * nbBPdef);

    allocate_default_listOfPointsAndVectors(nbIPdef, nbBPdef, nbIPHZdef, nbBPHZdef);
}

//---------------------------------------------------------------------------
bool DLMFD_3Dcylinder::isIn(const geomVector &point) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dcylinder:: isIn");

    // return (ptr_FSrigidbody->isIn(point));

    // TODO: there is something to do with this function, see constructor.

    bool b_isIn = false;

    geomVector BottomToPoint(point - BottomCenter);
    double dot = (BottomToPoint, BottomToTopVec) / cylinder_height;

    if (dot < cylinder_height && dot > 0.)
        if (BottomToPoint.calcNormSquare() - dot * dot < cylinder_radius * cylinder_radius)
            b_isIn = true;

    if (periodic_directions)
    {
        for (size_t i = 0; (i < periodic_directions->size()) && !b_isIn; ++i)
        {
            BottomToPoint = point - (BottomCenter + (*periodic_directions)[i]);
            dot = (BottomToPoint, BottomToTopVec) / cylinder_height;

            if (dot <= cylinder_height && dot > 0)
                if (BottomToPoint.calcNormSquare() - dot * dot < cylinder_radius * cylinder_radius)
                    b_isIn = true;
        }
    }

    return b_isIn;
}

//---------------------------------------------------------------------------
void DLMFD_3Dcylinder::translateGeometricFeatures(geomVector const &newg)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dcylinder:: translateGeometricFeatures");

    if (gravity_center.calcDist(newg) > 1.e-2 * cylinder_radius)
    {
        geomVector translation_vec = newg - gravity_center;
        gravity_center = newg;
        BottomCenter += translation_vec;
        TopCenter += translation_vec;
    }
}
