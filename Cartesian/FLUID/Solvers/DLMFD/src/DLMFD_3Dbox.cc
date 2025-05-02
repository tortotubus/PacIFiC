#include <DLMFD_3Dbox.hh>
#include <DLMFD_RigidBody.hh>
#include <DLMFD_InteriorMultiplierPoint.hh>
#include <DLMFD_BoundaryMultiplierPoint.hh>
#include <FS_3Dbox.hh>
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
DLMFD_3Dbox::DLMFD_3Dbox() : DLMFD_RigidBody()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox:: DLMFD_3Dbox");
}

//---------------------------------------------------------------------------
DLMFD_3Dbox::DLMFD_3Dbox(FS_RigidBody *pgrb,
                         const bool &are_particles_fixed,
                         FV_DiscreteField *pField_,
                         double const critical_distance_) : DLMFD_RigidBody(pgrb, are_particles_fixed, pField_),
                                                            g2(NULL)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox:: DLMFD_3Dbox");

    // Set the FS 3Dbox additional parameters object
    set_ptr_FS_3Dbox_Additional_Param();

    gravity_center = *get_ptr_to_gravity_centre();
    radius = get_circumscribed_radius();
    corners = pagp->corners;
    facesVec = pagp->facesVec;

    coor_min.resize(3);
    coor_max.resize(3);

    set_all_MAC(pField_, critical_distance_);
}

//---------------------------------------------------------------------------
DLMFD_3Dbox::~DLMFD_3Dbox()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox:: ~DLMFD_3Dbox");

    corners.clear();

    for (size_t i = 0; i < facesVec.size(); ++i)
        facesVec[i].clear();
    facesVec.clear();

    if (g2)
        delete g2;
}

//---------------------------------------------------------------------------
void DLMFD_3Dbox::set_ptr_FS_3Dbox_Additional_Param()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox::set_ptr_FS_3Dbox_Additional_Param");

    pagp = get_ptr_FS_3Dbox_Additional_Param();
}

//---------------------------------------------------------------------------
void DLMFD_3Dbox::set_all_MAC(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox::set_all_MAC");

    FV_Mesh const *primary_grid = pField->primary_grid();

    size_t i, nper = periodic_directions ? periodic_directions->size() : 0;

    if (component_type == "O" || component_type == "PO")
    {
        if (g2 == NULL)
            g2 = new geomVector(3);
        (*g2)(0) = gravity_center(0) + 1e-4 * radius * random() / double(INT_MAX);
        (*g2)(1) = gravity_center(1) + 1e-4 * radius * random() / double(INT_MAX);
        (*g2)(2) = gravity_center(2) + 1e-4 * radius * random() / double(INT_MAX);
    }

    ndof = 0;
    nIP = 0;
    nBP = 0;
    nIPHZ = 0;
    nBPHZ = 0;

    bool in = false;
    if (primary_grid->is_in_domain_with_halozone_plus_ext(gravity_center(0),
                                                          gravity_center(1),
                                                          gravity_center(2),
                                                          radius))
    {
        in = true;
        if (interior_points.empty())
            allocate_default_listOfPointsAndVectors_3Dbox(critical_distance, pField);

        set_all_points(pField, critical_distance);
    }
}

//---------------------------------------------------------------------------
void DLMFD_3Dbox::set_all_points(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox:: set_all_points");

    //-- Boundary points
    set_boundary_points_list(pField, critical_distance);

    //-- Interior points
    set_interior_points_list(pField, critical_distance);
}

//---------------------------------------------------------------------------
void DLMFD_3Dbox::set_boundary_points_list(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox::set_boundary_points_list");

    FV_Mesh const *primary_grid = pField->primary_grid();

    list<DLMFD_BoundaryMultiplierPoint *>::iterator bp = boundary_points.begin();
    list<DLMFD_BoundaryMultiplierPoint *>::iterator bphz = halozone_boundary_points.begin();
    for (size_t i = 0; i < nBP; ++i)
        bp++;
    for (size_t i = 0; i < nBPHZ; ++i)
        bphz++;

    vector<vector<size_t>>::iterator faceIter;
    vector<geomVector>::iterator cornIter;

    set<pair<size_t, size_t>> facesBnd;
    pair<size_t, size_t> pairToAdd;
    double edge_length = 0.;
    bool b_isometric = true;

    // Generate nodes on the polyhedron edges
    for (faceIter = facesVec.begin(); faceIter != facesVec.end(); ++faceIter)
    {
        vector<size_t>::iterator itrLocNodes;
        for (itrLocNodes = faceIter->begin() + 1; itrLocNodes != faceIter->end(); ++itrLocNodes)
        {
            if (*(itrLocNodes - 1) > *itrLocNodes)
                pairToAdd = pair<size_t, size_t>(*itrLocNodes, *(itrLocNodes - 1));
            else
                pairToAdd = pair<size_t, size_t>(*(itrLocNodes - 1), *itrLocNodes);

            if (facesBnd.insert(pairToAdd).second)
            {
                setPtsOnEdge(corners[*(itrLocNodes - 1)],
                             corners[*itrLocNodes],
                             critical_distance,
                             false,
                             primary_grid,
                             bp,
                             bphz);

                double el = corners[*(itrLocNodes - 1)].calcDist(corners[*itrLocNodes]);

                if (!edge_length)
                    edge_length = el;
                else if (fabs(el - edge_length) > 0.01 * radius)
                    b_isometric = false;
            }
        }

        // add the last edge of the face to the list
        if (faceIter->front() > faceIter->back())
            pairToAdd = pair<size_t, size_t>(faceIter->back(), faceIter->front());
        else
            pairToAdd = pair<size_t, size_t>(faceIter->front(), faceIter->back());

        if (facesBnd.insert(pairToAdd).second)
        {
            setPtsOnEdge(corners[faceIter->front()],
                         corners[faceIter->back()],
                         critical_distance, false,
                         primary_grid,
                         bp,
                         bphz);

            double el = corners[faceIter->front()].calcDist(corners[faceIter->back()]);

            if (fabs(el - edge_length) > 0.01 * radius)
                b_isometric = false;
        }
    }

    // Generate the interior nodes on each face
    geomVector const *RefCorner = NULL;
    geomVector dir0(3), dir1(3), newPoint(3);
    size_t nbptsLocal0, nbptsLocal1;

    for (faceIter = facesVec.begin(); faceIter != facesVec.end(); ++faceIter)
    {
        RefCorner = &corners[faceIter->front()];
        dir0 = corners[*(faceIter->begin() + 1)] - *RefCorner;
        nbptsLocal0 = size_t(dir0.calcNorm() / critical_distance);
        dir0 /= double(nbptsLocal0);
        dir1 = corners[faceIter->back()] - *RefCorner;
        nbptsLocal1 = size_t(dir1.calcNorm() / critical_distance);
        dir1 /= double(nbptsLocal1);

        for (size_t i0 = 1; i0 <= nbptsLocal0 - 1; ++i0)
            for (size_t i1 = 1; i1 <= nbptsLocal1 - 1; ++i1)
            {
                newPoint = *RefCorner + dir0 * double(i0) + dir1 * double(i1);
                setBndPoint(newPoint, primary_grid, bp, bphz);
            }
    }

    // Add the corners of the particle
    for (cornIter = corners.begin(); cornIter != corners.end(); ++cornIter)
    {
        setBndPoint(*cornIter, primary_grid, bp, bphz);
    }
}

//---------------------------------------------------------------------------
void DLMFD_3Dbox::set_interior_points_list(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox:: set_interior_points_list");

    size_t dim = 3;
    size_t ncomps = pField->nb_components();
    vector<doubleVector const *> work(dim, NULL);
    vector<vector<doubleVector const *>> field_mesh(ncomps, work);

    double x, y, z;
    geomVector node(3);

    vector<geomVector>::iterator cornIter;

    // Cuboid sub-domain occupied by the component
    // Correction by 0.1*critical_distance avoids rounding errors associated
    // to comparison of doubles and truncation errors from Grains3D on corners
    // coordinates
    coor_min.set(1.e20);
    coor_max.set(-1.e20);
    for (size_t m = 0; m < dim; ++m)
        for (cornIter = corners.begin(); cornIter != corners.end(); cornIter++)
        {
            coor_min(m) = (*cornIter)(m) < coor_min(m) ? (*cornIter)(m) : coor_min(m);
            coor_max(m) = (*cornIter)(m) > coor_max(m) ? (*cornIter)(m) : coor_max(m);
        }

    for (size_t m = 0; m < dim; ++m)
        for (size_t comp = 0; comp < ncomps; ++comp)
        {
            field_mesh[comp][m] = pField->get_DOF_coordinates_vector(comp, m);
            (*index_min)(comp, m) = FV_Mesh::min_index(field_mesh[comp][m], coor_min(m) - 0.1 * critical_distance);
            (*index_max)(comp, m) = FV_Mesh::max_index(field_mesh[comp][m], coor_max(m) + 0.1 * critical_distance);
        }

    // Interior points
    list<DLMFD_InteriorMultiplierPoint *>::iterator ip = interior_points.begin();
    list<DLMFD_InteriorMultiplierPoint *>::iterator iphz = halozone_interior_points.begin();
    for (size_t i = 0; i < nIP; ++i)
        ip++;
    for (size_t i = 0; i < nIPHZ; ++i)
        iphz++;

    size_t __nip = interior_points.size(),
           __niphz = halozone_interior_points.size();

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
                                set_interior_point(comp, node, i, j, k, ip, __nip);
                            }
                            else
                            {
                                set_halozone_interior_point(comp, node, i, j, k, iphz, __niphz);
                            }
                        }
                    }
                }
            }
        }
}

//---------------------------------------------------------------------------
void DLMFD_3Dbox::update()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox::update");

    gravity_center = *get_ptr_to_gravity_centre();
    radius = get_circumscribed_radius();
    translational_velocity = get_rigid_body_translational_velocity();
    angular_velocity_3D = get_rigid_body_angular_velocity();
    inertia_3D = get_rigid_body_inertia();

    if (periodic_directions)
        periodic_directions = get_ptr_to_periodic_directions();

    corners = pagp->corners;
    facesVec = pagp->facesVec;
}

//---------------------------------------------------------------------------
FS_3Dbox_Additional_Param const *DLMFD_3Dbox::get_ptr_FS_3Dbox_Additional_Param()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox::get_ptr_FS_3Dbox_Additional_Param");

    return (dynamic_cast<FS_3Dbox *>(ptr_FSrigidbody)->get_ptr_FS_3Dbox_Additional_Param());
}

//---------------------------------------------------------------------------
void DLMFD_3Dbox::allocate_default_listOfPointsAndVectors_3Dbox(const double &critical_distance, FV_DiscreteField *pField)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DS_3Dbox::allocate_default_listOfPointsAndVectors_3Dbox()");

    coor_min.set(1.e20);
    coor_max.set(-1.e20);
    for (size_t m = 0; m < 3; ++m)
        for (vector<geomVector>::iterator cornIter = corners.begin();
             cornIter != corners.end(); cornIter++)
        {
            coor_min(m) = (*cornIter)(m) < coor_min(m) ? (*cornIter)(m) : coor_min(m);
            coor_max(m) = (*cornIter)(m) > coor_max(m) ? (*cornIter)(m) : coor_max(m);
        }

    double mesh_size = critical_distance / sqrt(3.);
    size_t security_bandwidth = pField->primary_grid()->get_security_bandwidth();

    size_t nxI = size_t((coor_max(0) - coor_min(0)) / mesh_size) + 1;
    size_t nyI = size_t((coor_max(1) - coor_min(1)) / mesh_size) + 1;
    size_t nzI = size_t((coor_max(2) - coor_min(2)) / mesh_size) + 1;

    size_t nxB = size_t((coor_max(0) - coor_min(0)) / critical_distance) + 1;
    size_t nyB = size_t((coor_max(1) - coor_min(1)) / critical_distance) + 1;
    size_t nzB = size_t((coor_max(2) - coor_min(2)) / critical_distance) + 1;

    size_t nbIPdef = 3 * nxI * nyI * nzI;
    size_t nbBPdef = 2 * (nxB * nyB + nxB * nzB + nyB * nzB);
    size_t nbIPHZdef = security_bandwidth * (nxI * nyI + nxI * nzI + nyI * nzI);
    size_t nbBPHZdef = size_t(0.5 * nbBPdef);

    allocate_default_listOfPointsAndVectors(nbIPdef, nbBPdef, nbIPHZdef, nbBPHZdef);
}

//---------------------------------------------------------------------------
double DLMFD_3Dbox::calcPointDeterm4by4(const geomVector &pointOne,
                                        const geomVector &pointTwo,
                                        const geomVector &pointThree,
                                        const geomVector &pointFour) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox::calcPointDeterm4by4");

    double x2 = pointTwo(0);
    double y2 = pointTwo(1);
    double z2 = pointTwo(2);
    double x3 = pointThree(0);
    double y3 = pointThree(1);
    double z3 = pointThree(2);
    double x4 = pointFour(0);
    double y4 = pointFour(1);
    double z4 = pointFour(2);

    double retVal =
        pointOne(0) * (y2 * z3 + y3 * z4 + y4 * z2 - y4 * z3 - y3 * z2 - y2 * z4) -
        pointOne(1) * (x2 * z3 + x3 * z4 + x4 * z2 - x4 * z3 - x3 * z2 - x2 * z4) +
        pointOne(2) * (x2 * y3 + x3 * y4 + x4 * y2 - x4 * y3 - x3 * y2 - x2 * y4) -
        (x2 * y3 * z4 + x3 * y4 * z2 + x4 * y2 * z3 - x4 * y3 * z2 - x3 * y2 * z4 - x2 * y4 * z3);

    return retVal;
}

//---------------------------------------------------------------------------
bool DLMFD_3Dbox::checkPointInTetrahedron(const geomVector &pointOne,
                                          const geomVector &pointTwo,
                                          const geomVector &pointThree,
                                          const geomVector &pointFour,
                                          const geomVector &pointToCheck) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox::checkPointInTetrahedron");

    // Link: http://steve.hollasch.net/cgindex/geometry/ptintet.html

    double detTot = calcPointDeterm4by4(pointOne, pointTwo, pointThree, pointFour);

    double detOne = calcPointDeterm4by4(pointToCheck, pointTwo, pointThree, pointFour);

    double detTwo = calcPointDeterm4by4(pointOne, pointToCheck, pointThree, pointFour);

    double detThree = calcPointDeterm4by4(pointOne, pointTwo, pointToCheck, pointFour);

    double detFour = calcPointDeterm4by4(pointOne, pointTwo, pointThree, pointToCheck);

    bool in = false;

    double sumSubElem = detOne + detTwo + detThree + detFour;

    if (fabs(detTot - sumSubElem) > THRESHOLD)
        std::cout << "ERROR: summation error in determinat 3D : "
                  << detTot - sumSubElem << endl;

    if (detTot == 0)
    {
        std::cout << "degenerated tetrahedron: det == 0 " << endl;
        abort();
    }

    if ((detTot * detOne >= -0.) && (detTot * detTwo >= -0.) && (detTot * detThree >= -0.) && (detTot * detFour >= -0.))
        in = true;

    return in;
}

//---------------------------------------------------------------------------
bool DLMFD_3Dbox::isIn(const geomVector &point) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox:: isIn");

    // return (ptr_FSrigidbody->isIn(point));

    // TODO: there is something to do with this function, see constructor.

    bool b_isIn = false;

    bool is_an_obstacle = component_type == "O" || component_type == "PO";

    for (vector<vector<size_t>>::const_iterator faceIter = facesVec.begin();
         faceIter != facesVec.end() && !b_isIn;
         ++faceIter)
        for (size_t idxPts = 2; idxPts < 4 && !b_isIn; ++idxPts)
            if (checkPointInTetrahedron(corners[(*faceIter)[0]],
                                        corners[(*faceIter)[idxPts - 1]],
                                        corners[(*faceIter)[idxPts]],
                                        is_an_obstacle ? *g2 : gravity_center,
                                        point))

                b_isIn = true;

    return (b_isIn);
}

//---------------------------------------------------------------------------
bool DLMFD_3Dbox::proximityQuery(DLMFD_RigidBody const *second_component, const double &distance) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox::proximityQuery");

    geomVector const *gc2 = second_component->get_ptr_to_gravity_centre();
    double radius2 = second_component->get_circumscribed_radius();

    bool condition_one = fabs(0.5 * (coor_min(0) + coor_max(0)) - (*gc2)(0)) <= 0.5 * (coor_max(0) - coor_min(0)) + radius2;
    bool condition_two = fabs(0.5 * (coor_min(1) + coor_max(1)) - (*gc2)(1)) <= 0.5 * (coor_max(1) - coor_min(1)) + radius2;
    bool condition_three = fabs(0.5 * (coor_min(2) + coor_max(2)) - (*gc2)(2)) <= 0.5 * (coor_max(2) - coor_min(2)) + radius2;

    return (condition_one && condition_two && condition_three);
}

//---------------------------------------------------------------------------
void DLMFD_3Dbox::translateGeometricFeatures(geomVector const &newg)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_3Dbox:: translateGeometricFeatures");

    gravity_center = newg;
}
