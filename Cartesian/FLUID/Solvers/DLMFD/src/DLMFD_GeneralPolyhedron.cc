#include <DLMFD_GeneralPolyhedron.hh>
#include <DLMFD_RigidBody.hh>
#include <DLMFD_InteriorMultiplierPoint.hh>
#include <DLMFD_BoundaryMultiplierPoint.hh>
#include <DLMFD_FictitiousDomain.hh>
#include <FS_GeneralPolyhedron.hh>
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
DLMFD_GeneralPolyhedron::DLMFD_GeneralPolyhedron() : DLMFD_RigidBody()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron:: DLMFD_GeneralPolyhedron");
}

//---------------------------------------------------------------------------
DLMFD_GeneralPolyhedron::DLMFD_GeneralPolyhedron(FS_RigidBody *pgrb,
                                                 const bool &are_particles_fixed,
                                                 FV_DiscreteField *pField_,
                                                 double const critical_distance_) : DLMFD_RigidBody(pgrb, are_particles_fixed, pField_),
                                                                                    g2(NULL)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron:: DLMFD_GeneralPolyhedron");

    // Set the FS GeneralPolyhedron additional parameters object
    set_ptr_FS_GeneralPolyhedron_Additional_Param();

    gravity_center = *get_ptr_to_gravity_centre();
    radius = get_circumscribed_radius();
    corners = pagp->corners;
    facesVec = pagp->facesVec;

    set_all_MAC(pField_, critical_distance_);
}

//---------------------------------------------------------------------------
DLMFD_GeneralPolyhedron::~DLMFD_GeneralPolyhedron()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron:: ~DLMFD_GeneralPolyhedron");

    corners.clear();

    for (size_t i = 0; i < facesVec.size(); ++i)
        facesVec[i].clear();
    facesVec.clear();

    if (g2)
        delete g2;
}

//---------------------------------------------------------------------------
void DLMFD_GeneralPolyhedron::set_ptr_FS_GeneralPolyhedron_Additional_Param()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron::set_ptr_FS_GeneralPolyhedron_Additional_Param");

    pagp = get_ptr_FS_GeneralPolyhedron_Additional_Param();
}

//---------------------------------------------------------------------------
void DLMFD_GeneralPolyhedron::set_all_MAC(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron::set_all_MAC");

    FV_Mesh const *primary_grid = pField->primary_grid();

    if (component_type == "O" || component_type == "PO")
    {
        if (g2 == NULL)
            g2 = new geomVector(3);
        (*g2)(0) = gravity_center(0) + 1e-4 * radius * random() / double(INT_MAX);
        (*g2)(1) = gravity_center(1) + 1e-4 * radius * random() / double(INT_MAX);
        (*g2)(2) = gravity_center(2) + 1e-4 * radius * random() / double(INT_MAX);
    }
    geomVector *g2ref = NULL;

    ndof;
    nIP = 0;
    nBP = 0;
    nIPHZ = 0;
    nBPHZ = 0;

    geomVector node(3);

    bool in = false;
    if (periodic_directions)
    {
        size_t i, j, ncorners = corners.size(), nper = periodic_directions->size();

        geomVector gvref(gravity_center);
        if (g2)
            g2ref = new geomVector(*g2);

        vector<geomVector> cornersref(ncorners, node);

        for (j = 0; j < ncorners; ++j)
            cornersref[j] = corners[j];

        for (i = 0; i < nper; ++i)
        {
            gravity_center = gvref + (*periodic_directions)[i];

            if (g2)
                *g2 = *g2ref + (*periodic_directions)[i];

            if (primary_grid->is_in_domain_with_halozone_plus_ext(gravity_center(0),
                                                                  gravity_center(1),
                                                                  gravity_center(2),
                                                                  radius))
            {
                in = true;
                for (j = 0; j < ncorners; ++j)
                    corners[j] = cornersref[j] + (*periodic_directions)[i];

                if (interior_points.empty())
                    allocate_default_listOfPointsAndVectors_GeneralPolyhedron(critical_distance, pField);

                set_all_points(pField, critical_distance);
            }
        }
        gravity_center = gvref;
        if (g2)
            *g2 = *g2ref;
        for (j = 0; j < ncorners; ++j)
            corners[j] = cornersref[j];
    }
    if (g2)
        delete g2ref;

    // !!! IMPORTANT !!!
    // DLM/FD points for the primary particle position must always be constructed
    // AFTER those for its periodic clones in case of a periodic particle
    if (primary_grid->is_in_domain_with_halozone_plus_ext(gravity_center(0),
                                                          gravity_center(1),
                                                          gravity_center(2),
                                                          radius))
    {
        in = true;
        if (interior_points.empty())
            allocate_default_listOfPointsAndVectors_GeneralPolyhedron(critical_distance, pField);

        set_all_points(pField, critical_distance);
    }

    // This method must be called once, and only once.
    // When temperature and particle_as_fixed_obtacle are combined, set_all_MAC
    // is called at each time step while allocate_exact_listOfPointInfosAndVectors
    // must be called only once. The boolean b_exactAllocation_done prevents from
    // performing the allocation twice as the return value of the method
    // allocate_exact_listOfPointInfosAndVectors is true
    if (is_particle_fixed && !b_exactAllocation_done)
        b_exactAllocation_done = allocate_exact_listOfPointInfosAndVectors();
}

//---------------------------------------------------------------------------
void DLMFD_GeneralPolyhedron::set_all_points(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron:: set_all_points");

    //-- Boundary points
    set_boundary_points_list(pField, critical_distance);

    //-- Interior points
    set_interior_points_list(pField, critical_distance);
}

//---------------------------------------------------------------------------
void DLMFD_GeneralPolyhedron::set_boundary_points_list(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron::set_boundary_points_list");

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
    double edge_length = 0., spacing = DLMFD_FictitiousDomain::BoundaryPointsSpacing_coef * critical_distance;
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
                             spacing,
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
                         spacing,
                         false,
                         primary_grid,
                         bp,
                         bphz);

            double el = corners[faceIter->front()].calcDist(corners[faceIter->back()]);
            if (fabs(el - edge_length) > 0.01 * radius)
                b_isometric = false;
        }
    }

    // Generate the interior nodes on each face
    if (b_isometric && facesVec.size() == 4)
        setBndPts_IsometricTetrahedron(corners, spacing, primary_grid, bp, bphz);
    else if (b_isometric && facesVec.size() == 6)
        setBndPts_IsometricCube(corners, spacing, primary_grid, bp, bphz);

    // Add the corners of the particle
    for (cornIter = corners.begin(); cornIter != corners.end(); ++cornIter)
        setBndPoint(*cornIter, primary_grid, bp, bphz);
}

//---------------------------------------------------------------------------
void DLMFD_GeneralPolyhedron::set_interior_points_list(FV_DiscreteField *pField, double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron:: set_interior_points_list");

    size_t dim = 3;
    doubleVector coor_min(dim, 1.e20), coor_max(dim, -1.e20);
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
                                (*ip)->set(comp, node, i, j, k, gravity_center);
                                ++nIP;
                                if (nIP == __nip)
                                {
                                    extend_ip_list(DLMFD_RigidBody::BlockSize_InteriorPoints);
                                    __nip = interior_points.size();
                                }
                                ip++;
                            }
                            else
                            {
                                (*iphz)->set(comp, node, i, j, k, gravity_center);
                                ++nIPHZ;
                                if (nIPHZ == __niphz)
                                {
                                    extend_iphz_list(DLMFD_RigidBody::BlockSize_HZ_InteriorPoints);
                                    __niphz = halozone_interior_points.size();
                                }
                                iphz++;
                            }
                        }
                    }
                }
            }
        }
}

//---------------------------------------------------------------------------
void DLMFD_GeneralPolyhedron::setBndPts_IsometricCube(vector<geomVector> const &corners_,
                                                      double spacing,
                                                      FV_Mesh const *primary_grid,
                                                      list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp,
                                                      list<DLMFD_BoundaryMultiplierPoint *>::iterator &bphz)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron::setBndPts_IsometricCube");

    vector<vector<size_t>>::iterator faceIter;
    geomVector RefCorner(3), dir0(3), dir1(3), newPoint(3);

    // Number of points
    faceIter = facesVec.begin();
    size_t nbptsLocal = int(corners_[faceIter->front()].calcDist(corners_[faceIter->back()]) / spacing);

    if (nbptsLocal > 1)
    {
        for (faceIter = facesVec.begin(); faceIter != facesVec.end(); ++faceIter)
        {
            RefCorner = corners_[faceIter->front()];
            dir0 = (corners_[*(faceIter->begin() + 1)] - RefCorner) * (1. / nbptsLocal);
            dir1 = (corners_[faceIter->back()] - RefCorner) * (1. / nbptsLocal);

            for (size_t i0 = 1; i0 <= nbptsLocal - 1; ++i0)
                for (size_t i1 = 1; i1 <= nbptsLocal - 1; ++i1)
                {
                    newPoint = RefCorner + dir0 * double(i0) + dir1 * double(i1);
                    setBndPoint(newPoint, primary_grid, bp, bphz);
                }
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_GeneralPolyhedron::setBndPts_IsometricTetrahedron(vector<geomVector> const &corners_,
                                                             double spacing,
                                                             FV_Mesh const *primary_grid,
                                                             list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp,
                                                             list<DLMFD_BoundaryMultiplierPoint *>::iterator &bphz)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron::setBndPts_IsometricTetrahedron");

    vector<vector<size_t>>::iterator faceIter;
    geomVector RefCorner(3), dir0(3), dir1(3), newPoint(3);

    // Number of points
    faceIter = facesVec.begin();
    size_t nbptsLocal = int(corners_[faceIter->front()].calcDist(corners_[faceIter->back()]) / spacing);

    if (nbptsLocal > 2)
    {
        for (faceIter = facesVec.begin(); faceIter != facesVec.end(); ++faceIter)
        {
            RefCorner = corners_[faceIter->front()];
            dir0 = (corners_[*(faceIter->begin() + 1)] - RefCorner) * (1. / nbptsLocal);
            dir1 = (corners_[faceIter->back()] - RefCorner) * (1. / nbptsLocal);

            for (size_t i0 = 1; i0 <= nbptsLocal - 2; ++i0)
                for (size_t i1 = 1; i1 <= nbptsLocal - (i0 + 1); ++i1)
                {
                    newPoint = RefCorner + dir0 * double(i0) + dir1 * double(i1);
                    setBndPoint(newPoint, primary_grid, bp, bphz);
                }
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_GeneralPolyhedron::update()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron::update");

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
FS_GeneralPolyhedron_Additional_Param const *DLMFD_GeneralPolyhedron::get_ptr_FS_GeneralPolyhedron_Additional_Param()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron::get_ptr_FS_GeneralPolyhedron_Additional_Param");

    return (dynamic_cast<FS_GeneralPolyhedron *>(ptr_FSrigidbody)->get_ptr_FS_GeneralPolyhedron_Additional_Param());
}

//---------------------------------------------------------------------------
void DLMFD_GeneralPolyhedron::allocate_default_listOfPointsAndVectors_GeneralPolyhedron(const double &critical_distance, FV_DiscreteField *pField)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DS_GeneralPolyhedron::allocate_default_listOfPointsAndVectors_GeneralPolyhedron()");

    if (is_particle_fixed)
        initialize_listOfDLMFDPoints();
    else
    {

        double pi = acos(-1.);
        double mesh_size = critical_distance / sqrt(3.);
        double spacing = DLMFD_FictitiousDomain::BoundaryPointsSpacing_coef * critical_distance;

        size_t security_bandwidth = pField->primary_grid()->get_security_bandwidth();
        size_t nb = size_t(pi * radius / spacing);
        size_t n = size_t(2. * radius / mesh_size) + 1;
        size_t nbIPdef = 3. * size_t(pi * n * n * n / 6.);
        size_t nbBPdef = size_t(1.5 * nb * nb);
        size_t nbIPHZdef = security_bandwidth * n * n * 3 * 2;
        size_t nbBPHZdef = size_t(0.25 * nbBPdef);

        allocate_default_listOfPointsAndVectors(nbIPdef, nbBPdef, nbIPHZdef, nbBPHZdef);
    }
}

//---------------------------------------------------------------------------
void DLMFD_GeneralPolyhedron::erase_critical_interior_points_PerProc(double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DS_GeneralPolyhedron::erase_critical_interior_points_PerProc");

    list<DLMFD_BoundaryMultiplierPoint *>::iterator imp, ihp;
    list<DLMFD_InteriorMultiplierPoint *>::iterator ivi;

    // Compute unit normals to faces
    size_t nfaces = facesVec.size();
    geomVector nullvec(3);
    vector<geomVector> face_normals(nfaces, nullvec);
    bool test = true;
    for (size_t i = 0; i < nfaces; ++i)
    {
        face_normals[i] = (corners[facesVec[i][1]] - corners[facesVec[i][0]]) ^
                          (corners[facesVec[i][2]] - corners[facesVec[i][0]]);
        face_normals[i] /= face_normals[i].calcNorm();
    }

    if (!interior_points.empty())
    {
        ivi = interior_points.begin();
        for (size_t j = 0; j < nIP; ++j, ivi++)
            if ((*ivi)->isValid())
            {
                // Test normal distance to faces
                // If normal distance <  critical_distance, test = true
                test = false;
                for (size_t m = 0; m < nfaces && !test; ++m)
                    if (fabs((face_normals[m], (*ivi)->get_coordinates() - corners[facesVec[m][0]])) < critical_distance)
                        test = true;

                if (periodic_directions && !test)
                {
                    size_t nper = periodic_directions->size();
                    for (size_t i = 0; i < nper && !test; ++i)
                        for (size_t m = 0; m < nfaces && !test; ++m)
                            if (fabs((face_normals[m], (*ivi)->get_coordinates() - corners[facesVec[m][0]] - (*periodic_directions)[i])) < critical_distance)
                                test = true;
                }

                // If test, compare to boundary points
                if (test)
                {
                    // Compare to boundary points on this process
                    if (!boundary_points.empty())
                    {
                        imp = boundary_points.begin();
                        for (size_t i = 0; i < nBP && (*ivi)->isValid(); ++i, imp++)
                            if ((*imp)->isValid())
                                if (((*imp)->get_coordinates()).calcDist((*ivi)->get_coordinates()) < critical_distance)
                                    (*ivi)->set_validity(false);
                    }

                    // Compare to boundary points in halo zone
                    if (!halozone_boundary_points.empty() && (*ivi)->isValid())
                    {
                        ihp = halozone_boundary_points.begin();
                        for (size_t i = 0; i < nBPHZ && (*ivi)->isValid(); ++i, ihp++)
                            if ((*ihp)->isValid())
                                if (((*ihp)->get_coordinates()).calcDist((*ivi)->get_coordinates()) < critical_distance)
                                    (*ivi)->set_validity(false);
                    }
                }
            }
    }
}

//---------------------------------------------------------------------------
double DLMFD_GeneralPolyhedron::calcPointDeterm4by4(const geomVector &pointOne,
                                                    const geomVector &pointTwo,
                                                    const geomVector &pointThree,
                                                    const geomVector &pointFour) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron::calcPointDeterm4by4");

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
bool DLMFD_GeneralPolyhedron::checkPointInTetrahedron(const geomVector &pointOne,
                                                      const geomVector &pointTwo,
                                                      const geomVector &pointThree,
                                                      const geomVector &pointFour,
                                                      const geomVector &pointToCheck) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron::checkPointInTetrahedron");

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
bool DLMFD_GeneralPolyhedron::isIn(const geomVector &point) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron:: isIn");

    // return (ptr_FSrigidbody->isIn(point));

    // TODO: there is something to do with this function, see constructor.

    bool inSolid = false;

    bool is_an_obstacle = component_type == "O" || component_type == "PO";

    // The particle is a tetrahedron
    if (corners.size() == 4 && facesVec.size() == 4)
        inSolid = checkPointInTetrahedron(corners[0], corners[1], corners[2], corners[3], point);

    // The particle is a polyhedron with more than 4 faces
    else
    {
        for (vector<vector<size_t>>::const_iterator faceIter = facesVec.begin();
             faceIter != facesVec.end() && !inSolid; ++faceIter)
        {
            size_t noPntsFace = faceIter->size();
            if (noPntsFace == 3) // face is triangle
            {
                if (checkPointInTetrahedron(corners[(*faceIter)[0]],
                                            corners[(*faceIter)[1]],
                                            corners[(*faceIter)[2]],
                                            is_an_obstacle ? *g2 : gravity_center,
                                            point))
                    inSolid = true;
            }
            else
            {
                for (size_t idxPts = 2; idxPts < noPntsFace && !inSolid; ++idxPts)
                    if (checkPointInTetrahedron(corners[(*faceIter)[0]],
                                                corners[(*faceIter)[idxPts - 1]],
                                                corners[(*faceIter)[idxPts]],
                                                is_an_obstacle ? *g2 : gravity_center,
                                                point))
                        inSolid = true;
            }
        }
    }

    return inSolid;
}

//---------------------------------------------------------------------------
void DLMFD_GeneralPolyhedron::translateGeometricFeatures(geomVector const &newg)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeneralPolyhedron:: translateGeometricFeatures");

    if (gravity_center.calcDist(newg) > 1.e-2 * radius)
    {
        geomVector translation_vec = newg - gravity_center;
        gravity_center = newg;

        if (g2)
            *g2 += translation_vec;

        size_t j, ncorners = corners.size();
        for (j = 0; j < ncorners; ++j)
            corners[j] += translation_vec;
    }
}
