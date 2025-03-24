#include <DLMFD_RigidBody.hh>
#include <DLMFD_FictitiousDomain.hh>
#include <FV_Mesh.hh>
#include <FS_RigidBody.hh>
#include <geomVector.hh>
#include <doubleVector.hh>
#include <math.h>
#include <fstream>
using namespace std;

size_t_array3D *DLMFD_RigidBody::Q2numb = NULL;

//---------------------------------------------------------------------------
DLMFD_RigidBody::DLMFD_RigidBody() : VEC_r(*DLMFD_FictitiousDomain::dbnull),
                                     VEC_x(*DLMFD_FictitiousDomain::dbnull),
                                     VEC_lambda(*DLMFD_FictitiousDomain::dbnull),
                                     VEC_w(*DLMFD_FictitiousDomain::dbnull)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: DLMFD_RigidBody");
}

//---------------------------------------------------------------------------
DLMFD_RigidBody::DLMFD_RigidBody(FS_RigidBody *pgrb,
                                 const bool &are_particles_fixed,
                                 FV_DiscreteField *pField_) : ptr_FSrigidbody(pgrb),
                                                              pField_(pField_),
                                                              VEC_r(*DLMFD_FictitiousDomain::dbnull),
                                                              VEC_x(*DLMFD_FictitiousDomain::dbnull),
                                                              VEC_lambda(*DLMFD_FictitiousDomain::dbnull),
                                                              VEC_w(*DLMFD_FictitiousDomain::dbnull),
                                                              periodic_directions(NULL),
                                                              ndof(0)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: DLMFD_RigidBody");

    dim = 3;
    index_min = new size_t_array2D(dim, dim, 0);
    index_max = new size_t_array2D(dim, dim, 0);

    NDOF_comp = NULL;
    NDOF_leverage = NULL;
    NDOF_nfieldUNK = NULL;
    NDOF_globalpos = NULL;
    NDOF_deltaOmega = NULL;
    NDOF_FVTriplet = NULL;

    fluidsolid_coupling_factor = 1.;

    // Resize vectors
    translational_velocity.resize(3);
    angular_velocity_3D.resize(3);

    gravity_center.resize(3);

    t_tran.resize(3);
    q_tran.resize(3);
    t_rot_3D.resize(3);
    q_rot_3D.resize(3);

    // Set the finite element counting
    if (Q2numb == NULL)
    {
        Q2numb = new size_t_array3D(3, 3, 3);
        (*Q2numb)(0, 0, 0) = 0;
        (*Q2numb)(1, 0, 0) = 1;
        (*Q2numb)(2, 0, 0) = 2;
        (*Q2numb)(0, 1, 0) = 5;
        (*Q2numb)(1, 1, 0) = 4;
        (*Q2numb)(2, 1, 0) = 3;
        (*Q2numb)(0, 2, 0) = 6;
        (*Q2numb)(1, 2, 0) = 7;
        (*Q2numb)(2, 2, 0) = 8;

        (*Q2numb)(0, 0, 1) = 9;
        (*Q2numb)(1, 0, 1) = 10;
        (*Q2numb)(2, 0, 1) = 11;
        (*Q2numb)(0, 1, 1) = 14;
        (*Q2numb)(1, 1, 1) = 13;
        (*Q2numb)(2, 1, 1) = 12;
        (*Q2numb)(0, 2, 1) = 15;
        (*Q2numb)(1, 2, 1) = 16;
        (*Q2numb)(2, 2, 1) = 17;

        (*Q2numb)(0, 0, 2) = 18;
        (*Q2numb)(1, 0, 2) = 19;
        (*Q2numb)(2, 0, 2) = 20;
        (*Q2numb)(0, 1, 2) = 23;
        (*Q2numb)(1, 1, 2) = 22;
        (*Q2numb)(2, 1, 2) = 21;
        (*Q2numb)(0, 2, 2) = 24;
        (*Q2numb)(1, 2, 2) = 25;
        (*Q2numb)(2, 2, 2) = 26;
    }

    // Set component type
    set_component_type();

    // Set the translational velocity
    set_translational_velocity();

    // Set the angular velocity
    set_angular_velocity_3D();
}

//---------------------------------------------------------------------------
DLMFD_RigidBody::~DLMFD_RigidBody()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: ~DLMFD_RigidBody");

    if (Q2numb)
    {
        delete Q2numb;
        Q2numb = NULL;
    }

    list<DLMFD_InteriorMultiplierPoint *>::iterator iip;
    list<DLMFD_BoundaryMultiplierPoint *>::iterator ibp;

    for (iip = interior_points.begin(); iip != interior_points.end(); iip++)
        delete *iip;
    interior_points.clear();

    for (ibp = boundary_points.begin(); ibp != boundary_points.end(); ibp++)
        delete *ibp;
    boundary_points.clear();

    list<struct ULBD_RHSInfos>::iterator il;
    for (il = points_infos.begin(); il != points_infos.end(); il++)
    {
        il->NodeNoU.clear();
        il->positionU.clear();
        il->omega_delta.clear();
    }
    points_infos.clear();

    for (iip = halozone_interior_points.begin(); iip != halozone_interior_points.end();
         iip++)
        delete *iip;
    halozone_interior_points.clear();

    for (ibp = halozone_boundary_points.begin(); ibp != halozone_boundary_points.end();
         ibp++)
        delete *ibp;
    halozone_boundary_points.clear();

    if (NDOF_comp)
        delete[] NDOF_comp;
    if (NDOF_leverage)
        delete[] NDOF_leverage;
    if (NDOF_nfieldUNK)
        delete[] NDOF_nfieldUNK;

    if (NDOF_globalpos)
        delete[] NDOF_globalpos;
    if (NDOF_deltaOmega)
        delete[] NDOF_deltaOmega;
    if (NDOF_FVTriplet)
        delete[] NDOF_FVTriplet;

    if (periodic_directions)
    {
        periodic_directions->clear();
        delete periodic_directions;
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_component_type()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_component_type");

    component_type = ptr_FSrigidbody->get_type();
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_ptr_constrained_field(FV_DiscreteField *pField__)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::set_ptr_constrained_field");

    pField_ = pField__;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::initialize_listOfDLMFDPoints()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::initialize_listOfDLMFDPoints");

    DLMFD_InteriorMultiplierPoint *imp = NULL;
    DLMFD_BoundaryMultiplierPoint *bmp = NULL;

    geomVector virtual_point = geomVector(0., 0., 0.);

    if (interior_points.empty())
    {
        interior_points.push_back(imp);
        interior_points.back() = new DLMFD_InteriorMultiplierPoint(0, virtual_point, 0, 0, 0, gravity_center);
    }

    if (halozone_interior_points.empty())
    {
        halozone_interior_points.push_back(imp);
        halozone_interior_points.back() = new DLMFD_InteriorMultiplierPoint(0, virtual_point, 0, 0, 0, gravity_center);
    }

    if (boundary_points.empty())
    {
        boundary_points.push_back(bmp);
        boundary_points.back() = new DLMFD_BoundaryMultiplierPoint(virtual_point, gravity_center);
    }

    if (halozone_boundary_points.empty())
    {
        halozone_boundary_points.push_back(bmp);
        halozone_boundary_points.back() = new DLMFD_BoundaryMultiplierPoint(virtual_point, gravity_center);
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_boundary_point(const geomVector &point, list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_boundary_point");

    (*bp)->set(0, point, gravity_center);
    ++nBP;
    if (nBP == boundary_points.size())
        extend_bp_list(nBP);
    bp++;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_halozone_boundary_point(const geomVector &point, list<DLMFD_BoundaryMultiplierPoint *>::iterator &bphz)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::set_halozone_boundary_point");

    (*bphz)->set(0, point, gravity_center);
    ++nBPHZ;
    if (nBPHZ == halozone_boundary_points.size())
        extend_bphz_list(nBPHZ);
    bphz++;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_interior_point(const size_t &comp,
                                         const geomVector &point,
                                         size_t i, size_t j, size_t k,
                                         list<DLMFD_InteriorMultiplierPoint *>::iterator &ip)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_interior_point");

    (*ip)->set(comp, point, i, j, k, gravity_center);
    ++nIP;
    if (nIP == interior_points.size())
        extend_ip_list(nIP);
    ip++;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_halozone_interior_point(const size_t &comp,
                                                  const geomVector &point,
                                                  size_t i, size_t j, size_t k,
                                                  list<DLMFD_InteriorMultiplierPoint *>::iterator &iphz)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_interior_point");

    (*iphz)->set(comp, point, i, j, k, gravity_center);
    ++nIPHZ;
    if (nIPHZ == halozone_interior_points.size())
        extend_iphz_list(nIPHZ);
    iphz++;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::setBndPoint(geomVector const &point,
                                  FV_Mesh const *primary_grid,
                                  list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp,
                                  list<DLMFD_BoundaryMultiplierPoint *>::iterator &bphz)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::setBndPoint");

    if (primary_grid->is_in_domain_with_halozone(point(0), point(1), point(2)))
    {
        if (primary_grid->is_in_domain_on_current_processor(point(0), point(1), point(2)))
        {
            set_boundary_point(point, bp);
        }
        else
        {
            set_halozone_boundary_point(point, bphz);
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::setPtsOnEdge(const geomVector &firstPoint,
                                   const geomVector &secondPoint,
                                   const double &critical_distance, const bool &addBeginPnt,
                                   FV_Mesh const *primary_grid,
                                   list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp,
                                   list<DLMFD_BoundaryMultiplierPoint *>::iterator &bphz)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::setPtsOnEdge");

    double tol = 0.01;
    size_t nbptsLocal = size_t(firstPoint.calcDist(secondPoint) / critical_distance + tol);

    geomVector point;

    if (nbptsLocal)
    {
        double dx = (firstPoint(0) - secondPoint(0)) / nbptsLocal,
               dy = (firstPoint(1) - secondPoint(1)) / nbptsLocal,
               dz = (firstPoint(2) - secondPoint(2)) / nbptsLocal, x, y, z;
        size_t beginIdx = addBeginPnt ? 0 : 1;

        for (size_t idxSub = beginIdx; idxSub <= nbptsLocal - 1; ++idxSub)
        {
            x = secondPoint(0) + idxSub * dx;
            y = secondPoint(1) + idxSub * dy;
            z = secondPoint(2) + idxSub * dz;
            point = geomVector(x, y, z);

            setBndPoint(point, primary_grid, bp, bphz);
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::erase_critical_interior_points_PerProc(double critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::erase_critical_interior_points_PerProc");

    list<DLMFD_BoundaryMultiplierPoint *>::iterator imp, ihp;
    list<DLMFD_InteriorMultiplierPoint *>::iterator ivi;

    if (!interior_points.empty())
    {
        // Compare to boundary points on this process
        if (!boundary_points.empty())
        {
            imp = boundary_points.begin();
            for (size_t i = 0; i < nBP; ++i, imp++)
            {
                if ((*imp)->isValid())
                {
                    ivi = interior_points.begin();
                    for (size_t j = 0; j < nIP; ++j, ivi++)
                        if ((*ivi)->isValid())
                        {
                            if (((*imp)->get_coordinates()).calcDist((*ivi)->get_coordinates()) < critical_distance)
                            {
                                (*ivi)->set_validity(false);
                            }
                        }
                }
            }
        }

        // Compare to boundary points in halo zone
        if (!halozone_boundary_points.empty())
        {
            ihp = halozone_boundary_points.begin();
            for (size_t i = 0; i < nBPHZ; ++i, ihp++)
                if ((*ihp)->isValid())
                {
                    ivi = interior_points.begin();
                    for (size_t j = 0; j < nIP; ++j, ivi++)
                        if ((*ivi)->isValid())
                            if (((*ihp)->get_coordinates()).calcDist((*ivi)->get_coordinates()) < critical_distance)
                                (*ivi)->set_validity(false);
                }
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::erase_critical_boundary_points_ptp_PerProc(DLMFD_RigidBody *second_component, const double &critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::erase_critical_boundary_points_ptp_PerProc");

    size_t nclones2 = second_component->get_number_periodicClones();

    erase_critical_boundary_points_ptp_PerProc_oneGC(second_component, critical_distance);

    if (periodic_directions)
    {
        geomVector gcref = gravity_center;
        for (size_t i = 0; i < periodic_directions->size(); ++i)
        {
            // gravity_center = gcref + (*periodic_directions)[i];
            translateGeometricFeatures(gcref + (*periodic_directions)[i]);
            erase_critical_boundary_points_ptp_PerProc_oneGC(second_component, critical_distance);
        }
        //    gravity_center = gcref ;
        translateGeometricFeatures(gcref);
    }

    if (nclones2)
    {
        geomVector gcref2 = *(second_component->get_ptr_to_gravity_centre());
        vector<geomVector> *pd2 = second_component->get_periodicClones_DirectionsVector();
        for (size_t j = 0; j < nclones2; ++j)
        {
            second_component->translateGeometricFeatures(gcref2 + (*pd2)[j]);
            erase_critical_boundary_points_ptp_PerProc_oneGC(second_component, critical_distance);
        }
        second_component->translateGeometricFeatures(gcref2);

        if (periodic_directions)
        {
            geomVector gcref = gravity_center;
            for (size_t i = 0; i < periodic_directions->size(); ++i)
            {
                // gravity_center = gcref + (*periodic_directions)[i];
                translateGeometricFeatures(gcref + (*periodic_directions)[i]);
                for (size_t j = 0; j < nclones2; ++j)
                {
                    second_component->translateGeometricFeatures(gcref2 + (*pd2)[j]);
                    erase_critical_boundary_points_ptp_PerProc_oneGC(second_component,
                                                                     critical_distance);
                    // second_component->set_gravityCenter( gcref2 ) ;
                }
            }
            // gravity_center = gcref ;
            translateGeometricFeatures(gcref);
            second_component->translateGeometricFeatures(gcref2);
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::erase_critical_boundary_points_ptp_PerProc_oneGC(DLMFD_RigidBody *second_component, const double &critical_distance)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::erase_critical_boundary_points_ptp_PerProc_oneGC");

    // TO CLEAN : see the issue with pointer/object for isIn etc

    if (proximityQuery(second_component, critical_distance))
    {
        // Determine points in contact region
        list<DLMFD_BoundaryMultiplierPoint *> BP_contactRegion,
            BPHZ_contactRegion,
            C2_BP_contactRegion,
            C2_BPHZ_contactRegion;

        DLMFDPoints_in_ContactRegion(BP_contactRegion, BPHZ_contactRegion,
                                     *second_component->get_ptr_to_gravity_centre(),
                                     second_component->get_circumscribed_radius() + critical_distance);

        second_component->DLMFDPoints_in_ContactRegion(C2_BP_contactRegion,
                                                       C2_BPHZ_contactRegion,
                                                       gravity_center,
                                                       radius + critical_distance);

        // Procedure to tag boundary points (BP)
        // 1. find BP of the former component that lie inside the latter and
        // vice-versa
        // 2. compare BP & HP of the 2 particles
        list<DLMFD_BoundaryMultiplierPoint *>::iterator imp, imp2, ihp, ihp2;
        double dist = 0.;

        // Find BP & HP of the former component that lie inside the latter and
        // vice-versa
        for (imp = BP_contactRegion.begin(); imp != BP_contactRegion.end(); imp++)
            if ((*imp)->isValid())
                if (second_component->isIn((*imp)->get_coordinates()))
                    (*imp)->set_validity(false);

        for (ihp = BPHZ_contactRegion.begin(); ihp != BPHZ_contactRegion.end(); ihp++)
            if ((*ihp)->isValid())
                if (second_component->isIn((*ihp)->get_coordinates()))
                    (*ihp)->set_validity(false);

        for (imp2 = C2_BP_contactRegion.begin(); imp2 != C2_BP_contactRegion.end();
             imp2++)
            if ((*imp2)->isValid())
                if (isIn((*imp2)->get_coordinates()))
                    (*imp2)->set_validity(false);

        for (ihp2 = C2_BPHZ_contactRegion.begin(); ihp2 != C2_BPHZ_contactRegion.end();
             ihp2++)
            if ((*ihp2)->isValid())
                if (isIn((*ihp2)->get_coordinates()))
                    (*ihp2)->set_validity(false);

        // Compare BP & HP of the 2 particles
        for (imp = BP_contactRegion.begin(); imp != BP_contactRegion.end(); imp++)
            for (imp2 = C2_BP_contactRegion.begin();
                 imp2 != C2_BP_contactRegion.end(); imp2++)
            {
                dist = ((*imp2)->get_coordinates()).calcDist((*imp)->get_coordinates());
                if (dist < critical_distance)
                {
                    (*imp)->set_validity(false);
                    (*imp2)->set_validity(false);
                }
            }

        for (ihp = BPHZ_contactRegion.begin(); ihp != BPHZ_contactRegion.end(); ihp++)
            for (imp2 = C2_BP_contactRegion.begin();
                 imp2 != C2_BP_contactRegion.end(); imp2++)
            {
                dist = ((*imp2)->get_coordinates()).calcDist((*ihp)->get_coordinates());
                if (dist < critical_distance)
                {
                    (*ihp)->set_validity(false);
                    (*imp2)->set_validity(false);
                }
            }

        for (ihp2 = C2_BPHZ_contactRegion.begin(); ihp2 != C2_BPHZ_contactRegion.end();
             ihp2++)
            for (imp = BP_contactRegion.begin(); imp != BP_contactRegion.end(); imp++)
            {
                dist = ((*imp)->get_coordinates()).calcDist((*ihp2)->get_coordinates());
                if (dist < critical_distance)
                {
                    (*ihp2)->set_validity(false);
                    (*imp)->set_validity(false);
                }
            }

        for (ihp = BPHZ_contactRegion.begin(); ihp != BPHZ_contactRegion.end(); ihp++)
            for (ihp2 = C2_BPHZ_contactRegion.begin();
                 ihp2 != C2_BPHZ_contactRegion.end(); ihp2++)
            {
                dist = ((*ihp)->get_coordinates()).calcDist((*ihp2)->get_coordinates());
                if (dist < critical_distance)
                {
                    (*ihp)->set_validity(false);
                    (*ihp2)->set_validity(false);
                }
            }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_points_infos(FV_DiscreteField *pField)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::set_points_infos");

    list<struct ULBD_RHSInfos>::iterator pinfo = points_infos.begin();
    size_t globalDOFNo, compIdx;
    FV_TRIPLET const *ijk;
    FV_TRIPLET ijk_init;
    ijk_init.i = 0;
    ijk_init.j = 0;
    ijk_init.k = 0;
    ndof = 0;
    ntotal_fieldunk = 0;

    // Interior points
    if (nIP)
    {
        list<DLMFD_InteriorMultiplierPoint *>::const_iterator iterIntPnts = interior_points.begin();
        for (size_t i = 0; i < nIP; ++i, iterIntPnts++)
            if ((*iterIntPnts)->isValid())
            {
                // Pointer of the internal point
                pinfo->ptr_point = (const DLMFD_ParticlePoint *)(*iterIntPnts);

                // Set the field component
                compIdx = (*iterIntPnts)->get_compNumber();

                pinfo->compIdx = compIdx;

                // Global DOF number
                // A priori always >= 0 i.e. DOF is unknown because DOFs that are not
                // unknowns have already been discarded in the set_all_points method
                // of the component
                ijk = (*iterIntPnts)->get_localNodeTriplet();
                globalDOFNo = pField->DOF_global_number(ijk->i, ijk->j, ijk->k, compIdx);

                // Set infos
                if (pinfo->positionU.size() != 1)
                {
                    pinfo->MacTripletU.clear();
                    pinfo->positionU.clear();
                    pinfo->omega_delta.clear();
                    pinfo->MacTripletU.push_back(ijk_init);
                    pinfo->positionU.push_back(0);
                    pinfo->omega_delta.push_back(0.);
                }
                pinfo->positionU.front() = globalDOFNo;
                pinfo->MacTripletU.front() = *ijk;
                pinfo->omega_delta.front() = 1.;

                // Move the iterator and update ndof & ntotal_fieldunk
                ++ndof;
                ++ntotal_fieldunk;
                pinfo++;
            }
    }

    // Boundary points
    if (nBP)
    {
        list<DLMFD_BoundaryMultiplierPoint *>::const_iterator iterBndPnts = boundary_points.begin();
        size_t ncomps = pField->nb_components();
        for (size_t i = 0; i < nBP; ++i, iterBndPnts++)
            if ((*iterBndPnts)->isValid())
                for (compIdx = 0; compIdx < ncomps; ++compIdx)
                // Loop on all components
                {
                    // get the mesh ID of the element containing the boundary point
                    pinfo->ptr_point = (const DLMFD_ParticlePoint *)(*iterBndPnts);

                    // Set the field component
                    pinfo->compIdx = compIdx;

                    // Clear the lists of OnePointInfos
                    pinfo->MacTripletU.clear();
                    pinfo->positionU.clear();
                    pinfo->omega_delta.clear();

                    // Get surrouding dofs
                    fill_DLMFD_pointInfos(pField, *pinfo, DLMFD_FictitiousDomain::b_SecondOrderInterpol);

                    // Move the iterator and update ndof & ntotal_fieldunk
                    ntotal_fieldunk += pinfo->positionU.size();

                    ++ndof;
                    pinfo++;
                }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::check_allocation_DLMFD_Cvectors()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::check_allocation_DLMFD_Cvectors");

    if (ndof)
    {
        if (NDOF_comp)
            delete[] NDOF_comp;
        if (NDOF_leverage)
            delete[] NDOF_leverage;
        if (NDOF_nfieldUNK)
            delete[] NDOF_nfieldUNK;

        NDOF_comp = new size_t[ndof];
        NDOF_leverage = new double[(gravity_center.getVecSize() - 1) * ndof];
        NDOF_nfieldUNK = new size_t[ndof];

        if (NDOF_globalpos)
            delete[] NDOF_globalpos;
        if (NDOF_deltaOmega)
            delete[] NDOF_deltaOmega;
        if (NDOF_FVTriplet)
            delete[] NDOF_FVTriplet;

        NDOF_globalpos = new size_t[ntotal_fieldunk];
        NDOF_deltaOmega = new double[ntotal_fieldunk];
        NDOF_FVTriplet = new size_t[3 * ntotal_fieldunk];
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::fill_DLMFD_Cvectors()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::fill_DLMFD_Cvectors");

    if (ndof)
    {
        list<struct ULBD_RHSInfos>::iterator il = points_infos.begin();
        size_t comp = 0;
        for (size_t i = 0; i < ndof; ++i, il++)
        {
            comp = il->compIdx;
            NDOF_comp[i] = comp;
            switch (comp)
            {
            case 0:
                NDOF_leverage[2 * i] =
                    il->ptr_point->get_oneCoordinate_GCPointVector(2);
                NDOF_leverage[2 * i + 1] =
                    il->ptr_point->get_oneCoordinate_GCPointVector(1);
                break;
            case 1:
                NDOF_leverage[2 * i] =
                    il->ptr_point->get_oneCoordinate_GCPointVector(2);
                NDOF_leverage[2 * i + 1] =
                    il->ptr_point->get_oneCoordinate_GCPointVector(0);
                break;
            default:
                NDOF_leverage[2 * i] =
                    il->ptr_point->get_oneCoordinate_GCPointVector(1);
                NDOF_leverage[2 * i + 1] =
                    il->ptr_point->get_oneCoordinate_GCPointVector(0);
                break;
            }

            NDOF_nfieldUNK[i] = il->positionU.size();
        }

        // Same allocation strategy for ntotal_fieldunk
        list<size_t>::const_iterator ilpositionU;
        list<double>::const_iterator ildom;
        list<FV_TRIPLET>::iterator ilFVTripletU;
        il = points_infos.begin();
        size_t l = 0, m = 0;
        for (size_t i = 0; i < ndof; ++i, il++)
        {
            ildom = il->omega_delta.begin();
            ilFVTripletU = il->MacTripletU.begin();
            for (ilpositionU = il->positionU.begin(); ilpositionU != il->positionU.end();
                 ilpositionU++, ildom++, ilFVTripletU++)
            {
                NDOF_globalpos[l] = *ilpositionU;
                NDOF_deltaOmega[l] = *ildom;
                NDOF_FVTriplet[m] = ilFVTripletU->i;
                NDOF_FVTriplet[m + 1] = ilFVTripletU->j;
                NDOF_FVTriplet[m + 2] = ilFVTripletU->k;
                ++l;
                m += 3;
            }
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_coupling_factor(double const &rho_f, bool const &explicit_treatment)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_coupling_factor");

    fluidsolid_coupling_factor = 1.;
    if (!explicit_treatment)
        fluidsolid_coupling_factor -= rho_f / rho_s;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_mass_and_density_and_inertia()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_mass_and_density_and_inertia");

    mass = get_rigid_body_mass();
    rho_s = get_rigid_body_density();
    inertia_3D = get_rigid_body_inertia();
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_volume()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::set_volume");

    volume = mass / rho_s;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_translational_velocity()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_translational_velocity");

    translational_velocity = get_rigid_body_translational_velocity();
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_translational_velocity(const geomVector &vtran)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_translational_velocity");

    translational_velocity = vtran;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_angular_velocity_3D()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_translational_velocity");

    angular_velocity_3D = get_rigid_body_angular_velocity();
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_angular_velocity_3D(const geomVector &vrot)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_translational_velocity");

    angular_velocity_3D = vrot;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_Qu(const geomVector &qtran)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_Qu");

    q_tran = qtran;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_Qrot(const geomVector &qrot)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_Qrot");

    q_rot_3D = qrot;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_Tu(const geomVector &ttran)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_Tu");

    t_tran = ttran;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_Trot(const geomVector &trot)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: set_Trot");

    t_rot_3D = trot;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::set_ttran_ncomp(const size_t &ncomp_)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::set_ttran_ncomp");

    t_tran.resize(ncomp_);
}

//---------------------------------------------------------------------------
int DLMFD_RigidBody::get_number_periodicClones() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_number_periodicClones");

    return (periodic_directions ? periodic_directions->size() : 0);
}

//---------------------------------------------------------------------------
vector<geomVector> *DLMFD_RigidBody::get_periodicClones_DirectionsVector() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_periodicClones_DirectionsVector");

    return periodic_directions;
}

//---------------------------------------------------------------------------
string DLMFD_RigidBody::get_component_type() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: get_component_type");

    return component_type;
}

//---------------------------------------------------------------------------
GEOMETRICSHAPE DLMFD_RigidBody::get_GeometricObjectType() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_GeometricObjectType");

    return (ptr_FSrigidbody->get_shape_type());
}

//---------------------------------------------------------------------------
double DLMFD_RigidBody::get_circumscribed_radius() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: get_circumscribed_radius");

    return (ptr_FSrigidbody->get_circumscribed_radius());
}

//---------------------------------------------------------------------------
geomVector const *DLMFD_RigidBody::get_ptr_to_gravity_centre() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: get_ptr_to_gravity_centre");

    return (dynamic_cast<FS_RigidBody *>(ptr_FSrigidbody)->get_ptr_to_gravity_centre());
}

//---------------------------------------------------------------------------
size_t DLMFD_RigidBody::get_npts_output(bool const &withIntPts)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: get_npts_output");

    list<DLMFD_BoundaryMultiplierPoint *>::const_iterator imp = boundary_points.begin();
    list<DLMFD_InteriorMultiplierPoint *>::const_iterator ivi = interior_points.begin();
    size_t output_npts = 0;

    for (size_t i = 0; i < nBP; ++i, imp++)
        if ((*imp)->isValid())
            ++output_npts;

    if (withIntPts)
        for (size_t i = 0; i < nIP; ++i, ivi++)
            if ((*ivi)->isValid())
                ++output_npts;

    return output_npts;
}

//---------------------------------------------------------------------------
geomVector DLMFD_RigidBody::get_rigid_body_velocity(geomVector const &point) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: get_rigid_body_velocity");

    return (ptr_FSrigidbody->rigid_body_velocity(point));
}

//---------------------------------------------------------------------------
geomVector DLMFD_RigidBody::get_rigid_body_angular_velocity() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: get_rigid_body_angular_velocity");

    return (ptr_FSrigidbody->rigid_body_angular_velocity());
}

//---------------------------------------------------------------------------
geomVector DLMFD_RigidBody::get_angular_velocity_3D() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_angular_velocity_3D");

    return angular_velocity_3D;
}

//---------------------------------------------------------------------------
geomVector DLMFD_RigidBody::get_rigid_body_translational_velocity() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_rigid_body_translational_velocity");

    return (ptr_FSrigidbody->rigid_body_translational_velocity());
}

//---------------------------------------------------------------------------
geomVector DLMFD_RigidBody::get_translational_velocity() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_translational_velocity");

    return translational_velocity;
}

//---------------------------------------------------------------------------
double DLMFD_RigidBody::get_rigid_body_mass() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_mass");

    double mass_ = get<0>(ptr_FSrigidbody->get_mass_and_density_and_moi());

    return mass_;
}

//---------------------------------------------------------------------------
double DLMFD_RigidBody::get_rigid_body_density() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_rigid_body_density");

    double density_ = get<1>(ptr_FSrigidbody->get_mass_and_density_and_moi());

    return density_;
}

//---------------------------------------------------------------------------
double DLMFD_RigidBody::get_density() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_density");

    return rho_s;
}

//---------------------------------------------------------------------------
vector<vector<double>> DLMFD_RigidBody::get_rigid_body_inertia() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_mass_and_density_and_inertia");

    return ptr_FSrigidbody->get_inertia();
}

//---------------------------------------------------------------------------
geomVector const DLMFD_RigidBody::get_Qu() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_Qu");

    return q_tran;
}

//---------------------------------------------------------------------------
geomVector const DLMFD_RigidBody::get_Qrot() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_Qrot");

    return q_rot_3D;
}

//---------------------------------------------------------------------------
geomVector const DLMFD_RigidBody::get_Tu() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_Tu");

    return t_tran;
}

//---------------------------------------------------------------------------
geomVector const DLMFD_RigidBody::get_Trot() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_Trot");

    return t_rot_3D;
}

//---------------------------------------------------------------------------
double DLMFD_RigidBody::get_volume() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::get_volume");

    return volume;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::extend_bp_list(size_t const &np)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::extent_bp_list");

    DLMFD_BoundaryMultiplierPoint *bmp = NULL;
    geomVector virtual_point = geomVector(0., 0., 0.);

    for (size_t i = 0; i < np; ++i)
    {
        boundary_points.push_back(bmp);
        boundary_points.back() =
            new DLMFD_BoundaryMultiplierPoint(virtual_point, gravity_center);
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::extend_bphz_list(size_t const &np)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::extend_bphz_list");

    DLMFD_BoundaryMultiplierPoint *bmphz = NULL;
    geomVector virtual_point = geomVector(0., 0., 0.);

    for (size_t i = 0; i < np; ++i)
    {
        halozone_boundary_points.push_back(bmphz);
        halozone_boundary_points.back() = new DLMFD_BoundaryMultiplierPoint(virtual_point, gravity_center);
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::extend_ip_list(size_t const &np)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::extent_ip_list");

    DLMFD_InteriorMultiplierPoint *imp = NULL;
    geomVector virtual_point = geomVector(0., 0., 0.);

    for (size_t i = 0; i < np; ++i)
    {
        interior_points.push_back(imp);
        interior_points.back() = new DLMFD_InteriorMultiplierPoint(0, virtual_point, 0, 0, 0, gravity_center);
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::extend_iphz_list(size_t const &np)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::extend_iphz_list");

    DLMFD_InteriorMultiplierPoint *imphz = NULL;
    geomVector virtual_point = geomVector(0., 0., 0.);

    for (size_t i = 0; i < np; ++i)
    {
        halozone_interior_points.push_back(imphz);
        halozone_interior_points.back() = new DLMFD_InteriorMultiplierPoint(0, virtual_point, 0, 0, 0, gravity_center);
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::allocate_default_listOfPointsAndVectors(size_t const &nbIPdef,
                                                              size_t const &nbBPdef,
                                                              size_t const &nbIPHZdef,
                                                              size_t const &nbBPHZdef)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_Sphere:: allocate_default_listOfPointsAndVectors");

    DLMFD_InteriorMultiplierPoint *imp = NULL;
    DLMFD_BoundaryMultiplierPoint *bmp = NULL;

    geomVector virtual_point = geomVector(0., 0., 0.);

    for (size_t i = 0; i < nbIPdef; ++i)
    {
        interior_points.push_back(imp);
        interior_points.back() = new DLMFD_InteriorMultiplierPoint(0, virtual_point, 0, 0, 0, gravity_center);
    }

    for (size_t i = 0; i < nbIPHZdef; ++i)
    {
        halozone_interior_points.push_back(imp);
        halozone_interior_points.back() = new DLMFD_InteriorMultiplierPoint(0, virtual_point, 0, 0, 0, gravity_center);
    }

    for (size_t i = 0; i < nbBPdef; ++i)
    {
        boundary_points.push_back(bmp);
        boundary_points.back() = new DLMFD_BoundaryMultiplierPoint(virtual_point, gravity_center);
    }

    for (size_t i = 0; i < nbBPHZdef; ++i)
    {
        halozone_boundary_points.push_back(bmp);
        halozone_boundary_points.back() = new DLMFD_BoundaryMultiplierPoint(virtual_point, gravity_center);
    }

    size_t ndofdef = nbIPdef + nbBPdef * gravity_center.getVecSize();
    struct ULBD_RHSInfos OnePointInfos;
    OnePointInfos.ptr_point = NULL;
    OnePointInfos.compIdx = 0;
    for (size_t i = 0; i < ndofdef; ++i)
        points_infos.push_back(OnePointInfos);

    VEC_r.re_initialize(ndofdef, 0.);
    VEC_x.re_initialize(ndofdef, 0.);
    VEC_lambda.re_initialize(ndofdef, 0.);
    VEC_w.re_initialize(ndofdef, 0.);
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::update()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: update");

    gravity_center = *get_ptr_to_gravity_centre();
    radius = get_circumscribed_radius();
    translational_velocity = get_rigid_body_translational_velocity();
    angular_velocity_3D = get_rigid_body_angular_velocity();
    inertia_3D = get_rigid_body_inertia();
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::update_RB_position_and_velocity(geomVector const &pos,
                                                      geomVector const &vel,
                                                      geomVector const &ang_vel,
                                                      vector<geomVector> const &periodic_directions, double const &time_step)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: update_RB_position_and_velocity");

    return (ptr_FSrigidbody->update_RB_position_and_velocity(pos, vel, ang_vel, periodic_directions, time_step));
}

//---------------------------------------------------------------------------
bool DLMFD_RigidBody::proximityQuery(DLMFD_RigidBody const *second_component,
                                     const double &distance) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::proximityQuery");

    bool close = false;

    switch (second_component->get_GeometricObjectType())
    {
    case GEOM_3DBOX:
        close = second_component->proximityQuery(this, distance);
        break;

    default:
        if (gravity_center.calcDist(*second_component->get_ptr_to_gravity_centre()) <= radius + second_component->get_circumscribed_radius() + distance)
            close = true;
    }

    return close;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::DLMFDPoints_in_ContactRegion(list<DLMFD_BoundaryMultiplierPoint *> &BP_contactRegion,
                                                   list<DLMFD_BoundaryMultiplierPoint *> &BPHZ_contactRegion,
                                                   geomVector const &refPoint,
                                                   double const &distance_contactRegion)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::DLMFDPoints_in_ContactRegion");

    list<DLMFD_BoundaryMultiplierPoint *>::iterator imp;

    // Boundary points in the contact region
    imp = boundary_points.begin();
    for (size_t i = 0; i < nBP; ++i, imp++)
        if (((*imp)->get_coordinates()).calcDist(refPoint) < distance_contactRegion)
            BP_contactRegion.push_back(*imp);

    // Boundary points in halo zone in the contact region
    imp = halozone_boundary_points.begin();
    for (size_t i = 0; i < nBPHZ; ++i, imp++)
        if (((*imp)->get_coordinates()).calcDist(refPoint) < distance_contactRegion)
            BPHZ_contactRegion.push_back(*imp);
}

//---------------------------------------------------------------------------
bool DLMFD_RigidBody::hasDLMFDPointsOnProc() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::hasDLMFDPointsOnProc");

    bool has_valid_point = false;

    // Interior points
    list<DLMFD_InteriorMultiplierPoint *>::const_iterator ivi = interior_points.begin();

    for (size_t i = 0; i < nIP && !has_valid_point; ++i, ivi++)
        if ((*ivi)->isValid())
            has_valid_point = true;

    // Boundary points
    if (!has_valid_point)
    {
        list<DLMFD_BoundaryMultiplierPoint *>::const_iterator imp = boundary_points.begin();

        for (size_t i = 0; i < nBP && !has_valid_point; ++i, imp++)
            if ((*imp)->isValid())
                has_valid_point = true;
    }

    return has_valid_point;
}

//---------------------------------------------------------------------------
bool DLMFD_RigidBody::hasDLMFDPoints_inHalozone_OnProc() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::hasDLMFDPointsOnProc");

    bool has_halozone_point = false;

    // Interior points
    list<DLMFD_InteriorMultiplierPoint *>::const_iterator ivi = halozone_interior_points.begin();

    for (size_t i = 0; i < nIPHZ && !has_halozone_point; ++i, ivi++)
        if ((*ivi)->isValid())
            has_halozone_point = true;

    // Boundary points
    if (!has_halozone_point)
    {
        list<DLMFD_BoundaryMultiplierPoint *>::const_iterator imp = halozone_boundary_points.begin();

        for (size_t i = 0; i < nBPHZ && !has_halozone_point; ++i, imp++)
            if ((*imp)->isValid())
                has_halozone_point = true;
    }

    return has_halozone_point;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::print_partPointsCoordinates(ofstream &f,
                                                  geomVector const *translated_distance_vector,
                                                  const string &text2write_before,
                                                  const string &text2write_after,
                                                  bool const &withIntPts) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: print_partPointsCoordinates");

    list<DLMFD_BoundaryMultiplierPoint *>::const_iterator bmp;
    list<DLMFD_InteriorMultiplierPoint *>::const_iterator imp;
    size_t dim = 3;

    // Interior points
    imp = interior_points.begin();
    if (withIntPts)
        for (size_t i = 0; i < nIP; ++i, imp++)
            if ((*imp)->isValid())
            {
                f << text2write_before;
                for (size_t index = 0; index < dim; index++)
                    f << (*imp)->get_oneCoordinate(index) << "\t";
                f << text2write_after << endl;
            }

    // Boundary points in the contact region
    bmp = boundary_points.begin();
    for (size_t i = 0; i < nBP; ++i, bmp++)
        if ((*bmp)->isValid())
        {
            f << text2write_before;
            for (size_t index = 0; index < dim; index++)
                f << (*bmp)->get_oneCoordinate(index) << "\t";
            f << text2write_after << endl;
        }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::clear_listOfPointsAndVectors()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: clear_listOfPointsAndVectors");

    list<DLMFD_InteriorMultiplierPoint *>::iterator iip;
    list<DLMFD_BoundaryMultiplierPoint *>::iterator ibp;
    for (iip = interior_points.begin(); iip != interior_points.end(); iip++)
        delete *iip;
    interior_points.clear();
    for (ibp = boundary_points.begin(); ibp != boundary_points.end(); ibp++)
        delete *ibp;
    boundary_points.clear();
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::nullify_Uzawa_vectors()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: nullify_Uzawa_vectors");

    VEC_r.set(0.);
    VEC_x.set(0.);
    VEC_lambda.set(0.);
    VEC_w.set(0.);
}

//---------------------------------------------------------------------------
bool DLMFD_RigidBody::is_interior_points_empty()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: is_interior_points_empty");

    return interior_points.empty();
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::compute_Qu(bool init)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::compute_Qu");

    q_tran.setVecZero();

    if (ndof)
    {
        list<struct ULBD_RHSInfos>::iterator il = points_infos.begin();
        doubleVector const *work = &VEC_w;

        if (init)
            work = &VEC_lambda;

        for (size_t i = 0; i < ndof; ++i, il++)
        {
            q_tran.addOneComp(NDOF_comp[i], -(*work)(i));
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::compute_Qrot(bool init)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::compute_Qrot");

    q_rot_3D.setVecZero();

    if (ndof)
    {
        doubleVector const *work = &VEC_w;
        if (init)
            work = &VEC_lambda;
        double lambda_value = 0.;
        size_t comp = 0;
        for (size_t i = 0; i < ndof; ++i)
        {
            comp = NDOF_comp[i];
            lambda_value = (*work)(i);
            switch (comp)
            {
            case 0:
                q_rot_3D.addOneComp(1, -lambda_value * NDOF_leverage[2 * i]);
                q_rot_3D.addOneComp(2, lambda_value * NDOF_leverage[2 * i + 1]);
                break;
            case 1:
                q_rot_3D.addOneComp(0, lambda_value * NDOF_leverage[2 * i]);
                q_rot_3D.addOneComp(2, -lambda_value * NDOF_leverage[2 * i + 1]);
                break;
            default:
                q_rot_3D.addOneComp(0, -lambda_value * NDOF_leverage[2 * i]);
                q_rot_3D.addOneComp(1, lambda_value * NDOF_leverage[2 * i + 1]);
                break;
            }
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::add_to_Qu(const geomVector &qtran)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::add_to_Qu");

    q_tran += qtran;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::add_to_Qrot(const geomVector &qrot)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::add_to_Qrot");

    q_rot_3D += qrot;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::correctQvectorsAndInitUzawa_Velocity(const double &rho_f,
                                                           const double &timestep,
                                                           const geomVector &gravity_vector_split,
                                                           const geomVector &gravity_vector)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::correctQvectorsAndInitUzawa_Velocity");

    // inverse the inertia matrix
    inversedInertia_3D = calcInvers3by3Matrix(inertia_3D);

    double coeffF = 1. - rho_f / rho_s;

    // correct Q vectors
    q_tran = translational_velocity * (fluidsolid_coupling_factor * mass / timestep) + gravity_vector_split * (coeffF * mass) - q_tran;

    MAC_ASSERT(inversedInertia_3D.size() == dim);

    geomVector Fomega(3);
    for (size_t i = 0; i < dim; ++i)
    {
        double tempVal = 0.;
        for (size_t j = 0; j < dim; ++j)
            tempVal += inertia_3D[i][j] * angular_velocity_3D(j);
        
        Fomega(i) = tempVal * fluidsolid_coupling_factor / timestep;
    }

    q_rot_3D = Fomega - q_rot_3D;
}

//---------------------------------------------------------------------------
vector<vector<double>> DLMFD_RigidBody::calcInvers3by3Matrix(const vector<vector<double>> &oldMatrix)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::calcInvers3by3Matrix");

    MAC_ASSERT(oldMatrix.size() == dim);

    vector<vector<double>> inversedMatrix(dim, vector<double>(dim, 0.));

    double block1 = oldMatrix[1][1] * oldMatrix[2][2] - oldMatrix[1][2] * oldMatrix[2][1];
    double block2 = oldMatrix[1][2] * oldMatrix[2][0] - oldMatrix[1][0] * oldMatrix[2][2];
    double block3 = oldMatrix[1][0] * oldMatrix[2][1] - oldMatrix[1][1] * oldMatrix[2][0];

    double determinat = oldMatrix[0][0] * block1 +
                        oldMatrix[0][1] * block2 +
                        oldMatrix[0][2] * block3;

    inversedMatrix[0][0] = block1 / determinat;
    inversedMatrix[0][1] = (oldMatrix[0][2] * oldMatrix[2][1] - oldMatrix[0][1] * oldMatrix[2][2]) / determinat;
    inversedMatrix[0][2] = (oldMatrix[0][1] * oldMatrix[1][2] - oldMatrix[0][2] * oldMatrix[1][1]) / determinat;
    inversedMatrix[1][0] = block2 / determinat;
    inversedMatrix[1][1] = (oldMatrix[0][0] * oldMatrix[2][2] - oldMatrix[0][2] * oldMatrix[2][0]) / determinat;
    inversedMatrix[1][2] = (oldMatrix[0][2] * oldMatrix[1][0] - oldMatrix[0][0] * oldMatrix[1][2]) / determinat;
    inversedMatrix[2][0] = block3 / determinat;
    inversedMatrix[2][1] = (oldMatrix[0][1] * oldMatrix[2][0] - oldMatrix[0][0] * oldMatrix[2][1]) / determinat;
    inversedMatrix[2][2] = (oldMatrix[0][0] * oldMatrix[1][1] - oldMatrix[0][1] * oldMatrix[1][0]) / determinat;

    return inversedMatrix;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::calculateParticleVelocities(double const &rho_f, double const &timestep)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::calculateParticleVelocities");

    double alpha = fluidsolid_coupling_factor;
    double coeff = timestep / alpha;

    if (component_type != "O")
    {
        t_tran = q_tran * (coeff / mass);

        MAC_ASSERT(inversedInertia_3D.size() == dim);

        t_rot_3D.setVecZero();
        for (size_t i = 0; i < dim; ++i)
        {
            double tempVal = 0.;
            for (size_t j = 0; j < dim; ++j)
                tempVal += inversedInertia_3D[i][j] * q_rot_3D(j);

            t_rot_3D(i) = tempVal * coeff;
        }
    }
    else
    {
        t_tran.setVecZero();
        t_rot_3D.setVecZero();
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::updateSolutionVectors_Velocity(const double &alpha, const double &beta)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::updateSolutionVectors_Velocity");

    translational_velocity = t_tran * alpha + translational_velocity * beta;
    angular_velocity_3D = t_rot_3D * alpha + angular_velocity_3D * beta;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::setTVectors_constant()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::updateSolutionVectors_Velocity");

    t_tran = translational_velocity;
    t_rot_3D = angular_velocity_3D;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::compute_fluid_rhs(DLMFD_ProjectionNavierStokesSystem *GLOBAL_EQ, bool init)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::compute_fluid_rhs");

    if (ndof)
    {
        doubleVector const *work = &VEC_w;
        double lambda_value = 0., coef = 1.;
        if (init)
        {
            work = &VEC_lambda;
            coef = -1.;
        }

        // Optimized version for multi-core architecture using classical
        // aligned-data C vectors

        size_t j, start = 0;

        for (size_t i = 0; i < ndof; ++i)
        {
            lambda_value = (*work)(i);
            for (j = 0; j < NDOF_nfieldUNK[i]; ++j)
                GLOBAL_EQ->assemble_inQUvector(NDOF_deltaOmega[j + start] * lambda_value, NDOF_globalpos[j + start], coef);
            start += NDOF_nfieldUNK[i];
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::compute_x_residuals_Velocity(FV_DiscreteField *pField)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::compute_x_residuals_Velocity");

    double solid_velocity_component = 0.;
    size_t comp = 0, i;
    size_t j, l = 0, m = 0;

    for (i = 0; i < ndof; ++i)
    {
        // Solid part
        comp = NDOF_comp[i];
        solid_velocity_component = t_tran(comp);

        switch (comp)
        {
        case 0:
            solid_velocity_component += t_rot_3D(1) * NDOF_leverage[2 * i] - t_rot_3D(2) * NDOF_leverage[2 * i + 1];
            break;
        case 1:
            solid_velocity_component += t_rot_3D(2) * NDOF_leverage[2 * i + 1] - t_rot_3D(0) * NDOF_leverage[2 * i];
            break;
        default:
            solid_velocity_component += t_rot_3D(0) * NDOF_leverage[2 * i] - t_rot_3D(1) * NDOF_leverage[2 * i + 1];
            break;
        }

        VEC_x(i) = -solid_velocity_component;

        // Fluid part
        for (j = 0; j < NDOF_nfieldUNK[i]; ++j)
        {
            VEC_x(i) += NDOF_deltaOmega[l] * pField->DOF_value(NDOF_FVTriplet[m],
                                                               NDOF_FVTriplet[m + 1],
                                                               NDOF_FVTriplet[m + 2],
                                                               comp,
                                                               0);

            ++l;
            m += 3;
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::compute_r_and_w_FirstUzawaIteration()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::compute_r_and_w_FirstUzawaIteration");

    if (ndof)
    {
        for (size_t i = 0; i < ndof; ++i)
        {
            VEC_r(i) = -VEC_x(i);
            VEC_w(i) = VEC_r(i);
        }
    }
}

//---------------------------------------------------------------------------
double DLMFD_RigidBody::compute_r_dot_r() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::compute_r_dot_r");

    double rdotr = 0.;

    if (ndof)
        for (size_t i = 0; i < ndof; ++i)
            rdotr += VEC_r(i) * VEC_r(i);

    return rdotr;
}

//---------------------------------------------------------------------------
double DLMFD_RigidBody::compute_w_dot_x() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::compute_w_dot_x");

    double wdotx = 0.;

    if (ndof)
        for (size_t i = 0; i < ndof; ++i)
            wdotx += VEC_w(i) * VEC_x(i);

    return wdotx;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::update_lambda_and_r(const double &alpha)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::update_lambda_and_r");

    if (ndof)
    {
        for (size_t i = 0; i < ndof; ++i)
        {
            VEC_lambda(i) -= alpha * VEC_w(i);
            VEC_r(i) -= alpha * VEC_x(i);
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::update_w(const double &beta)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::update_lambda_and_r");

    if (ndof)
        for (size_t i = 0; i < ndof; ++i)
            VEC_w(i) = VEC_r(i) + beta * VEC_w(i);
}
//---------------------------------------------------------------------------
void DLMFD_RigidBody::compute_fluid_DLMFD_explicit(DLMFD_ProjectionNavierStokesSystem *GLOBAL_EQ,
                                                   bool bulk,
                                                   FV_DiscreteField *pField)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::compute_fluid_DLMFD_explicit");

    if (ndof)
    {
        list<struct ULBD_RHSInfos>::const_iterator il = points_infos.begin();
        double lambda_value = 0.;
        pair<size_t, double> Onepair;
        list<size_t>::const_iterator ilpositionU;
        list<double>::const_iterator ildom;
        doubleVector const *work = &VEC_lambda;
        double coef = -1.;

        for (size_t i = 0; i < ndof; ++i, il++)
        {
            lambda_value = (*work)(i);
            ildom = il->omega_delta.begin();

            if (il->omega_delta.size() != 1)
                for (ilpositionU = il->positionU.begin(); ilpositionU != il->positionU.end(); ilpositionU++, ildom++)
                    GLOBAL_EQ->assemble_inQUvector(*ildom * lambda_value, *ilpositionU, coef);

            else if (bulk)
            {
                // Coef 1 / 8
                size_t globalDOFNo;
                FV_TRIPLET ijk = il->MacTripletU.front();
                globalDOFNo = pField->DOF_global_number(ijk.i, ijk.j, ijk.k,
                                                        il->compIdx);
                GLOBAL_EQ->assemble_inQUvector(0.125 * lambda_value, globalDOFNo,
                                               coef);

                // Coef 1 / 16
                if (pField->DOF_is_unknown(ijk.i + 1, ijk.j, ijk.k, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i + 1, ijk.j, ijk.k,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.0625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i - 1, ijk.j, ijk.k, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i - 1, ijk.j, ijk.k,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.0625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i, ijk.j + 1, ijk.k, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i, ijk.j + 1, ijk.k,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.0625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i, ijk.j - 1, ijk.k, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i, ijk.j - 1, ijk.k,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.0625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i, ijk.j, ijk.k + 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i, ijk.j, ijk.k + 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.0625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i, ijk.j, ijk.k - 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i, ijk.j, ijk.k - 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.0625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                // Coef 1 / 32
                if (pField->DOF_is_unknown(ijk.i + 1, ijk.j, ijk.k - 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i + 1, ijk.j, ijk.k - 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.03125 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i, ijk.j + 1, ijk.k - 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i, ijk.j + 1, ijk.k - 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.03125 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i, ijk.j - 1, ijk.k - 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i, ijk.j - 1, ijk.k - 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.03125 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i - 1, ijk.j, ijk.k - 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i - 1, ijk.j, ijk.k - 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.03125 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i + 1, ijk.j, ijk.k + 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i + 1, ijk.j, ijk.k + 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.03125 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i, ijk.j + 1, ijk.k + 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i, ijk.j + 1, ijk.k + 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.03125 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i, ijk.j - 1, ijk.k + 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i, ijk.j - 1, ijk.k + 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.03125 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i - 1, ijk.j, ijk.k + 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i - 1, ijk.j, ijk.k + 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.03125 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i + 1, ijk.j + 1, ijk.k, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i + 1, ijk.j + 1, ijk.k,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.03125 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i + 1, ijk.j - 1, ijk.k, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i + 1, ijk.j - 1, ijk.k,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.03125 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i - 1, ijk.j + 1, ijk.k, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i - 1, ijk.j + 1, ijk.k,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.03125 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i - 1, ijk.j - 1, ijk.k, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i - 1, ijk.j - 1, ijk.k,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.03125 * lambda_value, globalDOFNo,
                                                   coef);
                }

                // Coef 1 / 64
                if (pField->DOF_is_unknown(ijk.i + 1, ijk.j + 1, ijk.k - 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i + 1, ijk.j + 1, ijk.k - 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.015625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i + 1, ijk.j - 1, ijk.k - 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i + 1, ijk.j - 1, ijk.k - 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.015625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i - 1, ijk.j + 1, ijk.k - 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i - 1, ijk.j + 1, ijk.k - 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.015625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i - 1, ijk.j - 1, ijk.k - 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i - 1, ijk.j - 1, ijk.k - 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.015625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i + 1, ijk.j + 1, ijk.k + 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i + 1, ijk.j + 1, ijk.k + 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.015625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i + 1, ijk.j - 1, ijk.k + 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i + 1, ijk.j - 1, ijk.k + 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.015625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i - 1, ijk.j + 1, ijk.k + 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i - 1, ijk.j + 1, ijk.k + 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.015625 * lambda_value, globalDOFNo,
                                                   coef);
                }

                if (pField->DOF_is_unknown(ijk.i - 1, ijk.j - 1, ijk.k + 1, il->compIdx))
                {
                    globalDOFNo = pField->DOF_global_number(ijk.i - 1, ijk.j - 1, ijk.k + 1,
                                                            il->compIdx);
                    GLOBAL_EQ->assemble_inQUvector(0.015625 * lambda_value, globalDOFNo,
                                                   coef);
                }
            }
        }
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::fill_DLMFD_pointInfos(FV_DiscreteField *pField,
                                            ULBD_RHSInfos &OnePointInfos,
                                            bool const &SecondOrderBP)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::fill_DLMFD_pointInfos");

    size_t i, j, k, globalDOFNo;
    vector<doubleVector const *> mesh(dim, NULL);
    doubleVector DOFcoor(dim, 0.), cellsize(dim, 0.);
    geomVector point(dim);
    size_t_vector indices(dim);
    double weightFunct = 1.;
    FV_TRIPLET ijk;
    ijk.i = 0;
    ijk.j = 0;
    ijk.k = 0;

    // Localize in structured mesh
    for (i = 0; i < dim; ++i)
    {
        mesh[i] = pField->get_DOF_coordinates_vector(OnePointInfos.compIdx, i);
        point(i) = OnePointInfos.ptr_point->get_oneCoordinate(i);

        if (component_type != "P")
            FV_Mesh::between(mesh[i], point(i), indices(i));
        else
            FV_Mesh::between_subinterval(mesh[i],
                                         point(i),
                                         indices(i),
                                         (*index_min)(OnePointInfos.compIdx, i),
                                         (*index_max)(OnePointInfos.compIdx, i));

        cellsize(i) = (*mesh[i])(indices(i) + 1) - (*mesh[i])(indices(i));
    }

    // If the particle is periodic, some boundary points correspond to the primary
    // particle position and others to its periodic clones.
    // Hence, it is mandatory to set the geometric features corresponding
    // to this DLM/FD point using GC = point - point_to_GC
    // where point_to_GC is ParticlePoint::GCPointVector
    // and is available from OnePointInfos.ptr_point->get_GCPointVector()
    // !!! The reason is: the method isIn(x,y,z) uses the gravity center
    // position to function !!!
    if (component_type == "PP" && SecondOrderBP)
        translateGeometricFeatures(point - OnePointInfos.ptr_point->get_GCPointVector());

    // Check if second order interpolation is available
    size_t iref, jref, kref;
    bool Q2_interpol = false;
    // In case of Lower set interpolation if the excluding points
    // are at the left or bottom or behind of stencil, we need
    // to iterate inversely (reducing) in x, y or z directions
    bool right = true;
    bool top = true;
    bool front = true;

    geomVector correction_point;

    if (SecondOrderBP)
    {
        correction_point = geomVector(0.01 * radius, 0., 0.);
        if (isIn(point + correction_point))
        {
            if (indices(0) >= pField->get_min_index_unknown_on_proc(
                                  OnePointInfos.compIdx, 0) +
                                  1 &&
                indices(0) + 1 <= pField->get_max_index_unknown_on_proc(
                                      OnePointInfos.compIdx, 0))
                Q2_interpol = true;
            iref = indices(0) - 1;
        }
        else
        {
            if (indices(0) >= pField->get_min_index_unknown_on_proc(
                                  OnePointInfos.compIdx, 0) &&
                indices(0) + 2 <= pField->get_max_index_unknown_on_proc(
                                      OnePointInfos.compIdx, 0))
                Q2_interpol = true;
            iref = indices(0);
        }

        if (Q2_interpol)
        {
            Q2_interpol = false;

            correction_point = geomVector(0., 0.01 * radius, 0.);
            if (isIn(point + correction_point))
            {
                if (indices(1) >= pField->get_min_index_unknown_on_proc(
                                      OnePointInfos.compIdx, 1) +
                                      1 &&
                    indices(1) + 1 <= pField->get_max_index_unknown_on_proc(
                                          OnePointInfos.compIdx, 1))
                    Q2_interpol = true;
                jref = indices(1) - 1;
            }
            else
            {
                if (indices(1) >= pField->get_min_index_unknown_on_proc(
                                      OnePointInfos.compIdx, 1) &&
                    indices(1) + 2 <= pField->get_max_index_unknown_on_proc(
                                          OnePointInfos.compIdx, 1))
                    Q2_interpol = true;
                jref = indices(1);
            }
        }

        if (Q2_interpol)
        {
            Q2_interpol = false;

            correction_point = geomVector(0., 0., 0.01 * radius);
            if (isIn(point + correction_point))
            {
                if (indices(2) >= pField->get_min_index_unknown_on_proc(
                                      OnePointInfos.compIdx, 2) +
                                      1 &&
                    indices(2) + 1 <= pField->get_max_index_unknown_on_proc(
                                          OnePointInfos.compIdx, 2))
                    Q2_interpol = true;
                kref = indices(2) - 1;
            }
            else
            {
                if (indices(2) >= pField->get_min_index_unknown_on_proc(
                                      OnePointInfos.compIdx, 2) &&
                    indices(2) + 2 <= pField->get_max_index_unknown_on_proc(
                                          OnePointInfos.compIdx, 2))
                    Q2_interpol = true;
                kref = indices(2);
            }
        }
    }

    // // Computation of weight of nodes that contribute to the field constraint
    // // * If Q2_interpol is true, use Q2 interpolation with 3 nodes in each
    // // direction, i.e., a 27-nodes stencil,
    // // * Otherwise, use a multi-linear interpolation
    if (Q2_interpol && SecondOrderBP)
    {
        size_t ii, jj, kk;
        double xtransf, ytransf, ztransf, xref, yref, zref;
        // Please refer to Dyn and Floater for definition of variables,
        // and to MAC_solidcomponent2D for more detailed comments on the
        // algorithm

        xref = (*mesh[0])(iref);
        xtransf = (point(0) - xref) / (2. * cellsize(0));

        yref = (*mesh[1])(jref);
        ytransf = (point(1) - yref) / (2. * cellsize(1));

        zref = (*mesh[2])(kref);
        ztransf = (point(2) - zref) / (2. * cellsize(2));

        for (ii = 0; ii < 3; ++ii)
        {
            i = iref + ii;
            for (jj = 0; jj < 3; ++jj)
            {
                j = jref + jj;
                for (kk = 0; kk < 3; ++kk)
                {
                    k = kref + kk;
                    globalDOFNo = pField->DOF_global_number(i, j, k, OnePointInfos.compIdx);
                    OnePointInfos.positionU.push_back(globalDOFNo);

                    ijk.i = i;
                    ijk.j = j;
                    ijk.k = k;
                    OnePointInfos.MacTripletU.push_back(ijk);

                    weightFunct = Q2weighting((*Q2numb)(ii, jj, kk), xtransf, ytransf, ztransf);
                    OnePointInfos.omega_delta.push_back(weightFunct);
                }
            }
        }
    }
    else
    {
        //    cout << "Degenerate as Linear Interpolation" << endl;
        for (i = indices(0); i < indices(0) + 2; ++i)
        {
            DOFcoor(0) = (*mesh[0])(i);
            for (j = indices(1); j < indices(1) + 2; ++j)
            {
                DOFcoor(1) = (*mesh[1])(j);
                for (k = indices(2); k < indices(2) + 2; ++k)
                {
                    DOFcoor(2) = (*mesh[2])(k);
                    globalDOFNo = pField->DOF_global_number(i, j, k,
                                                            OnePointInfos.compIdx);
                    OnePointInfos.positionU.push_back(globalDOFNo);
                    ijk.i = i;
                    ijk.j = j;
                    ijk.k = k;
                    OnePointInfos.MacTripletU.push_back(ijk);
                    weightFunct = (1. - fabs(point(0) - DOFcoor(0)) / cellsize(0)) * (1. - fabs(point(1) - DOFcoor(1)) / cellsize(1)) * (1. - fabs(point(2) - DOFcoor(2)) / cellsize(2));
                    OnePointInfos.omega_delta.push_back(weightFunct);
                }
            }
        }
    }

    // !!! IMPORTANT !!!
    // Remark : In case of a periodic particle, there is no need to reset the
    // gravity center to the primary particle position since the list of
    // boundary points is constructed such that it ends with the primary
    // particle position and this method
    // DLMFD_RigidBody::fill_DLMFD_pointInfos is called when looping on the
    // full list of boundary points
}

//---------------------------------------------------------------------------
double DLMFD_RigidBody::compute_weight_generic(double const &x,
                                               size_t const &i, size_t const &order)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::compute_weight_generic");

    double result = 0.;
    switch (order)
    {
    case 0:
        if (i == 0)
            result = 1.;
        else
            cout << "!!order 0 interpole i is 0!!" << endl;
        break;
    case 1:
        if (i == 0)
            result = 1. - 2. * x;
        else if (i == 1)
            result = 2. * x;
        else
            cout << "!!order 1 interpole i is 0 or 1!!" << endl;
        break;
    case 2:
        if (i == 0)
            result = (x - 0.5) * (x - 1.) / 0.5;
        else if (i == 1)
            result = x * (1. - x) / 0.25;
        else if (i == 2)
            result = x * (x - 0.5) / 0.5;
        else
            cout << "!!order 2 interpole i is 0, 1 or 2!!" << endl;
        break;
    }
    return (result);
}

//---------------------------------------------------------------------------
double DLMFD_RigidBody::Q2weighting(size_t const &i, double const &x,
                                    double const &y, double const &z)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody::Q2weighting");

    double result = 0.;

    switch (i)
    {
    case 0:
        result = 8.0 * (1.0 - x) * (0.5 - x) * (1.0 - y) * (0.5 - y) * (1.0 - z) * (0.5 - z);
        break;
    case 1:
        result = 16.0 * x * (1.0 - x) * (1.0 - y) * (0.5 - y) * (1.0 - z) * (0.5 - z);
        break;
    case 2:
        result = -8.0 * x * (0.5 - x) * (1.0 - y) * (0.5 - y) * (1.0 - z) * (0.5 - z);
        break;
    case 3:
        result = -16.0 * x * (0.5 - x) * y * (1.0 - y) * (1.0 - z) * (0.5 - z);
        break;
    case 4:
        result = 32.0 * x * (1.0 - x) * y * (1.0 - y) * (1.0 - z) * (0.5 - z);
        break;
    case 5:
        result = 16.0 * (1.0 - x) * (0.5 - x) * y * (1.0 - y) * (1.0 - z) * (0.5 - z);
        break;
    case 6:
        result = -8.0 * (1.0 - x) * (0.5 - x) * y * (0.5 - y) * (1.0 - z) * (0.5 - z);
        break;
    case 7:
        result = -16.0 * x * (1.0 - x) * y * (0.5 - y) * (1.0 - z) * (0.5 - z);
        break;
    case 8:
        result = 8.0 * x * (0.5 - x) * y * (0.5 - y) * (1.0 - z) * (0.5 - z);
        break;
    case 9:
        result = 16.0 * (1.0 - x) * (0.5 - x) * (1.0 - y) * (0.5 - y) * z * (1.0 - z);
        break;
    case 10:
        result = 32.0 * x * (1.0 - x) * (1.0 - y) * (0.5 - y) * z * (1.0 - z);
        break;
    case 11:
        result = -16.0 * x * (0.5 - x) * (1.0 - y) * (0.5 - y) * z * (1.0 - z);
        break;
    case 12:
        result = -32.0 * x * (0.5 - x) * y * (1.0 - y) * z * (1.0 - z);
        break;
    case 13:
        result = 64.0 * x * (1.0 - x) * y * (1.0 - y) * z * (1.0 - z);
        break;
    case 14:
        result = 32.0 * (1.0 - x) * (0.5 - x) * y * (1.0 - y) * z * (1.0 - z);
        break;
    case 15:
        result = -16.0 * (1.0 - x) * (0.5 - x) * y * (0.5 - y) * z * (1.0 - z);
        break;
    case 16:
        result = -32.0 * x * (1.0 - x) * y * (0.5 - y) * z * (1.0 - z);
        break;
    case 17:
        result = 16.0 * x * (0.5 - x) * y * (0.5 - y) * z * (1.0 - z);
        break;
    case 18:
        result = -8.0 * (1.0 - x) * (0.5 - x) * (1.0 - y) * (0.5 - y) * z * (0.5 - z);
        break;
    case 19:
        result = -16.0 * x * (1.0 - x) * (1.0 - y) * (0.5 - y) * z * (0.5 - z);
        break;
    case 20:
        result = 8.0 * x * (0.5 - x) * (1.0 - y) * (0.5 - y) * z * (0.5 - z);
        break;
    case 21:
        result = 16.0 * x * (0.5 - x) * y * (1.0 - y) * z * (0.5 - z);
        break;
    case 22:
        result = -32.0 * x * (1.0 - x) * y * (1.0 - y) * z * (0.5 - z);
        break;
    case 23:
        result = -16.0 * (1.0 - x) * (0.5 - x) * y * (1.0 - y) * z * (0.5 - z);
        break;
    case 24:
        result = 8.0 * (1.0 - x) * (0.5 - x) * y * (0.5 - y) * z * (0.5 - z);
        break;
    case 25:
        result = 16.0 * x * (1.0 - x) * y * (0.5 - y) * z * (0.5 - z);
        break;
    case 26:
        result = -8.0 * x * (0.5 - x) * y * (0.5 - y) * z * (0.5 - z);
        break;
    }

    return (result);
}
