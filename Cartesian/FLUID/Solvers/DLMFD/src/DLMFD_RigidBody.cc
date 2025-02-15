#include <DLMFD_RigidBody.hh>
#include <FS_RigidBody.hh>
#include <geomVector.hh>
#include <fstream>
using namespace std;

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

    size_t dim = 3;
    index_min = new size_t_array2D(dim, dim, 0);
    index_max = new size_t_array2D(dim, dim, 0);
}

//---------------------------------------------------------------------------
DLMFD_RigidBody::~DLMFD_RigidBody()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: ~DLMFD_RigidBody");
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::setBndPoint(double const &x, double const &y,
                                  double const &z, FV_Mesh const *primary_grid,
                                  list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: setBndPoint");

    (*bp)->set(0, x, y, z, gravity_center);
    ++nBP;
    if (nBP == boundary_points.size())
        extent_bp_list(nBP);
    bp++;
}

//---------------------------------------------------------------------------
size_t DLMFD_RigidBody::get_npts_output(bool const &withIntPts)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: get_npts_output");

    size_t output_npts = nBP;

    if (withIntPts)
        output_npts += nIP;

    return output_npts;
}

// //---------------------------------------------------------------------------
// void DLMFD_RigidBody::add_boundary_point(const geomVector &point, double critical_distance)
// //---------------------------------------------------------------------------
// {
//     MAC_LABEL("DLMFD_RigidBody:: add_boundary_point");

//     DLMFD_BoundaryMultiplierPoint *bmp = new DLMFD_BoundaryMultiplierPoint(point, gravity_center);
//     boundarypoints.push_back(bmp);
// }

// //---------------------------------------------------------------------------
// void DLMFD_RigidBody::add_interior_point(const geomVector &point, double critical_distance)
// //---------------------------------------------------------------------------
// {
//     MAC_LABEL("DLMFD_RigidBody:: add_interior_point");

//     DLMFD_InteriorMultiplierPoint *imp = new DLMFD_InteriorMultiplierPoint(point, gravity_center);
//     interiorpoints.push_back(imp);
// }

//---------------------------------------------------------------------------
void DLMFD_RigidBody::extent_ip_list(size_t const &np)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: extent_ip_list");

    DLMFD_InteriorMultiplierPoint *imp = NULL;
    geomVector virtual_point = geomVector(0., 0., 0.);

    for (size_t i = 0; i < np; ++i)
    {
        interior_points.push_back(imp);
        interior_points.back() =
            new DLMFD_InteriorMultiplierPoint(0, virtual_point, 0, 0, 0, gravity_center);
    }
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::extent_bp_list(size_t const &np)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: extent_bp_list");

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
    size_t dim = translated_distance_vector->getVecSize();

    // Interior points
    imp = interior_points.begin();
    if (withIntPts)
        for (size_t i = 0; i < nIP; ++i, imp++)
        {
            f << text2write_before;
            for (size_t index = 0; index < dim; index++)
                f << (*imp)->get_oneCoordinate(index) - (*translated_distance_vector)(index) << "\t";
            f << text2write_after << endl;
        }

    // Boundary points in the contact region
    bmp = boundary_points.begin();
    for (size_t i = 0; i < nBP; ++i, bmp++)
    {
        f << text2write_before;
        for (size_t index = 0; index < dim; index++)
            f << (*bmp)->get_oneCoordinate(index) - (*translated_distance_vector)(index) << "\t";
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