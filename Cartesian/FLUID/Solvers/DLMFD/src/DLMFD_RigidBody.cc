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
void DLMFD_RigidBody::setBndPoint(const geomVector &point, list<DLMFD_BoundaryMultiplierPoint *>::iterator &bp)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: setBndPoint");

    (*bp)->set(0, point, gravity_center);
    ++nBP;
    if (nBP == boundary_points.size())
        extend_bp_list(nBP);
    bp++;
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::setIntPoint(const size_t &comp,
                                  const geomVector &point,
                                  size_t i, size_t j, size_t k,
                                  list<DLMFD_InteriorMultiplierPoint *>::iterator &ip)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: setBndPoint");

    (*ip)->set(comp, point, i, j, k, gravity_center);
    ++nIP;
    if (nIP == interior_points.size())
        extend_ip_list(nIP);
    ip++;
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
tuple<double, double, double> DLMFD_RigidBody::get_mass_and_density_and_moi() const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: get_mass_and_density_and_moi");

    return (ptr_FSrigidbody->get_mass_and_density_and_moi());
}

//---------------------------------------------------------------------------
void DLMFD_RigidBody::extend_ip_list(size_t const &np)
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
void DLMFD_RigidBody::extend_bp_list(size_t const &np)
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
    size_t dim = 3;

    // Interior points
    imp = interior_points.begin();
    if (withIntPts)
        for (size_t i = 0; i < nIP; ++i, imp++)
        {
            f << text2write_before;
            for (size_t index = 0; index < dim; index++)
                f << (*imp)->get_oneCoordinate(index) << "\t";
            // cout << (*imp)->get_oneCoordinate(dim - 3)
            //      << " " << (*imp)->get_oneCoordinate(dim - 2)
            //      << " " << (*imp)->get_oneCoordinate(dim - 1)
            //      << endl;
            f << text2write_after << endl;
        }

    // Boundary points in the contact region
    bmp = boundary_points.begin();
    for (size_t i = 0; i < nBP; ++i, bmp++)
    {
        f << text2write_before;
        for (size_t index = 0; index < dim; index++)
            f << (*bmp)->get_oneCoordinate(index) << "\t";
        // cout << (*bmp)->get_oneCoordinate(dim - 3)
        //      << " " << (*bmp)->get_oneCoordinate(dim - 2)
        //      << " " << (*bmp)->get_oneCoordinate(dim - 1)
        //      << endl;
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
bool DLMFD_RigidBody::is_interior_points_empty()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_RigidBody:: is_interior_points_empty");

    return interior_points.empty();
}