#include <DLMFD_GeomBoundary.hh>

//---------------------------------------------------------------------------
DLMFD_GeomBoundary::DLMFD_GeomBoundary(const string &geomtype_)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeomBoundary::DLMFD_GeomBoundary");

    geomtype = geomtype_;
}




//---------------------------------------------------------------------------
DLMFD_GeomBoundary::~DLMFD_GeomBoundary()
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeomBoundary::~DLMFD_GeomBoundary");
}




//---------------------------------------------------------------------------
DLMFD_GeomBoundary::DLMFD_GeomBoundary(const DLMFD_GeomBoundary &M)
//---------------------------------------------------------------------------
{
    MAC_LABEL(
        "DLMFD_GeomBoundary::DLMFD_GeomBoundary(const DLMFD_GeomBoundary &M)");

    geomtype = M.geomtype;
}




//---------------------------------------------------------------------------
ostream &operator<<(ostream &f, const DLMFD_GeomBoundary &G)
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeomBoundary::operator <<");

    f << "Type = " << G.geomtype;

    return f;
}




//---------------------------------------------------------------------------
void DLMFD_GeomBoundary::display(ostream &f) const
//---------------------------------------------------------------------------
{
    MAC_LABEL("DLMFD_GeomBoundary::display");

    f << "Type = " << geomtype;
}
