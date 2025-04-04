#include <DLMFD_AllGeomBoundaries.hh>
#include <DLMFD_GeomBoundary.hh>
#include <DLMFD_XYZBoxPlan.hh>
#include <FV_DiscreteField.hh>

/* Constructor with arguments
-----------------------------*/
DLMFD_AllGeomBoundaries::DLMFD_AllGeomBoundaries(FV_DiscreteField *pField)
{
    MAC_LABEL("DLMFD_AllGeomBoundaries::DLMFD_AllGeomBoundaries(pField)");

    list<pair<size_t, double>> gbFromField = pField->main_geometric_boundaries_on_proc();
    list<pair<size_t, double>>::iterator il;
    DLMFD_GeomBoundary *pGB = NULL;

    for (il = gbFromField.begin(); il != gbFromField.end(); il++)
    {
        all_geomboundaries.push_back(pGB);
        all_geomboundaries.back() = new DLMFD_XYZBoxPlan("DLMFD_XYZBoxPlan", il->first, il->second);
    }
}

/* Destructor
-------------*/
DLMFD_AllGeomBoundaries::~DLMFD_AllGeomBoundaries()
{
    MAC_LABEL("DLMFD_AllGeomBoundaries::~DLMFD_AllGeomBoundaries");

    list<DLMFD_GeomBoundary *>::iterator il;
    for (il = all_geomboundaries.begin(); il != all_geomboundaries.end(); il++)
        delete *il;
    all_geomboundaries.clear();
}

/* Return a pointer to the list of geometric boundaries
-------------------------------------------------------*/
const list<DLMFD_GeomBoundary *> *DLMFD_AllGeomBoundaries::get_ptr_list_geomboundaries() const
{
    MAC_LABEL("DLMFD_AllGeomBoundaries::get_ptr_list_geomboundaries");

    return &all_geomboundaries;
}

/* Operator <<
--------------*/
ostream &operator<<(ostream &f, const DLMFD_AllGeomBoundaries &G)
{
    MAC_LABEL("DLMFD_AllGeomBoundaries::operator <<");

    f << "Number of geometric boundaries = " << G.all_geomboundaries.size() << endl;
    list<DLMFD_GeomBoundary *>::const_iterator il;
    for (il = G.all_geomboundaries.begin(); il != G.all_geomboundaries.end(); il++)
        f << *(*il) << endl;

    return f;
}

/* Display
----------*/
void DLMFD_AllGeomBoundaries::display(ostream &f) const
{
    MAC_LABEL("DLMFD_AllGeomBoundaries::display <<");

    f << "Number of geometric boundaries = " << all_geomboundaries.size() << endl;
    list<DLMFD_GeomBoundary *>::const_iterator il;
    for (il = all_geomboundaries.begin(); il != all_geomboundaries.end(); il++)
    {
        (*il)->display(f);
        f << endl;
    }
}

/* Translate geometric boundaries
---------------------------------*/
void DLMFD_AllGeomBoundaries::translate(const geomVector &translation_vector,
                                  const size_t &translation_direction)
{
    MAC_LABEL("DLMFD_AllGeomBoundaries::translate");

    list<DLMFD_GeomBoundary *>::iterator il;
    for (il = all_geomboundaries.begin(); il != all_geomboundaries.end(); il++)
        (*il)->translate(translation_vector, translation_direction);
}
