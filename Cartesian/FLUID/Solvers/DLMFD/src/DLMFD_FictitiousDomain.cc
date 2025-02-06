#include <DLMFD_FictitiousDomain.hh>

//---------------------------------------------------------------------------
DLMFD_FictitiousDomain::DLMFD_FictitiousDomain(MAC_Object *a_owner) : MAC_Object(a_owner)
//--------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: DLMFD_FictitiousDomain");
}

//---------------------------------------------------------------------------
DLMFD_FictitiousDomain::~DLMFD_FictitiousDomain(void)
//--------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: ~DLMFD_FictitiousDomain");
}

//---------------------------------------------------------------------------
DLMFD_FictitiousDomain *DLMFD_FictitiousDomain::create(MAC_Object *a_owner)
//---------------------------------------------------------------------------
{
   MAC_LABEL("DLMFD_FictitiousDomain:: create");
   DLMFD_FictitiousDomain *result = new DLMFD_FictitiousDomain(a_owner);
   return (result);
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::do_one_inner_iteration()
//---------------------------------------------------------------------------
{
   run_DLMFD_UzawaSolver();
}

//---------------------------------------------------------------------------
void DLMFD_FictitiousDomain::run_DLMFD_UzawaSolver()
//---------------------------------------------------------------------------
{
   printf("Hello World from Uzawa algorithm \n");
}
