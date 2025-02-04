#ifndef DLMFD_FictitiousDomain_HH
#define DLMFD_FictitiousDomain_HH

#include <MAC_DoubleVector.hh>

/** @brief The Class DLMFD_FictitiousDomain.

Solver for the coupling with particles using a Distributed Lagrange
Multiplier/Fictitious Domain method.

@author A. Wachs - Pacific project 2024-2025 */

class DLMFD_FictitiousDomain: public MAC_Object
{
   public: //-----------------------------------------------------------------

      static DLMFD_FictitiousDomain* create(MAC_Object* a_owner);

      void do_one_inner_iteration();

      void run_DLMFD_UzawaSolver();

   
   protected: //--------------------------------------------------------------
   
   private: //----------------------------------------------------------------

      DLMFD_FictitiousDomain(MAC_Object* a_owner);

      ~DLMFD_FictitiousDomain(void);
};

#endif
