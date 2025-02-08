#ifndef DLMFD_FictitiousDomain_HH
#define DLMFD_FictitiousDomain_HH

#include <MAC_DoubleVector.hh>
#include <FS_SolidPlugIn.hh>
#include <DLMFD_AllRigidBodies.hh>
#include <FV_DomainAndFields.hh>
#include <MAC_ModuleExplorer.hh>
#include <FV_TimeIterator.hh>

/** @brief The Class DLMFD_FictitiousDomain.

Solver for the coupling with particles using a Distributed Lagrange
Multiplier/Fictitious Domain method.

@author A. Wachs & M. Houlette - Pacific project 2024-2025 */

class DLMFD_FictitiousDomain : public MAC_Object
{
public: //-----------------------------------------------------------------
   static DLMFD_FictitiousDomain *create(MAC_Object *a_owner, FV_DomainAndFields const *dom,
                                         MAC_ModuleExplorer const *exp);

   void do_one_inner_iteration(FV_TimeIterator const *t_it);

   void update_rigid_bodies(FV_TimeIterator const *t_it);

   void run_DLMFD_UzawaSolver(FV_TimeIterator const *t_it);

   void initialize_DLMFD_problem(FV_TimeIterator const *t_it);

protected: //--------------------------------------------------------------
private:   //----------------------------------------------------------------
   DLMFD_FictitiousDomain(MAC_Object *a_owner, FV_DomainAndFields const *dom,
                          MAC_ModuleExplorer const *exp);

   ~DLMFD_FictitiousDomain();

   FV_DiscreteField *UU;
   FV_DiscreteField *PP;

   // Physical Parameters
   size_t dim;
   double density;
   double viscosity;

   // Numerical parameters
   size_t sub_prob_number;

   // Grains3D variables
   string solidSolverType;
   FS_SolidPlugIn *solidSolver;
   bool b_solidSolver_parallel;
   string solidSolver_insertionFile;
   string solidSolver_simulationFile;
   istringstream *solidFluid_transferStream;
   DLMFD_AllRigidBodies *allrigidbodies;
   bool are_particles_fixed;

   // Restart
   bool b_restart;
};

#endif
