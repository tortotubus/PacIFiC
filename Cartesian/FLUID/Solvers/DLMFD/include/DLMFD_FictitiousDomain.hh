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
   //-- Substeps of the step by step progression

   /** @brief Create and initialize an instance of DLMFD_FictitiousDomain
   @param a_owner The MAC-based object
   @param dom Domain
   @param exp To read the data file */
   static DLMFD_FictitiousDomain *create(MAC_Object *a_owner, FV_DomainAndFields const *dom,
                                         MAC_ModuleExplorer const *exp);

   /** @name Substeps of the step by step progression */
   //@{

   /** @brief Tasks performed at the main loop
   @param t_it Time iterator */
   void do_one_inner_iteration(FV_TimeIterator const *t_it);

   //@}

   //-- DLMFD solver methods
   /** @name DLMFD solver methods */
   //@{

   /** @brief Prediction of the rigid body attributes (Newton's law)
   @param t_it Time iterator */
   void update_rigid_bodies(FV_TimeIterator const *t_it);

   /** @brief Correction of the rigid body attributes with the Fictitious Domain
   method
   @param t_it Time iterator */
   void run_DLMFD_UzawaSolver(FV_TimeIterator const *t_it);

   /** @brief Initialization of the DLMFD problem to solve
   @param t_it Time iterator */
   void DLMFD_construction(FV_TimeIterator const *t_it);

   //@}

protected: //--------------------------------------------------------------
private:   //----------------------------------------------------------------
   //-- Substeps of the step by step progression

   /** @name Constructors & Destructor */
   //@{

   /** @brief Constructor with arguments
   @param a_owner The MAC-based object
   @param dom Domain
   @param exp To read the data file */
   DLMFD_FictitiousDomain(MAC_Object *a_owner, FV_DomainAndFields const *dom,
                          MAC_ModuleExplorer const *exp);

   /** @brief Constructor with arguments */
   ~DLMFD_FictitiousDomain();

   //@}

   //-- Attributes

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
