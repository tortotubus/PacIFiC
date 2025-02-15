#ifndef DLMFD_FictitiousDomain_HH
#define DLMFD_FictitiousDomain_HH

#include <MAC_DoubleVector.hh>
#include <FS_SolidPlugIn.hh>
#include <DLMFD_AllRigidBodies.hh>
#include <FV_DomainAndFields.hh>
#include <MAC_ModuleExplorer.hh>
#include <FV_TimeIterator.hh>
#include <MAC_Communicator.hh>
using namespace std;

/** @brief The Class DLMFD_FictitiousDomain.

Solver for the coupling with particles using a Distributed Lagrange
Multiplier/Fictitious Domain method.

@author A. Wachs & M. Houlette - Pacific project 2024-2025 */

struct NavierStokes2FluidSolid
{
   string solid_resDir;
};

class DLMFD_FictitiousDomain : public MAC_Object
{
public: //-----------------------------------------------------------------
   //-- Substeps of the step by step progression

   /** @brief Create and initialize an instance of DLMFD_FictitiousDomain
   @param a_owner The MAC-based object
   @param dom Domain
   @param exp To read the data file */
   static DLMFD_FictitiousDomain *create(MAC_Object *a_owner, FV_DomainAndFields const *dom,
                                         MAC_ModuleExplorer const *exp, NavierStokes2FluidSolid transfert);

   /** @name Substeps of the step by step progression */
   //@{

   /** @brief Tasks performed at the main loop
   @param t_it Time iterator */
   void do_one_inner_iteration(FV_TimeIterator const *t_it);

   /** @brief Tasks performed for additional savings
   @param t_it Time iterator */
   void do_additional_savings(int const &cycleNumber,
                              FV_TimeIterator const *t_it);
   // const double &translated_distance,
   // const size_t &translation_direction);

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

   //-- Output methods
   /** @name Output methods */
   //@{

   /** @brief Writing PVTU
   @param filename File name */
   void write_PVTU_multiplier_file(string const &filename) const;

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
                          MAC_ModuleExplorer const *exp, NavierStokes2FluidSolid transfert);

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

   // MPI data
   MAC_Communicator const *pelCOMM;
   size_t size_proc;
   size_t rank;

   // Output
   string SolidSolverResultsDirectory;
   ostringstream Paraview_saveMultipliers_pvd;
   geomVector Paraview_translated_distance_vector;
};

#endif
