#ifndef DLMFD_FictitiousDomain_HH
#define DLMFD_FictitiousDomain_HH

#include <MAC_DoubleVector.hh>
#include <FS_SolidPlugIn.hh>
#include <DLMFD_AllRigidBodies.hh>
#include <DLMFD_ProjectionNavierStokesSystem.hh>
#include <FV_DomainAndFields.hh>
#include <MAC_ModuleExplorer.hh>
#include <FV_TimeIterator.hh>
#include <MAC_Communicator.hh>
#include <MAC_DoubleVector.hh>
using namespace std;

/** @brief The Class DLMFD_FictitiousDomain.

Solver for the coupling with particles using a Distributed Lagrange
Multiplier/Fictitious Domain method.

@author A. Wachs & M. Houlette - Pacific project 2024-2025 */

struct NavierStokes2FluidSolid
{
   // Output
   string solid_resDir;

   // Parameters
   double rho_f;
   geomVector gravity_vector;
   geomVector split_gravity_vector;

   // Linear resolution
   DLMFD_ProjectionNavierStokesSystem *GLOBAL_EQ;
};

class DLMFD_FictitiousDomain : public MAC_Object
{
public: //-----------------------------------------------------------------
        //-- Public class attributes
   static bool b_SecondOrderInterpol;
   static bool b_LowerSetInterpol;

   //-- Substeps of the step by step progression

   /** @brief Create and initialize an instance of DLMFD_FictitiousDomain
   @param a_owner The MAC-based object
   @param dom Domain
   @param exp To read the data file */
   static DLMFD_FictitiousDomain *create(MAC_Object *a_owner, FV_DomainAndFields const *dom,
                                         MAC_ModuleExplorer const *exp, NavierStokes2FluidSolid transfert);

   /** @name Substeps of the step by step progression */
   //@{

   /** @brief Tasks performed before the main loop
   @param t_it Time iterator */
   void do_before_time_stepping(FV_TimeIterator const *t_it);

   /** @brief Tasks performed at the main loop
   @param t_it Time iterator */
   void do_one_inner_iteration(FV_TimeIterator const *t_it);

   /** @brief Tasks performed for additional savings
   @param t_it Time iterator */
   void do_additional_savings(int const &cycleNumber,
                              FV_TimeIterator const *t_it);

   //@}

   //-- Set methods
   /** @name Set methods */
   //@{

   /** @brief Setting the critical distance attribute
   @param t_icritical_distance_ Critical distance to set */
   void set_critical_distance(double critical_distance_);

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

   /** @brief Solve the DLMFD problem
   @param t_it Time iterator */
   void DLMFD_solving(FV_TimeIterator const *t_it);

   /** @brief Compute sum of q_tran=-<lambda,V>_P and q_rot = -<lambda,xi^GM>_P
   for shared solid components of each process on master process
   @param t_it Time iterator */
   void calculate_ReductionFor_qtranAndqrot(FV_TimeIterator const *t_it);

   /** @brief Complete the initialization of DLM/FD-Uzawa solving algorithm
   @param t_it Time iterator */
   void velocity_broadcast_andUpdate_First(FV_TimeIterator const *t_it);

   /** @brief Solve one DLM/FD-Uzawa iteration for the particles system on
   master process, then for each particle, broadcast t_tran and t_rot
   from master process on each process.
   @param t_it Time iterator */
   void velocity_broadcast_andUpdateInOneIt(FV_TimeIterator const *t_it);

   /** @brief  Broadcast t_tran and t_rot of all particles shared by processes
   from master to processes that own these particles
   @param t_it Time iterator */
   void Broadcast_tVectors_sharedParticles_MasterToAll(FV_TimeIterator const *t_it);

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

   DLMFD_ProjectionNavierStokesSystem *GLOBAL_EQ;

   // Physical Parameters
   size_t dim;
   double rho_f;
   geomVector gravity_vector;
   geomVector split_gravity_vector;

   // Numerical parameters
   size_t sub_prob_number;
   double critical_distance;

   // Grains3D variables
   string solidSolverType;
   FS_SolidPlugIn *solidSolver;
   string solidSolver_insertionFile;
   string solidSolver_simulationFile;
   istringstream *solidFluid_transferStream;
   DLMFD_AllRigidBodies *allrigidbodies;

   // Booleans
   bool b_restart;
   bool b_explicit_added_mass;
   bool are_particles_fixed;
   bool b_solidSolver_parallel;

   // MPI data
   MAC_Communicator const *pelCOMM;
   size_t size_proc;
   size_t rank;
   size_t master;

   // Output
   string SolidSolverResultsDirectory;
   ostringstream Paraview_saveMultipliers_pvd;
   geomVector Paraview_translated_distance_vector;
};

#endif
