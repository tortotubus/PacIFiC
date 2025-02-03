#ifndef DLMFD_FictitiousDomain_HH
#define DLMFD_FictitiousDomain_HH

#include <mpi.h>
#include <FV_OneStepIteration.hh>
#include <PAC_computingtime.hh>
#include <PAC_solvercomputingtime.hh>
#include <DLMFD_NavierStokes.hh>
#include <MAC_DoubleVector.hh>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using namespace std;

class MAC_Communicator ;
class FV_DiscreteField ;
class FV_Mesh ;
class FS_SolidPlugIn ;
class DLMFD_AllRigidBodies ;

/** @brief The Class DLMFD_FictitiousDomain.

Server for the intiating of the NavierStokes class.

@author A. Goyal - Pacific project 2022 */

class DLMFD_FictitiousDomain : public FV_OneStepIteration,
	public PAC_ComputingTime, public PAC_SolverComputingTime
{
   public: //-----------------------------------------------------------------

   //-- Substeps of the step by step progression

      /** @name Substeps of the step by step progression */
      //@{
      /** @brief Tasks performed at initialization of the algorithm, before
      starting the time stepping loop
      @param t_it time iterator */
      virtual void do_before_time_stepping( FV_TimeIterator const* t_it,
      	std::string const& basename ) ;

      /** @brief Perform one time step
      @param t_it time iterator */
      virtual void do_one_inner_iteration( FV_TimeIterator const* t_it ) ;

      /** @brief Tasks performed at initialization of each time step
      @param t_it time iterator */
      virtual void do_before_inner_iterations_stage(
      	FV_TimeIterator const* t_it );

      /** @brief Tasks performed after of each time step
      @param t_it time iterator */
      virtual void do_after_inner_iterations_stage(
      	FV_TimeIterator const* t_it );

      /** @brief Tasks performed at the end of the time stepping loop */
      virtual void do_after_time_stepping( void );

      /** @brief Save additional data than fields
      @param t_it time iterator
      @param cycleNumber cycle number */
      virtual void do_additional_savings( FV_TimeIterator const* t_it,
      	int const& cycleNumber  );

      /** @brief Inbuilt function for post processing data
      @param dom domain and fields
      @param exp module explorer */
      virtual void do_more_post_processing( FV_DomainAndFields * dom,
                                           MAC_ModuleExplorer const* exp ) ;
      //@}

   protected: //--------------------------------------------------------------

   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /** @name Constructors & Destructor */
      //@{
      /** @brief Destructor */
      ~DLMFD_FictitiousDomain( void ) ;

      /** @brief Copy constructor */
      DLMFD_FictitiousDomain( DLMFD_FictitiousDomain const& other ) ;

      /** @brief Operator ==
      @param other the right hand side */
      DLMFD_FictitiousDomain& operator=( DLMFD_FictitiousDomain const& other ) ;

      /** @brief Constructor with arguments
      @param a_owner the MAC-based object
      @param exp to read the data file */
      DLMFD_FictitiousDomain( MAC_Object* a_owner,
      		              FV_DomainAndFields const* dom,
	                       MAC_ModuleExplorer const* exp ) ;

      /** @brief Constructor without argument */
      DLMFD_FictitiousDomain( void ) ;

      /** @brief Create a clone
      @param a_owner the MAC-based object
      @param dom mesh and fields
      @param prms set of parameters
      @param exp to read the data file */
      virtual DLMFD_FictitiousDomain* create_replica(
                     		MAC_Object* a_owner,
                     		FV_DomainAndFields const* dom,
                     		MAC_ModuleExplorer* exp ) const ;
      //@}


      //-- Projection-Translation methods

      /** @name Projection-Translation methods */
      //@{
      /** @brief Set translation vector and direction */
      void set_translation_vector();
      //@}
      

   private: //----------------------------------------------------------------

   //-- Class attributes

      static DLMFD_FictitiousDomain const* PROTOTYPE ;

   //-- Attributes

      double rho;
      double mu;
      double kai;
      string AdvectionScheme;
      string StencilCorrection;
      bool is_CConlyDivergence;
      double FluxRedistThres;
      size_t AdvectionTimeAccuracy;
      size_t space_dimensions;
      bool b_restart;
      bool is_solids;
      bool is_GRAINS;
      bool is_STL;
      string STL_file;
      bool is_HE, is_NS, is_NSwithHE;
      double RBTemp;

      string insertion_type;
      bool is_stressCal;
      string ViscousStressOrder;
      string PressureStressOrder;
      double surface_cell_scale;
      bool is_surfacestressOUT;
      size_t stressCalFreq;
      bool is_par_motion;
      MAC_DoubleVector* gravity_vector ;

      DLMFD_NavierStokes* FlowSolver ;

      // Grains3D variable
      string solidSolverType;
      FS_SolidPlugIn* solidSolver;
      bool b_solidSolver_parallel;
      string solidSolver_insertionFile;
      string solidSolver_simulationFile;
      istringstream* solidFluid_transferStream;
      DLMFD_AllRigidBodies* allrigidbodies;
      bool b_particles_as_fixed_obstacles;
      vector< vector<double> >* hydroFT;

      // Grid motion
      bool b_projection_translation;
      FV_Mesh const* primary_grid;
      double critical_distance_translation;
      geomVector MVQ_translation_vector;
      size_t translation_direction;
      double bottom_coordinate;
      double translated_distance;      

      MAC_Communicator const* macCOMM;
      size_t nb_procs;
      size_t my_rank;
      size_t is_master;      

} ;

#endif
