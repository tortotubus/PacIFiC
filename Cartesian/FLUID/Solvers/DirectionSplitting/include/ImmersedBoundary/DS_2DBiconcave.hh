#ifndef _DS_2DBICONCAVE__
#define _DS_2DBICONCAVE__

#include <DS_ImmersedBoundary.hh>
#include <string>
using std::string;


/** @brief The class DS_2DBiconcave.

A moving or stationary 2D Circular of axisymmetric cross-section in the
Direction Splitting solver.

@author A. Goyal - Pacific project 2024 */

class DS_2DBiconcave: public DS_ImmersedBoundary
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_2DBiconcave();

      /** @brief Constructor with arguments
      @param pgrb pointer to the corresponding geometric immersed boundary */
      DS_2DBiconcave( FS_RigidBody* pgrb );

      /** @brief Destructor */
      ~DS_2DBiconcave();
      //@}


   //-- Set methods

      /**@name Set methods */
      //@{

      //@}


   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Compute the surface points by discretizing the Disc
      surface in approximately equal areas (if possible) */
      void compute_surface_parameters();

      /** @brief Compute number of points on a sphere
      @param surface_cell_scale scale of surface cell compared with the grid
      @param dx grid size */
      void compute_number_of_surface_variables( double const& surface_cell_scale
                                              , double const& dx );

      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{

      //@}


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Copy constructor
      @param copy copied DS_2DBiconcave object */
      DS_2DBiconcave( DS_2DBiconcave const& copy );
      //@}
};

#endif
