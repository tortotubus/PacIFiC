#ifndef _DS_2DRBC__
#define _DS_2DRBC__

#include <DS_ImmersedBoundary.hh>
#include <string>
using std::string;


/** @brief The class DS_2DRBC.

A moving or stationary 2D RBC of axisymmetric cross-section in the
Direction Splitting solver.

@author A. Goyal - Pacific project 2024 */

class DS_2DRBC: public DS_ImmersedBoundary
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_2DRBC();

      /** @brief Destructor */
      ~DS_2DRBC();
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
      void compute_surface_points();

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
      @param copy copied DS_2DRBC object */
      DS_2DRBC( DS_2DRBC const& copy );
      //@}
};

#endif
