#ifndef _DS_3DRBC__
#define _DS_3DRBC__

#include <DS_ImmersedBoundary.hh>
#include <string>
using std::string;


/** @brief The class DS_3DRBC.

A moving or stationary rigid 3D RBC of axisymmetric cross-section in the
Direction Splitting solver.

@author A. Goyal - Pacific project 2024 */

class DS_3DRBC: public DS_ImmersedBoundary
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_3DRBC();

      /** @brief Destructor */
      ~DS_3DRBC();
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
      @param copy copied DS_3DRBC object */
      DS_3DRBC( DS_3DRBC const& copy );
      //@}
};

#endif
