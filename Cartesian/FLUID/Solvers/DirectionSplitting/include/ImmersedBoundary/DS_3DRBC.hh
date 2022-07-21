#ifndef _DS_3DRBC__
#define _DS_3DRBC__

#include <DS_ImmersedBoundary.hh>
#include <string>
using std::string;


/** @brief The class DS_3DRBC.

A moving or stationary rigid 3D RBC of axisymmetric cross-section in the
Direction Splitting solver.

@author A. Goyal - Pacific project 2022 */

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
      /** @brief writes one point of RBC mesh to .vtu file **/
      void write_one_point_to_VTK( double const& time
                                         , size_t const& cyclenum );
      //@}


   //-- Methods

      /**@name Methods */
      //@{


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
