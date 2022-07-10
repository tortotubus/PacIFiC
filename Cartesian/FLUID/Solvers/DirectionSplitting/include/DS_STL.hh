#ifndef _DS_STL__
#define _DS_STL__

#include <DS_RigidBody.hh>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
using std::string;
using std::ostream;

#define Nopx 10
#define Nopy 10
#define Nopz 1

/** @brief The class DS_STL.

A stationary rigid STL in the Direction Splitting solver.

@author A. Goyal - Pacific project 2022 */

class DS_STL: public DS_RigidBody
{
   public: //-----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Default constructor */
      DS_STL();

      /** @brief Destructor */
      ~DS_STL();
      //@}


   //-- Set methods

      /**@name Set methods */
      //@{
      /** @brief Updates the rigid body features from its corresponding
      geometric rigid body */
      void update();
      //@}

   //-- Get Methods
      /**@name Get methods */
      //@{
      /** @brief Return the halo zones of the rigid body */
      vector<geomVector*> get_rigid_body_haloZone( ) const;

      /** @brief Returns rigid body velocity including rotation speed at pt
      @param pt the point */
      geomVector get_rigid_body_velocity( geomVector const& pt ) const;

      /** @brief Returns rigid body angular velocity */
      geomVector get_rigid_body_angular_velocity( ) const;

      /** @brief Returns pointer to the rigid body gravity center */
      geomVector const* get_ptr_to_gravity_centre( ) const;

      /** @brief Returns circumscribed radius */
      double get_circumscribed_radius( ) const;

      /** @brief Returns a tuple of mass and density of RB */
      std::tuple<double,double> get_mass_and_density() const;

      /** @brief Returns the distance of a point with the rigid body
      with a given ray vector and source
      @param pt the point
      @param direction x, y or z (=0, 1 or 2)
      @param positive true if search in the positive direction of the coordinate
      axis and false otherwise */
      double get_distanceTo( geomVector const& source,
                             geomVector const& rayDir,
                             double const& delta ) const;

      /** @brief Return the surface points on the rigid body */
      vector<geomVector*> get_rigid_body_surface_points( ) const;

      /** @brief Return the area of the surface points on the rigid body */
      vector<geomVector*> get_rigid_body_surface_areas( ) const;

      /** @brief Return the normal of the surface points on the rigid body */
      vector<geomVector*> get_rigid_body_surface_normals( ) const;

      //@}

   //-- Methods

      /**@name Methods */
      //@{
      /** @brief Writes the attributes in a stream
      @param out output stream
      @param indent_width indentation width */
      void display( ostream& out, size_t const& indent_width ) const;

      /** @brief Compute the halozone of a rigid body
      @param out output stream
      @param indent_width indentation width */
      void compute_rigid_body_halozone( );

      /** @brief Returns whether a point is inside the rigid body
      @param pt the point */
      bool isIn( geomVector const& pt ) const;

      /** @brief Returns whether a point is inside the rigid body
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      bool isIn( double const& x, double const& y, double const& z ) const;

      /** @brief Returns the level set value of a point from the rigid body
      @param pt the point */
      double level_set_value( geomVector const& pt ) const;

      /** @brief Returns the level set value of a point from the rigid body
      @param x x-coordinate of the point
      @param y x-coordinate of the point
      @param z x-coordinate of the point */
      double level_set_value( double const& x
                            , double const& y
                            , double const& z ) const;

      /** @brief Compute the surface points by discretizing the rigid body
      surface in approximately equal areas (if possible) */
      void compute_surface_points( );


      /** @brief Compute number of points on a rigid body
      @param surface_cell_scale scale of surface cell compared with the grid
      @param dx grid size */
      void compute_number_of_surface_variables(
               double const& surface_cell_scale, double const& dx );

      //@}


   protected: //--------------------------------------------------------------

   //-- Attributes

      /**@name Parameters */
      //@{
      vector<geomVector*> m_surface_points; /**< vector of points distributed on
      	the surface of the particle to compute surface integrals */
      vector<geomVector*> m_surface_area; /**< vector of the area associated
         with the points distributed on the surface of the particle */
      vector<geomVector*> m_surface_normal; /**< vector of the normal associated
         with the points distributed on the surface of the particle */
      vector<geomVector> m_surface_Pforce; /**< vector of the pressure force
         on the points distributed on the surface of the particle */
      vector<geomVector> m_surface_Vforce; /**< vector of the viscous force
         on the points distributed on the surface of the particle */
      vector<double> m_surface_Tgrad; /**< vector of the temperature gradient
         on the points distributed on the surface of the particle */
      vector<geomVector*> m_halo_zone; /**< vector of min and max extents
         of rigid body halozone, required for void fraction detection */
      size_t Ntot; /** < Stores the total number of surface points on a rigid
         body */
      //@}

      vector<tuple<double,double,double>> Llvls; /**< vector containing the
    	vertices of the triangulation */
      vector<tuple<double,double,double>> Llvns; /**< vector containing the
    	normals associated to the triangles i.e. one per 3 vertices */

      vector<tuple<double,double,double>> v; /**< intermediate variable */
      vector<vector<tuple<double,double,double>>> vz =
                 vector<vector<tuple<double,double,double>>>(5,v);
      /**< intermediate variable */
      vector<vector<tuple<double,double,double>>> vy =
	         vector<vector<tuple<double,double,double>>>(5,v);
      /**< intermediate variable */

      vector<vector<vector<tuple<double,double,double>>>> ttrbox_xz =
	      vector<vector<vector<tuple<double,double,double>>>>(Nopx,vz);
      /**< vector containing the list of triangles belonging to halos
       * defined accross the xz plane */
      vector<vector<vector<tuple<double,double,double>>>> ttrbox_xy =
	      vector<vector<vector<tuple<double,double,double>>>>(Nopx,vy);
      /**< vector containing the list of triangles belonging to halos
       * defined accross the xz plane */
      vector<vector<vector<tuple<double,double,double>>>> ttrbox_yz =
	      vector<vector<vector<tuple<double,double,double>>>>(Nopy,vz);
      /**< vector containing the list of triangles belonging to halos
       * defined accross the yz plane */

      double tridx_xz[Nopx][Nopz]; /**< Array indicating the number of triangles
				     per halo accross the xz plane */
      double tridx_xy[Nopx][Nopy]; /**< Array indicating the number of triangles
				     per halo accross the xy plane */
      double tridx_yz[Nopy][Nopz]; /**< Array indicating the number of triangles
				     per halo accross the xz plane */

      int Npls; /**<Number of triangles */

   //-- Methods

      /**@name Methods */
      //@{

      //@}


   private: //----------------------------------------------------------------

   //-- Constructors & Destructor

      /**@name Constructors & Destructor */
      //@{
      /** @brief Copy constructor
      @param copy copied DS_STL object */
//      DS_STL( DS_RigidBody const& copy );
      //@}

      // AM 7.7.2022
      //-- Methods

      /**@name Methods */
      //@{

      /** @brief Returns the dot product of two vectors
      @param vect_A[] first vector
      @param vect_B[] second vector */
      double dotProduct(double vect_A[],
		        double vect_B[]) const;

      /** @brief Computes the cross product of two vectors
      @param vect_A[]  first vector
      @param vect_B[]  second vector
      @param cross_P[] result */
      void crossProduct(double vect_A[],
		        double vect_B[],
			double cross_P[]) const;

      /** @brief compute the difference between two vectors
      @param vect_A[] first vector
      @param vect_B[] second vector
      @param diff_P[] result */
      void diffProduct(double vect_A[],
		       double vect_B[],
		       double diff_P[]) const;

      /** @brief Computes the sign of the signed volume
      // of the tetrahedron (A,B,C,D)
      // 1 -> positive or 0 -> negative
      @param vect_A[] first vector
      @param vect_B[] second vector
      @param vect_C[] third vector
      @param vect_D[] fourth vector */
      int orient3d(double vect_A[],
		   double vect_B[],
		   double vect_C[],
		   double vect_D[]) const;

      /** @brief returns 1 if a segment [q1 q2]intersects
       * a triangle (tri1, tri2, tri3)
      @param q1[] first point
      @param q2[] second point
      @param tri1[] first vertex
      @param tri2[] second vertex
      @param tri3[] third vertex */
      int intersect3d(double q1[],
		      double q2[],
		      double tri1[],
		      double tri2[],
		      double tri3[]) const;
      /** @brief reads and STL file
       * and fills the variables Llvls LLvns */
      void readSTL();

      //@}
};

#endif
