#ifndef _DLMFD_ALLRIGIDBODIES__
#define _DLMFD_ALLRIGIDBODIES__

#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <FS_AllRigidBodies.hh>
#include <DLMFD_RigidBody.hh>
#include <FV_DiscreteField.hh>
using namespace std;

/** @brief The class DLMFD_AllRigidBodies.

The array of all rigid bodies in the Fictitious Domain solver.

@author M. Houlette - Pacific project 2025 */

class DLMFD_AllRigidBodies
{
public: //-----------------------------------------------------------------
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{

    /** @brief Constructor with arguments
    @param dimens Number of space dimensions
    @param in Input stream where features of rigid bodies are read
    @param are_particles_fixed True if the all the particles are
    obstacles
    @param UF Pointer to flow field UF
    @param PF Pointer to flow field PF */
    DLMFD_AllRigidBodies(size_t &dim,
                         istringstream &solidFluid_transferStream,
                         bool const &are_particles_fixed,
                         FV_DiscreteField const *UF,
                         FV_DiscreteField const *PF);

    /** @brief Destructor */
    ~DLMFD_AllRigidBodies();
    //@}

    // -- Set methods
    /** @name Set methods */
    //@{

    /** @brief Return the total number of rigid bodies
    @param critical_distance Critical distance */
    void set_all_points(double critical_distance) const;

    /** @brief Set the list of IDs of on proc */
    void set_listIdOnProc();

    //@}

    // -- Get methods
    /** @name Get methods */
    //@{

    /** @brief Return the total number of rigid bodies */
    size_t get_number_rigid_bodies() const;

    //@}

    // -- Update methods
    /** @name Get methods */
    //@{

    /** @brief Update method
    @param solidFluid_transferStream */
    void update(istringstream &solidFluid_transferStream);

    /** @brief Update the RB position and velocity
    @param pos updated position
    @param vel updated translation velocity */
    void update_RB_position_and_velocity(geomVector const &pos,
                                         geomVector const &vel,
                                         geomVector const &ang_vel,
                                         vector<geomVector> const &periodic_directions,
                                         double const &time_step);

    //@}

    // -- Output methods
    /** @name Output methods */
    //@{

    /** @brief Output DLMFD points in Paraview
    @param filename File name */
    void output_DLMFDPoints_PARAVIEW(const string &filename,
                                     geomVector const *translated_distance_vector,
                                     const bool &withIntPts,
                                     size_t rank) const;

    //@}

protected: //----------------------------------------------------------------
private:   //----------------------------------------------------------------
    //-- Attributes

    /** @name Parameters */
    //@{

    size_t RBs_number;                                     /**< Number of rigid bodies */
    FS_AllRigidBodies *ptr_FSallrigidbodies;               /**< Pointer to the geometric rigid bodies */
    vector<DLMFD_RigidBody *> vec_ptr_DLMFDallrigidbodies; /**<  Pointer to the vector of DLMFD rigid bodies */
    list<int> onProc;                                      /**< list of ids of solid component located on this process */

    //@}

    // Pointers to the constant fields and primary grid

    FV_DiscreteField const *pField; /**< Pointer to constrained field*/
    FV_Mesh const *MESH;
};

#endif