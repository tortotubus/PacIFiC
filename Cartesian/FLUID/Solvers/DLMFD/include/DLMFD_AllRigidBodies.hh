#ifndef _DLMFD_ALLRIGIDBODIES__
#define _DLMFD_ALLRIGIDBODIES__

#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <FS_AllRigidBodies.hh>
using std::istream;

/** @brief The class DS_AllRigidBodies.

The array of all rigid bodies in the Fictitious Domain solver.

@author M. Houlette - Pacific project 2025 */

class DLMFD_AllRigidBodies
{
public: //-----------------------------------------------------------------
    DLMFD_AllRigidBodies(size_t &dim, istringstream &in, bool const &are_particles_fixed);

    ~DLMFD_AllRigidBodies();

    size_t get_number_rigid_bodies() const;

    void update(istringstream &in);

protected: //--------------------------------------------------------------
private:   //----------------------------------------------------------------
    FS_AllRigidBodies *ptr_FSallrigidbodies;
    size_t RBs_number;
};

#endif