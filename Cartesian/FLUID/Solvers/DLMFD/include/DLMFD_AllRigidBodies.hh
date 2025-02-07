#ifndef _DLMFD_ALLRIGIDBODIES__
#define _DLMFD_ALLRIGIDBODIES__

#include <size_t_vector.hh>
#include <size_t_array2D.hh>
#include <FS_AllRigidBodies.hh>
#include <DLMFD_RigidBody.hh>
#include <FV_DiscreteField.hh>
using namespace std;

/** @brief The class DS_AllRigidBodies.

The array of all rigid bodies in the Fictitious Domain solver.

@author M. Houlette - Pacific project 2025 */

class DLMFD_AllRigidBodies
{
public: //-----------------------------------------------------------------
    DLMFD_AllRigidBodies(size_t &dim,
                         istringstream &in,
                         bool const &are_particles_fixed,
                         FV_DiscreteField const *UF,
                         FV_DiscreteField const *PF);

    ~DLMFD_AllRigidBodies();

    size_t get_number_rigid_bodies() const;

    void update(istringstream &in);

protected: //--------------------------------------------------------------
private:   //----------------------------------------------------------------
    FS_AllRigidBodies *ptr_FSallrigidbodies;
    vector<DLMFD_RigidBody *> vec_ptr_DLMFDallrigidbodies;
    size_t RBs_number;
};

#endif