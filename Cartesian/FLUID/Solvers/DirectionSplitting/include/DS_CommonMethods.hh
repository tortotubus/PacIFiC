#ifndef DS_CommonMethods_HH
#define DS_CommonMethods_HH

#include <mpi.h>
#include <FV_OneStepIteration.hh>
#include <PAC_computingtime.hh>
#include <PAC_solvercomputingtime.hh>
#include <DS_NavierStokes.hh>
#include <DS_HeatTransfer.hh>
#include <MAC_DoubleVector.hh>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
using namespace std;

class MAC_Communicator;
class FV_DiscreteField;
class FV_Mesh;
class FS_SolidPlugIn;
class DS_AllRigidBodies;
class DS_AllImmersedBoundaries;

/** @brief The Class DS_CommonMethods.

Contains essential functions for both momentum and energy conservation equations.

@author A. Goyal - Pacific project 2024 */

class DS_CommonMethods
{
public: //-----------------------------------------------------------------
        
    //-- Constructors & Destructor
    /** @name Constructors & Destructor */
    //@{
    /** @brief Destructor */
    ~DS_CommonMethods(void);

    /** @brief Constructor with arguments */
    DS_CommonMethods(size_t &dimens, MAC_Communicator const *arb_macCOMM
                                   , FV_DiscreteField *arb_UF
                                   , FV_DiscreteField *arb_PF);

    /** @brief Constructor with arguments */
    DS_CommonMethods(size_t &dimens, MAC_Communicator const *arb_macCOMM
                                   , FV_DiscreteField *arb_TF);

    /** @brief Constructor with arguments */
    DS_CommonMethods(size_t &dimens, DS_AllRigidBodies *arb_allRB
                                   , MAC_Communicator const *arb_macCOMM
                                   , FV_DiscreteField *arb_UF
                                   , FV_DiscreteField *arb_PF);

    /** @brief Constructor with arguments */
    DS_CommonMethods(size_t &dimens, DS_AllRigidBodies *arb_allRB
                                   , MAC_Communicator const *arb_macCOMM
                                   , FV_DiscreteField *arb_TF);

    /** @brief Constructor without argument */
    DS_CommonMethods(void);

    //@}

    //-- Methods
    /** @name Methods */
    //@{

    /** @brief Initialize the values on the field nodes in MAC grid*/
    void initialize_grid_nodes_on_rigidbody(FV_DiscreteField *FF
                                          , vector<size_t> const &list);

    //@}

protected: //--------------------------------------------------------------
private:   //----------------------------------------------------------------


private: //----------------------------------------------------------------
         //-- Class attributes


    size_t dim;
    DS_AllRigidBodies *allrigidbodies;
    MAC_Communicator const *macCOMM;
    FV_Mesh const *MESH;
    FV_DiscreteField *UF;
    FV_DiscreteField *PF;
    FV_DiscreteField *TF;
    bool is_periodic[3];
};

#endif