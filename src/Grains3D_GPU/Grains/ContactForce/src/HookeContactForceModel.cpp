#include "HookeContactForceModel.hh"
#include "GrainsUtils.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ HookeContactForceModel<T>::HookeContactForceModel()
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node
template <typename T>
__HOST__ HookeContactForceModel<T>::HookeContactForceModel(DOMNode* root)
{
    GAssert(ReaderXML::hasNodeAttr(root, "kn"), "kn not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "en"), "en not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "etat"), "etat not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "muc"), "muc not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "kr"), "kr not defined! Aborting Grains!");

    m_kn   = T(ReaderXML::getNodeAttr_Double(root, "kn"));
    m_en   = T(ReaderXML::getNodeAttr_Double(root, "en"));
    m_etat = T(ReaderXML::getNodeAttr_Double(root, "etat"));
    m_muc  = T(ReaderXML::getNodeAttr_Double(root, "muc"));
    m_kr   = T(ReaderXML::getNodeAttr_Double(root, "kr"));

    m_muen = log(m_en) / sqrt(PI<T> * PI<T> + log(m_en) * log(m_en));
}

// -------------------------------------------------------------------------------------------------
// Constructor with five values as contact parameters
template <typename T>
__HOSTDEVICE__ HookeContactForceModel<T>::HookeContactForceModel(T kn, T en, T etat, T muc, T kr)
    : m_kn(kn)
    , m_en(en)
    , m_etat(etat)
    , m_muc(muc)
    , m_kr(kr)
{
    m_muen = log(m_en) / sqrt(PI<T> * PI<T> + log(m_en) * log(m_en));
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ HookeContactForceModel<T>::~HookeContactForceModel()
{
}

// -------------------------------------------------------------------------------------------------
// Gets the ContactForceModel type
template <typename T>
__HOSTDEVICE__ ContactForceModelType HookeContactForceModel<T>::getContactForceModelType() const
{
    return (HOOKE);
}

// -------------------------------------------------------------------------------------------------
// Gets the parameters of the Hooke contact force model
template <typename T>
__HOSTDEVICE__ void HookeContactForceModel<T>::getContactForceModelParameters(
    T& kn, T& en, T& etat, T& muc, T& kr) const
{
    kn   = m_kn;
    en   = m_en;
    etat = m_etat;
    muc  = m_muc;
    kr   = m_kr;
}

// -------------------------------------------------------------------------------------------------
// Performs forces & torques computation
template <typename T>
__HOSTDEVICE__ void
    HookeContactForceModel<T>::performForcesCalculus(const Vector3<T>& contactVector,
                                                     const Vector3<T>& relVelocityAtContact,
                                                     const Vector3<T>& relAngVelocity,
                                                     const T           overlapDistance,
                                                     const T           averageMass,
                                                     Vector3<T>&       delFN,
                                                     Vector3<T>&       delFT,
                                                     Vector3<T>&       delM) const
{
    // Notes:
    // - contactVector is a unit vector pointing from A to B
    // - overlapDistance is negative when there is penetration

    // Normal linear elastic force
    // We do this here as we want to modify the penetration vector later
    delFN = m_kn * overlapDistance * contactVector;

    // Unit normal vector at contact point
    Vector3<T> v_n = (relVelocityAtContact * contactVector) * contactVector;
    Vector3<T> v_t = relVelocityAtContact - v_n;

    // Unit tangential vector along relative velocity at contact point
    T          normv_t = norm(v_t);
    Vector3<T> tangent(0, 0, 0);
    if(normv_t > EPS<T>)
        tangent = v_t / normv_t;

    // Normal dissipative force
    T gamman = -T(2) * m_muen * sqrt(averageMass * m_kn);
    delFN -= gamman * v_n;
    T normFN = norm(delFN);

    // Tangential dissipative force
    // If m_etat = -1, we compute its value such that gamma_n = gamma_t, i.e., same damping in the
    // normal and tangential directions
    T etat   = (m_etat == T(-1)) ? (-m_muen * sqrt(m_kn / averageMass)) : m_etat;
    delFT    = (-T(2) * etat * averageMass) * v_t;
    T normFT = norm(delFT);

    // Tangential Coulomb saturation
    T fn = m_muc * normFN;
    if(fn < normFT)
        delFT = (-fn) * tangent;

    // Rolling resistance moment
    delM = Vector3<T>(0, 0, 0);
    if(m_kr)
    {
        // Relative angular velocity at contact point
        Vector3<T> wn     = (relAngVelocity * contactVector) * contactVector;
        Vector3<T> wt     = relAngVelocity - wn;
        T          normwt = norm(wt);

        // Anti-spinning effect along the normal wn
        delM = -m_kr * normFN * T(0.001) * wn;

        // Classical rolling resistance moment
        if(normwt > EPS<T>)
            delM -= m_kr * normFN * wt;
    }
}

// -------------------------------------------------------------------------------------------------
// Returns a torce based on the contact information
template <typename T>
__HOSTDEVICE__ void HookeContactForceModel<T>::computeForces(const ContactInfo<T>& contactInfos,
                                                             const Vector3<T>& relVelocityAtContact,
                                                             const Vector3<T>& relAngVelocity,
                                                             const Vector3<T>& vA,
                                                             const Vector3<T>& vB,
                                                             ContactHistory<T>* contactHistory,
                                                             Torce<T>&          torceA,
                                                             Torce<T>&          torceB) const
{
    // Note: contactHistory is unused for non-memory models
    (void)contactHistory;

    // Get snapshot with all contact information
    auto snapshot = contactInfos.getSnapshot();

    // Compute contact force and torque
    Vector3<T> delFN, delFT, delM;
    performForcesCalculus(snapshot.contactVector,
                          relVelocityAtContact,
                          relAngVelocity,
                          snapshot.overlapDistance,
                          snapshot.averageMass,
                          delFN,
                          delFT,
                          delM);

    delFN += delFT;
    torceA.addForce(delFN, snapshot.contactPoint - vA);
    torceB.addForce(-delFN, snapshot.contactPoint - vB);
    if(m_kr)
    {
        torceA.addTorque(delM);
        torceB.addTorque(-delM);
    }
}

// -------------------------------------------------------------------------------------------------
// Computes and prints contact parameter estimates for a head-on collision at velocity v0.
template <typename T>
__HOST__ void HookeContactForceModel<T>::computeEstimates(T             massA,
                                                          T             massB,
                                                          T             crustA,
                                                          T             crustB,
                                                          T             v0,
                                                          T             dt,
                                                          std::string   labelA,
                                                          std::string   labelB,
                                                          std::ostream& out) const
{
    // Reduced mass  m* = mA*mB / (mA+mB);  for particle-obstacle pairs massB = 1e20
    T avmass = T(1) / (T(1) / massA + T(1) / massB);

    // Normal damping coefficient eta_n
    T eta_n = -T(2) * m_muen * sqrt(m_kn * avmass);

    // Natural frequency, damped frequency, contact time
    T omega0 = sqrt(m_kn / avmass);
    T theta  = sqrt(omega0 * omega0 - eta_n * eta_n);
    T Tc     = PI<T> / theta;

    // Maximum penetration depth: Newton iteration on  d/dt[delta(t)] = 0
    // delta(t) = (v0/theta) * exp(-eta_n * t) * sin(theta * t)
    T t0 = (m_en > T(0.2)) ? (Tc / T(2)) : (Tc / T(100));
    for(int iter = 0; iter < 1000; ++iter)
    {
        T et = exp(-eta_n * t0);
        T st = sin(theta * t0);
        T ct = cos(theta * t0);
        T f  = (v0 / theta) * et * (-eta_n * st + theta * ct);
        T df = (v0 / theta) * et
               * (eta_n * eta_n * st - T(2) * eta_n * theta * ct - theta * theta * st);
        T dt_ = f / df;
        t0 -= dt_;
        if(fabs(dt_) < T(1e-10) * fabs(t0) + T(1e-14))
            break;
    }
    T delta_max = (v0 / theta) * exp(-eta_n * t0) * sin(theta * t0);

    T crustAllowed = crustA + crustB;
    T maxForce     = m_kn * delta_max;
    T weight       = avmass * T(9.81);

    out << "  Pair [" << labelA << " | " << labelB << "]  (kn=" << m_kn << ", en=" << m_en << ")\n";
    out << "    Crust budget          : " << crustAllowed << " m\n";
    out << "    v0                    : " << v0 << " m/s\n";
    out << "    v0 * dt / crust       : " << v0 * dt / crustAllowed << "\n";
    out << "    Contact time Tc       : " << Tc << " s  (" << Tc / dt << " time steps)\n";
    out << "    Max penetration       : " << delta_max << " m\n";
    out << "    eta_n                 : " << eta_n << "\n";
    out << "    Max elastic force     : " << maxForce << " N\n";
    out << "    Force / weight        : " << maxForce / weight << "\n";
    if(delta_max > crustAllowed)
        out << "    *** WARNING: delta_max (" << delta_max << ") exceeds crust budget ("
            << crustAllowed << "). Increase kn or crust. ***\n";
    out << "\n";
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class HookeContactForceModel<float>;
template class HookeContactForceModel<double>;