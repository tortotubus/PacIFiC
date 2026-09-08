#include "HookeMemoryContactForceModel.hh"
#include "GrainsParameters.hh"
#include "GrainsUtils.hh"
#include "QuaternionMath.hh"
#include "VectorMath.hh"

// -------------------------------------------------------------------------------------------------
// Default constructor
template <typename T>
__HOSTDEVICE__ HookeMemoryContactForceModel<T>::HookeMemoryContactForceModel()
{
}

// -------------------------------------------------------------------------------------------------
// Constructor with an XML node
template <typename T>
__HOST__ HookeMemoryContactForceModel<T>::HookeMemoryContactForceModel(DOMNode* root)
{
    GAssert(ReaderXML::hasNodeAttr(root, "kn"), "kn not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "en"), "en not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "kt"), "kt not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "etat"), "etat not defined! Aborting Grains!");
    GAssert(ReaderXML::hasNodeAttr(root, "muc"), "muc not defined! Aborting Grains!");
    GAssert(GrainsParameters<T>::m_dt > T(0), "Time step not defined! Aborting Grains!");

    m_kn   = T(ReaderXML::getNodeAttr_Double(root, "kn"));
    m_en   = T(ReaderXML::getNodeAttr_Double(root, "en"));
    m_kt   = T(ReaderXML::getNodeAttr_Double(root, "kt"));
    m_etat = T(ReaderXML::getNodeAttr_Double(root, "etat"));
    m_muc  = T(ReaderXML::getNodeAttr_Double(root, "muc"));
    if(ReaderXML::hasNodeAttr(root, "mur"))
        m_mur = T(ReaderXML::getNodeAttr_Double(root, "mur"));
    if(ReaderXML::hasNodeAttr(root, "etarpf"))
        m_etarpf = T(ReaderXML::getNodeAttr_Double(root, "etarpf"));

    m_muen = log(m_en) / sqrt(PI<T> * PI<T> + log(m_en) * log(m_en));
    m_dt   = GrainsParameters<T>::m_dt;
}

// -------------------------------------------------------------------------------------------------
// Constructor with eight values as contact parameters
template <typename T>
__HOSTDEVICE__ HookeMemoryContactForceModel<T>::HookeMemoryContactForceModel(
    T kn, T en, T kt, T etat, T muc, T mur, T etarpf, T dt)
    : m_kn(kn)
    , m_en(en)
    , m_kt(kt)
    , m_etat(etat)
    , m_muc(muc)
    , m_mur(mur)
    , m_etarpf(etarpf)
    , m_dt(dt)
{
    m_muen = log(m_en) / sqrt(PI<T> * PI<T> + log(m_en) * log(m_en));
}

// -------------------------------------------------------------------------------------------------
// Destructor
template <typename T>
__HOSTDEVICE__ HookeMemoryContactForceModel<T>::~HookeMemoryContactForceModel()
{
}

// -------------------------------------------------------------------------------------------------
// Gets the ContactForceModel type
template <typename T>
__HOSTDEVICE__ ContactForceModelType
    HookeMemoryContactForceModel<T>::getContactForceModelType() const
{
    return (HOOKEMEMORY);
}

// -------------------------------------------------------------------------------------------------
// Gets the parameters of the HookeMemory contact force model
template <typename T>
__HOSTDEVICE__ void HookeMemoryContactForceModel<T>::getContactForceModelParameters(
    T& kn, T& en, T& kt, T& etat, T& muc, T& mur, T& etarpf) const
{
    kn     = m_kn;
    en     = m_en;
    kt     = m_kt;
    etat   = m_etat;
    muc    = m_muc;
    mur    = m_mur;
    etarpf = m_etarpf;
}

// -------------------------------------------------------------------------------------------------
// Performs forces & torques computation with optional memory tracking
template <typename T>
__HOSTDEVICE__ void
    HookeMemoryContactForceModel<T>::performForcesCalculus(const Vector3<T>&  contactVector,
                                                           const Vector3<T>&  relVelocityAtContact,
                                                           const Vector3<T>&  relAngVelocity,
                                                           const T            overlapDistance,
                                                           const T            averageMass,
                                                           const T            averageRadius,
                                                           ContactHistory<T>* contactHistory,
                                                           Vector3<T>&        delFN,
                                                           Vector3<T>&        delFT,
                                                           Vector3<T>&        delM) const
{
    // Guard against a null history pointer. This can happen when the contact hash table is
    // full and findOrInsert() had to return false (no empty / tombstone slot available). In
    // that case we fall back to a zero-initialised local history so the contact is treated
    // as a new one for this timestep. The simulation stays numerically correct (normal force
    // is exact; tangential spring starts fresh) and no illegal memory access occurs.
    ContactHistory<T> fallbackHistory;
    if(contactHistory == nullptr)
        contactHistory = &fallbackHistory;

    // Notes:
    // - contactVector is a unit vector pointing from A to B (the normal)
    // - overlapDistance is negative when there is penetration

    // Penetration vector
    Vector3<T> normal = contactVector;

    // Relative velocity components
    Vector3<T> v_n = (relVelocityAtContact * normal) * normal;
    Vector3<T> v_t = relVelocityAtContact - v_n;

    // 1) Compute normal force
    // Normal linear elastic force
    delFN = m_kn * overlapDistance * normal;

    // Normal dissipative force
    T gamman = -T(2) * m_muen * sqrt(averageMass * m_kn);
    delFN -= gamman * v_n;
    T normFN = norm(delFN);

    // 2) Compute tangential force with memory
    // Check if this is a new contact (previousNormal is zero)
    Vector3<T> prevNormal     = contactHistory->m_previousNormal;
    bool       contactExisted = (norm(prevNormal) > EPS<T>);

    // Rotate both cumulative displacements to the current contact plane
    if(contactExisted)
    {
        Quaternion<T> qrot;
        qrot.setRotFromTwoVectors(prevNormal, normal);
        contactHistory->m_tangentialDisplacement = qrot >> contactHistory->m_tangentialDisplacement;
        contactHistory->m_rollingDisplacement    = qrot >> contactHistory->m_rollingDisplacement;
    }

    // Add contribution of current timestep
    contactHistory->m_tangentialDisplacement += m_dt * v_t;

    // Update the normal vector in history
    contactHistory->m_previousNormal = normal;

    // Compute tangential force direction
    // If m_etat = -1, we compute its value such that gamma_n = gamma_t, i.e., same damping in the
    // normal and tangential directions
    T          etat        = (m_etat == T(-1)) ? (-m_muen * sqrt(m_kn / averageMass)) : m_etat;
    Vector3<T> viscousFT   = (-T(2) * etat * averageMass) * v_t;
    Vector3<T> tentativeFT = -m_kt * contactHistory->m_tangentialDisplacement + viscousFT;
    T          normFT      = norm(tentativeFT);
    Vector3<T> tangentDir = (normFT > EPS<T>) ? tentativeFT / normFT : Vector3<T>(T(0), T(0), T(0));

    // Compute tangential force with Coulomb limit
    if(normFT <= m_muc * normFN)
    {
        // Below Coulomb limit
        delFT = normFT * tangentDir;
    }
    else
    {
        // Above Coulomb limit - apply saturation and adjust history
        delFT = m_muc * normFN * tangentDir;
        if(m_kt > EPS<T>)
        {
            contactHistory->m_tangentialDisplacement
                = (-m_muc * normFN * tangentDir + viscousFT) / m_kt;
        }
    }

    // 3) Compute rolling resistance torque with memory (if applicable)
    delM = Vector3<T>(T(0), T(0), T(0));
    if(m_mur > EPS<T>)
    {
        // Using Jiang et al (2005, 2015) formulation
        T Req       = averageRadius;
        T kr        = T(3) * m_kn * m_mur * m_mur * Req * Req;    // torque spring stiffness
        T etar      = T(3) * gamman * m_mur * m_mur * Req * Req;  // torque dashpot coefficient
        T maxNormMk = m_mur * Req * normFN;                       // saturation torque

        // Rolling component of relative angular velocity (exclude spinning about the normal)
        Vector3<T> wt_r = relAngVelocity - (relAngVelocity * normal) * normal;

        // Update rolling friction spring (rotation was applied above; only rolling component)
        contactHistory->m_rollingDisplacement -= kr * m_dt * wt_r;
        T normMk = norm(contactHistory->m_rollingDisplacement);

        // Apply saturation
        if(normMk > maxNormMk)
            contactHistory->m_rollingDisplacement *= maxNormMk / normMk;

        delM = contactHistory->m_rollingDisplacement - m_etarpf * etar * wt_r;
    }
}

// -------------------------------------------------------------------------------------------------
// Returns a force based on the contact information with memory tracking
template <typename T>
__HOSTDEVICE__ void
    HookeMemoryContactForceModel<T>::computeForces(const ContactInfo<T>& contactInfos,
                                                   const Vector3<T>&     relVelocityAtContact,
                                                   const Vector3<T>&     relAngVelocity,
                                                   const Vector3<T>&     vA,
                                                   const Vector3<T>&     vB,
                                                   ContactHistory<T>*    contactHistory,
                                                   Torce<T>&             torceA,
                                                   Torce<T>&             torceB) const
{
    // Get snapshot with all contact information
    auto snapshot = contactInfos.getSnapshot();

    // Compute contact forces and torques
    Vector3<T> delFN, delFT, delM;
    performForcesCalculus(snapshot.contactVector,
                          relVelocityAtContact,
                          relAngVelocity,
                          snapshot.overlapDistance,
                          snapshot.averageMass,
                          snapshot.averageRadius,
                          contactHistory,
                          delFN,
                          delFT,
                          delM);

    // Apply forces and torques
    Vector3<T> totalForce = delFN + delFT;
    torceA.addForce(totalForce, snapshot.contactPoint - vA);
    torceB.addForce(-totalForce, snapshot.contactPoint - vB);

    if(m_mur > EPS<T>)
    {
        torceA.addTorque(delM);
        torceB.addTorque(-delM);
    }
}

// -------------------------------------------------------------------------------------------------
// Computes and prints contact parameter estimates for a head-on collision at velocity v0.
template <typename T>
__HOST__ void HookeMemoryContactForceModel<T>::computeEstimates(T             massA,
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

    // Normal damping coefficient eta_n (uses m_muen computed at construction)
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
template class HookeMemoryContactForceModel<float>;
template class HookeMemoryContactForceModel<double>;
