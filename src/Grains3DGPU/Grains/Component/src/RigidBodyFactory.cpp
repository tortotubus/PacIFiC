#include "RigidBodyFactory.hh"
#include "Box.hh"
#include "Cone.hh"
#include "Convex.hh"
#include "ConvexFactory.hh"
#include "Cylinder.hh"
#include "GrainsParameters.hh"
#include "Matrix3.hh"
#include "QuaternionMath.hh"
#include "Rectangle.hh"
#include "RigidBody.hh"
#include "Sphere.hh"
#include "Superquadric.hh"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

/* ============================================================================================== */
/* High-Level Methods                                                                             */
/* ============================================================================================== */
// Creates and stores a RigidBody object in the host memory.
template <typename T>
__HOST__ void RigidBodyFactory<T>::create(grains::xml::OptionalElement obstacles,
                                          grains::xml::OptionalElement particles,
                                          ParticleData<T>& obstacleData,
                                          ParticleData<T>& particleData,
                                          uint&            numObstacles,
                                          uint&            numParticles)
{
    // Obstacles
    numObstacles = 0;
    const auto allObstacles = obstacles ? obstacles->getChildren()
                                        : std::vector<grains::xml::Element>{};
    uint numRefObs = static_cast<uint>(allObstacles.size());

    obstacleData.numTemplates = numRefObs;
    if(numRefObs > 0)
    {
        obstacleData.refRB.initialize(numRefObs);
        obstacleData.subBodiesPerTemplate.initialize(numRefObs);
        obstacleData.numEachRef.initialize(numRefObs);
        obstacleData.isComposite.initialize(numRefObs);
        obstacleData.refLocalPos.initialize(numRefObs);
        obstacleData.refLocalQuat.initialize(numRefObs);
        obstacleData.refInitialPos.initialize(numRefObs);
        obstacleData.refInitialOri.initialize(numRefObs);
        uint obsIdx = 0;
        for(const auto& nObs : allObstacles)
        {
            obstacleData.numEachRef[obsIdx]           = 1;
            obstacleData.subBodiesPerTemplate[obsIdx] = 1;
            obstacleData.isComposite[obsIdx]          = 0u;
            obstacleData.refRB[obsIdx]                = new RigidBody<T>(nObs);
            obstacleData.refLocalPos[obsIdx]          = Vector3<T>(T(0), T(0), T(0));
            obstacleData.refLocalQuat[obsIdx]         = Quaternion<T>(T(0), T(0), T(0), T(1));
            Vector3<T>    centre(T(0), T(0), T(0));
            Quaternion<T> rotation(T(0), T(0), T(0), T(1));
            auto nTransform = nObs.getChild("Transformation");
            if(nTransform)
            {
                auto nCentre = nTransform->getChild("Centre");
                if(nCentre)
                    centre = Vector3<T>(nCentre);
                auto nRotation = nTransform->getChild("AngularPosition");
                if(nRotation)
                    rotation = Quaternion<T>(nRotation);
            }
            obstacleData.refInitialPos[obsIdx] = centre;
            obstacleData.refInitialOri[obsIdx] = rotation;
            numObstacles += obstacleData.numEachRef[obsIdx];
            ++obsIdx;
        }
    }

    // Particles
    numParticles = 0;
    const auto allParticles = particles ? particles->getChildren()
                                        : std::vector<grains::xml::Element>{};
    uint numRefPar = static_cast<uint>(allParticles.size());
    particleData.numTemplates = 0;
    if(numRefPar > 0)
    {
        // Intermediate storage (vectors before copying into GrainsMemBuffers)
        std::vector<RigidBody<T>*> rbVec;         // per-sub-body prototype
        std::vector<Vector3<T>>    localPosVec;   // per-sub-body prototype
        std::vector<Quaternion<T>> localQuatVec;  // per-sub-body prototype
        std::vector<uint>          subBodiesVec;  // per-template
        std::vector<uint>          numEachVec;    // per-template
        std::vector<uint>          isCompVec;     // per-template
        std::vector<Vector3<T>>    initPosVec;    // per-template
        std::vector<Quaternion<T>> initOriVec;    // per-template
        uint                       numTemplates = 0;

        for(const auto& nParticle : allParticles)
        {
            std::string type = nParticle.attributeOr<std::string>("Type", "");
            uint M_i = nParticle.requireAttribute<uint>("Number");

            // World transform for this template
            Vector3<T>    centre(T(0), T(0), T(0));
            Quaternion<T> rotation(T(0), T(0), T(0), T(1));
            auto nTransform = nParticle.getChild("Transformation");
            if(nTransform)
            {
                auto nCentre = nTransform->getChild("Centre");
                if(nCentre)
                    centre = Vector3<T>(nCentre);
                auto nRotation = nTransform->getChild("AngularPosition");
                if(nRotation)
                    rotation = Quaternion<T>(nRotation);
            }

            if(type != "Composite")
            {
                // Standalone particle
                rbVec.push_back(new RigidBody<T>(nParticle));
                localPosVec.push_back(Vector3<T>(T(0), T(0), T(0)));
                localQuatVec.push_back(Quaternion<T>(T(0), T(0), T(0), T(1)));
                subBodiesVec.push_back(1u);
                numEachVec.push_back(M_i);
                isCompVec.push_back(0u);
                initPosVec.push_back(centre);
                initOriVec.push_back(rotation);
                numParticles += M_i;
                ++numTemplates;
            }
            else
            {
                // Composite particle -- collect and sort sub-bodies by LocalIdx
                std::vector<std::pair<uint, grains::xml::Element>> sortedSubs;
                for(const auto& child : nParticle.getChildren())
                {
                    if(child.name() == "SubBody")
                    {
                        uint localIdx = child.requireAttribute<uint>("LocalIdx");
                        sortedSubs.push_back({localIdx, child});
                    }
                }
                GAssert(!sortedSubs.empty(), "Composite particle has no SubBody children!");
                std::sort(sortedSubs.begin(), sortedSubs.end(), [](const auto& a, const auto& b) {
                    return a.first < b.first;
                });

                uint S_i = sortedSubs.size();
                subBodiesVec.push_back(S_i);
                numEachVec.push_back(M_i);
                isCompVec.push_back(1u);
                initPosVec.push_back(centre);
                initOriVec.push_back(rotation);

                // Density and Material are read from the parent <Particle> node;
                // individual <SubBody> nodes may override Density but share Material.
                T subDensity = T(nParticle.attributeOr<double>("Density", 0.0));
                std::string parentMat
                    = nParticle.requireAttribute<std::string>("Material");
                if(GrainsParameters<T>::m_materialMap.count(parentMat) == 0)
                    GrainsParameters<T>::m_materialMap.emplace(
                        parentMat,
                        (uint)GrainsParameters<T>::m_materialMap.size());
                uint matId = GrainsParameters<T>::m_materialMap[parentMat];

                for(auto& [localIdx, nSub] : sortedSubs)
                {
                    auto nConvex = nSub.getChild("Convex");
                    T ct = T(grains::xml::requireElement(nConvex)
                                 .requireAttribute<double>("CrustThickness"));
                    Convex<T>* convex = ConvexFactory<T>::create(nConvex);
                    T d = T(nSub.attributeOr<double>("Density", subDensity));
                    rbVec.push_back(new RigidBody<T>(convex, ct, d, matId));

                    Vector3<T>    localPos(T(0), T(0), T(0));
                    Quaternion<T> localQuat(T(0), T(0), T(0), T(1));
                    auto nLocalTf = nSub.getChild("LocalTransformation");
                    if(nLocalTf)
                    {
                        auto nLC = nLocalTf->getChild("Centre");
                        if(nLC)
                            localPos = Vector3<T>(nLC);
                        auto nLR = nLocalTf->getChild("AngularPosition");
                        if(nLR)
                            localQuat = Quaternion<T>(nLR);
                    }
                    localPosVec.push_back(localPos);
                    localQuatVec.push_back(localQuat);
                }

                // Compute composite mass + full 3x3 inertia in master frame, then
                // diagonalize via Jacobi so the master uses principal-axis inertia.
                // All local transforms are updated to reflect the new frame; the
                // initial orientation is adjusted by the inverse rotation so world
                // placement is unchanged.
                {
                    uint base = static_cast<uint>(rbVec.size()) - S_i;

                    // Accumulate mass and 3x3 inertia tensor in master body frame
                    T totalMass = T(0);
                    T I[3][3]   = {};
                    for(uint k = 0; k < S_i; ++k)
                    {
                        RigidBody<T>* rb = rbVec[base + k];
                        T             mk = rb->getMass();
                        totalMass += mk;

                        // Rotate sub-body's diagonal inertia into master frame via q_local
                        T Ik[3];
                        rb->getInertia(Ik);
                        Matrix3<T> Rk = localQuatVec[base + k].toMatrix();
                        // I_contribution[i][j] = sum_l R(i,l) * Ik[l] * R(j,l)
                        for(int ii = 0; ii < 3; ++ii)
                            for(int jj = 0; jj < 3; ++jj)
                                for(int l = 0; l < 3; ++l)
                                    I[ii][jj] += Rk(ii, l) * Ik[l] * Rk(jj, l);

                        // Parallel-axis theorem: I += mk*(|r|^2 * delta_ij - r_i*r_j)
                        Vector3<T> rk = localPosVec[base + k];
                        T          r2 = rk[X] * rk[X] + rk[Y] * rk[Y] + rk[Z] * rk[Z];
                        I[0][0] += mk * (r2 - rk[X] * rk[X]);
                        I[1][1] += mk * (r2 - rk[Y] * rk[Y]);
                        I[2][2] += mk * (r2 - rk[Z] * rk[Z]);
                        I[0][1] -= mk * rk[X] * rk[Y];
                        I[1][0] = I[0][1];
                        I[0][2] -= mk * rk[X] * rk[Z];
                        I[2][0] = I[0][2];
                        I[1][2] -= mk * rk[Y] * rk[Z];
                        I[2][1] = I[1][2];
                    }

                    // Jacobi eigendecomposition (cyclic sweeps, max 50 sweeps).
                    // V accumulates Givens rotations; columns of V are eigenvectors.
                    // After convergence, V^T * I * V = diag(eigenvalues).
                    T V[3][3] = {{T(1), T(0), T(0)}, {T(0), T(1), T(0)}, {T(0), T(0), T(1)}};
                    for(int sweep = 0; sweep < 50; ++sweep)
                    {
                        bool converged = true;
                        for(int p = 0; p < 2; ++p)
                        {
                            for(int q = p + 1; q < 3; ++q)
                            {
                                T Ipq = I[p][q];
                                if(std::abs(Ipq) < T(1e-14))
                                    continue;
                                converged = false;
                                T tau     = (I[p][p] - I[q][q]) / (T(2) * Ipq);
                                T t = (tau >= T(0)) ? T(1) / (tau + std::sqrt(T(1) + tau * tau))
                                                    : T(1) / (tau - std::sqrt(T(1) + tau * tau));
                                T c = T(1) / std::sqrt(T(1) + t * t);
                                T s = t * c;
                                // Update diagonal (Ipq -> 0)
                                I[p][p] -= t * Ipq;
                                I[q][q] += t * Ipq;
                                I[p][q] = I[q][p] = T(0);
                                // Update remaining off-diagonal row/col
                                for(int r = 0; r < 3; ++r)
                                {
                                    if(r == p || r == q)
                                        continue;
                                    T Irp   = c * I[r][p] - s * I[r][q];
                                    T Irq   = s * I[r][p] + c * I[r][q];
                                    I[r][p] = I[p][r] = Irp;
                                    I[r][q] = I[q][r] = Irq;
                                }
                                // Accumulate eigenvector matrix: V_new = V_old * G
                                for(int r = 0; r < 3; ++r)
                                {
                                    T Vrp   = c * V[r][p] - s * V[r][q];
                                    T Vrq   = s * V[r][p] + c * V[r][q];
                                    V[r][p] = Vrp;
                                    V[r][q] = Vrq;
                                }
                            }
                        }
                        if(converged)
                            break;
                    }

                    // q_p rotates from old master frame to principal frame (= V^T).
                    // Matrix3 layout is row-major: buf[3*i+j] = M(i,j).
                    T             Rp_buf[9] = {V[0][0],
                                               V[1][0],
                                               V[2][0],
                                               V[0][1],
                                               V[1][1],
                                               V[2][1],
                                               V[0][2],
                                               V[1][2],
                                               V[2][2]};
                    Matrix3<T>    Rp(Rp_buf);
                    Quaternion<T> q_p(Rp);
                    T             qn = norm(q_p);
                    if(qn > T(1e-12))
                        q_p *= (T(1) / qn);

                    // Rotate local positions and local quaternions into principal frame
                    for(uint k = 0; k < S_i; ++k)
                    {
                        localPosVec[base + k]  = q_p >> localPosVec[base + k];
                        localQuatVec[base + k] = q_p * localQuatVec[base + k];
                        T qnk                  = norm(localQuatVec[base + k]);
                        if(qnk > T(1e-12))
                            localQuatVec[base + k] *= (T(1) / qnk);
                    }

                    // Compensate initial orientation: q_init_new = q_init_old * conj(q_p)
                    // so the composite's world placement is unchanged.
                    // Use a const reference to select the value-returning conjugate() overload
                    // (the non-const overload is in-place void and would cause a compile error).
                    const Quaternion<T>& cq_p = q_p;
                    initOriVec.back()         = initOriVec.back() * conjugate(cq_p);
                    T qni                     = norm(initOriVec.back());
                    if(qni > T(1e-12))
                        initOriVec.back() *= (T(1) / qni);

                    // Store composite mass + principal inertia in the master prototype (k=0).
                    // All M_i instances are copied from this prototype, so they all
                    // get the correct properties without repeating the computation.
                    rbVec[base]->setCompositeProperties(totalMass, I[0][0], I[1][1], I[2][2]);
                }

                numParticles += M_i * S_i;
                ++numTemplates;
            }
        }

        // Copy from vectors into GrainsMemBuffers
        uint totalSubBodyProtos   = rbVec.size();
        particleData.numTemplates = numTemplates;
        if(numTemplates > 0)
        {
            particleData.subBodiesPerTemplate.initialize(numTemplates);
            particleData.numEachRef.initialize(numTemplates);
            particleData.isComposite.initialize(numTemplates);
            particleData.refInitialPos.initialize(numTemplates);
            particleData.refInitialOri.initialize(numTemplates);
        }
        if(totalSubBodyProtos > 0)
        {
            particleData.refRB.initialize(totalSubBodyProtos);
            particleData.refLocalPos.initialize(totalSubBodyProtos);
            particleData.refLocalQuat.initialize(totalSubBodyProtos);
        }
        for(uint t = 0; t < numTemplates; ++t)
        {
            particleData.subBodiesPerTemplate[t] = subBodiesVec[t];
            particleData.numEachRef[t]           = numEachVec[t];
            particleData.isComposite[t]          = isCompVec[t];
            particleData.refInitialPos[t]        = initPosVec[t];
            particleData.refInitialOri[t]        = initOriVec[t];
        }
        for(uint r = 0; r < totalSubBodyProtos; ++r)
        {
            particleData.refRB[r]        = rbVec[r];
            particleData.refLocalPos[r]  = localPosVec[r];
            particleData.refLocalQuat[r] = localQuatVec[r];
        }
    }
}

// -------------------------------------------------------------------------------------------------
// Explicit instantiation
template class RigidBodyFactory<float>;
template class RigidBodyFactory<double>;
