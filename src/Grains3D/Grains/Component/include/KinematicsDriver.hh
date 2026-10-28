#ifndef GRAINS_KINEMATICS_DRIVER_HH
#define GRAINS_KINEMATICS_DRIVER_HH

#include <list>
class Particle;

/** Group-level motion dispatch. Drivers do not own particles or the simulation
 * clock. The stepper supplies the integration stage's separate time increments.
 * Concrete drivers must select only the particles they are authorized to move.
 */
class KinematicsDriver {
public:
  virtual ~KinematicsDriver() = default;
  virtual void move(std::list<Particle*> const& particles, double time,
                    double dtVelocity, double dtPosition) const = 0;
  virtual void advanceVelocity(std::list<Particle*> const& particles,
                               double time, double dtVelocity) const = 0;
};

/** Uses existing Particle/CompositeParticle motion and m_kinematics integration.
 * Skips externally driven particles and MPI/periodic ghosts. Stateless: mutable
 * integrator history stays with each particle and its existing time integrator.
 */
class RigidBodyKinematicsDriver final : public KinematicsDriver {
public:
  void move(std::list<Particle*> const& particles, double time,
            double dtVelocity, double dtPosition) const override;
  void advanceVelocity(std::list<Particle*> const& particles,
                       double time, double dtVelocity) const override;
};

#endif
