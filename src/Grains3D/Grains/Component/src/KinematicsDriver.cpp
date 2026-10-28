#include "KinematicsDriver.hh"
#include "Particle.hh"

void RigidBodyKinematicsDriver::move(std::list<Particle*> const& particles,
    double time, double dtVelocity, double dtPosition) const
{
  for (auto* particle : particles)
    if (particle->getTag() != 2 && particle->isIntegratedByGrains())
      particle->Move(time, dtVelocity, dtPosition);
}

void RigidBodyKinematicsDriver::advanceVelocity(std::list<Particle*> const& particles,
    double time, double dtVelocity) const
{
  for (auto* particle : particles)
    if (particle->getTag() != 2 && particle->isIntegratedByGrains())
      particle->advanceVelocity(time, dtVelocity);
}
