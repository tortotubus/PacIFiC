#ifndef GRAINS_EXTERNAL_PARTICLE_STATE_HH
#define GRAINS_EXTERNAL_PARTICLE_STATE_HH

#include "Point3.hh"
#include "Quaternion.hh"
 
struct ExternalParticleState {
  int id;
  Point3 position;
  Quaternion orientation;
  Vector3 translationalVelocity;
  Vector3 angularVelocity;
};

#endif
