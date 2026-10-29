#ifndef GRAINS_EXTERNAL_MODEL_HH
#define GRAINS_EXTERNAL_MODEL_HH

#include "ExternalParticleState.hh"
#include "Torsor.hh"
#include <memory>
#include <string>
#include <vector>

class Particle;

/** One external model and its independently moving contact particles. */
class ExternalModel {
public:
  struct ContactLoad {
    int particleId;
    Torsor wrench; // contact action on particle; world axes, unwrapped centre
  };

  /** Valid only during initialize(); Grains3D owns adopted particles. */
  class Registration {
  public:
    virtual ~Registration() = default;
    /** Geometry must be initialized; initial centres use the primary domain.
     * Returns the stable assigned physical ID.
     */
    virtual int add(std::string const& attachment, std::unique_ptr<Particle> particle) = 0;
    virtual void bind(std::string const& attachment, int existingParticleId) = 0;
    virtual void excludeContacts(std::string const& first, std::string const& second) = 0;
  };

  virtual ~ExternalModel() = default;
  virtual void initialize(Registration& registration) = 0;
  /** Advance this model once over [time, time+dt], using only its own contact
   * loads, held constant over that interval. Return exactly one state per bound
   * particle at time+dt, using continuous unwrapped positions and world velocities.
   * Loads act on the particle in world axes about its unwrapped centre at time.
   * Only contact loads are supplied. The model defines its own inertial properties
   * and applies gravity or other body loads; Grains does not integrate those for it.
   * A failed coupled advance cannot be retried: the model may already have moved.
   */
  virtual std::vector<ExternalParticleState> advance(double time, double dt,
      std::vector<ContactLoad> const& loads) = 0;
};

#endif
