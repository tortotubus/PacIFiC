#ifndef GRAINS_PENDULUM_HH
#define GRAINS_PENDULUM_HH

#include "ExternalModel.hh"
#include "ReaderXML.hh"
#include "Particle.hh"

/** Pinned, inextensible bending rod with one contact spherocylinder per link.
 * The rod owns its implicit dynamics; Grains owns its collision particles.
 */
class Pendulum final : public ExternalModel {
public:
  explicit Pendulum(DOMNode* node);
  ~Pendulum() override;
  void initialize(Registration& registration) override;
  std::vector<ExternalParticleState> advance(double time, double dt,
      std::vector<ContactLoad> const& loads) override;
  double maximumRadius() const { return m_radius+0.5*m_length/m_count; }
private:
  Point3 m_anchor;
  Vector3 m_direction;
  size_t m_count=0;
  double m_length=0., m_radius=0., m_density=0., m_crustThickness=0.;
  double m_bendingRigidity=0., m_dampingRigidity=0.;
  std::string m_material;
  std::vector<ExternalParticleState> m_states;
  bool m_initialized=false;
  struct Dynamics;
  std::unique_ptr<Dynamics> m_dynamics;
};

#endif
