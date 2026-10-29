#pragma once
#include <vector>
#include <Eigen/Dense>
namespace GrainsRod {
using real_t = double;
using Vec3 = Eigen::Vector3d;
enum class Boundary { Pinned, Free, PinnedRight };
struct BendingRodImplicitConfig {
  int links = 14;
  real_t length = 1.0;
  real_t total_mass = 1.0;
  real_t ei = 0.0;
  Boundary boundary = Boundary::Pinned;
  Vec3 pin = Vec3::Zero();
  std::vector<real_t> masses;
  std::vector<real_t> inertias;
  std::vector<real_t> bending_stiffness;
  real_t damping_rigidity = 0.0;
  std::vector<real_t> bending_damping;
  real_t spectral_radius = 0.8;
  real_t nonlinear_tolerance = 1e-9;
  int max_iterations = 20;
};
/** Adapted from the supplied ELFF BendingRodImplicit implementation.
 * Generalized-alpha integration on unit-tangent charts, with AutoDiff Newton
 * Jacobians. Forces are integrated nodal forces; torques are world-frame link
 * couples held constant during a step. Axial spin is constrained. No ELFF
 * runtime, mesh, checkpoint or force-density interface is required here.
 */
class BendingRodImplicit {
public:
  explicit BendingRodImplicit(const BendingRodImplicitConfig &config);
  std::vector<real_t> nodal_masses() const;
  real_t energy() const { return energy(state_); }
  struct StateParts {
    Vec3 center = Vec3::Zero();
    Vec3 center_velocity = Vec3::Zero();
    std::vector<Vec3> tangents;
    std::vector<Vec3> omega;
  };
  void advance(real_t dt, const std::vector<Vec3>& nodal_forces,
               const std::vector<Vec3>& link_torques);
  const Eigen::VectorXd &state() const { return state_; }
  void setState(const Eigen::VectorXd &state);
  Eigen::MatrixXd nodes(const Eigen::VectorXd &state) const;
  Eigen::MatrixXd velocities(const Eigen::VectorXd &state) const;
  StateParts unpack(const Eigen::VectorXd &state) const;
  Eigen::VectorXd initialState(const std::vector<Vec3> &tangents = {},
                              const std::vector<Vec3> &angular_velocity = {},
                              Vec3 center_of_mass = Vec3::Zero(),
                              Vec3 center_velocity = Vec3::Zero()) const;
  Eigen::VectorXd rhs(real_t time, const Eigen::VectorXd &state,
                      const std::vector<Vec3> &nodal_forces = {}) const;
private:
  Eigen::VectorXd integrateState(const Eigen::VectorXd &state, real_t dt,
                                const std::vector<Vec3> &loads) const;
  int stateSize() const;
  bool isFree() const;
  real_t energy(const Eigen::VectorXd &state) const;
  void normalizeState(Eigen::VectorXd &state) const;
  static void validateNonnegative(real_t value, const char *name);
  std::vector<Vec3> broadcastLinkValues(const std::vector<Vec3> &supplied,
                                      const Vec3 &fallback,
                                      const char *name) const;
  int tangentOffset() const;
  int velocityOffset() const;
  int omegaOffset() const;
  std::vector<Vec3> jointTorques(const std::vector<Vec3> &tangents,
                               const std::vector<Vec3> &omega) const;
  std::vector<Vec3> angularAcceleration(const std::vector<Vec3> &tangents,
                                      const std::vector<Vec3> &omega,
                                      const std::vector<Vec3> &nodal_forces) const;
  std::vector<Vec3> applied_torques_;
  int n_;
  real_t length_;
  real_t h_ = 0.0;
  real_t total_mass_;
  Boundary boundary_;
  Vec3 pin_;
  std::vector<real_t> masses_;
  std::vector<real_t> inertias_;
  std::vector<real_t> bending_stiffness_;
  std::vector<real_t> bending_damping_;
  Eigen::VectorXd com_weights_;
  Eigen::VectorXd state_;
  Eigen::MatrixXd tangent_mass_;
  real_t spectral_radius_, tolerance_;
  int max_iterations_;
  mutable Eigen::VectorXd acceleration_;
  mutable real_t previous_dt_ = 0.0;
};
} // namespace GrainsRod
