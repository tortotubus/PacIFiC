#include "Grains.hh"
#include <algorithm>
#include <cmath>
#include <climits>
#include <stdexcept>

namespace {
// Every participant throws the same error before any live state is changed.
void validateCollectively(std::string const& error, MPI_Comm comm, int rank)
{
  if (comm == MPI_COMM_NULL) {
    if (!error.empty()) throw std::invalid_argument(error);
    return;
  }
  int source = error.empty() ? INT_MAX : rank;
  MPI_Allreduce(MPI_IN_PLACE, &source, 1, MPI_INT, MPI_MIN, comm);
  if (source == INT_MAX) return;
  int size = rank == source ? static_cast<int>(error.size()) : 0;
  MPI_Bcast(&size, 1, MPI_INT, source, comm);
  std::string message = rank == source ? error : std::string(size, '\0');
  MPI_Bcast(&message[0], size, MPI_CHAR, source, comm);
  throw std::invalid_argument(message);
}

double wrap(double value, double origin, double length)
{
  double offset = std::fmod(value - origin, length);
  if (offset < 0.) offset += length;
  return origin + (offset < length ? offset : 0.);
}
} // namespace

void Grains::importExternalParticleStates(
    std::vector<ExternalParticleState> const& input)
{
  constexpr int loadWidth = 6 * Component::loadSourceCount;
  if (m_wrapper && !m_processorIsActive)
    throw std::invalid_argument("External states: inactive MPI rank");
  MPI_Comm comm = m_wrapper ? m_wrapper->get_active_procs_comm() : MPI_COMM_NULL;
  int rank = 0;
  if (comm != MPI_COMM_NULL) MPI_Comm_rank(comm, &rank);
  std::string error;
  if (!m_collision || m_dimension != 3 || m_periodicity.size() != 3)
    error = "External states require an initialized 3D simulation";
  if (rank != 0 && !input.empty())
    error = "External states must be submitted only by active rank 0";
  if (input.size() > static_cast<size_t>(INT_MAX / std::max(13, loadWidth)))
    error = "External state snapshot is too large";
  validateCollectively(error, comm, rank);

  int count = rank == 0 ? static_cast<int>(input.size()) : 0;
  if (comm != MPI_COMM_NULL) MPI_Bcast(&count, 1, MPI_INT, 0, comm);
  std::vector<ExternalParticleState> states(count);
  std::vector<int> ids(count);
  // Explicit scalar wire format: position, quaternion (x,y,z,w), v, omega.
  std::vector<double> values(static_cast<size_t>(count) * 13);
  if (rank == 0) {
    for (int i = 0; i < count; ++i) {
      ids[i] = input[i].id;
      for (int d = 0; d < 3; ++d) {
        values[13*i+d] = input[i].position[d];
        values[13*i+7+d] = input[i].translationalVelocity[d];
        values[13*i+10+d] = input[i].angularVelocity[d];
      }
      for (int d = 0; d < 4; ++d) values[13*i+3+d] = input[i].orientation[d];
    }
  }
  if (comm != MPI_COMM_NULL) {
    MPI_Bcast(ids.data(), count, MPI_INT, 0, comm);
    MPI_Bcast(values.data(), count*13, MPI_DOUBLE, 0, comm);
  }
  std::map<int, size_t> index;
  error.clear();
  for (int i = 0; i < count; ++i) {
    auto& s = states[i];
    s.id = ids[i];
    if (!index.emplace(s.id, i).second) error = "Duplicate external particle ID";
    for (int d = 0; d < 13; ++d)
      if (!std::isfinite(values[13*i+d])) error = "Nonfinite external state";
    s.position = Point3(values[13*i], values[13*i+1], values[13*i+2]);
    s.orientation = Quaternion(values[13*i+3], values[13*i+4],
                               values[13*i+5], values[13*i+6]);
    s.translationalVelocity = Vector3(values[13*i+7], values[13*i+8], values[13*i+9]);
    s.angularVelocity = Vector3(values[13*i+10], values[13*i+11], values[13*i+12]);
    double qnorm = std::hypot(std::hypot(s.orientation[0], s.orientation[1]),
                              std::hypot(s.orientation[2], s.orientation[3]));
    if (!(qnorm > 0.) || !std::isfinite(qnorm)) error = "Invalid external quaternion";
    else for (int d = 0; d < 4; ++d) s.orientation[d] /= qnorm;
  }
  validateCollectively(error, comm, rank);

  Vector3 lengths, origin;
  App::get_size(lengths[X], lengths[Y], lengths[Z]);
  App::get_origin(origin[X], origin[Y], origin[Z]);
  Vector3 widths(m_collision->getCellSize(X), m_collision->getCellSize(Y),
                 m_collision->getCellSize(Z));
  for (int d = 0; d < 3; ++d)
    if (!(lengths[d] > 0.) || !(widths[d] > 0.)) error = "Invalid collision domain";
  validateCollectively(error, comm, rank);

  std::vector<int> owners(count, 0);
  std::vector<double> loads(static_cast<size_t>(count) * loadWidth, 0.);
  std::vector<std::pair<Particle*, Transform>> pending;
  for (Particle* p : *m_allcomponents.getActiveParticles()) {
    if (p->getTag() == 2) continue;
    auto found = index.find(p->getID());
    if (p->isIntegratedByGrains()) {
      if (found != index.end()) error = "External state targets an internally driven particle";
      continue;
    }
    if (found == index.end()) {
      error = "External snapshot omits an active externally driven particle";
      continue;
    }
    size_t i = found->second;
    ++owners[i];
    auto breakdown = p->getLoads();
    for (size_t source = 0; source < breakdown.size(); ++source)
      breakdown[source].copyForceTorque(loads.data(), loadWidth*i + 6*source);
    auto const& s = states[i];
    Point3 previous = *p->getPosition();
    auto history = m_externalStateHistory.find(s.id);
    if (history != m_externalStateHistory.end()) previous = history->second.position;
    else {
      // First import uses the primary-domain branch, not a ghost image.
      for (int d = 0; d < 3; ++d)
        if (m_periodicity[d]) previous[d] = wrap(previous[d], origin[d], lengths[d]);
    }
    Vector3 delta = s.position - previous;
    Quaternion old = *p->getQuaternionRotation();
    double norm = Norm(old);
    old *= 1. / norm;
    Quaternion relative = s.orientation * old.Conjugate();
    double angle = 2. * std::atan2(Norm(*relative.getVector3()),
                                   std::abs(relative.getdouble()));
    double motion = Norm(delta) + p->getCircumscribedRadius() * angle;
    if (!std::isfinite(motion) || motion > p->getCrustThickness())
      error = "External particle " + std::to_string(s.id) + " exceeds its motion/crust bound";
    for (int d = 0; d < 3; ++d) {
      if (!std::isfinite(delta[d]) || std::abs(delta[d]) >= widths[d])
        error = "External motion exceeds the linked-cell halo";
      if (!m_periodicity[d] &&
          (s.position[d] < origin[d] || s.position[d] >= origin[d] + lengths[d]))
        error = "External particle exits a nonperiodic domain boundary";
    }
    Transform transform = *p->getRigidBody()->getTransform();
    transform.setBasis(Matrix(s.orientation));
    // Preserve a small displacement in the current image until clone promotion.
    transform.setOrigin(*p->getPosition() + delta);
    pending.emplace_back(p, transform);
  }
  if (comm != MPI_COMM_NULL)
    MPI_Allreduce(MPI_IN_PLACE, owners.data(), count, MPI_INT, MPI_SUM, comm);
  for (int n : owners)
    if (n != 1) error = "External particle ID has no unique physical owner";
  validateCollectively(error, comm, rank);
  if (states.empty()) {
    m_externalStateHistory.clear();
    return;
  }
  if (comm != MPI_COMM_NULL)
    MPI_Allreduce(MPI_IN_PLACE, loads.data(), count*loadWidth, MPI_DOUBLE, MPI_SUM, comm);

  // Allocate the replacement history before committing live state.
  std::map<int, ExternalParticleState> nextHistory;
  for (auto const& s : states) nextHistory.emplace(s.id, s);
  if (GrainsExec::m_partialPer_is_active)
    m_allcomponents.setPositionDir_nm1(GrainsExec::getPartialPeriodicity()->dir);
  try {
    for (auto const& entry : pending) {
      auto const& s = states[index.at(entry.first->getID())];
      entry.first->applyExternalState(entry.second, s.orientation,
                                     s.translationalVelocity, s.angularVelocity);
    }
    SCT_set_start("Move");
    synchronizeParticleStates();
    // Reference-based clone creation can otherwise restore a reference's driver.
    // Reapply metadata and consistent child kinematics to every resulting image.
    auto updateImage = [&](Particle* p) {
      auto it = index.find(p->getID());
      if (it == index.end()) return;
      auto const& s = states[it->second];
      p->setKinematicDriver(Particle::KinematicDriver::External);
      Transform pose = *p->getRigidBody()->getTransform();
      p->applyExternalState(pose, s.orientation,
                            s.translationalVelocity, s.angularVelocity);
      size_t i = it->second;
      Component::LoadBreakdown breakdown;
      for (size_t source = 0; source < breakdown.size(); ++source) {
        size_t j = loadWidth*i + 6*source;
        breakdown[source] = Torsor(*p->getPosition(),
            Vector3(loads[j], loads[j+1], loads[j+2]),
            Vector3(loads[j+3], loads[j+4], loads[j+5]));
      }
      p->setLoads(breakdown);
    };
    for (auto* p : *m_allcomponents.getActiveParticles()) updateImage(p);
    for (auto* p : *m_allcomponents.getCloneParticles()) updateImage(p);
    for (auto const& p : *m_allcomponents.getPeriodicCloneParticles()) updateImage(p.second);
  } catch (...) {
    if (comm != MPI_COMM_NULL) MPI_Abort(comm, 1);
    throw;
  }
  m_externalStateHistory.swap(nextHistory);
}
