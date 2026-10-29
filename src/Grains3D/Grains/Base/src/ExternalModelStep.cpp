/** @file ExternalModelStep.cpp
 * @brief Batched model callbacks within the normal Grains integration stages.
 * Contacts and integration history stay on the live distributed particles.
 */
#include "ExternalModelInternal.hh"
#include <algorithm>
#include <cmath>
#include <climits>
#include <stdexcept>

using namespace GrainsExternalModelDetail;

/** Gather only contact wrenches from physical owners, then call each model once.
 * No particle moves until every callback has returned and the joint batch has
 * passed importParticleStates() validation. Callback order cannot change loads.
 */
std::vector<ExternalParticleState> Grains::advanceModelBatches(double time, double dt)
{
  MPI_Comm comm = m_wrapper ? m_wrapper->get_active_procs_comm() : MPI_COMM_NULL;
  int rank = 0;
  if (comm != MPI_COMM_NULL) MPI_Comm_rank(comm, &rank);
  std::map<int, size_t> indices;
  for (auto const& model : m_externalModels->models)
    for (auto const& binding : model.bindings)
      indices.emplace(binding.particleId, indices.size());
  check(indices.size() > static_cast<size_t>(INT_MAX / 9)
      ? "External load batch too large" : "", comm);
  std::vector<int> owners(indices.size(), 0);
  // Centre, force, torque. Moments stay about the particle centre when relabeling
  // a periodic image with its continuous unwrapped position.
  std::vector<double> values(9 * indices.size(), 0.);
  Vector3 origin, length;
  App::get_origin(origin[X], origin[Y], origin[Z]);
  App::get_size(length[X], length[Y], length[Z]);
  for (auto* p : *m_allcomponents.getActiveParticles()) {
    auto it = indices.find(p->getID());
    if (p->getTag() == 2 || it == indices.end()) continue;
    size_t i = it->second;
    ++owners[i];
    Point3 centre = *p->getPosition();
    auto history = m_externalStateHistory.find(p->getID());
    if (history != m_externalStateHistory.end()) centre = history->second.position;
    else for (int d = 0; d < 3; ++d) if (m_periodicity[d]) {
      double offset = std::fmod(centre[d] - origin[d], length[d]);
      if (offset < 0.) offset += length[d];
      centre[d] = origin[d] + offset;
    }
    for (int d = 0; d < 3; ++d) values[9*i+d] = centre[d];
    p->getContactLoad().copyForceTorque(values.data(), 9*i+3);
  }
  if (comm != MPI_COMM_NULL) {
    MPI_Allreduce(MPI_IN_PLACE, owners.data(), owners.size(), MPI_INT, MPI_SUM, comm);
    if (rank == 0)
      MPI_Reduce(MPI_IN_PLACE, values.data(), values.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
    else
      MPI_Reduce(values.data(), nullptr, values.size(), MPI_DOUBLE, MPI_SUM, 0, comm);
  }
  check(std::any_of(owners.begin(), owners.end(), [](int n) { return n != 1; })
      ? "External particle has no unique physical owner" : "", comm);

  std::vector<ExternalParticleState> result;
  std::string error;
  if (rank == 0) try {
    for (auto& model : m_externalModels->models) {
      std::set<int> expected;
      std::vector<ExternalModel::ContactLoad> loads;
      for (auto const& binding : model.bindings) {
        int id = binding.particleId;
        size_t j = 9 * indices.at(id);
        expected.insert(id);
        loads.push_back({id, Torsor(Point3(values[j], values[j+1], values[j+2]),
            Vector3(values[j+3], values[j+4], values[j+5]),
            Vector3(values[j+6], values[j+7], values[j+8]))});
      }
      auto states = model.instance->advance(time, dt, loads);
      if (states.size() != expected.size())
        throw std::invalid_argument("External model must return exactly its bound particle states");
      for (auto const& state : states)
        if (!expected.erase(state.id))
          throw std::invalid_argument("External model returned an unknown or duplicate particle ID");
      result.insert(result.end(), states.begin(), states.end());
    }
  } catch (std::exception const& e) { error = e.what(); }
    catch (...) { error = "External model advance failed"; }
  status(error, comm);
  return result;
}

/** Shared serial/MPI step. The caller has already advanced the clock to t+dt.
 * Leapfrog reuses endpoint forces, initializes them at t if necessary, and
 * evaluates forces again after drift. Other schemes evaluate once before motion.
 * An external solver may mutate its own state before throwing, so any failed
 * model step poisons the simulation; no transactional rollback is promised.
 */
void Grains::advanceParticleStep()
{
  bool external = hasExternalModels();
  MPI_Comm comm = m_wrapper ? m_wrapper->get_active_procs_comm() : MPI_COMM_NULL;
  if (external) {
    if (m_externalModels->busy) throw std::logic_error("Nested coupled steps are forbidden");
    if (m_externalModels->poisoned)
      throw std::logic_error("A failed coupled advance requires a new simulation");
    m_externalModels->busy = true;
  }
  // Native contact/motion errors must reach all ranks before model collectives.
  auto collectively = [&](auto operation) {
    if (!external) { operation(); return; }
    std::string error;
    try { operation(); }
    catch (std::exception const& e) { error = e.what(); }
    catch (...) { error = "Grains contact or motion evaluation failed"; }
    check(error, comm);
  };
  auto forces = [&]() {
    if (GrainsExec::m_postprocess_forces_at_this_time) m_collision->resetPPForceIndex();
    SCT_set_start("ComputeForces");
    collectively([&] { computeParticlesForceAndTorque(); });
    SCT_get_elapsed_time("ComputeForces");
  };
  try {
    bool leapfrog = GrainsExec::m_TIScheme == "SecondOrderLeapFrog";
    if (external && leapfrog && !m_externalModels->contactsReady) {
      // Initialize acceleration without advancing tangential contact history.
      double endpoint = m_time, dt = m_dt;
      m_time -= dt;
      m_dt = 0.;
      try { forces(); }
      catch (...) { m_time = endpoint; m_dt = dt; throw; }
      m_time = endpoint;
      m_dt = dt;
    } else if (!leapfrog) forces();

    std::vector<ExternalParticleState> states;
    if (external) {
      states = advanceModelBatches(m_time - m_dt, m_dt);
      importParticleStates(states, false);
    }
    SCT_set_start("Move");
    collectively([&] {
      m_allcomponents.Move(m_time, leapfrog ? 0.5*m_dt : m_dt,
                           m_dt, m_dt, m_collision);
    });
    // Import applies the complete external batch and synchronizes all particles,
    // including the native particles just moved above, using normal MPI migration.
    if (external) importParticleStates(states, true);
    else synchronizeParticleStates();
    if (leapfrog) {
      forces();
      collectively([&] { m_allcomponents.advanceParticlesVelocity(m_time, 0.5*m_dt); });
    }
    if (external) {
      m_externalModels->contactsReady = true;
      m_externalModels->busy = false;
    }
  } catch (...) {
    if (external) {
      m_externalModels->busy = false;
      m_externalModels->poisoned = true;
    }
    throw;
  }
}

/** Low-level collective stepping entry point, sharing Simulation() stages.
 * Repeated calls provide fixed explicit subcycling. Scheduled output belongs to
 * Simulation(); this method advances physics and the simulation clock only.
 */
void Grains::advanceExternalModels(double dt)
{
  if (!hasExternalModels()) throw std::logic_error("Register external models before advancing");
  if (m_externalModels->busy) throw std::logic_error("Nested coupled steps are forbidden");
  if (m_externalModels->poisoned)
    throw std::logic_error("A failed coupled advance requires a new simulation");
  if (m_wrapper && !m_processorIsActive) throw std::logic_error("Inactive model rank");
  MPI_Comm comm = m_wrapper ? m_wrapper->get_active_procs_comm() : MPI_COMM_NULL;
  check(!std::isfinite(dt) || dt <= 0. ? "Invalid external-model timestep" : "", comm);
  double low[] = {m_time, dt}, high[] = {m_time, dt};
  if (comm != MPI_COMM_NULL) {
    MPI_Allreduce(MPI_IN_PLACE, low, 2, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(MPI_IN_PLACE, high, 2, MPI_DOUBLE, MPI_MAX, comm);
  }
  check(low[0] != high[0] || low[1] != high[1] ? "MPI clocks or timesteps differ" : "", comm);
  m_dt = dt;
  m_time += dt;
  ++GrainsExec::m_time_counter;
  SCT_insert_app("ComputeForces");
  SCT_insert_app("Move");
  SCT_insert_app("UpdateParticleActivity");
  SCT_insert_app("LinkUpdate");
  advanceParticleStep();
  m_allcomponents.computeNumberParticles(m_wrapper);
}
