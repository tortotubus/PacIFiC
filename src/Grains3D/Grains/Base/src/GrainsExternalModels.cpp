/** @file GrainsExternalModels.cpp
 * @brief XML creation and Grains ownership of external contact particles.
 */
#include "ExternalModelInternal.hh"
#include "Pendulum.hh"
#include <algorithm>
#include <climits>
#include <cmath>
#include <iomanip>
#include <stdexcept>

using namespace GrainsExternalModelDetail;

/** @brief Return whether at least one external model has completed
 * registration. */
bool Grains::hasExternalModels() const {
  return m_externalModels && !m_externalModels->models.empty();
}

/** @brief Return whether external-model registration or stepping is currently
 * in progress. */
bool Grains::isExternalModelStepActive() const {
  return m_externalModels && m_externalModels->busy;
}

/**
 * @brief Transfer a prepared contact particle into local distributed ownership.
 *
 * Force insertion bypasses overlap rejection while retaining normal
 * linked-cell, ghost and periodic-clone handling. Ranks whose grid does not
 * contain the particle discard their temporary copy; retained bodies become
 * Grains-owned.
 */
void Grains::adoptExternalParticle(std::unique_ptr<Particle> particle) {
  bool keep = true;
  if (m_wrapper) {
    auto result = m_collision->insertParticleParallel(
        m_time, particle.get(), m_allcomponents.getActiveParticles(),
        m_allcomponents.getCloneParticles(),
        m_allcomponents.getReferenceParticles(), m_periodic, true, m_wrapper);
    keep = result.first;
  } else
    m_collision->insertParticleSerial(
        particle.get(), m_allcomponents.getActiveParticles(),
        m_allcomponents.getPeriodicCloneParticles(),
        m_allcomponents.getReferenceParticles(), m_periodic, true);
  if (keep)
    m_allcomponents.adoptActiveParticle(std::move(particle),
                                        m_wrapper != nullptr);
}

/**
 * @brief Parse XML model definitions before collision-grid and contact-law
 * setup.
 * @param construction Construction element containing optional ExternalModels.
 * @return Maximum contact-particle circumscribed radius, or zero if absent.
 *
 * Every active rank constructs pending models so materials and geometry bounds
 * are available locally. The current XML factory accepts Pendulum in fresh 3D
 * Standard/MPI runs. Output uses the normal simulation writers.
 * Particle adoption is deferred to registerPendingExternalModels().
 */
double Grains::readExternalModels(DOMNode *construction) {
  auto *models = ReaderXML::getNode(construction, "ExternalModels");
  if (!models)
    return 0.;
  if (m_dimension != 3 || m_restart)
    throw std::invalid_argument(
        "XML external models require a fresh 3D simulation");
  auto type =
      ReaderXML::getNodeAttr_String(construction->getParentNode(), "Type");
  if (type != "Standard" && type != "MPI")
    throw std::invalid_argument(
        "XML external models support Standard and MPI simulations");
  std::set<std::string> names;
  double radius = 0.;
  for (auto *node = models->getFirstChild(); node;
       node = node->getNextSibling()) {
    if (node->getNodeType() != DOMNode::ELEMENT_NODE)
      continue;
    if (ReaderXML::getNodeName(node) != "Pendulum")
      throw std::invalid_argument("Unknown external model type");
    auto name = ReaderXML::getNodeAttr_String(node, "Name");
    if (name.empty() || !names.insert(name).second)
      throw std::invalid_argument(
          "External model names must be unique and nonempty");
    std::unique_ptr<Pendulum> model(new Pendulum(node));
    radius = std::max(radius, model->maximumRadius());
    m_pendingExternalModels.emplace_back(name, std::move(model));
  }
  if (names.empty())
    throw std::invalid_argument("ExternalModels requires at least one model");
  return radius;
}

/**
 * @brief Collectively register parsed models after normal initial insertion.
 *
 * Only root transfers its pending solver instance into the registry. Other
 * ranks participate with null instances, then discard their parsing-only
 * models.
 */
void Grains::registerPendingExternalModels() {
  for (auto &pending : m_pendingExternalModels) {
    std::unique_ptr<ExternalModel> model;
    if (m_rank == 0)
      model = std::move(pending.second);
    registerExternalModel(pending.first, std::move(model));
  }
  m_pendingExternalModels.clear();
}

