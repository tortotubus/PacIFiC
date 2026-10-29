/** @file ExternalModelRegistration.cpp
 * @brief Model initialization, attachment validation and distributed particle adoption.
 * Registration stages new objects before publishing their ownership to Grains.
 */
#include "ExternalModelInternal.hh"
#include "ContactBuilderFactory.hh"
#include <algorithm>
#include <climits>
#include <cmath>
#include <iomanip>
#include <stdexcept>

using namespace GrainsExternalModelDetail;

/** Temporary implementation of the initialize() registration callback. */
class Grains::ModelRegistration final : public ExternalModel::Registration {
  struct NewParticle {
    std::string attachment;
    std::unique_ptr<Particle> body;
  };
  Grains &g;
  ExternalModelData &registry;
  MPI_Comm comm = MPI_COMM_NULL;
  int rank = 0, ranks = 1;
  std::vector<NewParticle> additions;
  std::vector<ParticleBinding> pendingBindings;
  std::vector<std::pair<std::string, std::string>> pendingExclusions;
  std::string registering;
  int nextId = 0;
  void check(std::string const& error) const { GrainsExternalModelDetail::check(error, comm); }
  void status(std::string error) const { GrainsExternalModelDetail::status(error, comm); }
  void broadcast(std::string& wire) const { GrainsExternalModelDetail::broadcast(wire, comm); }
  std::vector<std::string> gather(std::string const& wire) const {
    return GrainsExternalModelDetail::gather(wire, comm);
  }
public:
  ModelRegistration(Grains&, ExternalModelData&);
  int add(std::string const&, std::unique_ptr<Particle>) override;
  void bind(std::string const&, int) override;
  void excludeContacts(std::string const&, std::string const&) override;
  void registerModel(std::string const&, std::unique_ptr<ExternalModel>);
};

Grains::ModelRegistration::ModelRegistration(Grains &simulation, ExternalModelData &data)
    : g(simulation), registry(data) {
  if (g.m_wrapper && !g.m_processorIsActive)
    throw std::logic_error("Inactive rank cannot participate in coupling");
  if (g.m_wrapper) {
    comm = g.m_wrapper->get_active_procs_comm();
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &ranks);
  }
  check(!g.m_collision || g.m_dimension != 3 || g.m_periodicity.size() != 3
            ? "External models require an initialized 3D simulation"
            : "");
}

int Grains::ModelRegistration::add(std::string const &attachment, std::unique_ptr<Particle> p) {
  if (registering.empty() || !p || !p->getRigidBody())
    throw std::invalid_argument("Invalid particle registration");
  if (attachment.empty() ||
      std::any_of(pendingBindings.begin(), pendingBindings.end(),
                  [&](ParticleBinding const &b) {
                    return b.attachment == attachment;
                  }))
    throw std::invalid_argument("Duplicate or empty attachment name");
  if (nextId == INT_MAX)
    throw std::overflow_error("Particle ID space exhausted");
  int id = ++nextId;
  p->setID(id);
  p->setKinematicDriver(Driver::External);
  pendingBindings.push_back({registering, attachment, id});
  additions.push_back({attachment, std::move(p)});
  return id;
}

void Grains::ModelRegistration::bind(std::string const &attachment, int id) {
  if (registering.empty() || attachment.empty() ||
      std::any_of(pendingBindings.begin(), pendingBindings.end(),
                  [&](ParticleBinding const &b) {
                    return b.attachment == attachment || b.particleId == id;
                  }))
    throw std::invalid_argument(
        "Invalid or duplicate existing particle binding");
  pendingBindings.push_back({registering, attachment, id});
}

void Grains::ModelRegistration::excludeContacts(std::string const &a, std::string const &b) {
  if (a == b)
    throw std::invalid_argument("A particle cannot exclude itself");
  pendingExclusions.emplace_back(a, b);
}

/**
 * @brief Collectively initialize a model and publish its particles and
 * bindings.
 * @param name Unique, nonempty name, supplied on root.
 * @param instance Model transferred to the registry on root; null on other
 * ranks.
 *
 * Root calls initialize(), validates pending attachments and geometry, and
 * broadcasts particle definitions, materials and exclusions. Every rank
 * prepares reference geometry before adoption through the normal insertion
 * rules. Material changes are restored on preparation failure. Publication
 * failures abort MPI because distributed ownership may already have changed.
 */
void Grains::ModelRegistration::registerModel(
    std::string const &name, std::unique_ptr<ExternalModel> instance) {
  check(rank != 0 && instance ? "Only root may supply an external model" : "");
  int maximum = std::max(g.m_allcomponents.getMaxParticleIDnumber(),
                         Particle::getMaxIDnumber());
  if (comm != MPI_COMM_NULL)
    MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_INT, MPI_MAX, comm);
  nextId = maximum;
  auto oldMaterials = ContactBuilderFactory::captureMaterials();
  auto rollbackMaterials = [&]() {
    ContactBuilderFactory::restoreMaterials(oldMaterials);
  };
  int types =
          static_cast<int>(g.m_allcomponents.getReferenceParticles()->size()),
      low = types, high = types;
  if (comm != MPI_COMM_NULL) {
    MPI_Allreduce(MPI_IN_PLACE, &low, 1, MPI_INT, MPI_MIN, comm);
    MPI_Allreduce(MPI_IN_PLACE, &high, 1, MPI_INT, MPI_MAX, comm);
  }
  check(low != high ? "MPI reference particle tables differ" : "");
  std::ostringstream local = output();
  for (auto *p : *g.m_allcomponents.getActiveParticles())
    if (p->getTag() != 2)
      local << p->getID() << ' ';
  auto physical = gather(local.str());
  std::string error, wire;
  if (rank == 0)
    try {
      if (!instance || name.empty())
        throw std::invalid_argument(
            "Registration requires a name and a model on root");
      for (auto const &model : registry.models)
        if (model.name == name)
          throw std::invalid_argument("Duplicate model name");
      registering = name;
      instance->initialize(*this);
      registering.clear();
      if (pendingBindings.empty())
        throw std::invalid_argument(
            "Model must register at least one particle");
      std::set<int> existing;
      for (auto const &packet : physical) {
        std::istringstream in(packet);
        int id;
        while (in >> id)
          existing.insert(id);
      }
      for (auto const &added : additions)
        existing.insert(added.body->getID());
      std::set<int> bound;
      for (auto const &model : registry.models)
        for (auto const &b : model.bindings)
          bound.insert(b.particleId);
      for (auto const &b : pendingBindings)
        if (!existing.count(b.particleId) || !bound.insert(b.particleId).second)
          throw std::invalid_argument("Unknown or already bound particle ID");
      Vector3 origin, length;
      App::get_origin(origin[X], origin[Y], origin[Z]);
      App::get_size(length[X], length[Y], length[Z]);
      std::ostringstream out = output();
      ContactBuilderFactory::writeMaterials(out);
      out << std::quoted(name) << ' ' << nextId << ' ' << additions.size()
          << ' ';
      for (size_t i = 0; i < additions.size(); ++i) {
        auto *p = additions[i].body.get();
        if (!(p->getCrustThickness() > 0.) ||
            !std::isfinite(p->getCircumscribedRadius()) ||
            !p->getRigidBody()->getTransform()->getBasis().isRotation())
          throw std::invalid_argument("Invalid registered particle geometry");
        for (int d = 0; d < 3; ++d)
          if (!std::isfinite((*p->getPosition())[d]) ||
              (*p->getPosition())[d] < origin[d] ||
              (*p->getPosition())[d] >= origin[d] + length[d] ||
              !std::isfinite((*p->getTranslationalVelocity())[d]) ||
              !std::isfinite((*p->getAngularVelocity())[d]))
            throw std::invalid_argument(
                "Registered particle must start inside the primary domain with "
                "finite velocities");
        p->setGeometricType(types + i);
        p->setActivity(COMPUTE);
        Quaternion q;
        q.setQuaternion(p->getRigidBody()->getTransform()->getBasis());
        p->applyExternalState(*p->getRigidBody()->getTransform(), q,
                              *p->getTranslationalVelocity(),
                              *p->getAngularVelocity());
        p->InitializeForce(false);
        std::ostringstream definition = output();
        p->write(definition);
        out << std::quoted(definition.str()) << ' ';
      }
      out << pendingBindings.size() << ' ';
      std::map<std::string, int> attachments;
      for (auto const &b : pendingBindings) {
        attachments[b.attachment] = b.particleId;
        out << std::quoted(b.attachment) << ' ' << b.particleId << ' ';
      }
      out << pendingExclusions.size() << ' ';
      for (auto const &pair : pendingExclusions) {
        if (!attachments.count(pair.first) || !attachments.count(pair.second))
          throw std::invalid_argument("Unknown contact-exclusion attachment");
        out << attachments.at(pair.first) << ' ' << attachments.at(pair.second)
            << ' ';
      }
      wire = out.str();
      if (wire.size() > static_cast<size_t>(INT_MAX))
        throw std::overflow_error("External-model payload too large");
    } catch (std::exception const &e) {
      error = e.what();
    } catch (...) {
      error = "Model particle registration failed";
    }
  try {
    status(error);
  } catch (...) {
    rollbackMaterials();
    throw;
  }
  broadcast(wire);
  std::vector<std::unique_ptr<Particle>> bodies, refs;
  ExternalModelData::Model model;
  std::set<std::pair<int, int>> exclusions;
  error.clear();
  try {
    std::istringstream in(wire);
    size_t count;
    ContactBuilderFactory::readMaterials(in);
    in >> std::quoted(model.name) >> nextId >> count;
    for (size_t i = 0; i < count; ++i) {
      std::string definition;
      in >> std::quoted(definition);
      std::unique_ptr<Particle> body;
      if (rank == 0)
        body = std::move(additions[i].body);
      else {
        std::istringstream data(definition);
        std::string tag, shape = "none";
        data >> tag;
        if (tag == "<CompositeParticle>") {
          std::string marker;
          data >> marker >> shape;
        }
        body = makeParticle(tag == "<CompositeParticle>", shape);
        body->read(data);
        body->setKinematicDriver(Driver::External);
        body->setActivity(COMPUTE);
        if (!data)
          throw std::runtime_error("Cannot deserialize registered particle");
      }
      ContactBuilderFactory::defineMaterial(body->getMaterial(), false);
      std::unique_ptr<Particle> ref(body->createCloneCopy(false));
      ref->setID(GrainsExec::m_ReferenceParticleDefaultID);
      refs.push_back(std::move(ref));
      bodies.push_back(std::move(body));
    }
    in >> count;
    for (size_t i = 0; i < count; ++i) {
      ParticleBinding b;
      b.model = model.name;
      in >> std::quoted(b.attachment) >> b.particleId;
      model.bindings.push_back(b);
    }
    in >> count;
    for (size_t i = 0; i < count; ++i) {
      int a, b;
      in >> a >> b;
      exclusions.emplace(std::min(a, b), std::max(a, b));
    }
    if (!in)
      throw std::runtime_error("Invalid model registration payload");
    auto &all = g.m_allcomponents;
    all.reserveReferenceParticles(refs.size());
    registry.models.reserve(registry.models.size() + 1);
  } catch (std::exception const &e) {
    error = e.what();
  }
  try {
    check(error);
  } catch (...) {
    rollbackMaterials();
    throw;
  }
  // Publication follows the same force-insertion and ownership rules as normal
  // insertion. No caller retains ownership of adopted component objects.
  try {
    for (auto &ref : refs)
      g.m_allcomponents.adoptReferenceParticle(std::move(ref));
    for (auto &body : bodies)
      g.adoptExternalParticle(std::move(body));
    std::set<int> bound;
    for (auto const &b : model.bindings)
      bound.insert(b.particleId);
    for (auto *p : *g.m_allcomponents.getActiveParticles())
      if (bound.count(p->getID()))
        p->setKinematicDriver(Driver::External);
    Particle::setMaxIDnumber(std::max(Particle::getMaxIDnumber(), nextId));
    for (auto const& pair : exclusions)
      g.m_collision->excludeContactPair(pair.first, pair.second);
    if (rank == 0)
      model.instance = std::move(instance);
    registry.models.push_back(std::move(model));
    registry.contactsReady = false;
    g.m_allcomponents.computeNumberParticles(g.m_wrapper);
  } catch (...) {
    if (comm != MPI_COMM_NULL)
      MPI_Abort(comm, 1);
    throw;
  }
}

/**
 * @brief Collectively register a model after simulation geometry is
 * initialized.
 * @param name Model name, validated on root for uniqueness.
 * @param instance Owned model on root; other active ranks supply null.
 *
 * Creates the registry lazily, rejects reentrant or poisoned use, and guards
 * registration with the busy flag. The flag is cleared on success or exception.
 */
void Grains::registerExternalModel(std::string const &name,
                                   std::unique_ptr<ExternalModel> instance) {
  if (m_externalModels && m_externalModels->busy)
    throw std::logic_error("Cannot register a model during a coupled step");
  if (!m_externalModels)
    m_externalModels = std::make_shared<ExternalModelData>();
  if (m_externalModels->poisoned)
    throw std::logic_error(
        "A failed coupled advance requires a new simulation");
  ModelRegistration registration(*this, *m_externalModels);
  m_externalModels->busy = true;
  try {
    registration.registerModel(name, std::move(instance));
    m_externalModels->busy = false;
  } catch (...) {
    m_externalModels->busy = false;
    throw;
  }
}

