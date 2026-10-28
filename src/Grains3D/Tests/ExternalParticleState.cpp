#include "GrainsBuilderFactory.hh"
#include "GrainsMPI.hh"
#include "Sphere.hh"
#include "CompositeParticle.hh"
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/util/XMLString.hpp>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>

namespace {
void require(bool value, const char* message)
{
  if (!value) throw std::runtime_error(message);
}
bool near(double a, double b) { return std::abs(a-b) < 1.e-10; }

std::unique_ptr<Grains> context()
{
  ReaderXML::initialize();
  auto* core = xercesc::XMLString::transcode("Core");
  auto* root = xercesc::XMLString::transcode("Grains3D");
  auto* type = xercesc::XMLString::transcode("Type");
  auto* standard = xercesc::XMLString::transcode("Standard");
  auto* impl = xercesc::DOMImplementationRegistry::getDOMImplementation(core);
  auto* doc = impl->createDocument(nullptr, root, nullptr);
  doc->getDocumentElement()->setAttribute(type, standard);
  std::unique_ptr<Grains> result(GrainsBuilderFactory::create(doc->getDocumentElement()));
  doc->release();
  xercesc::XMLString::release(&core);
  xercesc::XMLString::release(&root);
  xercesc::XMLString::release(&type);
  xercesc::XMLString::release(&standard);
  ReaderXML::terminate();
  return result;
}

struct Domain : App {
  static void configure(bool periodic) {
    m_domain_global_periodic_vectors.clear();
    m_periodic_vector_indices.clear();
    set_dimensions(10., 10., 10., 0., 0., 0.);
    set_periodicity({periodic, false, false});
  }
};

class InspectableComposite : public CompositeParticle {
public:
  InspectableComposite(DOMElement* root) : CompositeParticle(root, 1) {}
  InspectableComposite(InspectableComposite const& other) : CompositeParticle(other, false) {}
  Particle* child(size_t i) { return m_elementaryParticles.at(i); }
};

class Simulation : public GrainsMPI {
public:
  Simulation(bool parallel, bool periodic) {
    Domain::configure(periodic);
    m_periodicity.assign(3, false);
    m_periodicity[X] = periodic;
    m_periodic = periodic;
    m_dt = 0.001;
    m_time = 0.;
    m_collision = new LinkedCell;
    m_allApp.push_back(m_collision);
    if (parallel) {
      m_wrapper = new GrainsMPIWrapper(2, 1, 1, periodic, 0, 0, "");
      GrainsExec::m_MPI = true;
      GrainsExec::setComm(m_wrapper);
      m_rank = m_wrapper->get_rank();
      m_nprocs = 2;
      if (periodic) m_wrapper->setMPIperiodicVectors(10., 10., 10.);
      App::set_local_domain_size(5., 10., 10.);
      App::set_local_domain_origin(m_wrapper->get_nb_procs_direction(),
                                   m_wrapper->get_MPI_coordinates());
      require(!m_collision->set(1., m_wrapper->get_nb_procs_direction(),
                  m_wrapper->get_MPI_coordinates(), m_wrapper->get_MPI_neighbors(), ""),
              "MPI grid initialization");
    } else {
      require(!m_collision->set(1., ""), "serial grid initialization");
    }
    m_collision->Link(m_allcomponents.getObstacles());
    SCT_insert_app("Move");
    SCT_insert_app("LinkUpdate");
    SCT_insert_app("UpdateParticleActivity");
    auto* ref = new Particle(new RigidBodyWithCrust(new Sphere(0.1),
                    Transform(), false, 0.02), 2500., "test", 0);
    ref->setKinematicDriver(Particle::KinematicDriver::External);
    m_allcomponents.AddReferenceParticle(ref, 0);
  }
  ~Simulation() { GrainsExec::m_MPI = false; GrainsExec::setComm(nullptr); }

  void add(int id, double x) {
    auto* ref = m_allcomponents.getReferenceParticles()->front();
    Particle* p = ref->createCloneCopy(false);
    p->setID(id);
    p->setPosition(Point3(x, 5., 5.));
    p->setActivity(COMPUTE);
    p->InitializeForce(false);
    p->addBodyForce(Vector3(0.25, 0.5, 0.75), Component::LoadSource::Contact);
    p->addBodyForce(Vector3(0.25, 0.5, 0.75));
    p->addBodyForce(Vector3(0.5, 1., 1.5), Component::LoadSource::Hydrodynamic);
    p->addTorque(Vector3(1., 2., 3.), 0, Component::LoadSource::Contact);
    p->addTorque(Vector3(3., 3., 3.), 0);
    std::map<std::tuple<int,int,int>,
        std::tuple<bool, Vector3, Vector3, Vector3>> history;
    history.emplace(std::make_tuple(42, 0, 0), std::make_tuple(true,
        Vector3(0.1, 0., 0.), Vector3(0., 1., 0.), Vector3(0., 0., 0.2)));
    p->setContactMap(history);
    bool keep = true;
    if (m_wrapper) {
      auto inserted = m_collision->insertParticleParallel(m_time, p,
          m_allcomponents.getActiveParticles(), m_allcomponents.getCloneParticles(),
          m_allcomponents.getReferenceParticles(), m_periodic, true, m_wrapper);
      keep = inserted.first;
    } else {
      m_collision->insertParticleSerial(p, m_allcomponents.getActiveParticles(),
          m_allcomponents.getPeriodicCloneParticles(),
          m_allcomponents.getReferenceParticles(), m_periodic, true);
    }
    if (!keep) { delete p; return; }
    m_allcomponents.getActiveParticles()->push_back(p);
    if (m_wrapper && p->getTag() == 1)
      m_allcomponents.getParticlesInBufferzone()->push_back(p);
    if (m_wrapper && p->getTag() == 2)
      m_allcomponents.getCloneParticles()->push_back(p);
  }
  void submit(std::vector<ExternalParticleState> states) {
    if (m_rank != 0) states.clear();
    importExternalParticleStates(states);
  }
  InspectableComposite* addComposite(DOMElement* root) {
    auto* ref = new InspectableComposite(root);
    ref->setKinematicDriver(Particle::KinematicDriver::External);
    m_allcomponents.AddReferenceParticle(ref, 0);
    auto* p = new InspectableComposite(*ref);
    p->setID(3);
    p->setPosition(Point3(5., 5., 5.));
    p->setActivity(COMPUTE);
    p->InitializeForce(false);
    m_collision->Link(p);
    m_allcomponents.getActiveParticles()->push_back(p);
    return p;
  }
  void reject(std::vector<ExternalParticleState> states) {
    bool rejected = false;
    try { submit(states); } catch (std::invalid_argument const&) { rejected = true; }
    require(rejected, "invalid batch must be rejected on every rank");
  }
  bool contactAt(double x) {
    Particle probe(new RigidBodyWithCrust(new Sphere(0.1), Transform(), false, 0.02),
                   2500., "test", 0);
    probe.setPosition(Point3(x, 5., 5.));
    return m_collision->isContactWithCrust(&probe);
  }
  void driver(int id, Particle::KinematicDriver driver) {
    for (auto* p : *m_allcomponents.getActiveParticles())
      if (p->getID() == id) p->setKinematicDriver(driver);
  }
  void checkInternalMotion() {
    auto* p = m_allcomponents.getActiveParticles()->front();
    p->setKinematicDriver(Particle::KinematicDriver::Grains);
    Point3 before = *p->getPosition();
    moveParticlesAndObstacles(0.0005, 0.001, 0.001);
    m_allcomponents.advanceParticlesVelocity(0., 0.0005);
    require((*p->getPosition())[X] > before[X] &&
            (*p->getTranslationalVelocity())[X] > 0.,
            "internal motion and leapfrog kick survive bookkeeping extraction");
  }
  void rejectNonRootInput() {
    bool rejected = false;
    std::vector<ExternalParticleState> batch;
    if (m_rank != 0) batch.push_back({1, Point3(), Quaternion(), Vector3(), Vector3()});
    try { importExternalParticleStates(batch); }
    catch (std::invalid_argument const&) { rejected = true; }
    require(rejected, "nonroot submission must fail collectively");
  }
  void checkObstacleLoadReduction() {
    SimpleObstacle obstacle("test", new RigidBodyWithCrust(new Sphere(0.1),
        Transform(), false, 0.02), "test", false, false);
    obstacle.InitializeForce(false);
    obstacle.addForce(Point3(1., 0., 0.), Vector3(0., m_rank+1., 0.),
                      0, Component::LoadSource::Contact);
    obstacle.addBodyForce(Vector3(1., 0., 0.), Component::LoadSource::Hydrodynamic);
    obstacle.addTorque(Vector3(0., 0., 1.), 0);
    m_wrapper->sumObstaclesLoad({&obstacle});
    auto contact = obstacle.getContactLoad();
    require(near((*contact.getForce())[Y], 3.) && near((*contact.getTorque())[Z], 3.),
            "reduce obstacle contact wrench independently");
    require(near((*obstacle.getLoad(Component::LoadSource::Hydrodynamic).getForce())[X], 2.) &&
            near((*obstacle.getLoad(Component::LoadSource::Applied).getTorque())[Z], 2.),
            "reduce other obstacle sources independently");
    require(near((*obstacle.getForce())[Y], 3.) && near((*obstacle.getTorque())[Z], 5.),
            "obstacle source reduction preserves original totals");
  }
  void verify(int id, double x, double vx = 2.) {
    int owners = 0;
    for (auto* p : *m_allcomponents.getActiveParticles()) {
      if (p->getID() != id) continue;
      require(!p->isIntegratedByGrains(), "external driver survives clone reconstruction");
      require(near((*p->getTranslationalVelocity())[X], vx), "image velocity refresh");
      require(near((*p->getAngularVelocity())[Z], 3.), "image angular velocity refresh");
      require(near((*p->getForce())[Y], 2.) && near((*p->getTorque())[Z], 6.),
              "preserve resultant loads through image/owner changes");
      auto contactLoad = p->getContactLoad();
      require(near((*contactLoad.getForce())[Y], 0.5) &&
              near((*contactLoad.getTorque())[Z], 3.),
              "preserve contact-only wrench through image/owner changes");
      require(near((*p->getLoad(Component::LoadSource::Hydrodynamic).getForce())[Y], 1.),
              "preserve hydrodynamic source through image/owner changes");
      require(p->getContactMap()->count(std::make_tuple(42, 0, 0)) == 1,
              "preserve contact history through import");
      if (p->getTag() != 2) {
        ++owners;
        require(near((*p->getPosition())[X], x), "physical owner position");
        require(p->getCell() == m_collision->getCell(*p->getPosition()), "cell relinking");
      }
    }
    if (m_wrapper) MPI_Allreduce(MPI_IN_PLACE, &owners, 1, MPI_INT, MPI_SUM,
                                 m_wrapper->get_active_procs_comm());
    if (owners != 1) {
      std::cerr << "Expected x=" << x << " owners=" << owners << '\n';
      for (auto* p : *m_allcomponents.getActiveParticles())
        std::cerr << "id=" << p->getID() << " x=" << (*p->getPosition())[X]
                  << " tag=" << p->getTag() << '\n';
    }
    require(owners == 1, "one physical owner after synchronization");
    require(near(m_time, 0.), "import must not advance simulation time");
  }
protected:
  void synchronizeParticleStates() override {
    if (m_wrapper) GrainsMPI::synchronizeParticleStates();
    else Grains::synchronizeParticleStates();
  }
};

ExternalParticleState state(int id, double x)
{
  return {id, Point3(x, 5., 5.), Quaternion(), Vector3(2., 0., 0.), Vector3(0., 0., 3.)};
}

void validation()
{
  Simulation sim(false, false);
  sim.add(1, 4.995);
  sim.add(2, 7.);
  require(!sim.contactAt(5.2), "initial probe must be outside contact");
  auto a = state(1, 5.005), b = state(2, 7.);
  sim.reject({a});
  sim.reject({a, a, b});
  sim.reject({state(99, 5.), b});
  auto bad = a;
  bad.orientation = Quaternion(0., 0., 0., 0.);
  sim.reject({bad, b});
  bad = a;
  bad.translationalVelocity[X] = std::numeric_limits<double>::infinity();
  sim.reject({bad, b});
  sim.reject({state(1, 6.), b});
  bad = a;
  bad.orientation = Quaternion(0., 0., std::sin(0.5), std::cos(0.5));
  sim.reject({bad, b});
  // A failed batch must not change the history used by the next valid import.
  a.orientation = Quaternion(0., 0., 0., -2.);
  sim.submit({a, b});
  sim.verify(1, 5.005);
  sim.verify(2, 7.);
  require(sim.contactAt(5.2), "contact detection must use imported cell/geometry");
  sim.driver(2, Particle::KinematicDriver::Grains);
  sim.reject({a, b});
}

void crossing(bool parallel, bool periodic)
{
  Simulation sim(parallel, periodic);
  double start = periodic ? 9.995 : 4.995;
  sim.add(1, start);
  if (parallel) sim.rejectNonRootInput();
  sim.submit({state(1, start)});
  sim.submit({state(1, start + 0.01)});
  sim.verify(1, periodic ? 0.005 : 5.005);
  sim.reject({state(1, start + 10.01)}); // An unwrapped full lap must not alias.
  sim.submit({state(1, start)});
  sim.verify(1, start);
  if (parallel) sim.checkObstacleLoadReduction();
}

void composite()
{
  Simulation sim(false, false);
  ReaderXML::initialize();
  auto* p = sim.addComposite(ReaderXML::getRoot(GRAINS_COMPOSITE_FIXTURE));
  ReaderXML::terminate();
  auto s = state(3, 5.001);
  s.orientation = Quaternion(0., 0., std::sin(0.005), std::cos(0.005));
  // Force construction of the old transformed geometry before invalidation.
  Point3 oldOrigin = *p->child(0)->getRigidBody()->getTransformWithCrust()->getOrigin();
  sim.submit({s});
  Matrix rotation(s.orientation);
  for (size_t i = 0; i < 2; ++i) {
    auto* child = p->child(i);
    Vector3 offset = rotation * Vector3(i == 0 ? -0.2 : 0.2, 0., 0.);
    require(near(child->getPosition()->DistanceTo(s.position + offset), 0.),
            "composite child position");
    require(near(child->getPosition()->DistanceTo(
            *child->getRigidBody()->getTransformWithCrust()->getOrigin()), 0.),
            "invalidate the child's previously computed crust transform");
    Vector3 v = s.translationalVelocity + (s.angularVelocity ^ offset);
    require(near(Norm(*child->getTranslationalVelocity() - v), 0.),
            "composite child rigid-body velocity");
    Matrix qbasis(*child->getQuaternionRotation());
    Vector3 axis(1., 0., 0.);
    require(near(Norm(qbasis * axis -
            child->getRigidBody()->getTransform()->getBasis() * axis), 0.),
            "composite child orientation agrees with its own transform");
  }
  require(near((*p->getAngularVelocity())[Z], 3.), "parent angular velocity");
  require(!near(oldOrigin[X], (*p->child(0)->getPosition())[X]), "test moved cached geometry");
  sim.checkInternalMotion();
  Vector3 offset = p->getRigidBody()->getTransform()->getBasis() * Vector3(-0.2, 0., 0.);
  require(near(p->child(0)->getPosition()->DistanceTo(*p->getPosition() + offset), 0.),
          "rigid-body driver preserves composite virtual motion dispatch");
}

void internalMotionAndBoundary()
{
  {
    Simulation sim(false, false);
    sim.add(1, 4.);
    sim.checkInternalMotion();
  }
  {
    Simulation sim(false, false);
    sim.add(1, 0.005);
    sim.reject({state(1, -0.005)});
    sim.submit({state(1, 0.005)});
    sim.verify(1, 0.005);
  }
}
} // namespace

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  try {
    auto ctx = context();
    bool parallel = argc > 1;
    if (parallel) {
      // The existing MPI wrapper has process-global logging lifetime: use one
      // simulation per MPI test process, matching the executable's lifecycle.
      crossing(true, std::string(argv[1]) == "--mpi-periodic");
    } else {
      validation();
      composite();
      internalMotionAndBoundary();
      crossing(false, false);
      crossing(false, true);
    }
  } catch (std::exception const& error) {
    std::cerr << error.what() << '\n';
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  MPI_Finalize();
  return 0;
}
