#include "GrainsBuilderFactory.hh"
#include "KinematicsDriver.hh"
#include "AppPRSHydroFT.hh"
#include "HookeContactForceModel.hh"
#include "HookeMemoryContactForceModel.hh"
#include "HertzContactForceModel.hh"
#include "HertzMemoryContactForceModel.hh"
#include "SimpleObstacle.hh"
#include "CompositeObstacle.hh"
#include "Sphere.hh"
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationRegistry.hpp>
#include <xercesc/util/XMLString.hpp>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {
using Source = Component::LoadSource;
void require(bool value, const char* message)
{
  if (!value) throw std::runtime_error(message);
}
void equal(Vector3 const& a, Vector3 const& b, const char* message)
{
  require(Norm(a-b) < 1.e-9 * (1. + Norm(a) + Norm(b)), message);
}

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

RigidBodyWithCrust* sphere()
{
  return new RigidBodyWithCrust(new Sphere(1.), Transform(), false, 0.1);
}

void checkSum(Component& p)
{
  auto loads = p.getLoads();
  Torsor sum(*p.getTorsor()->getPoint());
  for (auto const& load : loads) sum += load;
  equal(*sum.getForce(), *p.getForce(), "source forces must sum to integration total");
  equal(*sum.getTorque(), *p.getTorque(), "source moments must sum to integration total");
}

void sourceLoads()
{
  GrainsExec::m_vgravity = Vector3(0., -9.81, 0.);
  Particle p(sphere(), 2500., "test", 0);
  p.setPosition(Point3(1., 2., 3.));
  p.InitializeForce(true);
  Vector3 gravity = *p.getForce();
  p.addBodyForce(Vector3(1., 2., 3.));
  p.addForce(Point3(2., 2., 3.), Vector3(0., 4., 0.), 0, Source::Contact);
  p.addTorque(Vector3(0., 0., 5.), 0, Source::Contact);
  p.addTorque(Vector3(1., 0., 0.), 0); // Legacy calls remain applied loads.
  AppPRSHydroFT hydro;
  hydro.allocateHydroFT(1);
  std::vector<std::vector<double>> hydroData{{7., 8., 9., 1., 2., 3.}};
  hydro.setHydroFT(&hydroData);
  std::list<Particle*> particles{&p};
  hydro.ComputeForces(0., 0.001, &particles);

  auto contact = p.getContactLoad();
  equal(*contact.getForce(), Vector3(0., 4., 0.), "contact excludes gravity and hydro");
  equal(*contact.getTorque(), Vector3(0., 0., 9.), "contact moment includes lever arm and couple");
  equal(*p.getLoad(Source::Gravity).getForce(), gravity, "gravity contribution");
  equal(*p.getLoad(Source::Hydrodynamic).getForce(), Vector3(7., 8., 9.), "hydro force tagging");
  equal(*p.getLoad(Source::Hydrodynamic).getTorque(), Vector3(1., 2., 3.), "hydro torque tagging");
  equal(*p.getLoad(Source::Applied).getForce(), Vector3(1., 2., 3.), "applied contribution");
  checkSum(p);

  contact.ChangeReferencePoint(Point3(0., 2., 3.));
  equal(*contact.getTorque(), Vector3(0., 0., 13.), "wrench transfer preserves force moment");
  p.setForce(Vector3(10., 20., 30.));
  p.setTorque(Vector3(3., 4., 5.));
  checkSum(p);
  equal(*p.getContactLoad().getForce(), Vector3(0., 4., 0.), "legacy setters preserve contact source");
  auto loads = p.getLoads();
  p.InitializeForce(false);
  p.setLoads(loads);
  checkSum(p);
  equal(*p.getForce(), Vector3(10., 20., 30.), "load snapshot restores total");
  p.InitializeForce(false);
  equal(*p.getContactLoad().getForce(), Vector3Null, "contact reset each force pass");
  equal(*p.getLoad(Source::Hydrodynamic).getTorque(), Vector3Null, "hydro reset each force pass");
  checkSum(p);
  GrainsExec::m_vgravity = Vector3Null;
}

template<class Law> void contactLaw()
{
  Particle a(sphere(), 2500., "test", 0), b(sphere(), 2500., "test", 0);
  a.setID(1); b.setID(2);
  a.setPosition(Point3(0., 0., 0.));
  b.setPosition(Point3(1.99, 0., 0.));
  a.InitializeForce(false); b.InitializeForce(false);
  a.addBodyForce(Vector3(0., 2., 0.));
  a.setTranslationalVelocity(Vector3(0., 0.1, 0.));
  a.setAngularVelocity(Vector3(0., 0., 0.1));
  std::map<std::string, double> parameters{{"kn", 1.e5}, {"en", 0.8},
      {"etat", 1.}, {"muc", 0.3}, {"kr", 0.01}, {"kt", 1.e4},
      {"mur", 0.01}, {"etarpf", 0.1}, {"epst", 1.e-10},
      {"Es", 1.e5}, {"Gs", 2.5e4}};
  Law law(parameters);
  PointContact point(Point3(0.995, 0., 0.), Vector3(-0.01, 0., 0.), -0.01, 1);
  law.computeForces(&a, &b, point, nullptr, 0.001, 1);
  auto ca = a.getContactLoad(), cb = b.getContactLoad();
  require(Norm(*ca.getForce()) > 0., "contact law must populate contact source");
  equal(*ca.getForce(), -*cb.getForce(), "equal and opposite contact reactions");
  equal(*a.getForce() - *ca.getForce(), Vector3(0., 2., 0.), "body force excluded from contact");
  equal(*a.getTorque(), *ca.getTorque(), "contact couples are tagged");
  checkSum(a); checkSum(b);
}

void obstacleLoads()
{
  CompositeObstacle root("root");
  auto* a = new SimpleObstacle("a", sphere(), "test", false, false);
  auto* b = new SimpleObstacle("b", sphere(), "test", false, false);
  root.append(a); root.append(b);
  a->InitializeForce(false); b->InitializeForce(false);
  a->addForce(Point3(1., 0., 0.), Vector3(0., 3., 0.), 2, Source::Contact);
  a->addTorque(Vector3(0., 0., 4.), 2, Source::Contact);
  equal(*a->getContactLoad().getForce(), Vector3Null, "ignore ghost reaction on obstacle");
  a->addForce(Point3(1., 0., 0.), Vector3(0., 3., 0.), 0, Source::Contact);
  b->addTorque(Vector3(0., 0., 4.), 0, Source::Contact);
  Torsor expected(*root.getPosition());
  expected += a->getContactLoad(); expected += b->getContactLoad();
  auto actual = root.getContactLoad();
  equal(*actual.getForce(), *expected.getForce(), "composite obstacle source aggregation");
  equal(*actual.getTorque(), *expected.getTorque(), "composite obstacle moment aggregation");
  root.InitializeForce(false);
  equal(*root.getContactLoad().getTorque(), Vector3Null, "reset composite obstacle sources");
}

void driverParity(const char* scheme)
{
  GrainsExec::m_TIScheme = scheme;
  GrainsExec::m_time_counter = 1;
  Particle direct(sphere(), 2500., "test", 0);
  direct.setTranslationalVelocity(Vector3(0.1, 0.2, 0.3));
  direct.setAngularVelocity(Vector3(0.01, 0.02, 0.03));
  direct.InitializeForce(false);
  direct.addBodyForce(Vector3(1., 2., 3.));
  direct.addTorque(Vector3(0.1, 0.2, 0.3), 0);
  Particle driven(direct, false), external(direct, false), ghost(direct, false);
  driven.InitializeForce(false); external.InitializeForce(false); ghost.InitializeForce(false);
  driven.setLoads(direct.getLoads());
  external.setKinematicDriver(Particle::KinematicDriver::External);
  ghost.setTag(2);
  Point3 externalPosition = *external.getPosition();
  Vector3 externalVelocity = *external.getTranslationalVelocity();
  RigidBodyKinematicsDriver rigid;
  KinematicsDriver const& driver = rigid;
  std::list<Particle*> particles{&driven, &external, &ghost};
  for (int step = 0; step < 2; ++step) {
    GrainsExec::m_time_counter = step + 1;
    double dt = 0.001;
    double dtVelocity = std::string(scheme) == "SecondOrderLeapFrog" ? dt/2. : dt;
    direct.Move(step*dt, dtVelocity, dt);
    driver.move(particles, step*dt, dtVelocity, dt);
    if (std::string(scheme) == "SecondOrderLeapFrog") {
      direct.advanceVelocity(step*dt, dt/2.);
      driver.advanceVelocity(particles, step*dt, dt/2.);
    }
    require(direct.getPosition()->DistanceTo(*driven.getPosition()) < 1.e-12,
            "driver must preserve existing integration positions");
    equal(*direct.getTranslationalVelocity(), *driven.getTranslationalVelocity(), "driver velocity parity");
    equal(*direct.getAngularVelocity(), *driven.getAngularVelocity(), "driver angular velocity parity");
    for (int d = 0; d < 4; ++d)
      require(std::abs((*direct.getQuaternionRotation())[d] -
                       (*driven.getQuaternionRotation())[d]) < 1.e-12, "driver orientation parity");
  }
  require(externalPosition.DistanceTo(*external.getPosition()) == 0. &&
          externalPosition.DistanceTo(*ghost.getPosition()) == 0., "driver skips external and ghost poses");
  equal(*external.getTranslationalVelocity(), externalVelocity, "driver skips external velocity kick");
  equal(*ghost.getTranslationalVelocity(), externalVelocity, "driver skips ghost velocity kick");
}
} // namespace

int main()
{
  try {
    auto ctx = context();
    sourceLoads();
    contactLaw<HookeContactForceModel>();
    contactLaw<HookeMemoryContactForceModel>();
    contactLaw<HertzContactForceModel>();
    contactLaw<HertzMemoryContactForceModel>();
    obstacleLoads();
    for (const char* scheme : {"FirstOrderExplicit", "SecondOrderExplicit",
                              "SecondOrderLeapFrog", "SecondOrderAdamsBashforth"})
      driverParity(scheme);
  } catch (std::exception const& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
  return 0;
}
