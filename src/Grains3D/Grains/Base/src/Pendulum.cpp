#include "Pendulum.hh"
#include "SpheroCylinder.hh"
#include "ContactBuilderFactory.hh"
#include "GrainsExec.hh"
#include "BendingRodImplicit.hh"
#include <algorithm>
#include <climits>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <xercesc/util/XMLString.hpp>

namespace {
std::string trim(std::string value)
{
  auto first=value.find_first_not_of(" \t\r\n");
  if (first==std::string::npos) return {};
  return value.substr(first,value.find_last_not_of(" \t\r\n")-first+1);
}
std::string required(DOMNode* node,char const* attribute)
{
  if (!ReaderXML::hasNodeAttr(node,attribute))
    throw std::invalid_argument(std::string("Pendulum requires attribute ")+attribute);
  auto value=trim(ReaderXML::getNodeAttr_String(node,attribute));
  if (value.empty()) throw std::invalid_argument("Pendulum attributes cannot be empty");
  return value;
}
double number(std::string const& value,bool positive)
{
  try {
    size_t consumed=0;
    double result=std::stod(value,&consumed);
    if (consumed==value.size() && std::isfinite(result) && (!positive || result>0.)) return result;
  } catch (std::exception const&) {}
  throw std::invalid_argument("Pendulum requires finite numbers and positive physical properties");
}
double attribute(DOMNode* node,char const* name)
{ return number(required(node,name),false); }
std::vector<DOMNode*> elements(DOMNode* node)
{
  std::vector<DOMNode*> result;
  for (auto* child=node->getFirstChild();child;child=child->getNextSibling())
    if (child->getNodeType()==DOMNode::ELEMENT_NODE) result.push_back(child);
  return result;
}
std::string text(DOMNode* node)
{
  if (!elements(node).empty()) throw std::invalid_argument("Pendulum properties must contain text only");
  auto* value=xercesc::XMLString::transcode(node->getTextContent());
  auto result=trim(value?value:""); xercesc::XMLString::release(&value);
  if (result.empty()) throw std::invalid_argument("Pendulum properties cannot be empty");
  return result;
}
}

struct Pendulum::Dynamics {
  explicit Dynamics(GrainsRod::BendingRodImplicitConfig const& config) : rod(config) {}
  GrainsRod::BendingRodImplicit rod;
  std::unordered_map<int,size_t> sites;
};

Pendulum::~Pendulum() = default;

Pendulum::Pendulum(DOMNode* node)
{
  if (!node || ReaderXML::getNodeName(node)!="Pendulum")
    throw std::invalid_argument("Expected Pendulum XML element");
  required(node,"Name");
  auto count=required(node,"N");
  if (count.find_first_not_of("0123456789")!=std::string::npos)
    throw std::invalid_argument("Pendulum N must be a positive integer");
  try {
    auto value=std::stoull(count);
    if (value==0 || value>INT_MAX/6-1) throw std::invalid_argument("Pendulum N exceeds solver capacity");
    m_count=static_cast<size_t>(value);
  } catch (std::exception const&) { throw std::invalid_argument("Pendulum N must be in [1, INT_MAX/6-1]"); }
  m_length=number(required(node,"L"),true);
  auto fields=elements(node);
  std::vector<std::string> names{"Anchor","Direction","Density","CrustThickness","Radius","Material"};
  if (fields.size()<names.size()) throw std::invalid_argument("Pendulum requires Anchor, Direction, Density, CrustThickness, Radius, Material");
  for (size_t i=0;i<names.size();++i)
    if (ReaderXML::getNodeName(fields[i])!=names[i]) throw std::invalid_argument("Missing, duplicate, or out-of-order Pendulum property");
  m_anchor=Point3(attribute(fields[0],"X"),attribute(fields[0],"Y"),attribute(fields[0],"Z"));
  m_direction=Vector3(attribute(fields[1],"X"),attribute(fields[1],"Y"),attribute(fields[1],"Z"));
  double scale=std::max({std::abs(m_direction[X]),std::abs(m_direction[Y]),std::abs(m_direction[Z])});
  if (scale==0.) throw std::invalid_argument("Pendulum Direction cannot be zero");
  m_direction/=scale; m_direction/=Norm(m_direction);
  m_density=number(text(fields[2]),true); m_crustThickness=number(text(fields[3]),true);
  m_radius=number(text(fields[4]),true); m_material=text(fields[5]);
  size_t field=names.size();
  for (auto property : {std::make_pair("BendingRigidity",&m_bendingRigidity),
                        std::make_pair("DampingRigidity",&m_dampingRigidity)}) {
    if (field<fields.size() && ReaderXML::getNodeName(fields[field])==property.first) {
      *property.second=number(text(fields[field++]),false);
      if (*property.second<0.) throw std::invalid_argument("Pendulum rigidities must be nonnegative");
    }
  }
  if (field!=fields.size()) throw std::invalid_argument("Unknown, duplicate, or out-of-order Pendulum property");
  Point3 endpoint=m_anchor+m_direction*m_length;
  for (int d=0;d<3;++d) if (!std::isfinite(endpoint[d])) throw std::invalid_argument("Pendulum endpoint is not finite");
  // Construction needs material IDs and geometry bounds before the grid exists.
  // Actual particles are generated only during model registration.
  ContactBuilderFactory::defineMaterial(m_material,false);
  GrainsExec::setMinCrustThickness(m_crustThickness);
}

void Pendulum::initialize(Registration& registration)
{
  if (m_initialized) throw std::logic_error("Pendulum is already initialized");
  GrainsRod::BendingRodImplicitConfig config;
  config.links=static_cast<int>(m_count); config.length=m_length;
  // Uniform slender rod mass/inertia, independent of contact-capsule overlap.
  config.total_mass=m_density*std::acos(-1.)*m_radius*m_radius*m_length;
  config.pin=GrainsRod::Vec3(m_anchor[X],m_anchor[Y],m_anchor[Z]);
  config.ei=m_bendingRigidity; config.damping_rigidity=m_dampingRigidity;
  m_dynamics.reset(new Dynamics(config));
  auto& rod=m_dynamics->rod;
  rod.setState(rod.initialState({GrainsRod::Vec3(m_direction[X],m_direction[Y],m_direction[Z])}));
  std::vector<std::string> attachments;
  for (size_t i=0;i<m_count;++i) {
    auto name="site-"+std::to_string(i+1);
    Point3 position=m_anchor+m_direction*(m_length*((static_cast<double>(i)+0.5)/m_count));
    std::unique_ptr<RigidBodyWithCrust> body(new RigidBodyWithCrust(new SpheroCylinder(m_radius,m_length/m_count),Transform(),false,m_crustThickness));
    std::unique_ptr<Particle> particle(new Particle(body.get(),m_density,m_material,0)); body.release();
    // Grains spherocylinders have their cylinder axis along local +Y.
    Eigen::Quaterniond rotation=Eigen::Quaterniond::FromTwoVectors(
        GrainsRod::Vec3::UnitY(),GrainsRod::Vec3(m_direction[X],m_direction[Y],m_direction[Z]));
    Quaternion orientation(rotation.x(),rotation.y(),rotation.z(),rotation.w());
    Transform pose; pose.setOrigin(position); pose.setBasis(Matrix(orientation));
    particle->applyExternalState(pose,orientation,Vector3Null,Vector3Null);
    particle->setKinematicDriver(Particle::KinematicDriver::External);
    particle->setTranslationalVelocity(Vector3Null); particle->setAngularVelocity(Vector3Null);
    ExternalParticleState state{0,position,*particle->getQuaternionRotation(),Vector3Null,Vector3Null};
    state.id=registration.add(name,std::move(particle));
    m_dynamics->sites.emplace(state.id,i);
    m_states.push_back(state); attachments.push_back(name);
  }
  for (size_t i=0;i<attachments.size();++i)
    for (size_t j=i+1;j<attachments.size();++j) registration.excludeContacts(attachments[i],attachments[j]);
  m_initialized=true;
}

std::vector<ExternalParticleState> Pendulum::advance(double time,double dt,std::vector<ContactLoad> const& loads)
{
  if (!m_initialized) throw std::logic_error("Initialize Pendulum before advancing it");
  if (!std::isfinite(time)) throw std::invalid_argument("Pendulum time must be finite");
  auto& rod=m_dynamics->rod;
  using GrainsRod::Vec3;
  auto masses=rod.nodal_masses();
  Vec3 gravity(GrainsExec::m_vgravity[X],GrainsExec::m_vgravity[Y],GrainsExec::m_vgravity[Z]);
  std::vector<Vec3> forces(m_count+1),torques(m_count,Vec3::Zero());
  for (size_t i=0;i<forces.size();++i) forces[i]=masses[i]*gravity;
  std::vector<bool> seen(m_count,false);
  for (auto const& load:loads) {
    auto site=m_dynamics->sites.find(load.particleId);
    if (site==m_dynamics->sites.end() || seen[site->second])
      throw std::invalid_argument("Pendulum received an unknown or duplicate particle load");
    size_t i=site->second; seen[i]=true;
    auto const& f=*load.wrench.getForce(); auto const& t=*load.wrench.getTorque();
    // The wrench is reduced at the link midpoint. Equal endpoint forces
    // preserve its resultant, moment and virtual work at every configuration.
    Vec3 halfForce=0.5*Vec3(f[X],f[Y],f[Z]);
    forces[i]+=halfForce; forces[i+1]+=halfForce;
    torques[i]=Vec3(t[X],t[Y],t[Z]);
  }
  auto before=rod.unpack(rod.state());
  rod.advance(dt,forces,torques);
  if (dt==0.) return m_states;
  auto after=rod.unpack(rod.state());
  auto positions=rod.nodes(rod.state()),velocities=rod.velocities(rod.state());
  for (size_t i=0;i<m_count;++i) {
    auto& state=m_states[i];
    for (int d=0;d<3;++d) {
      state.position[d]=0.5*(positions(i,d)+positions(i+1,d));
      state.translationalVelocity[d]=0.5*(velocities(i,d)+velocities(i+1,d));
      state.angularVelocity[d]=after.omega[i][d];
    }
    // Transport the capsule frame with its link, without introducing axial spin.
    auto const& q=state.orientation;
    Eigen::Quaterniond rotation=Eigen::Quaterniond::FromTwoVectors(before.tangents[i],after.tangents[i])
        * Eigen::Quaterniond(q[3],q[0],q[1],q[2]);
    rotation.normalize();
    state.orientation.setQuaternion(rotation.x(),rotation.y(),rotation.z(),rotation.w());
  }
  return m_states;
}
