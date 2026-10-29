#include "ExternalModel.hh"
#include "Pendulum.hh"
#include "Grains.hh"
#include "CompositeObstacle.hh"
#include "TrilobeCylinder.hh"
#include "QuadrilobeCylinder.hh"
#include "Dendrite.hh"
#include "ContactBuilderFactory.hh"
#include "KinematicsDriver.hh"
#include <algorithm>
#include <cmath>
#include <climits>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace {
using State = ExternalParticleState;
struct StepData {
  std::vector<State> particles;
  std::vector<ExternalModel::ContactLoad> contacts;
};
struct StepTiming { double time,dt; };
struct ParticleBinding {
  std::string model;
  std::string attachment;
  int particleId;
};
using Driver = Particle::KinematicDriver;
using ContactMap = std::map<std::tuple<int,int,int>,
    std::tuple<bool, Vector3, Vector3, Vector3>>;

double wrapped(double x, double origin, double length)
{
  double d = std::fmod(x-origin, length);
  if (d < 0.) d += length;
  return origin + (d < length ? d : 0.);
}
void vectorOut(std::ostream& out, solid::Group3 const& v)
{ for (int d=0; d<3; ++d) out << v[d] << ' '; }
Vector3 vectorIn(std::istream& in)
{ double x,y,z; in >> x >> y >> z; return Vector3(x,y,z); }
void stateOut(std::ostream& out, State const& s)
{
  out << s.id << ' '; vectorOut(out,s.position);
  for (int d=0; d<4; ++d) out << s.orientation[d] << ' ';
  vectorOut(out,s.translationalVelocity); vectorOut(out,s.angularVelocity);
}
State stateIn(std::istream& in)
{
  State s; in >> s.id; s.position=Point3(vectorIn(in));
  double q[4]; for (double& v:q) in >> v;
  s.orientation=Quaternion(q[0],q[1],q[2],q[3]);
  s.translationalVelocity=vectorIn(in); s.angularVelocity=vectorIn(in); return s;
}
std::unique_ptr<Particle> makeParticle(bool composite, std::string const& shape)
{
  if (!composite) return std::unique_ptr<Particle>(new Particle(false));
  if (shape=="TrilobeCylinder") return std::unique_ptr<Particle>(new TrilobeCylinder(false));
  if (shape=="QuadrilobeCylinder") return std::unique_ptr<Particle>(new QuadrilobeCylinder(false));
  if (shape=="Dendrite") return std::unique_ptr<Particle>(new Dendrite(false));
  if (shape!="none") throw std::invalid_argument("Unsupported serialized composite type");
  return std::unique_ptr<Particle>(new CompositeParticle(false));
}
class WorkingDiagnostics : public LinkedCell {
public:
  void addToContactsFeatures(double, PointContact const&) override {}
};
// Working calculations run on root without publishing diagnostics or output.
struct WorkingGlobals {
  bool mpi=GrainsExec::m_MPI, pp=GrainsExec::m_postprocess_forces_at_this_time;
  bool contact=GrainsExec::m_exception_Contact, motion=GrainsExec::m_exception_Motion;
  bool simulation=GrainsExec::m_exception_Simulation;
  unsigned long long gjk=GrainsExec::m_nb_GJK_calls;
  unsigned long long narrow=GrainsExec::m_nb_GJK_narrow_collision_detections;
  size_t counter=GrainsExec::m_time_counter;
  WorkingGlobals() { GrainsExec::m_MPI=false; GrainsExec::m_postprocess_forces_at_this_time=false; }
  ~WorkingGlobals() {
    GrainsExec::m_MPI=mpi; GrainsExec::m_postprocess_forces_at_this_time=pp;
    GrainsExec::m_exception_Contact=contact; GrainsExec::m_exception_Motion=motion;
    GrainsExec::m_exception_Simulation=simulation;
    GrainsExec::m_nb_GJK_calls=gjk; GrainsExec::m_nb_GJK_narrow_collision_detections=narrow;
    GrainsExec::m_time_counter=counter;
  }
};
} // namespace

class ExternalModelRegistry {
public:
  struct Model {
    std::string name;
    std::unique_ptr<ExternalModel> instance;
    std::vector<ParticleBinding> bindings;
  };
  std::vector<Model> models;
  std::set<std::pair<int,int>> excluded;
  std::map<int,State> unwrapped;
  bool busy=false;
  bool poisoned=false;
  bool contactsReady=false;
};

class ExternalModelSession final : public ExternalModel::Registration {
  using Snapshot = Component::InteractionState;
  struct StepState {
    StepData state;
    std::vector<std::unique_ptr<Particle>> particles;
    Snapshot obstacles;
  };
  struct NewParticle { std::string attachment; std::unique_ptr<Particle> body; };
  Grains& g;
  ExternalModelRegistry& registry;
  MPI_Comm comm=MPI_COMM_NULL;
  int rank=0, ranks=1;
  StepTiming context;
  std::unique_ptr<StepState> working;
  std::unique_ptr<StepState> accepted;
  std::map<int,int> owners;
  std::vector<NewParticle> additions;
  std::vector<ParticleBinding> pendingBindings;
  std::vector<std::pair<std::string,std::string>> pendingExclusions;
  std::string registering;
  int nextId=0;

  static Snapshot capture(Component* c) { return c->captureInteractionState(); }
  static void obstacleHistory(Component* c, bool finish) {
    if (dynamic_cast<CompositeObstacle*>(c)) {
      for (auto* child:static_cast<CompositeObstacle*>(c)->getObstacles()) obstacleHistory(child,finish);
    } else if (finish) c->updateContactMap();
    else c->setContactMapToFalse();
  }
  static void initializeObstacleLoads(Component* c) {
    if (dynamic_cast<CompositeObstacle*>(c)) {
      c->Component::InitializeForce(false);
      for (auto* child:static_cast<CompositeObstacle*>(c)->getObstacles()) initializeObstacleLoads(child);
    } else {
      auto loads=c->getLoads();
      loads[static_cast<size_t>(Component::LoadSource::Contact)]=Torsor(*c->getPosition());
      c->InitializeForce(false); c->setLoads(loads);
    }
  }
  static void restore(Component* c, Snapshot const& s, bool reanchor=false) {
    c->restoreInteractionState(s,reanchor);
  }
  static void snapshotOut(std::ostream& out, Snapshot const& s) {
    s.total.write(out); for (auto const& v:s.sources) v.write(out);
    out << s.contacts.size() << ' ';
    for (auto const& item:s.contacts) {
      out << std::get<0>(item.first) << ' ' << std::get<1>(item.first) << ' '
          << std::get<2>(item.first) << ' ' << std::get<0>(item.second) << ' ';
      vectorOut(out,std::get<1>(item.second)); vectorOut(out,std::get<2>(item.second));
      vectorOut(out,std::get<3>(item.second));
    }
    out << s.collisions.size() << ' ';
    for (auto const& v:s.collisions) { out << v.first << ' '; vectorOut(out,v.second); }
    out << s.coordination << ' ' << s.touching.size() << ' ';
    for (int id:s.touching) out << id << ' ';
    out << s.children.size() << ' ';
    for (auto const& child:s.children) snapshotOut(out,child);
  }
  static Snapshot snapshotIn(std::istream& in) {
    Snapshot s; s.total.read(in); for (auto& v:s.sources) v.read(in);
    size_t count; in >> count;
    for (size_t i=0;i<count;++i) {
      int a,b,c; bool active; in >> a >> b >> c >> active;
      Vector3 delta=vectorIn(in), normal=vectorIn(in), rolling=vectorIn(in);
      s.contacts.emplace(std::make_tuple(a,b,c),std::make_tuple(active,delta,normal,rolling));
    }
    in >> count;
    for (size_t i=0;i<count;++i) { int id; in >> id; s.collisions[id]=vectorIn(in); }
    in >> s.coordination >> count;
    for (size_t i=0;i<count;++i) { int id; in >> id; s.touching.insert(id); }
    in >> count; for (size_t i=0;i<count;++i) s.children.push_back(snapshotIn(in));
    if (!in) throw std::runtime_error("Invalid contact snapshot");
    return s;
  }
  static void mergeObstacleHistory(Snapshot& target, Snapshot const& other,
      int source, std::map<int,int> const& ownerMap) {
    for (auto const& item:other.contacts) {
      auto owner=ownerMap.find(std::get<1>(item.first));
      if (owner!=ownerMap.end() && owner->second==source) target.contacts[item.first]=item.second;
    }
    for (auto const& item:other.collisions) {
      auto owner=ownerMap.find(item.first);
      if (owner!=ownerMap.end() && owner->second==source) target.collisions[item.first]=item.second;
    }
    if (target.children.size()!=other.children.size()) throw std::runtime_error("MPI obstacle topology differs");
    for (size_t i=0;i<target.children.size();++i)
      mergeObstacleHistory(target.children[i],other.children[i],source,ownerMap);
  }
  static void clearHistory(Snapshot& s) {
    s.contacts.clear(); s.collisions.clear();
    for (auto& child:s.children) clearHistory(child);
  }
  static std::ostringstream output() {
    std::ostringstream out; out << std::setprecision(17); return out;
  }
  void broadcast(std::string& value) const {
    int size=rank==0 ? static_cast<int>(value.size()) : 0;
    if (comm!=MPI_COMM_NULL) {
      MPI_Bcast(&size,1,MPI_INT,0,comm);
      if (rank!=0) value.resize(size);
      MPI_Bcast(value.data(),size,MPI_CHAR,0,comm);
    }
  }
  void status(std::string error) const {
    broadcast(error);
    if (!error.empty()) throw std::runtime_error(error);
  }
  void check(std::string const& error) const {
    if (comm==MPI_COMM_NULL) { if (!error.empty()) throw std::runtime_error(error); return; }
    int source=error.empty()?INT_MAX:rank;
    MPI_Allreduce(MPI_IN_PLACE,&source,1,MPI_INT,MPI_MIN,comm);
    if (source==INT_MAX) return;
    int size=rank==source?static_cast<int>(error.size()):0;
    MPI_Bcast(&size,1,MPI_INT,source,comm);
    std::string message=rank==source?error:std::string(size,'\0');
    MPI_Bcast(message.data(),size,MPI_CHAR,source,comm);
    throw std::runtime_error(message);
  }
  std::vector<std::string> gather(std::string const& local) const {
    check(local.size()>static_cast<size_t>(INT_MAX/ranks)?"External-model payload too large":"");
    if (comm==MPI_COMM_NULL) return {local};
    int size=static_cast<int>(local.size());
    std::vector<int> sizes(ranks), offsets(ranks);
    MPI_Allgather(&size,1,MPI_INT,sizes.data(),1,MPI_INT,comm);
    int total=0; for (int i=0;i<ranks;++i) { offsets[i]=total; total+=sizes[i]; }
    std::string all(total,'\0');
    MPI_Allgatherv(local.data(),size,MPI_CHAR,all.data(),sizes.data(),offsets.data(),MPI_CHAR,comm);
    std::vector<std::string> result;
    for (int i=0;i<ranks;++i) result.push_back(all.substr(offsets[i],sizes[i]));
    return result;
  }
  State state(Particle* p) const {
    State s{p->getID(),*p->getPosition(),*p->getQuaternionRotation(),
            *p->getTranslationalVelocity(),*p->getAngularVelocity()};
    auto h=registry.unwrapped.find(s.id);
    auto e=g.m_externalStateHistory.find(s.id);
    if (e!=g.m_externalStateHistory.end()) s.position=e->second.position;
    else if (h!=registry.unwrapped.end()) {
      Vector3 origin,length; App::get_origin(origin[X],origin[Y],origin[Z]);
      App::get_size(length[X],length[Y],length[Z]);
      for (int d=0;d<3;++d) if (g.m_periodicity[d])
        s.position[d]+=std::round((h->second.position[d]-s.position[d])/length[d])*length[d];
    }
    return s;
  }
  void particleOut(std::ostream& out, Particle* p) const {
    std::ostringstream data=output(); p->write2014(data);
    out << static_cast<int>(p->getKinematicDriver()) << ' ' << std::quoted(data.str()) << ' ';
    snapshotOut(out,capture(p));
  }
  std::unique_ptr<Particle> particleIn(std::istream& in) const {
    int driver,type; std::string data; in >> driver >> std::quoted(data);
    std::istringstream input(data); input >> type;
    auto const& refs=*g.m_allcomponents.getReferenceParticles();
    if (type<0 || static_cast<size_t>(type)>=refs.size()) throw std::runtime_error("Unknown coupled particle type");
    auto* ref=refs[type];
    auto p=makeParticle(ref->isCompositeParticle(),ref->getSpecificCompositeShapeName());
    p->setGeometricType(type); p->read2014(input,&refs);
    p->setKinematicDriver(static_cast<Driver>(driver)); p->setActivity(COMPUTE); p->setTag(0);
    restore(p.get(),snapshotIn(in));
    if (!input || !in) throw std::runtime_error("Invalid coupled particle data");
    return p;
  }
  std::unique_ptr<StepState> cloneState(StepState const& base) const {
    auto t=std::unique_ptr<StepState>(new StepState);
    t->state=base.state; t->obstacles=base.obstacles;
    for (auto const& p:base.particles) {
      std::unique_ptr<Particle> copy(p->createCloneCopy(false));
      copy->setID(p->getID()); copy->setActivity(COMPUTE); copy->setTag(0);
      restore(copy.get(),capture(p.get())); t->particles.push_back(std::move(copy));
    }
    return t;
  }
  void validateStates(std::vector<State>& states, double elapsed) const {
    if (!std::isfinite(elapsed) || elapsed<0. || elapsed>context.dt)
      throw std::invalid_argument("Substep time is outside the coupled interval");
    if (states.size()!=accepted->particles.size()) throw std::invalid_argument("Substep must contain the complete joint state");
    std::map<int,State const*> initial;
    for (auto const& s:accepted->state.particles) initial[s.id]=&s;
    std::set<int> seen;
    Vector3 origin,length; App::get_origin(origin[X],origin[Y],origin[Z]); App::get_size(length[X],length[Y],length[Z]);
    for (auto& s:states) {
      auto it=initial.find(s.id);
      if (it==initial.end() || !seen.insert(s.id).second) throw std::invalid_argument("Unknown or duplicate substep particle ID");
      for (int d=0;d<3;++d)
        if (!std::isfinite(s.position[d]) || !std::isfinite(s.translationalVelocity[d]) || !std::isfinite(s.angularVelocity[d]))
          throw std::invalid_argument("Nonfinite substep state");
      double norm=0.; for (int d=0;d<4;++d) norm=std::hypot(norm,s.orientation[d]);
      if (!(norm>0.) || !std::isfinite(norm)) throw std::invalid_argument("Invalid substep quaternion");
      for (int d=0;d<4;++d) s.orientation[d]/=norm;
      auto* p=accepted->particles.at(std::distance(initial.begin(),it)).get();
      // Accepted particles and states are sorted by physical ID.
      Quaternion old=it->second->orientation;
      Quaternion rel=s.orientation*old.Conjugate();
      double angle=2.*std::atan2(Norm(*rel.getVector3()),std::abs(rel.getdouble()));
      Vector3 delta=s.position-it->second->position;
      double motion=Norm(delta)+p->getCircumscribedRadius()*angle;
      if (!std::isfinite(motion) || motion>p->getCrustThickness()) throw std::invalid_argument("Substep exceeds the accepted-boundary motion/crust bound");
      for (int d=0;d<3;++d) {
        if (std::abs(delta[d])>=g.m_collision->getCellSize(d)) throw std::invalid_argument("Substep exceeds the collision halo");
        if (!g.m_periodicity[d] && (s.position[d]<origin[d] || s.position[d]>=origin[d]+length[d]))
          throw std::invalid_argument("Substep exits a nonperiodic boundary");
      }
    }
    std::sort(states.begin(),states.end(),[](State const& a,State const& b){return a.id<b.id;});
  }
  void apply(Particle* p, State const& s) const {
    Vector3 origin,length; App::get_origin(origin[X],origin[Y],origin[Z]); App::get_size(length[X],length[Y],length[Z]);
    Point3 position=s.position;
    for (int d=0;d<3;++d) if (g.m_periodicity[d]) position[d]=wrapped(position[d],origin[d],length[d]);
    Transform transform=*p->getRigidBody()->getTransform();
    transform.setOrigin(position); transform.setBasis(Matrix(s.orientation));
    p->applyExternalState(transform,s.orientation,s.translationalVelocity,s.angularVelocity);
  }
  void refreshLoads() {
    auto& t=working;
    t->state.contacts.clear();
    for (size_t i=0;i<t->particles.size();++i) {
      Torsor load=t->particles[i]->getContactLoad();
      load.setPoint(t->state.particles[i].position);
      t->state.contacts.push_back({t->particles[i]->getID(),load});
    }
  }
  void gatherAccepted() {
    std::ostringstream out=output();
    size_t n=0; for (auto* p:*g.m_allcomponents.getActiveParticles()) if (p->getTag()!=2) ++n;
    out << n << ' ';
    for (auto* p:*g.m_allcomponents.getActiveParticles()) if (p->getTag()!=2) {
      stateOut(out,state(p)); particleOut(out,p);
    }
    snapshotOut(out,capture(g.m_allcomponents.getObstacles()));
    auto packets=gather(out.str());
    std::string error;
    if (rank==0) try {
      accepted.reset(new StepState);
      std::map<int,std::pair<State,std::unique_ptr<Particle>>> physical;
      std::vector<Snapshot> obstacleStates;
      for (int r=0;r<ranks;++r) {
        std::istringstream in(packets[r]); size_t count; in >> count;
        for (size_t i=0;i<count;++i) {
          State s=stateIn(in); auto p=particleIn(in);
          if (s.id!=p->getID() || !physical.emplace(s.id,std::make_pair(s,std::move(p))).second)
            throw std::runtime_error("Physical particle has multiple owners");
          owners[s.id]=r;
        }
        obstacleStates.push_back(snapshotIn(in));
      }
      for (auto& item:physical) {
        accepted->state.particles.push_back(item.second.first);
        accepted->particles.push_back(std::move(item.second.second));
      }
      accepted->obstacles=obstacleStates.front(); clearHistory(accepted->obstacles);
      for (int r=0;r<ranks;++r) mergeObstacleHistory(accepted->obstacles,obstacleStates[r],r,owners);
      std::set<int> bound;
      for (auto const& model:registry.models) for (auto const& binding:model.bindings) {
        if (!bound.insert(binding.particleId).second || !physical.count(binding.particleId))
          throw std::runtime_error("Missing or duplicate model particle binding");
      }
      for (auto const& p:accepted->particles)
        if (!p->isIntegratedByGrains() && !bound.count(p->getID()))
          throw std::runtime_error("Every external particle must have a registered model binding");
    } catch (std::exception const& e) { error=e.what(); }
    status(error);
  }

public:
  ExternalModelSession(Grains& simulation, ExternalModelRegistry& data):g(simulation),registry(data) {
    if (g.m_wrapper && !g.m_processorIsActive) throw std::logic_error("Inactive rank cannot participate in coupling");
    if (g.m_wrapper) { comm=g.m_wrapper->get_active_procs_comm(); MPI_Comm_rank(comm,&rank); MPI_Comm_size(comm,&ranks); }
    check(!g.m_collision || g.m_dimension!=3 || g.m_periodicity.size()!=3 ? "External models require an initialized 3D simulation":"");
  }
  int add(std::string const& attachment,std::unique_ptr<Particle> p) override {
    if (registering.empty() || !p || !p->getRigidBody()) throw std::invalid_argument("Invalid particle registration");
    if (attachment.empty() || std::any_of(pendingBindings.begin(),pendingBindings.end(),[&](ParticleBinding const& b){return b.attachment==attachment;}))
      throw std::invalid_argument("Duplicate or empty attachment name");
    if (nextId==INT_MAX) throw std::overflow_error("Particle ID space exhausted");
    int id=++nextId; p->setID(id); p->setKinematicDriver(Driver::External);
    pendingBindings.push_back({registering,attachment,id}); additions.push_back({attachment,std::move(p)}); return id;
  }
  void bind(std::string const& attachment,int id) override {
    if (registering.empty() || attachment.empty() || std::any_of(pendingBindings.begin(),pendingBindings.end(),[&](ParticleBinding const& b){return b.attachment==attachment || b.particleId==id;}))
      throw std::invalid_argument("Invalid or duplicate existing particle binding");
    pendingBindings.push_back({registering,attachment,id});
  }
  void excludeContacts(std::string const& a,std::string const& b) override {
    if (a==b) throw std::invalid_argument("A particle cannot exclude itself");
    pendingExclusions.emplace_back(a,b);
  }
  void registerModel(std::string const& name,std::unique_ptr<ExternalModel> instance);
  void run(double dt);

  void evaluate(std::vector<State> const& input,double elapsed);
  void advanceRigidBodies(double dtVelocity,double dtPosition,double elapsed);
  void advanceModels();
};

void ExternalModelSession::registerModel(std::string const& name,
    std::unique_ptr<ExternalModel> instance)
{
  check(rank!=0 && instance ? "Only root may supply an external model":"");
  int maximum=std::max(g.m_allcomponents.getMaxParticleIDnumber(), Particle::getMaxIDnumber());
  if (comm!=MPI_COMM_NULL) MPI_Allreduce(MPI_IN_PLACE,&maximum,1,MPI_INT,MPI_MAX,comm);
  nextId=maximum;
  auto oldMaterials=ContactBuilderFactory::captureMaterials();
  auto rollbackMaterials=[&]() { ContactBuilderFactory::restoreMaterials(oldMaterials); };
  int types=static_cast<int>(g.m_allcomponents.getReferenceParticles()->size()), low=types,high=types;
  if (comm!=MPI_COMM_NULL) {
    MPI_Allreduce(MPI_IN_PLACE,&low,1,MPI_INT,MPI_MIN,comm);
    MPI_Allreduce(MPI_IN_PLACE,&high,1,MPI_INT,MPI_MAX,comm);
  }
  check(low!=high?"MPI reference particle tables differ":"");
  std::ostringstream local=output();
  for (auto* p:*g.m_allcomponents.getActiveParticles()) if (p->getTag()!=2) local << p->getID() << ' ';
  auto physical=gather(local.str());
  std::string error,wire;
  if (rank==0) try {
    if (!instance || name.empty()) throw std::invalid_argument("Registration requires a name and a model on root");
    for (auto const& model:registry.models) if (model.name==name) throw std::invalid_argument("Duplicate model name");
    registering=name; instance->initialize(*this); registering.clear();
    if (pendingBindings.empty()) throw std::invalid_argument("Model must register at least one particle");
    std::set<int> existing;
    for (auto const& packet:physical) {
      std::istringstream in(packet); int id; while (in>>id) existing.insert(id);
    }
    for (auto const& added:additions) existing.insert(added.body->getID());
    std::set<int> bound;
    for (auto const& model:registry.models) for (auto const& b:model.bindings) bound.insert(b.particleId);
    for (auto const& b:pendingBindings)
      if (!existing.count(b.particleId) || !bound.insert(b.particleId).second)
        throw std::invalid_argument("Unknown or already bound particle ID");
    Vector3 origin,length; App::get_origin(origin[X],origin[Y],origin[Z]); App::get_size(length[X],length[Y],length[Z]);
    std::ostringstream out=output();
    ContactBuilderFactory::writeMaterials(out);
    out << std::quoted(name) << ' ' << nextId << ' ' << additions.size() << ' ';
    for (size_t i=0;i<additions.size();++i) {
      auto* p=additions[i].body.get();
      if (!(p->getCrustThickness()>0.) || !std::isfinite(p->getCircumscribedRadius()) ||
          !p->getRigidBody()->getTransform()->getBasis().isRotation())
        throw std::invalid_argument("Invalid registered particle geometry");
      for (int d=0;d<3;++d)
        if (!std::isfinite((*p->getPosition())[d]) || (*p->getPosition())[d]<origin[d] ||
            (*p->getPosition())[d]>=origin[d]+length[d] ||
            !std::isfinite((*p->getTranslationalVelocity())[d]) || !std::isfinite((*p->getAngularVelocity())[d]))
          throw std::invalid_argument("Registered particle must start inside the primary domain with finite velocities");
      p->setGeometricType(types+i); p->setActivity(COMPUTE);
      Quaternion q; q.setQuaternion(p->getRigidBody()->getTransform()->getBasis());
      p->applyExternalState(*p->getRigidBody()->getTransform(),q,
                           *p->getTranslationalVelocity(),*p->getAngularVelocity());
      p->InitializeForce(false);
      std::ostringstream definition=output(); p->write(definition);
      out << std::quoted(definition.str()) << ' ';
    }
    out << pendingBindings.size() << ' ';
    std::map<std::string,int> attachments;
    for (auto const& b:pendingBindings) {
      attachments[b.attachment]=b.particleId;
      out << std::quoted(b.attachment) << ' ' << b.particleId << ' ';
    }
    out << pendingExclusions.size() << ' ';
    for (auto const& pair:pendingExclusions) {
      if (!attachments.count(pair.first) || !attachments.count(pair.second)) throw std::invalid_argument("Unknown contact-exclusion attachment");
      out << attachments.at(pair.first) << ' ' << attachments.at(pair.second) << ' ';
    }
    wire=out.str();
    if (wire.size()>static_cast<size_t>(INT_MAX)) throw std::overflow_error("External-model payload too large");
  } catch (std::exception const& e) { error=e.what(); }
    catch (...) { error="Model particle registration failed"; }
  try { status(error); } catch (...) { rollbackMaterials(); throw; }
  broadcast(wire);
  std::vector<std::unique_ptr<Particle>> bodies,refs;
  ExternalModelRegistry::Model model;
  std::set<std::pair<int,int>> exclusions;
  error.clear();
  try {
    std::istringstream in(wire); size_t count;
    ContactBuilderFactory::readMaterials(in);
    in >> std::quoted(model.name) >> nextId >> count;
    for (size_t i=0;i<count;++i) {
      std::string definition; in >> std::quoted(definition);
      std::unique_ptr<Particle> body;
      if (rank==0) body=std::move(additions[i].body);
      else {
        std::istringstream data(definition); std::string tag,shape="none";
        data >> tag;
        if (tag=="<CompositeParticle>") { std::string marker; data >> marker >> shape; }
        body=makeParticle(tag=="<CompositeParticle>",shape);
        body->read(data); body->setKinematicDriver(Driver::External); body->setActivity(COMPUTE);
        if (!data) throw std::runtime_error("Cannot deserialize registered particle");
      }
      ContactBuilderFactory::defineMaterial(body->getMaterial(),false);
      std::unique_ptr<Particle> ref(body->createCloneCopy(false));
      ref->setID(GrainsExec::m_ReferenceParticleDefaultID);
      refs.push_back(std::move(ref)); bodies.push_back(std::move(body));
    }
    in >> count;
    for (size_t i=0;i<count;++i) {
      ParticleBinding b; b.model=model.name; in >> std::quoted(b.attachment) >> b.particleId;
      model.bindings.push_back(b);
    }
    in >> count;
    for (size_t i=0;i<count;++i) { int a,b; in >> a >> b; exclusions.emplace(std::min(a,b),std::max(a,b)); }
    if (!in) throw std::runtime_error("Invalid model registration payload");
    auto& all=g.m_allcomponents;
    all.reserveReferenceParticles(refs.size());
    registry.models.reserve(registry.models.size()+1);
  } catch (std::exception const& e) { error=e.what(); }
  try { check(error); }
  catch (...) { rollbackMaterials(); throw; }
  // Publication follows the same force-insertion and ownership rules as normal
  // insertion. No caller retains ownership of adopted component objects.
  try {
    for (auto& ref:refs) g.m_allcomponents.adoptReferenceParticle(std::move(ref));
    for (auto& body:bodies) g.adoptExternalParticle(std::move(body));
    std::set<int> bound;
    for (auto const& b:model.bindings) bound.insert(b.particleId);
    for (auto* p:*g.m_allcomponents.getActiveParticles()) if (bound.count(p->getID())) p->setKinematicDriver(Driver::External);
    Particle::setMaxIDnumber(std::max(Particle::getMaxIDnumber(),nextId));
    registry.excluded.insert(exclusions.begin(),exclusions.end());
    if (rank==0) model.instance=std::move(instance);
    registry.models.push_back(std::move(model));
    registry.contactsReady=false;
    g.m_allcomponents.computeNumberParticles(g.m_wrapper);
  } catch (...) { if (comm!=MPI_COMM_NULL) MPI_Abort(comm,1); throw; }
}

void ExternalModelSession::evaluate(std::vector<State> const& input,double elapsed)
{
  if (rank!=0) throw std::logic_error("Only root evaluates coupled contacts");
  std::vector<State> states=input; validateStates(states,elapsed);
  auto& t=working;
  t->state.particles=states;
  Obstacle* obstacles=g.m_allcomponents.getObstacles();
  Snapshot liveObstacles=capture(obstacles);
  WorkingGlobals globals;
  GrainsExec::m_time_counter=globals.counter+1;
  try {
    restore(obstacles,accepted->obstacles);
    initializeObstacleLoads(obstacles); obstacleHistory(obstacles,false);
    std::vector<std::unique_ptr<Particle>> images;
    Vector3 lengths; App::get_size(lengths[X],lengths[Y],lengths[Z]);
    for (size_t i=0;i<t->particles.size();++i) {
      auto* p=t->particles[i].get();
      restore(p,capture(accepted->particles[i].get()));
      apply(p,states[i]);
      auto noncontact=accepted->particles[i]->getLoads();
      noncontact[static_cast<size_t>(Component::LoadSource::Contact)]=Torsor(*p->getPosition());
      // Forces other than contact are frozen over this macrostep. Reanchor the
      // snapshot without inventing moments due to displacement between stages.
      for (auto& load:noncontact) load.setPoint(*p->getPosition());
      p->InitializeForce(false); p->setLoads(noncontact); p->setContactMapToFalse();
      for (int x=(g.m_periodicity[X]?-1:0);x<=(g.m_periodicity[X]?1:0);++x)
        for (int y=(g.m_periodicity[Y]?-1:0);y<=(g.m_periodicity[Y]?1:0);++y)
          for (int z=(g.m_periodicity[Z]?-1:0);z<=(g.m_periodicity[Z]?1:0);++z) {
            if (!x&&!y&&!z) continue;
            std::unique_ptr<Particle> image(p->createCloneCopy(false));
            image->setID(p->getID()); image->setTag(2); image->setActivity(COMPUTE);
            restore(image.get(),capture(p)); image->Translate(Vector3(x*lengths[X],y*lengths[Y],z*lengths[Z]));
            image->initialize_transformWithCrust_to_notComputed(); images.push_back(std::move(image));
          }
    }
    WorkingDiagnostics diagnostics;
    auto requireLaw=[&](Component* a,Component* b) {
      if (!ContactBuilderFactory::hasContactForceModel(a->getMaterial(),b->getMaterial()))
        throw std::runtime_error("No contact law configured for "+a->getMaterial()+" and "+b->getMaterial());
    };
    auto interact=[&](Particle* a,Particle* b) {
      if (a->getID()==b->getID() || registry.excluded.count({std::min(a->getID(),b->getID()),std::max(a->getID(),b->getID())})) return;
      if (!a->doBVolumeOverlap(b)) return;
      requireLaw(a,b);
      if (b->isCompositeParticle()) b->InterAction(a,elapsed,context.time+elapsed,&diagnostics);
      else a->InterAction(b,elapsed,context.time+elapsed,&diagnostics);
    };
    for (size_t i=0;i<t->particles.size();++i) {
      Particle* p=t->particles[i].get();
      for (size_t j=i+1;j<t->particles.size();++j) interact(p,t->particles[j].get());
      for (auto const& image:images) interact(p,image.get());
      for (auto* obstacle:obstacles->getObstacles()) {
        if (!p->doBVolumeOverlap(obstacle)) continue;
        requireLaw(p,obstacle);
        if (p->isCompositeParticle()) p->InterAction(obstacle,elapsed,context.time+elapsed,&diagnostics);
        else obstacle->InterAction(p,elapsed,context.time+elapsed,&diagnostics);
      }
      p->updateContactMap();
    }
    obstacleHistory(obstacles,true); t->obstacles=capture(obstacles);
    restore(obstacles,liveObstacles);
  } catch (std::exception const&) { restore(obstacles,liveObstacles); throw; }
    catch (...) { restore(obstacles,liveObstacles); throw std::runtime_error("Coupled contact evaluation failed"); }
  refreshLoads();
}

void ExternalModelSession::advanceRigidBodies(double dtVelocity,double dtPosition,double elapsed)
{
  if (!std::isfinite(dtVelocity) || !std::isfinite(dtPosition) || dtVelocity<0. || dtPosition<0. ||
      dtVelocity>context.dt || dtPosition>context.dt) throw std::invalid_argument("Invalid integration-stage increments");
  auto& t=working;
  WorkingGlobals globals; GrainsExec::m_time_counter=globals.counter+1;
  std::list<Particle*> particles;
  for (auto const& p:t->particles) particles.push_back(p.get());
  RigidBodyKinematicsDriver driver;
  if (dtPosition==0.) driver.advanceVelocity(particles,context.time+elapsed,dtVelocity);
  else driver.move(particles,context.time+elapsed,dtVelocity,dtPosition);
  std::vector<State> states;
  Vector3 lengths; App::get_size(lengths[X],lengths[Y],lengths[Z]);
  for (size_t i=0;i<t->particles.size();++i) {
    auto* p=t->particles[i].get(); State s=t->state.particles[i];
    if (p->isIntegratedByGrains()) {
      Point3 old=s.position;
      s.position=*p->getPosition();
      for (int d=0;d<3;++d) if (g.m_periodicity[d]) s.position[d]+=std::round((old[d]-s.position[d])/lengths[d])*lengths[d];
      s.orientation=*p->getQuaternionRotation(); s.translationalVelocity=*p->getTranslationalVelocity(); s.angularVelocity=*p->getAngularVelocity();
    }
    states.push_back(s);
  }
  validateStates(states,elapsed); t->state.particles=states;
  refreshLoads();
}

void ExternalModelSession::run(double dt)
{
  check(registry.models.empty()?"No external models are registered":"");
  check(!std::isfinite(dt) || dt<=0.?"Invalid coupled step or nonroot scheme":"");
  double minDt=dt,maxDt=dt,minTime=g.m_time,maxTime=g.m_time;
  if (comm!=MPI_COMM_NULL) {
    MPI_Allreduce(MPI_IN_PLACE,&minDt,1,MPI_DOUBLE,MPI_MIN,comm);
    MPI_Allreduce(MPI_IN_PLACE,&maxDt,1,MPI_DOUBLE,MPI_MAX,comm);
    MPI_Allreduce(MPI_IN_PLACE,&minTime,1,MPI_DOUBLE,MPI_MIN,comm);
    MPI_Allreduce(MPI_IN_PLACE,&maxTime,1,MPI_DOUBLE,MPI_MAX,comm);
  }
  check(minDt!=maxDt || minTime!=maxTime?"MPI coupling clocks or increments differ":"");
  check(g.m_allcomponents.hasImposedObstacleMotion() ?
        "Explicit coupling currently requires fixed obstacles":"");
  context.time=g.m_time; context.dt=dt;
  gatherAccepted();
  std::string error,wire;
  if (rank==0) try {
    working=cloneState(*accepted);
    refreshLoads();
    bool leapfrog=GrainsExec::m_TIScheme=="SecondOrderLeapFrog";
    // All models receive loads from the same starting configuration. Ordinary
    // explicit schemes evaluate at the start; leapfrog reuses the last endpoint
    // force, initializing it without advancing friction history on the first step.
    if (!leapfrog || !registry.contactsReady)
      evaluate(working->state.particles,leapfrog?0.:dt);
    advanceModels();
    advanceRigidBodies(leapfrog?dt/2.:dt,dt,dt);
    if (leapfrog) {
      evaluate(working->state.particles,dt);
      advanceRigidBodies(dt/2.,0.,dt);
    }
    StepState const& chosen=*working;
    StepData const& selected=chosen.state;
    std::ostringstream out=output();
    out << accepted->state.particles.size() << ' ';
    for (auto const& s:accepted->state.particles) stateOut(out,s);
    out << selected.particles.size() << ' ';
    for (size_t i=0;i<chosen.particles.size();++i) {
      stateOut(out,selected.particles[i]); particleOut(out,chosen.particles[i].get());
    }
    snapshotOut(out,chosen.obstacles); wire=out.str();
    if (wire.size()>static_cast<size_t>(INT_MAX)) throw std::overflow_error("External-model payload too large");
  } catch (std::exception const& e) { error=e.what(); }
    catch (...) { error="Explicit coupled advance failed"; }
  status(error); broadcast(wire);

  // Parse and allocate all candidate state before validation or publication.
  std::map<int,std::unique_ptr<Particle>> prepared;
  std::map<int,Driver> drivers;
  std::map<int,State> seed,nextHistory,nextUnwrapped;
  Snapshot obstacleState;
  StepData result;
  error.clear();
  try {
    std::istringstream in(wire); size_t count;
    in >> count;
    for (size_t i=0;i<count;++i) { State s=stateIn(in); seed[s.id]=s; }
    in >> count;
    for (size_t i=0;i<count;++i) {
      State s=stateIn(in); auto p=particleIn(in);
      drivers[s.id]=p->getKinematicDriver(); nextUnwrapped[s.id]=s;
      if (!p->isIntegratedByGrains()) nextHistory[s.id]=s;
      Torsor load=p->getContactLoad(); load.setPoint(s.position);
      result.particles.push_back(s); result.contacts.push_back({s.id,load});
      prepared.emplace(s.id,std::move(p));
    }
    obstacleState=snapshotIn(in);
  } catch (std::exception const& e) { error=e.what(); }
  check(error);
  auto oldHistory=g.m_externalStateHistory;
  std::vector<std::pair<Particle*,Driver>> originalDrivers;
  for (auto* p:*g.m_allcomponents.getActiveParticles()) {
    originalDrivers.emplace_back(p,p->getKinematicDriver()); p->setKinematicDriver(Driver::External);
  }
  g.m_externalStateHistory=seed;
  try {
    g.importParticleStates(rank==0?result.particles:std::vector<State>(),false);
  } catch (...) {
    g.m_externalStateHistory.swap(oldHistory);
    for (auto const& item:originalDrivers) item.first->setKinematicDriver(item.second);
    throw;
  }
  // This is the only mutation boundary. Distributed synchronization failures
  // after this point are fatal, matching the existing import API contract.
  try {
    g.m_time=context.time+dt; g.m_dt=dt;
    g.importParticleStates(rank==0?result.particles:std::vector<State>(),true);
    for (auto* p:*g.m_allcomponents.getActiveParticles()) {
      auto const& source=prepared.at(p->getID());
      p->copyIntegrationState(*source);
      Transform pose=*p->getRigidBody()->getTransform();
      p->applyExternalState(pose,*source->getQuaternionRotation(),
                            *source->getTranslationalVelocity(),*source->getAngularVelocity());
      p->setKinematicDriver(drivers.at(p->getID()));
      restore(p,capture(source.get()),true);
    }
    restore(g.m_allcomponents.getObstacles(),obstacleState);
    g.m_externalStateHistory.swap(nextHistory); registry.unwrapped.swap(nextUnwrapped);
    g.m_allcomponents.computeNumberParticles(g.m_wrapper);
    if (comm!=MPI_COMM_NULL) MPI_Barrier(comm);
    registry.contactsReady=true;
    ++GrainsExec::m_time_counter;
  } catch (...) { registry.poisoned=true; if (comm!=MPI_COMM_NULL) MPI_Abort(comm,1); throw; }
  return;
}

void ExternalModelSession::advanceModels()
{
  auto joint=working->state.particles;
  std::map<int,size_t> indices;
  for (size_t i=0;i<joint.size();++i) indices[joint[i].id]=i;
  for (auto& model:registry.models) {
    std::set<int> expected;
    for (auto const& binding:model.bindings) expected.insert(binding.particleId);
    std::vector<ExternalModel::ContactLoad> loads;
    for (auto const& load:working->state.contacts)
      if (expected.count(load.particleId)) loads.push_back(load);
    auto states=model.instance->advance(context.time,context.dt,loads);
    if (states.size()!=expected.size())
      throw std::invalid_argument("External model must return exactly its bound particle states");
    for (auto const& state:states) {
      if (!expected.erase(state.id))
        throw std::invalid_argument("External model returned an unknown or duplicate particle ID");
      joint.at(indices.at(state.id))=state;
    }
  }
  validateStates(joint,context.dt);
  for (size_t i=0;i<joint.size();++i) {
    apply(working->particles[i].get(),joint[i]);

  }
  working->state.particles=std::move(joint);
}

void Grains::registerExternalModel(std::string const& name,
    std::unique_ptr<ExternalModel> instance)
{
  if (m_externalModels && m_externalModels->busy) throw std::logic_error("Cannot register a model during a coupled step");
  if (!m_externalModels) m_externalModels=std::make_shared<ExternalModelRegistry>();
  if (m_externalModels->poisoned) throw std::logic_error("A failed coupled advance requires a new simulation");
  ExternalModelSession session(*this,*m_externalModels);
  m_externalModels->busy=true;
  try {
    session.registerModel(name,std::move(instance)); m_externalModels->busy=false;
  } catch (...) { m_externalModels->busy=false; throw; }
}

bool Grains::hasExternalModels() const
{ return m_externalModels && !m_externalModels->models.empty(); }

bool Grains::isExternalModelStepActive() const
{ return m_externalModels && m_externalModels->busy; }

void Grains::advanceExternalModels(double dt)
{
  if (!m_externalModels) throw std::logic_error("Register external models before advancing a coupled step");
  if (m_externalModels->busy) throw std::logic_error("Nested coupled steps are forbidden");
  if (m_externalModels->poisoned) throw std::logic_error("A failed coupled advance requires a new simulation");
  ExternalModelSession session(*this,*m_externalModels);
  m_externalModels->busy=true;
  try {
    session.run(dt);
    m_externalModels->busy=false;
  } catch (...) {
    m_externalModels->poisoned=true;
    m_externalModels->busy=false; throw;
  }
}

void Grains::adoptExternalParticle(std::unique_ptr<Particle> particle)
{
  bool keep=true;
  if (m_wrapper) {
    auto result=m_collision->insertParticleParallel(m_time,particle.get(),
        m_allcomponents.getActiveParticles(),m_allcomponents.getCloneParticles(),
        m_allcomponents.getReferenceParticles(),m_periodic,true,m_wrapper);
    keep=result.first;
  } else m_collision->insertParticleSerial(particle.get(),
      m_allcomponents.getActiveParticles(),m_allcomponents.getPeriodicCloneParticles(),
      m_allcomponents.getReferenceParticles(),m_periodic,true);
  if (keep) m_allcomponents.adoptActiveParticle(std::move(particle),m_wrapper!=nullptr);
}

// Read geometry on every active rank before contact-law setup. Only root keeps
// the solver instance after collective registration; all ranks need its materials
// and maximum radius while constructing the collision grid.
double Grains::readExternalModels(DOMNode* construction)
{
  auto* models=ReaderXML::getNode(construction,"ExternalModels");
  if (!models) return 0.;
  if (m_dimension!=3 || m_restart)
    throw std::invalid_argument("XML external models require a fresh 3D simulation");
  auto type=ReaderXML::getNodeAttr_String(construction->getParentNode(),"Type");
  if (type!="Standard" && type!="MPI")
    throw std::invalid_argument("XML external models support Standard and MPI simulations");
  auto* simulation=ReaderXML::getNode(construction->getParentNode(),"Simulation");
  auto* output=simulation?ReaderXML::getNode(simulation,"PostProcessing"):nullptr;
  if (output && ReaderXML::getNode(output,"ForceStats"))
    throw std::invalid_argument("External model contact statistics output is not implemented");
  auto* writers=output?ReaderXML::getNode(output,"Writers"):nullptr;
  if (writers) for (auto* writer=writers->getFirstChild();writer;writer=writer->getNextSibling())
    if (writer->getNodeType()==DOMNode::ELEMENT_NODE && ReaderXML::hasNodeAttr(writer,"ForceNetwork") &&
        ReaderXML::getNodeAttr_String(writer,"ForceNetwork")=="True")
      throw std::invalid_argument("External model force-network output is not implemented");
  std::set<std::string> names;
  double radius=0.;
  for (auto* node=models->getFirstChild();node;node=node->getNextSibling()) {
    if (node->getNodeType()!=DOMNode::ELEMENT_NODE) continue;
    if (ReaderXML::getNodeName(node)!="Pendulum") throw std::invalid_argument("Unknown external model type");
    auto name=ReaderXML::getNodeAttr_String(node,"Name");
    if (name.empty() || !names.insert(name).second) throw std::invalid_argument("External model names must be unique and nonempty");
    std::unique_ptr<Pendulum> model(new Pendulum(node));
    radius=std::max(radius,model->maximumRadius());
    m_pendingExternalModels.emplace_back(name,std::move(model));
  }
  if (names.empty()) throw std::invalid_argument("ExternalModels requires at least one model");
  return radius;
}

void Grains::registerPendingExternalModels()
{
  for (auto& pending:m_pendingExternalModels) {
    std::unique_ptr<ExternalModel> model;
    if (m_rank==0) model=std::move(pending.second);
    registerExternalModel(pending.first,std::move(model));
  }
  m_pendingExternalModels.clear();
}

void Grains::simulateExternalModels(double timeInterval)
{
  if (!m_processorIsActive) return;
  if (!std::isfinite(m_dt) || m_dt<=0. || !std::isfinite(timeInterval) || timeInterval<0.)
    throw std::invalid_argument("External model simulation requires a positive fixed timestep");
  if (m_insertion_mode==IM_OVERTIME || GrainsExec::m_partialPer_is_active || m_allApp.size()!=1)
    throw std::invalid_argument("External model runs currently support initial insertion, full periodicity, and constant gravity only");
  double end=timeInterval>0.?std::min(m_tend,m_time+timeInterval):m_tend;
  SCT_insert_app("Move"); SCT_insert_app("UpdateParticleActivity");
  SCT_insert_app("LinkUpdate"); SCT_insert_app("OutputResults");
  while (end-m_time>0.01*m_dt) {
    double next=m_time+m_dt;
    m_lastTime_save=false;
    GrainsExec::m_output_data_at_this_time=false;
    GrainsExec::m_postprocess_forces_at_this_time=false;
    if (m_timeSave!=m_save.end() && *m_timeSave-next<0.01*m_dt) {
      GrainsExec::m_output_data_at_this_time=true;
      ++m_timeSave; m_lastTime_save=true;
    }
    advanceExternalModels(m_dt);
    m_allcomponents.computeObstaclesLoad(m_time,m_dt,m_wrapper);
    m_allcomponents.outputObstaclesLoad(m_time,m_dt,false,false,m_rank);
    if (GrainsExec::m_output_data_at_this_time) {
      SCT_set_start("OutputResults");
      double vmax=0.,vmean=0.;
      m_allcomponents.ComputeMaxMeanVelocity(vmax,vmean,m_wrapper);
      if (m_rank==0) fVitMax << GrainsExec::doubleToString(ios::scientific,6,m_time)
                           << '\t' << vmax << '\t' << vmean << '\n';
      m_allcomponents.PostProcessing(m_time,m_dt,m_collision,m_insertion_windows);
      saveReload(m_time);
      SCT_get_elapsed_time("OutputResults");
    }
  }
}
