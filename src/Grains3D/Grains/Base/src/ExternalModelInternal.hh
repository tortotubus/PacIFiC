/** Private model ownership and MPI helpers; concrete models use ExternalModel.hh. */
#pragma once
#include "Grains.hh"
#include <set>
#include <sstream>

namespace GrainsExternalModelDetail {
using Driver = Particle::KinematicDriver;
struct ParticleBinding {
  std::string model;
  std::string attachment;
  int particleId;
};
std::unique_ptr<Particle> makeParticle(bool composite, std::string const &shape);
std::ostringstream output();
void broadcast(std::string &value, MPI_Comm comm);
void status(std::string error, MPI_Comm comm);
void check(std::string const &error, MPI_Comm comm);
std::vector<std::string> gather(std::string const &local, MPI_Comm comm);
}

/** Grains-owned data only. Solver instances exist on active rank zero. */
struct Grains::ExternalModelData {
  struct Model {
    std::string name;
    std::unique_ptr<ExternalModel> instance;
    std::vector<GrainsExternalModelDetail::ParticleBinding> bindings;
  };
  std::vector<Model> models;
  bool busy = false;
  bool poisoned = false;
  bool contactsReady = false;
};
