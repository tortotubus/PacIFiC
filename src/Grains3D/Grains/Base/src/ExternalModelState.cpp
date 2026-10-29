/** @file ExternalModelState.cpp
 * @brief Registration wire helpers and collective error propagation.
 */
#include "ExternalModelInternal.hh"
#include "Dendrite.hh"
#include "QuadrilobeCylinder.hh"
#include "TrilobeCylinder.hh"
#include <climits>
#include <iomanip>
#include <stdexcept>
namespace GrainsExternalModelDetail {
std::unique_ptr<Particle> makeParticle(bool composite,
                                       std::string const &shape) {
  if (!composite)
    return std::unique_ptr<Particle>(new Particle(false));
  if (shape == "TrilobeCylinder")
    return std::unique_ptr<Particle>(new TrilobeCylinder(false));
  if (shape == "QuadrilobeCylinder")
    return std::unique_ptr<Particle>(new QuadrilobeCylinder(false));
  if (shape == "Dendrite")
    return std::unique_ptr<Particle>(new Dendrite(false));
  if (shape != "none")
    throw std::invalid_argument("Unsupported serialized composite type");
  return std::unique_ptr<Particle>(new CompositeParticle(false));
}
std::ostringstream output() {
  std::ostringstream out;
  out << std::setprecision(17);
  return out;
}

void broadcast(std::string &value, MPI_Comm comm) {
  int rank = 0;
  if (comm != MPI_COMM_NULL) MPI_Comm_rank(comm, &rank);
  int size = rank == 0 ? static_cast<int>(value.size()) : 0;
  if (comm != MPI_COMM_NULL) {
    MPI_Bcast(&size, 1, MPI_INT, 0, comm);
    if (rank != 0)
      value.resize(size);
    MPI_Bcast(value.data(), size, MPI_CHAR, 0, comm);
  }
}

void status(std::string error, MPI_Comm comm) {
  broadcast(error, comm);
  if (!error.empty())
    throw std::runtime_error(error);
}

void check(std::string const &error, MPI_Comm comm) {
  if (comm == MPI_COMM_NULL) {
    if (!error.empty())
      throw std::runtime_error(error);
    return;
  }
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  int source = error.empty() ? INT_MAX : rank;
  MPI_Allreduce(MPI_IN_PLACE, &source, 1, MPI_INT, MPI_MIN, comm);
  if (source == INT_MAX)
    return;
  int size = rank == source ? static_cast<int>(error.size()) : 0;
  MPI_Bcast(&size, 1, MPI_INT, source, comm);
  std::string message = rank == source ? error : std::string(size, '\0');
  MPI_Bcast(message.data(), size, MPI_CHAR, source, comm);
  throw std::runtime_error(message);
}

std::vector<std::string> gather(std::string const &local, MPI_Comm comm) {
  int ranks = 1;
  if (comm != MPI_COMM_NULL) MPI_Comm_size(comm, &ranks);
  check(local.size() > static_cast<size_t>(INT_MAX / ranks)
            ? "External-model payload too large"
            : "", comm);
  if (comm == MPI_COMM_NULL)
    return {local};
  int size = static_cast<int>(local.size());
  std::vector<int> sizes(ranks), offsets(ranks);
  MPI_Allgather(&size, 1, MPI_INT, sizes.data(), 1, MPI_INT, comm);
  int total = 0;
  for (int i = 0; i < ranks; ++i) {
    offsets[i] = total;
    total += sizes[i];
  }
  std::string all(total, '\0');
  MPI_Allgatherv(local.data(), size, MPI_CHAR, all.data(), sizes.data(),
                 offsets.data(), MPI_CHAR, comm);
  std::vector<std::string> result;
  for (int i = 0; i < ranks; ++i)
    result.push_back(all.substr(offsets[i], sizes[i]));
  return result;
}

}
