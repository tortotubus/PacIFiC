#pragma once

#include "common/DomainLocate.h"
#include "particle/ParticleConfig.h"

// ============================================================================
// Type definitions
// ============================================================================

/*!
 * @struct Particle
 *
 * @brief Particle state.
 */
typedef struct {
  uint64_t gid; /*!< Global ID */
  int pid;      /*!< MPI rank that owns this p: -1 if unknown */
} Particle;

#include "particle/ParticleFields.h"

// ============================================================================
// Function declarations
// ============================================================================

int particle_init(Particle *particle);
static int particle_set_gid(Particle *particle);
void particle_free(Particle *particle);
int particle_copy(Particle *dst, const Particle *src);
void particle_move(Particle *dst, Particle *src);
void particle_swap(Particle *a, Particle *b);
static inline bool particle_is_local(Particle *particle);

// ============================================================================
// Macros
// ============================================================================

macro PARTICLE_VARIABLES(Particle *particle = particle) {
  coord pos = {0};
  NOT_UNUSED(pos);
  foreach_dimension() { pos.x = pval(npos.x); }
  domain_wrap_coord(pos);
  double px = pos.x, py = pos.y, pz = pos.z;
}

// ============================================================================
// Function definitions
// ============================================================================

/*!
 * @brief Initialize an @c Particle.
 *
 * @param particle Node to initialize.
 * @return 0 on success, -1 on allocation failure or invalid input.
 *
 * @relates Particle
 */
int particle_init(Particle *particle) {
  if (!particle)
    return -1;

  foreach_dimension() {
    pval(npos.x) = 0.;
    pval(nvel.x) = 0.;
    pval(nforce.x) = 0.;
  }

#if _MPI
  particle->pid = -1;
#else
  particle->pid = 0;
#endif

  particle->gid = 0;

  return 0;
}

/*!
 * @brief Set a particle gid
 * @relates Particle
 */

static int particle_set_gid(Particle *particle) {
  static uint64_t counter = 1;

  if (!particle)
    return -1;

  if (particle->pid < 0)
    return -1;

  if (particle->gid)
    return particle->gid;

  if (counter >= (1ULL << 32))
    return -1;

  particle->gid = ((uint64_t)particle->pid << 32) | counter;

  counter++;

  return 0;
}

/*!
 * @brief Free resources owned by an @c Particle.
 *
 * @param particle Node to free.
 *
 * @relates Particle
 */
void particle_free(Particle *particle) {
  if (!particle)
    return;
}

/*!
 * @brief Deep-copy an @ref Particle.
 *
 * @param dst Destination particle (will be overwritten; any owned resources are
 * freed first).
 * @param src Source particle.
 * @return 0 on success, -1 on allocation failure or invalid input.
 *
 * @relates Particle
 */
int particle_copy(Particle *dst, const Particle *src) {
  if (!dst || !src)
    return -1;

  /* Free any resources currently owned by dst */
  particle_free(dst);

  /* Copy trivially-copiable fields first (but DO NOT keep
   * src->stencil.particle) */
  *dst = *src;

  return 0;
}

/*!
 * @brief Move an @ref Particle (transfer ownership).
 *
 * @param dst Destination particle (will be overwritten; any owned resources are
 * freed first).
 * @param src Source particle (left in a valid, empty state).
 *
 * @relates Particle
 */
void particle_move(Particle *dst, Particle *src) {
  if (!dst || !src)
    return;

  particle_free(dst);

  *dst = *src;
}

/*!
 * @brief Swap two @ref Particle values.
 *
 * @param a First particle.
 * @param b Second particle.
 *
 * @relates Particle
 */
void particle_swap(Particle *a, Particle *b) {
  if (!a || !b)
    return;
  Particle tmp = *a;
  *a = *b;
  *b = tmp;
}

/*!
 * @brief Returns true when @p particle is owned by this rank.
 *
 * @relates Particle
 */
static inline bool particle_is_local(Particle *particle) {
  if (!particle)
    return false;
#if _MPI
  return particle->pid == pid();
#else
  return true;
#endif
}

/*!
 * @brief Returns true when @p particle is a ghost
 *
 * @relates Particle
 */

static inline bool particle_is_ghost(Particle *particle) {
  return !particle_is_local(particle);
}
