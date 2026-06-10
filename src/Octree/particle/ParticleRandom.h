

/*!
 * @file ParticleRandom.h
 * @brief Counter-based deterministic random helpers for particle algorithms.
 *
 * These helpers avoid process-local global RNG state (`rand()`/`srand()`).
 * Every random value is a pure function of:
 *
 *   - a user-controlled seed,
 *   - a stream identifier,
 *   - an explicit counter/index.
 *
 * This is the intended pattern for particle code because serial and MPI runs
 * can visit particles in different orders and can migrate particles between
 * ranks. A call such as
 *
 *   particle_random_uniform01(seed, stream, index)
 *
 * returns the same value on every rank and in every run, provided the same
 * seed/stream/index tuple is used. Callers should therefore choose indices
 * from deterministic simulation identifiers rather than traversal order.
 *
 * These functions are for reproducible simulation noise and initialization,
 * not cryptographic randomness.
 * 
 * Reference:
 *   Guy L. Steele Jr., Doug Lea, and Christine H. Flood,
 *   "Fast splittable pseudorandom number generators", OOPSLA 2014,
 *   DOI: 10.1145/2714064.2660195.
 *
 * The 64-bit mixing constants match the widely used splitmix64 finalizer
 * popularized by Sebastiano Vigna's splitmix64 reference implementation.
 */

#include <stdint.h>
#include <math.h>

#include "common/DomainGeometry.h"

// ============================================================================
// Function declarations
// ============================================================================

static inline uint64_t particle_random_mix64 (uint64_t x);
static inline uint64_t particle_random_fold (uint64_t a, uint64_t b);
static inline uint64_t
particle_random_u64 (uint64_t seed, uint64_t stream, uint64_t index);
static inline double
particle_random_uniform01 (uint64_t seed, uint64_t stream, uint64_t index);
static inline double
particle_random_normal (uint64_t seed, uint64_t stream, uint64_t index);
static inline double
particle_random_uniform01_open (uint64_t seed, uint64_t stream, uint64_t index);
static inline coord
particle_random_gaussian_vector (uint64_t seed, uint64_t stream, uint64_t index);
static inline coord
particle_random_domain_position (uint64_t seed, uint64_t stream, uint64_t index);

/*!
 * Default seed for callers that do not expose their own seed parameter.
 * Production simulations should generally define their own seed in the case
 * setup so that random initialization and stochastic forcing are explicit.
 */
#ifndef PARTICLE_RANDOM_DEFAULT_SEED
#define PARTICLE_RANDOM_DEFAULT_SEED 1ULL
#endif

/*!
 * Built-in stream identifiers.
 *
 * Streams separate unrelated random choices. Reusing an index in two different
 * streams intentionally gives unrelated values.
 */

#define PARTICLE_RANDOM_STREAM_INIT 0x1111111111111111ULL
#define PARTICLE_RANDOM_STREAM_BROWNIAN 0x2222222222222222ULL
#define PARTICLE_RANDOM_STREAM_PLACEMENT 0x3333333333333333ULL

// ============================================================================
// Function definitions
// ============================================================================

/*!
 * @brief Mix a 64-bit integer into another 64-bit integer with good bit diffusion.
 *
 * This is the SplitMix64 finalizer. It is fast and deterministic across
 * platforms with standard unsigned 64-bit arithmetic.
 */
static inline uint64_t particle_random_mix64 (uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

/*!
 * @brief Combine two deterministic identifiers into one mixed counter.
 *
 * Use this to build indices from stable quantities, e.g. `(gid, iter)` or
 * `(particle_number, attempt)`.
 */
static inline uint64_t particle_random_fold (uint64_t a, uint64_t b) {
  return particle_random_mix64 (a ^ particle_random_mix64 (b));
}

/*!
 * @brief Return a deterministic 64-bit pseudorandom value.
 *
 * The output depends only on @p seed, @p stream, and @p index. It does not
 * depend on MPI rank, local particle order, or previous random draws.
 */
static inline uint64_t
particle_random_u64 (uint64_t seed, uint64_t stream, uint64_t index) {
  uint64_t x = particle_random_fold (seed, stream);
  return particle_random_fold (x, index);
}

/*!
 * @brief Return a deterministic uniform deviate in the half-open interval [0, 1).
 */
static inline double
particle_random_uniform01 (uint64_t seed, uint64_t stream, uint64_t index) {
  const uint64_t r = particle_random_u64 (seed, stream, index);
  return (double) (r >> 11) * 0x1.0p-53;
}

/*!
 * @brief Return a deterministic uniform deviate in the open interval (0, 1).
 *
 * This is useful for transforms such as Box-Muller that evaluate log(u).
 */
static inline double
particle_random_uniform01_open (uint64_t seed, uint64_t stream, uint64_t index) {
  double u = particle_random_uniform01 (seed, stream, index);
  return u > 0. ? u : 0x1.0p-53;
}

/*!
 * @brief Return a deterministic standard normal deviate.
 *
 * Uses Box-Muller from two counter-based uniform draws. The same
 * seed/stream/index tuple always gives the same normal deviate.
 */
static inline double
particle_random_normal (uint64_t seed, uint64_t stream, uint64_t index) {
  const double u1 = particle_random_uniform01_open (seed, stream, 2 * index);
  const double u2 = particle_random_uniform01 (seed, stream, 2 * index + 1);
  return sqrt (-2. * log (u1)) * cos (2. * pi * u2);
}

/*!
 * @brief Return a deterministic vector of independent standard normal deviates.
 *
 * The vector components are generated from component-specific counters derived
 * from @p index. This is suitable for per-particle Langevin/Brownian kicks when
 * @p index is built from stable identifiers such as `(particle->gid, iter)`.
 */
static inline coord
particle_random_gaussian_vector (uint64_t seed, uint64_t stream, uint64_t index) {
  coord g = {0., 0., 0.};
  g.x = particle_random_normal (seed, stream, 3 * index);
#if dimension >= 2
  g.y = particle_random_normal (seed, stream, 3 * index + 1);
#endif
#if dimension >= 3
  g.z = particle_random_normal (seed, stream, 3 * index + 2);
#endif
  return g;
}

/*!
 * @brief Return a deterministic random position in the physical domain.
 *
 * The position is uniform in each active coordinate direction. Rectangular
 * multigrid domains use the same length convention as the MPI particle-owner
 * code.
 *
 */
static inline coord
particle_random_domain_position (uint64_t seed, uint64_t stream, uint64_t index) {
  coord p = {0., 0., 0.};
  p.x = X0 + domain_extent_x () *
             particle_random_uniform01 (seed, stream, 3 * index);
#if dimension >= 2
  p.y = Y0 + domain_extent_y () *
             particle_random_uniform01 (seed, stream, 3 * index + 1);
#endif
#if dimension >= 3
  p.z = Z0 + domain_extent_z () *
             particle_random_uniform01 (seed, stream, 3 * index + 2);
#endif
  return p;
}
