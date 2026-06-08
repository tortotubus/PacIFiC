#pragma once

#include "grid/mempool.h"
#include "particle/ParticleFields.h"
#include "particle/ParticleList.h"

// ============================================================================
// Type definitions
// ============================================================================

/*!
 * @struct ParticleMempool
 *
 * @brief
 */
typedef struct {
  Mempool* pool;       /*!< Memory pool pointer */
  ParticleList active; /*!< list of pointers to active Particles */
  int len;             /*!< The (1D) size of the array */
  size_t datasize;     /*!< Extra bytes stored after each Particle */
} ParticleMempool;

// ============================================================================
// Function declarations
// ============================================================================

ParticleMempool particle_mempool_init (size_t pool_bytes, size_t datasize);
void particle_mempool_free (ParticleMempool* pmp);
Particle* particle_mempool_alloc_particle (ParticleMempool* pmp);
long particle_mempool_index_of (ParticleMempool* pmp, Particle* particle);
int particle_mempool_free_particle_ptr (ParticleMempool* pmp, Particle* particle);
void particle_mempool_free_particle (ParticleMempool* pmp, long i);
static inline size_t round_up_multiple (size_t n, size_t a);
static inline size_t particle_mempool_stride (const ParticleMempool* pmp);

// ============================================================================
// Macros
// ============================================================================

#define alignof(T)                                                                       \
  offsetof (                                                                             \
    struct {                                                                             \
      char c;                                                                            \
      T x;                                                                               \
    },                                                                                   \
    x)

// ============================================================================
// Function definitions
// ============================================================================

/*!
 * @brief
 *
 * @relates ParticleMempool
 */
static inline size_t round_up_multiple (size_t n, size_t a) {
  return a ? ((n + a - 1) / a) * a : n;
}

/*!
 * @brief
 *
 * @relates ParticleMempool
 */
static inline size_t particle_mempool_stride (const ParticleMempool* pmp) {
  assert (pmp);
  return round_up_multiple (sizeof (Particle) + pmp->datasize, 8);
}

/*!
 * @brief Initialize an ParticleMempool with a specified pool size and initial capacity.
 *
 * @param pool_bytes Total size in bytes for the memory pool.
 * @param datasize Total size in bytes of the data
 *
 * @return An initialized ParticleMempool structure.
 *
 * @details Creates a new memory pool for Particle objects and allocates an array to
 * track active nodes. The pool block size is rounded up to the system's max alignment.
 * The active list starts empty and grows dynamically as nodes are added.
 *
 * @relates ParticleMempool
 */
ParticleMempool particle_mempool_init (size_t pool_bytes, size_t datasize) {
  ParticleMempool pmp = {0};
  pmp.datasize = datasize;
  size_t block = particle_mempool_stride (&pmp);
  pmp.pool = mempool_new (pool_bytes, block);
  particle_list_init (&pmp.active, 0);
  return pmp;
}

/*!
 * @brief Destroy an ParticleMempool and free all associated resources.
 *
 * @param pmp Pointer to the ParticleMempool to destroy. Must not be NULL.
 *
 * @details Deallocates the active nodes list and destroys the underlying memory pool.
 * Resets all pointers to NULL and counts to 0. Does not free individual Particle objects
 * (they are managed by the memory pool).
 *
 * @relates ParticleMempool
 */
void particle_mempool_free (ParticleMempool* pmp) {
  assert (pmp);
  particle_list_free (&pmp->active);
  mempool_destroy (pmp->pool);
  pmp->pool = NULL;
}

/*!
 * @brief Create a new Particle from the memory pool and add it to the active list.
 *
 * @param pmp Pointer to the ParticleMempool. Must not be NULL.
 *
 * @return Pointer to the newly created Particle.
 *
 * @details Allocates a new Particle from the pool. The new node is added to the active
 * list. Program aborts if allocation fails.
 *
 * @relates ParticleMempool
 */
Particle* particle_mempool_alloc_particle (ParticleMempool* pmp) {
  Particle* node = (Particle*) mempool_alloc0 (pmp->pool);
  assert (node);
  particle_list_push (&pmp->active, node);
  return node;
}

/*!
 * @brief Find the index of an active Particle pointer in the mempool.
 *
 * @param pmp Pointer to the ParticleMempool. Must not be NULL.
 * @param node Pointer to the Particle to locate. Must not be NULL.
 *
 * @return Index of node in active array, or -1 if not found.
 *
 * @relates ParticleMempool
 */
long particle_mempool_index_of (ParticleMempool* pmp, Particle* particle) {
  if (!pmp || !particle)
    return -1;
  for (size_t i = 0; i < pmp->active.size; i++) {
    if (pmp->active.ptrs[i] == particle)
      return i;
  }
  return -1;
}

/*!
 * @brief Release and destroy an Particle by pointer.
 *
 * @param pmp Pointer to the ParticleMempool. Must not be NULL.
 * @param node Pointer to the Particle to release. Must not be NULL.
 *
 * @return 0 on success, -1 if the node is not found.
 *
 * @relates ParticleMempool
 */
int particle_mempool_free_particle_ptr (ParticleMempool* pmp, Particle* particle) {
  long idx = particle_mempool_index_of (pmp, particle);
  if (idx < 0)
    return -1;
  particle_mempool_free_particle (pmp, idx);
  return 0;
}

/*!
 * @brief Release and destroy an Particle at the given index.
 *
 * @param pmp Pointer to the ParticleMempool. Must not be NULL.
 * @param i Index of the node to release in the active array. Must be >= 0 and < size.
 *
 * @details Removes the node at index i from the active list using swap-with-last.
 * Calls particle_free() to clean up the node, then returns the memory to the pool
 * via mempool_free(). Index out of bounds causes program abort.
 *
 * @relates ParticleMempool
 */
void particle_mempool_free_particle (ParticleMempool* pmp, long i) {
  assert (i >= 0 && (size_t) i < pmp->active.size);
  Particle* dead = particle_list_remove_swap (&pmp->active, (size_t) i);
  particle_free (dead);
  mempool_free (pmp->pool, dead);
}
