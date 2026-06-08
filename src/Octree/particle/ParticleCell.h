#pragma once

#include "particle/Particle.h"
#include "particle/ParticleList.h"

// ============================================================================
// Type definitions
// ============================================================================

/*!
 * @struct ParticleCell
 */

typedef struct {
#if dimension == 1
  int i;
#elif dimension == 2
  int i, j;
#else
  int i, j, k;
#endif
  int level;
} ParticleCellPoint;

typedef struct {
  ParticleList particles;
  int pid;
} ParticleCell;

// ============================================================================
// Function declarations
// ============================================================================

int particle_cell_init(ParticleCell *pcell, size_t initial_capacity);
void particle_cell_free(ParticleCell *pcell);
int particle_cell_reserve(ParticleCell *pcell, size_t new_capacity);
int particle_cell_push(ParticleCell *pcell, Particle *particle);
Particle *particle_cell_remove_swap(ParticleCell *pcell, size_t i);
Particle *particle_cell_get(ParticleCell *pcell, size_t i);
void particle_cell_clear(ParticleCell *pcell);

ParticleCellPoint particle_cell_locate(coord c, int cell_level);
static inline bool particle_cell_point_wrap_or_reject(ParticleCellPoint *p);
static inline bool particle_cell_point_is_equal(ParticleCellPoint p1,
                                                ParticleCellPoint p2);
static inline bool particle_cell_point_is_after(ParticleCellPoint p1,
                                                ParticleCellPoint p2);

// ============================================================================
// Function definitions
// ============================================================================

/*!
 * @brief
 * @relates ParticleCell
 */
int particle_cell_init(ParticleCell *pcell, size_t initial_capacity) {
  return particle_list_init(&pcell->particles, initial_capacity);
}

/*!
 * @brief
 * @relates ParticleCell
 */
void particle_cell_free(ParticleCell *pcell) {
  particle_list_free(&pcell->particles);
}

/*!
 * @brief
 * @relates ParticleCell
 */
int particle_cell_reserve(ParticleCell *pcell, size_t new_capacity) {
  return particle_list_reserve(&pcell->particles, new_capacity);
}

/*!
 * @brief
 * @relates ParticleCell
 */
int particle_cell_push(ParticleCell *pcell, Particle *particle) {
  return particle_list_push(&pcell->particles, particle);
}

/*!
 * @brief
 * @relates ParticleCell
 */
Particle *particle_cell_remove_swap(ParticleCell *pcell, size_t i) {
  return particle_list_remove_swap(&pcell->particles, i);
}

/*!
 * @brief
 * @relates ParticleCell
 */
Particle *particle_cell_get(ParticleCell *pcell, size_t i) {
  return particle_list_get(&pcell->particles, i);
}

/*!
 * @brief
 * @relates ParticleCell
 */
void particle_cell_clear(ParticleCell *pcell) {
  particle_list_clear(&pcell->particles);
}

/*!
 * @brief
 * @relates ParticleCellPoint
 */
ParticleCellPoint particle_cell_locate(coord c, int level) {
  ParticleCellPoint pcell_point = {0};
  pcell_point.level = -1;

  const int n = 1 << level;
#if dimension == 1
  pcell_point.i = (c.x - X0) / domain_extent_x() * n + GHOSTS;
  if (pcell_point.i < GHOSTS || pcell_point.i >= n + GHOSTS)
    return pcell_point;
#elif dimension == 2
  pcell_point.i = (c.x - X0) / domain_extent_x() * n + GHOSTS;
  if (pcell_point.i < GHOSTS || pcell_point.i >= n + GHOSTS)
    return pcell_point;
  pcell_point.j = (c.y - Y0) / domain_extent_y() * n + GHOSTS;
  if (pcell_point.j < GHOSTS || pcell_point.j >= n + GHOSTS)
    return pcell_point;
#else
  pcell_point.i = (c.x - X0) / domain_extent_x() * n + GHOSTS;
  if (pcell_point.i < GHOSTS || pcell_point.i >= n + GHOSTS)
    return pcell_point;
  pcell_point.j = (c.y - Y0) / domain_extent_y() * n + GHOSTS;
  if (pcell_point.j < GHOSTS || pcell_point.j >= n + GHOSTS)
    return pcell_point;
  pcell_point.k = (c.z - Z0) / domain_extent_z() * n + GHOSTS;
  if (pcell_point.k < GHOSTS || pcell_point.k >= n + GHOSTS)
    return pcell_point;
#endif

  pcell_point.level = level;
  return pcell_point;
}

/*!
 * @brief Wraps a particle-cell point across periodic directions, or rejects it.
 *
 * @return true when the cell point is inside the valid real-cell range after
 * periodic wrapping; false when it lies outside a non-periodic direction.
 *
 * @relates ParticleCellPoint
 */
static inline bool particle_cell_point_wrap_or_reject(ParticleCellPoint *p) {
  if (!p || p->level < 0)
    return false;

  const int n = 1 << p->level;

#if dimension >= 1
  if (p->i < GHOSTS || p->i >= n + GHOSTS) {
    if (!Period.x)
      return false;
    int t = p->i - GHOSTS;
    t = (t % n + n) % n;
    p->i = t + GHOSTS;
  }
#endif
#if dimension >= 2
  if (p->j < GHOSTS || p->j >= n + GHOSTS) {
    if (!Period.y)
      return false;
    int t = p->j - GHOSTS;
    t = (t % n + n) % n;
    p->j = t + GHOSTS;
  }
#endif
#if dimension >= 3
  if (p->k < GHOSTS || p->k >= n + GHOSTS) {
    if (!Period.z)
      return false;
    int t = p->k - GHOSTS;
    t = (t % n + n) % n;
    p->k = t + GHOSTS;
  }
#endif

  return true;
}

/*!
 * @brief
 * @relates ParticleCellPoint
 */
static inline bool particle_cell_point_is_equal(ParticleCellPoint p1,
                                                ParticleCellPoint p2) {
#if dimension == 1
  if (p1.i != p2.i)
    return false;
#elif dimension == 2
  if (p1.i != p2.i)
    return false;
  if (p1.j != p2.j)
    return false;
#else
  if (p1.i != p2.i)
    return false;
  if (p1.j != p2.j)
    return false;
  if (p1.k != p2.k)
    return false;
#endif
  if (p1.level != p2.level)
    return false;

  return true;
}

/*!
 * @brief Tests whether @p p1 comes after @p p2 in cell-index order.
 *
 * @relates ParticleCellPoint
 */
static inline bool particle_cell_point_is_after(ParticleCellPoint p1,
                                                ParticleCellPoint p2) {
#if dimension >= 1
  if (p1.i != p2.i)
    return p1.i > p2.i;
#endif
#if dimension >= 2
  if (p1.j != p2.j)
    return p1.j > p2.j;
#endif
#if dimension >= 3
  if (p1.k != p2.k)
    return p1.k > p2.k;
#endif
  return p1.level > p2.level;
}
