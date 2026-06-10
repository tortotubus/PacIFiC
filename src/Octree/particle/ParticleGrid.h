

#include "particle/ParticleCell.h"
#include "particle/ParticleExchangeList.h"
#include "particle/ParticleMempool.h"
#if _MPI
#include <mpi.h>
#endif

// ============================================================================
// Type definitions
// ============================================================================

/**
 * @struct ParticleGrid
 */
typedef struct {
  ParticleMempool pool;
  bool dirty;
  int level;
#if dimension == 1
  ParticleCell *cells;
#elif dimension == 2
  ParticleCell **cells;
#else // dimension == 3
  ParticleCell ***cells;
#endif

#if _MPI

  ParticleExchangeList *snd_migrate;
  ParticleExchangeList *rcv_migrate;
  ParticleExchangeList *snd_boundary;
  ParticleExchangeList *rcv_boundary;

#if !TREE
  MPI_Comm cartcomm;
#endif
#endif
} ParticleGrid;

// ============================================================================
// Globals
// ============================================================================

ParticleGrid pg = {0};

// ============================================================================
// Function declarations
// ============================================================================

void particle_grid_init(int cell_level);
void particle_grid_free();
void particle_grid_update_cells();
Particle *particle_grid_add_particle(coord pos);
void particle_grid_delete_all_particles();
static inline int particle_grid_coord_pid(coord pos);

// ============================================================================
// Macros
// ============================================================================

#define PARTICLE_GRID_POOL_SIZE_BYTES (1 << 19)

/**
 * @def foreach_particle_cell
 * @brief Loops through all particle cells
 * @relates ParticleGrid
 */
macro foreach_particle_cell() {
  {
    const int n = (1 << pg.level) + 2 * GHOSTS;
#if dimension == 1
    for (int ci = GHOSTS - 1; ci < n - GHOSTS; ci++) {
      ParticleCellPoint pcell_point = {ci, pg.level};
      ParticleCell *pcell = &pg.cells[ci];
      // clang-format off
      {...}
      // clang-format on
    }
#elif dimension == 2
    for (int ci = GHOSTS; ci < n - GHOSTS; ci++) {
      for (int cj = GHOSTS; cj < n - GHOSTS; cj++) {
        ParticleCellPoint pcell_point = {ci, cj, pg.level};
        ParticleCell *pcell = &pg.cells[ci][cj];
        // clang-format off
        {...}
        // clang-format on
      }
    }
#else // dimension == 3
    for (int ci = GHOSTS; ci < n - GHOSTS; ci++) {
      for (int cj = GHOSTS; cj < n - GHOSTS; cj++) {
        for (int ck = GHOSTS; ck < n - GHOSTS; ck++) {
          ParticleCellPoint pcell_point = {ci, cj, ck, pg.level};
          ParticleCell *pcell = &pg.cells[ci][cj][ck];
          // clang-format off
          {...}
          // clang-format on
        }
      }
    }
#endif
  }
}

/**
 * @def foreach_particle
 * @brief Loops through all particles in all cells
 * @relates ParticleGrid
 */
macro foreach_particle() {
  {
    Point point = {0};
    int ig = 0, jg = 0, kg = 0;
    NOT_UNUSED(ig);
    NOT_UNUSED(jg);
    NOT_UNUSED(kg);

    foreach_particle_cell() {
      for (int cpid = 0; cpid < pcell->particles.size; cpid++) {
        Particle *particle = pcell->particles.ptrs[cpid];
        NOT_UNUSED(particle);
        PARTICLE_VARIABLES(particle);
        point = locate(pos.x, pos.y, pos.z);
        POINT_VARIABLES();

        // clang-format off
        {...}
        // clang-format on
      }
    }
  }
}

/**
 * @def foreach_particle_pair
 * @brief
 * @relates ParticleGrid
 */
macro foreach_particle_pair(int radius = 1) {
  {
    Point point = {0};
    int ig = 0, jg = 0, kg = 0;

    Particle *particle = NULL;
    Particle *particle_other = NULL;

    foreach_particle_cell() {
#if dimension == 1
#elif dimension == 2
      // In-cell pair traversal
      for (int p_idx = 0; p_idx < pcell->particles.size; p_idx++) {
        particle = pcell->particles.ptrs[p_idx];
        for (int p_other_idx = 0; p_other_idx < pcell->particles.size;
             p_other_idx++) {
          particle_other = pcell->particles.ptrs[p_other_idx];
          if (!particle->gid || !particle_other->gid ||
              particle->gid >= particle_other->gid)
            continue;
#if _MPI
          if (particle->pid != pid() && particle_other->pid != pid())
            continue;
#endif
          // clang-format off
          {...}
          // clang-format on
        }
      }

      // Neighbor-pair traversal
      for (int p_idx = 0; p_idx < pcell->particles.size; p_idx++) {
        particle = pcell->particles.ptrs[p_idx];
        for (int pci = -radius; pci <= radius; pci++) {
          for (int pcj = -radius; pcj <= radius; pcj++) {
            ParticleCellPoint pcell_other_point = {
                pcell_point.i + pci, pcell_point.j + pcj, pcell_point.level};
            if (pci == 0 && pcj == 0)
              continue;
            if (!particle_cell_point_wrap_or_reject(&pcell_other_point))
              continue;
            ParticleCell *pcell_other =
                &pg.cells[pcell_other_point.i][pcell_other_point.j];
            for (int p_other_idx = 0; p_other_idx < pcell_other->particles.size;
                 p_other_idx++) {
              particle_other = pcell_other->particles.ptrs[p_other_idx];
              if (!particle->gid || !particle_other->gid ||
                  particle->gid >= particle_other->gid)
                continue;
#if _MPI
              if (particle->pid != pid() && particle_other->pid != pid())
                continue;
#endif
              // clang-format off
              {...}
              // clang-format on
            }
          }
        }
      }
#else // dimension == 3
      // In-cell pair traversal
      for (int p_idx = 0; p_idx < pcell->particles.size; p_idx++) {
        particle = pcell->particles.ptrs[p_idx];
        for (int p_other_idx = 0; p_other_idx < pcell->particles.size;
             p_other_idx++) {
          particle_other = pcell->particles.ptrs[p_other_idx];
          if (!particle->gid || !particle_other->gid ||
              particle->gid >= particle_other->gid)
            continue;
#if _MPI
          if (particle->pid != pid() && particle_other->pid != pid())
            continue;
#endif
          // clang-format off
          {...}
          // clang-format on
        }
      }

      // Neighbor-pair traversal
      for (int p_idx = 0; p_idx < pcell->particles.size; p_idx++) {
        particle = pcell->particles.ptrs[p_idx];
        for (int pci = -radius; pci <= radius; pci++) {
          for (int pcj = -radius; pcj <= radius; pcj++) {
            for (int pck = -radius; pck <= radius; pck++) {
              if (pci == 0 && pcj == 0 && pck == 0)
                continue;

              ParticleCellPoint pcell_other_point = {
                  pcell_point.i + pci, pcell_point.j + pcj, pcell_point.k + pck,
                  pcell_point.level};
              if (!particle_cell_point_wrap_or_reject(&pcell_other_point))
                continue;

              ParticleCell *pcell_other =
                  &pg.cells[pcell_other_point.i][pcell_other_point.j]
                           [pcell_other_point.k];
              for (int p_other_idx = 0;
                   p_other_idx < pcell_other->particles.size; p_other_idx++) {
                particle_other = pcell_other->particles.ptrs[p_other_idx];
                if (!particle->gid || !particle_other->gid ||
                    particle->gid >= particle_other->gid)
                  continue;
#if _MPI
                if (particle->pid != pid() && particle_other->pid != pid())
                  continue;
#endif
                // clang-format off
                {...}
                // clang-format on
              }
            }
          }
        }
      }
#endif
    }
  }
}

// ============================================================================
// Function definitions
// ============================================================================

/**
 * @brief Return the MPI rank that owns @p pos, or -1 if @p pos is outside the
 * physical domain.
 *
 * @relates ParticleGrid
 */
static inline int particle_grid_coord_pid(coord pos) {
#if _MPI
  coord wpos = pos;
  domain_wrap_coord(wpos);

#if TREE
  {
    Point point = domain_locate_nonlocal(wpos.x, wpos.y, wpos.z);
    if (point.level < 0)
      return -1;

    int ig = 0, jg = 0, kg = 0;
    NOT_UNUSED(ig);
    NOT_UNUSED(jg);
    NOT_UNUSED(kg);
    POINT_VARIABLES();

    return allocated(0) ? cell.pid : -1;
  }
#else // !TREE
  {
    Point point = {0};
    point.level = depth();
    SET_DIMENSIONS();

    int coords[dimension];
    const int local_nx = point.n.x;
    const int global_nx = local_nx * Dimensions.x;
    const double Lx = domain_extent_x();

    double px = wpos.x;
    if (Period.x)
      px = X0 + fmod(fmod(px - X0, Lx) + Lx, Lx);

    int gi = (int)floor((px - X0) / Lx * global_nx);
    if (gi < 0 || gi >= global_nx) {
      if (!Period.x)
        return -1;
      gi = (gi % global_nx + global_nx) % global_nx;
    }
    if (gi == global_nx)
      gi = global_nx - 1;
    coords[0] = gi / local_nx;

#if dimension > 1
    const int local_ny = point.n.y;
    const int global_ny = local_ny * Dimensions.y;
    const double Ly = domain_extent_y();

    double py = wpos.y;
    if (Period.y)
      py = Y0 + fmod(fmod(py - Y0, Ly) + Ly, Ly);

    int gj = (int)floor((py - Y0) / Ly * global_ny);
    if (gj < 0 || gj >= global_ny) {
      if (!Period.y)
        return -1;
      gj = (gj % global_ny + global_ny) % global_ny;
    }
    if (gj == global_ny)
      gj = global_ny - 1;
    coords[1] = gj / local_ny;
#endif

#if dimension > 2
    const int local_nz = point.n.z;
    const int global_nz = local_nz * Dimensions.z;
    const double Lz = domain_extent_z();

    double pz = wpos.z;
    if (Period.z)
      pz = Z0 + fmod(fmod(pz - Z0, Lz) + Lz, Lz);

    int gk = (int)floor((pz - Z0) / Lz * global_nz);
    if (gk < 0 || gk >= global_nz) {
      if (!Period.z)
        return -1;
      gk = (gk % global_nz + global_nz) % global_nz;
    }
    if (gk == global_nz)
      gk = global_nz - 1;
    coords[2] = gk / local_nz;
#endif

    int owner = -1;
    MPI_Cart_rank(pg.cartcomm, coords, &owner);
    return owner;
  }
#endif // !TREE
#else  // !_MPI
  NOT_UNUSED(pos);
  return 0;
#endif // !_MPI
}

/**
 * @brief Initialize the particle grid
 * @param particle_count The number of particles you plan to have initially
 * @param level The basilisk grid depth where the particle cells live at
 * @relates ParticleGrid
 */
void particle_grid_init(int level) {
  if (npvar == 0)
    init_psolver();

#if _MPI
#if !TREE
  {
    int min_level = 0;
    int max_dim = Dimensions.x;
#if dimension > 1
    max_dim = max(max_dim, Dimensions.y);
#endif
#if dimension > 2
    max_dim = max(max_dim, Dimensions.z);
#endif

    while ((1 << min_level) < max_dim)
      min_level++;

    assert(level >= min_level);
  }
#endif
#endif

  pg.level = level;
  const int n = (1 << pg.level) + 2 * GHOSTS;
  pg.dirty = true;
  pg.pool = particle_mempool_init(PARTICLE_GRID_POOL_SIZE_BYTES,
                                  npvar * sizeof(double));

#if dimension == 1
  pg.cells = (ParticleCell *)calloc((size_t)n, sizeof(*pg.cells));
  assert(pg.cells);
#elif dimension == 2
  pg.cells = (ParticleCell **)malloc((size_t)n * sizeof(*pg.cells));
  assert(pg.cells);
  pg.cells[0] =
      (ParticleCell *)calloc((size_t)n * (size_t)n, sizeof(**pg.cells));
  assert(pg.cells[0]);
  for (int ci = 1; ci < n; ci++)
    pg.cells[ci] = pg.cells[0] + (size_t)ci * (size_t)n;
#else // dimension == 3
  pg.cells = (ParticleCell ***)malloc((size_t)n * sizeof(*pg.cells));
  assert(pg.cells);
  pg.cells[0] =
      (ParticleCell **)malloc((size_t)n * (size_t)n * sizeof(**pg.cells));
  assert(pg.cells[0]);
  pg.cells[0][0] = (ParticleCell *)calloc((size_t)n * (size_t)n * (size_t)n,
                                          sizeof(***pg.cells));
  assert(pg.cells[0][0]);

  for (int ci = 0; ci < n; ci++)
    pg.cells[ci] = pg.cells[0] + (size_t)ci * (size_t)n;
  for (int ci = 0; ci < n; ci++)
    for (int cj = 0; cj < n; cj++)
      pg.cells[ci][cj] =
          pg.cells[0][0] + ((size_t)ci * (size_t)n + (size_t)cj) * (size_t)n;
#endif

  foreach_particle_cell() { particle_cell_init(pcell, 0); }

#if _MPI
  pg.snd_boundary =
      (ParticleExchangeList *)calloc(npe(), sizeof(ParticleExchangeList));
  pg.rcv_boundary =
      (ParticleExchangeList *)calloc(npe(), sizeof(ParticleExchangeList));
  pg.snd_migrate =
      (ParticleExchangeList *)calloc(npe(), sizeof(ParticleExchangeList));
  pg.rcv_migrate =
      (ParticleExchangeList *)calloc(npe(), sizeof(ParticleExchangeList));

  for (int i = 0; i < npe(); i++) {
    particle_exchange_list_init(&pg.snd_boundary[i], 0);
    particle_exchange_list_set_pid(&pg.snd_boundary[i], i);
    particle_exchange_list_init(&pg.rcv_boundary[i], 0);
    particle_exchange_list_set_pid(&pg.rcv_boundary[i], i);
    particle_exchange_list_init(&pg.snd_migrate[i], 0);
    particle_exchange_list_set_pid(&pg.snd_migrate[i], i);
    particle_exchange_list_init(&pg.rcv_migrate[i], 0);
    particle_exchange_list_set_pid(&pg.rcv_migrate[i], i);
  }
#if !TREE

  int dims[dimension];
  int periods[dimension];
  dims[0] = Dimensions.x;
  periods[0] = Period.x;
#if dimension >= 2
  dims[1] = Dimensions.y;
  periods[1] = Period.y;
#endif
#if dimension >= 3
  dims[2] = Dimensions.z;
  periods[2] = Period.z;
#endif
  MPI_Cart_create(MPI_COMM_WORLD, dimension, dims, periods, 0, &pg.cartcomm);
#endif
#endif
}

/**
 * @brief Free the particle grid and its owned members.
 * @relates ParticleGrid
 */
void particle_grid_free() {
  if (pg.cells) {
    foreach_particle_cell() { particle_cell_free(pcell); }
  }

#if dimension == 1
  free(pg.cells);
  pg.cells = NULL;
#elif dimension == 2
  if (pg.cells) {
    free(pg.cells[0]);
    free(pg.cells);
    pg.cells = NULL;
  }
#else // dimension == 3
  if (pg.cells) {
    free(pg.cells[0][0]);
    free(pg.cells[0]);
    free(pg.cells);
    pg.cells = NULL;
  }
#endif

#if _MPI
  if (pg.snd_boundary) {
    for (int i = 0; i < npe(); i++)
      particle_exchange_list_free(&pg.snd_boundary[i]);
    free(pg.snd_boundary);
    pg.snd_boundary = NULL;
  }
  if (pg.rcv_boundary) {
    for (int i = 0; i < npe(); i++)
      particle_exchange_list_free(&pg.rcv_boundary[i]);
    free(pg.rcv_boundary);
    pg.rcv_boundary = NULL;
  }
  if (pg.snd_migrate) {
    for (int i = 0; i < npe(); i++)
      particle_exchange_list_free(&pg.snd_migrate[i]);
    free(pg.snd_migrate);
    pg.snd_migrate = NULL;
  }
  if (pg.rcv_migrate) {
    for (int i = 0; i < npe(); i++)
      particle_exchange_list_free(&pg.rcv_migrate[i]);
    free(pg.rcv_migrate);
    pg.rcv_migrate = NULL;
  }
#if !TREE
  if (pg.cartcomm != MPI_COMM_NULL)
    MPI_Comm_free(&pg.cartcomm);
#endif
#endif

  if (pg.pool.pool)
    particle_mempool_free(&pg.pool);

  pg.level = 0;
  pg.dirty = false;
}

/**
 * @brief Add a single particle to the grid
 * @relates ParticleGrid
 */
Particle *particle_grid_add_particle(coord pos) {

  // Check that the particle mempool actually exists first
  ParticleMempool *pmp = &pg.pool;
  if (!pmp->pool)
    return NULL;

  coord wpos = pos;
  domain_wrap_coord(wpos);

  ParticleCellPoint pcp = particle_cell_locate(wpos, pg.level);
  if (pcp.level < 0)
    return NULL;

  const int owner = particle_grid_coord_pid(wpos);
  if (owner < 0)
    return NULL;
#if _MPI
  if (owner != pid())
    return NULL;
#endif

  // Create the particle
  Particle *particle = particle_mempool_alloc_particle(pmp);
  assert(particle);
  particle_init(particle);
  particle->pid = owner;
  foreach_dimension() { pval(ppos.x) = wpos.x; }
  particle_set_gid(particle);

#if dimension == 1
  ParticleCell *pcell = &pg.cells[pcp.i];
#elif dimension == 2
  ParticleCell *pcell = &pg.cells[pcp.i][pcp.j];
#else
  ParticleCell *pcell = &pg.cells[pcp.i][pcp.j][pcp.k];
#endif
  assert(particle_cell_push(pcell, particle) == 0);
  return particle;
}

/**
 * @brief Deletes all particles from the grid
 * @relates ParticleGrid
 */
void particle_grid_delete_all_particles() {
  if (!pg.pool.pool)
    return;

  foreach_particle_cell() { particle_cell_clear(pcell); }

  while (pg.pool.active.size > 0)
    particle_mempool_free_particle(&pg.pool, (long)pg.pool.active.size - 1);
}

/**
 * @brief Rebuild cell lists from current particle positions
 * @relates ParticleGrid
 */
void particle_grid_update_cells() {
  if (!pg.cells)
    return;

  foreach_particle_cell() { particle_cell_clear(pcell); }

  for (size_t idx = 0; idx < pg.pool.active.size; idx++) {
    Particle *particle = pg.pool.active.ptrs[idx];
    coord pos = {0};
    foreach_dimension() { pos.x = pval(ppos.x); }
    domain_wrap_coord(pos);

    ParticleCellPoint pcp = particle_cell_locate(pos, pg.level);
    if (pcp.level < 0)
      continue;

#if dimension == 1
    ParticleCell *pcell = &pg.cells[pcp.i];
#elif dimension == 2
    ParticleCell *pcell = &pg.cells[pcp.i][pcp.j];
#else
    ParticleCell *pcell = &pg.cells[pcp.i][pcp.j][pcp.k];
#endif
    particle_cell_push(pcell, particle);
  }
}
