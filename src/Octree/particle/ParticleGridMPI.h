

#if _MPI
#include <mpi.h>
#endif
#include <limits.h>
#include <string.h>

#include "particle/ParticleGrid.h"

// ============================================================================
// Function declarations
// ============================================================================

void particle_grid_update_pid();
void particle_grid_clear_boundary();
void particle_grid_update_boundary_broadcast_cells(int radius);
void particle_boundary(Pscalar *fields = pall);

static inline size_t _particle_grid_migration_record_size();
static inline size_t _particle_boundary_record_size(int field_count);
static inline void _particle_grid_pack_particle(unsigned char *dst,
                                                Particle *particle);
static inline void _particle_grid_unpack_particle(Particle *particle,
                                                  const unsigned char *src);
static inline void _particle_boundary_pack_particle_fields(unsigned char *dst,
                                                           Particle *particle,
                                                           Pscalar *fields);
static inline void _particle_boundary_unpack_particle_fields(
    Particle *particle, const unsigned char *src, Pscalar *fields);
static inline Particle *_particle_grid_find_boundary_ghost(int owner_pid,
                                                           uint64_t gid);
static inline ParticleCell *
_particle_grid_cell_from_point(ParticleCellPoint pcp);
int _particle_grid_get_cell_pid(ParticleCellPoint pcp);
int _particle_grid_get_pid(Particle *particle);

// ============================================================================
// Type definitions
// ============================================================================

typedef struct {
  ParticleCellPoint *ptrs;
  size_t size;
  size_t capacity;
} ParticleCellPointList;

typedef struct {
  int requester_pid;
  ParticleCellPoint pcell_point;
} ParticleCellRequest;

// ============================================================================
// Private helper declarations
// ============================================================================

static inline void _particle_cell_point_list_init(ParticleCellPointList *list);
static inline void _particle_cell_point_list_free(ParticleCellPointList *list);
static inline int _particle_cell_point_list_reserve(ParticleCellPointList *list,
                                                    size_t capacity);
static inline int _particle_cell_point_list_push(ParticleCellPointList *list,
                                                 ParticleCellPoint point);
static inline int
_particle_cell_point_list_push_unique(ParticleCellPointList *list,
                                      ParticleCellPoint point);
static inline void
_particle_grid_build_cell_boundary_requests(ParticleCellPointList *requests,
                                            int radius);
static inline ParticleCellRequest *
_particle_grid_allgather_cell_requests(ParticleCellPointList *local_requests,
                                       int *global_request_count);
static inline void _particle_grid_build_snd_boundary_from_cell_requests(
    const ParticleCellRequest *requests, int request_count);
static inline void _particle_grid_exchange_boundary_particles();

// ============================================================================
// Function definitions
// ============================================================================

/*!
 * @brief
 */
static inline size_t _particle_grid_migration_record_size() {
  return sizeof(Particle) + pg.pool.datasize;
}

/*!
 * @brief Return the packed byte size for one particle boundary field record.
 */
static inline size_t _particle_boundary_record_size(int field_count) {
  assert(field_count >= 0);
  return sizeof(uint64_t) + (size_t)field_count * sizeof(double);
}

/*!
 * @brief
 */
static inline void _particle_grid_pack_particle(unsigned char *dst,
                                                Particle *particle) {
  assert(dst);
  assert(particle);
  memcpy(dst, particle, _particle_grid_migration_record_size());
}

/*!
 * @brief
 */
static inline void _particle_grid_unpack_particle(Particle *particle,
                                                  const unsigned char *src) {
  assert(particle);
  assert(src);
  memcpy(particle, src, _particle_grid_migration_record_size());
}

/*!
 * @brief Pack a particle gid and selected scalar values for ghost refresh.
 */
static inline void _particle_boundary_pack_particle_fields(unsigned char *dst,
                                                           Particle *particle,
                                                           Pscalar *fields) {
  assert(dst);
  assert(particle);
  assert(fields);

  unsigned char *cursor = dst;
  memcpy(cursor, &particle->gid, sizeof(particle->gid));
  cursor += sizeof(particle->gid);

  foreach_pscalar(fields) {
    const double value = *_pval(s, particle);
    memcpy(cursor, &value, sizeof(value));
    cursor += sizeof(value);
  }
}

/*!
 * @brief Unpack selected scalar values into an existing ghost particle.
 */
static inline void _particle_boundary_unpack_particle_fields(
    Particle *particle, const unsigned char *src, Pscalar *fields) {
  assert(particle);
  assert(src);
  assert(fields);

  const unsigned char *cursor = src + sizeof(uint64_t);
  foreach_pscalar(fields) {
    double value = 0.;
    memcpy(&value, cursor, sizeof(value));
    *_pval(s, particle) = value;
    cursor += sizeof(value);
  }
}

/*!
 * @brief Return the local ghost copy from owner rank @p owner_pid with @p gid.
 */
static inline Particle *_particle_grid_find_boundary_ghost(int owner_pid,
                                                           uint64_t gid) {
#if _MPI
  if (!pg.rcv_boundary || owner_pid < 0 || owner_pid >= npe())
    return NULL;

  ParticleExchangeList *list = &pg.rcv_boundary[owner_pid];
  for (size_t idx = 0; idx < list->particles.size; idx++) {
    Particle *particle = list->particles.ptrs[idx];
    if (particle && particle->gid == gid)
      return particle;
  }
#else
  NOT_UNUSED(owner_pid);
  NOT_UNUSED(gid);
#endif
  return NULL;
}

/*!
 * @brief Return the particle cell at @p pcp, or NULL if it is outside this
 * rank's replicated particle-cell array.
 */
static inline ParticleCell *
_particle_grid_cell_from_point(ParticleCellPoint pcp) {
  if (!pg.cells || pcp.level != pg.level)
    return NULL;

  const int n = (1 << pg.level) + 2 * GHOSTS;
#if dimension >= 1
  if (pcp.i < 0 || pcp.i >= n)
    return NULL;
#endif
#if dimension >= 2
  if (pcp.j < 0 || pcp.j >= n)
    return NULL;
#endif
#if dimension >= 3
  if (pcp.k < 0 || pcp.k >= n)
    return NULL;
#endif

#if dimension == 1
  return &pg.cells[pcp.i];
#elif dimension == 2
  return &pg.cells[pcp.i][pcp.j];
#else
  return &pg.cells[pcp.i][pcp.j][pcp.k];
#endif
}

static inline void _particle_cell_point_list_init(ParticleCellPointList *list) {
  assert(list);
  list->ptrs = NULL;
  list->size = 0;
  list->capacity = 0;
}

static inline void _particle_cell_point_list_free(ParticleCellPointList *list) {
  if (!list)
    return;
  free(list->ptrs);
  list->ptrs = NULL;
  list->size = 0;
  list->capacity = 0;
}

static inline int _particle_cell_point_list_reserve(ParticleCellPointList *list,
                                                    size_t capacity) {
  assert(list);
  if (capacity <= list->capacity)
    return 0;
  ParticleCellPoint *ptrs =
      (ParticleCellPoint *)realloc(list->ptrs, capacity * sizeof(*ptrs));
  if (!ptrs)
    return -1;
  list->ptrs = ptrs;
  list->capacity = capacity;
  return 0;
}

static inline int
_particle_cell_point_list_push_unique(ParticleCellPointList *list,
                                      ParticleCellPoint point) {
  assert(list);
  for (size_t idx = 0; idx < list->size; idx++) {
    if (particle_cell_point_is_equal(list->ptrs[idx], point))
      return 0;
  }
  if (list->size == list->capacity) {
    const size_t capacity = list->capacity ? 2 * list->capacity : 16;
    if (_particle_cell_point_list_reserve(list, capacity) != 0)
      return -1;
  }
  list->ptrs[list->size++] = point;
  return 0;
}

static inline int _particle_cell_point_list_push(ParticleCellPointList *list,
                                                 ParticleCellPoint point) {
  assert(list);
  if (list->size == list->capacity) {
    const size_t capacity = list->capacity ? 2 * list->capacity : 16;
    if (_particle_cell_point_list_reserve(list, capacity) != 0)
      return -1;
  }
  list->ptrs[list->size++] = point;
  return 0;
}

trace static inline void
_particle_grid_build_cell_boundary_requests(ParticleCellPointList *requests,
                                            int radius) {
  assert(requests);
  assert(radius >= 0);

  const int n_real = 1 << pg.level;
  size_t mark_count = (size_t)n_real;
#if dimension >= 2
  mark_count *= (size_t)n_real;
#endif
#if dimension >= 3
  mark_count *= (size_t)n_real;
#endif
  unsigned char *requested =
      mark_count > 0 ? (unsigned char *)calloc(mark_count, sizeof(*requested))
                     : NULL;
  assert(mark_count == 0 || requested);

  foreach_particle_cell() {
    bool has_owned_particle = false;
    for (size_t idx = 0; idx < pcell->particles.size; idx++) {
      Particle *particle = pcell->particles.ptrs[idx];
      if (particle && particle->pid == pid()) {
        has_owned_particle = true;
        break;
      }
    }
    if (!has_owned_particle)
      continue;

#if dimension == 1
    for (int di = -radius; di <= radius; di++) {
      ParticleCellPoint point = {pcell_point.i + di, pcell_point.level};
      if (!particle_cell_point_wrap_or_reject(&point))
        continue;
      size_t mark = (size_t)(point.i - GHOSTS);
      if (!requested[mark]) {
        requested[mark] = 1;
        assert(_particle_cell_point_list_push(requests, point) == 0);
      }
    }
#elif dimension == 2
    for (int di = -radius; di <= radius; di++) {
      for (int dj = -radius; dj <= radius; dj++) {
        ParticleCellPoint point = {pcell_point.i + di, pcell_point.j + dj,
                                   pcell_point.level};
        if (!particle_cell_point_wrap_or_reject(&point))
          continue;
        size_t mark = ((size_t)(point.i - GHOSTS) * (size_t)n_real +
                       (size_t)(point.j - GHOSTS));
        if (!requested[mark]) {
          requested[mark] = 1;
          assert(_particle_cell_point_list_push(requests, point) == 0);
        }
      }
    }
#else
    for (int di = -radius; di <= radius; di++) {
      for (int dj = -radius; dj <= radius; dj++) {
        for (int dk = -radius; dk <= radius; dk++) {
          ParticleCellPoint point = {pcell_point.i + di, pcell_point.j + dj,
                                     pcell_point.k + dk, pcell_point.level};
          if (!particle_cell_point_wrap_or_reject(&point))
            continue;
          size_t mark = (((size_t)(point.i - GHOSTS) * (size_t)n_real +
                          (size_t)(point.j - GHOSTS)) *
                             (size_t)n_real +
                         (size_t)(point.k - GHOSTS));
          if (!requested[mark]) {
            requested[mark] = 1;
            assert(_particle_cell_point_list_push(requests, point) == 0);
          }
        }
      }
    }
#endif
  }

  free(requested);
}

trace static inline ParticleCellRequest *
_particle_grid_allgather_cell_requests(ParticleCellPointList *local_requests,
                                       int *global_request_count) {
  assert(local_requests);
  assert(global_request_count);
  *global_request_count = 0;

#if _MPI
  const int nproc = npe();
  assert(local_requests->size <= (size_t)INT_MAX);
  const int local_count = (int)local_requests->size;

  int *recv_counts = (int *)calloc((size_t)nproc, sizeof(*recv_counts));
  int *recv_bytes = (int *)calloc((size_t)nproc, sizeof(*recv_bytes));
  int *recv_displs = (int *)calloc((size_t)nproc, sizeof(*recv_displs));
  assert(recv_counts);
  assert(recv_bytes);
  assert(recv_displs);

  MPI_Allgather(&local_count, 1, MPI_INT, recv_counts, 1, MPI_INT,
                MPI_COMM_WORLD);

  size_t total_bytes = 0;
  int total_count = 0;
  for (int peer = 0; peer < nproc; peer++) {
    const size_t peer_bytes =
        (size_t)recv_counts[peer] * sizeof(ParticleCellRequest);
    assert(peer_bytes <= (size_t)INT_MAX);
    assert(total_bytes <= (size_t)INT_MAX);
    recv_bytes[peer] = (int)peer_bytes;
    recv_displs[peer] = (int)total_bytes;
    total_bytes += peer_bytes;
    total_count += recv_counts[peer];
  }

  ParticleCellRequest *local_buffer =
      local_count > 0 ? (ParticleCellRequest *)malloc((size_t)local_count *
                                                      sizeof(*local_buffer))
                      : NULL;
  assert(local_count == 0 || local_buffer);
  for (int idx = 0; idx < local_count; idx++) {
    local_buffer[idx].requester_pid = pid();
    local_buffer[idx].pcell_point = local_requests->ptrs[idx];
  }

  ParticleCellRequest *global_buffer =
      total_count > 0 ? (ParticleCellRequest *)malloc(total_bytes) : NULL;
  assert(total_count == 0 || global_buffer);

  MPI_Allgatherv(local_buffer, local_count * (int)sizeof(ParticleCellRequest),
                 MPI_BYTE, global_buffer, recv_bytes, recv_displs, MPI_BYTE,
                 MPI_COMM_WORLD);

  free(local_buffer);
  free(recv_displs);
  free(recv_bytes);
  free(recv_counts);

  *global_request_count = total_count;
  return global_buffer;
#else
  return NULL;
#endif
}

trace static inline void _particle_grid_build_snd_boundary_from_cell_requests(
    const ParticleCellRequest *requests, int request_count) {
#if _MPI
  if (!requests || request_count <= 0)
    return;

  const int nproc = npe();
  for (int idx = 0; idx < request_count; idx++) {
    const int requester = requests[idx].requester_pid;
    if (requester < 0 || requester >= nproc || requester == pid())
      continue;

    ParticleCell *requested_pcell =
        _particle_grid_cell_from_point(requests[idx].pcell_point);
    if (!requested_pcell)
      continue;

    for (size_t pidx = 0; pidx < requested_pcell->particles.size; pidx++) {
      Particle *particle = requested_pcell->particles.ptrs[pidx];
      if (particle && particle->pid == pid())
        particle_exchange_list_push_unique(&pg.snd_boundary[requester],
                                           particle);
    }
  }
#endif
}

trace static inline void _particle_grid_exchange_boundary_particles() {
#if _MPI
  const int nproc = npe();
  const size_t record_size = _particle_grid_migration_record_size();
  assert(record_size <= (size_t)INT_MAX);

  int *send_counts = (int *)calloc((size_t)nproc, sizeof(*send_counts));
  int *recv_counts = (int *)calloc((size_t)nproc, sizeof(*recv_counts));
  int *send_bytes = (int *)calloc((size_t)nproc, sizeof(*send_bytes));
  int *recv_bytes = (int *)calloc((size_t)nproc, sizeof(*recv_bytes));
  int *send_displs = (int *)calloc((size_t)nproc, sizeof(*send_displs));
  int *recv_displs = (int *)calloc((size_t)nproc, sizeof(*recv_displs));
  assert(send_counts);
  assert(recv_counts);
  assert(send_bytes);
  assert(recv_bytes);
  assert(send_displs);
  assert(recv_displs);

  for (int peer = 0; peer < nproc; peer++) {
    assert(pg.snd_boundary[peer].particles.size <= (size_t)INT_MAX);
    send_counts[peer] = (int)pg.snd_boundary[peer].particles.size;
  }

  MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT,
               MPI_COMM_WORLD);

  size_t total_send_bytes = 0;
  size_t total_recv_bytes = 0;
  for (int peer = 0; peer < nproc; peer++) {
    const size_t peer_send_bytes = (size_t)send_counts[peer] * record_size;
    const size_t peer_recv_bytes = (size_t)recv_counts[peer] * record_size;
    assert(peer_send_bytes <= (size_t)INT_MAX);
    assert(peer_recv_bytes <= (size_t)INT_MAX);
    assert(total_send_bytes <= (size_t)INT_MAX);
    assert(total_recv_bytes <= (size_t)INT_MAX);

    send_bytes[peer] = (int)peer_send_bytes;
    recv_bytes[peer] = (int)peer_recv_bytes;
    send_displs[peer] = (int)total_send_bytes;
    recv_displs[peer] = (int)total_recv_bytes;

    total_send_bytes += peer_send_bytes;
    total_recv_bytes += peer_recv_bytes;
  }

  unsigned char *send_buffer =
      total_send_bytes > 0 ? (unsigned char *)malloc(total_send_bytes) : NULL;
  unsigned char *recv_buffer =
      total_recv_bytes > 0 ? (unsigned char *)malloc(total_recv_bytes) : NULL;
  assert(total_send_bytes == 0 || send_buffer);
  assert(total_recv_bytes == 0 || recv_buffer);

  for (int peer = 0; peer < nproc; peer++) {
    if (send_counts[peer] == 0)
      continue;
    unsigned char *dst = send_buffer + send_displs[peer];
    for (size_t idx = 0; idx < pg.snd_boundary[peer].particles.size; idx++) {
      Particle *particle = pg.snd_boundary[peer].particles.ptrs[idx];
      _particle_grid_pack_particle(dst + idx * record_size, particle);
    }
  }

  MPI_Alltoallv(send_buffer, send_bytes, send_displs, MPI_BYTE, recv_buffer,
                recv_bytes, recv_displs, MPI_BYTE, MPI_COMM_WORLD);

  for (int peer = 0; peer < nproc; peer++) {
    if (recv_counts[peer] == 0)
      continue;
    const unsigned char *src = recv_buffer + recv_displs[peer];
    for (int idx = 0; idx < recv_counts[peer]; idx++) {
      Particle *particle = particle_mempool_alloc_particle(&pg.pool);
      _particle_grid_unpack_particle(particle, src + (size_t)idx * record_size);
      particle_exchange_list_push(&pg.rcv_boundary[peer], particle);
    }
  }

  free(recv_buffer);
  free(send_buffer);
  free(recv_displs);
  free(send_displs);
  free(recv_bytes);
  free(send_bytes);
  free(recv_counts);
  free(send_counts);
#endif
}

/*!
 * @brief Remove ghost particles received by the previous boundary exchange and
 * clear boundary send/receive lists.
 */
trace void particle_grid_clear_boundary() {
#if _MPI
  if (!pg.pool.pool)
    return;

  for (int peer = 0; peer < npe(); peer++) {
    for (size_t idx = 0; idx < pg.rcv_boundary[peer].particles.size; idx++) {
      Particle *particle = pg.rcv_boundary[peer].particles.ptrs[idx];
      if (particle)
        particle_mempool_free_particle_ptr(&pg.pool, particle);
    }
    particle_exchange_list_clear(&pg.rcv_boundary[peer]);
    particle_exchange_list_clear(&pg.snd_boundary[peer]);
  }
#endif
}

/*!
 * @brief Build a request-driven ghost halo by broadcasting requested particle
 * cells and exchanging owned particles in those cells.
 *
 * Each rank requests all particle cells within @p radius of any cell containing
 * an owned particle. Requests are allgathered so tree cases where several ranks
 * contribute particles to the same particle cell do not need a unique cell
 * owner.
 */
trace void particle_grid_update_boundary_broadcast_cells(int radius) {
#if _MPI
  if (!pg.pool.pool)
    return;

  particle_grid_clear_boundary();
  particle_grid_update_cells();

  ParticleCellPointList local_requests = {0};
  _particle_cell_point_list_init(&local_requests);
  _particle_grid_build_cell_boundary_requests(&local_requests, radius);

  int global_request_count = 0;
  ParticleCellRequest *global_requests = _particle_grid_allgather_cell_requests(
      &local_requests, &global_request_count);

  _particle_grid_build_snd_boundary_from_cell_requests(global_requests,
                                                       global_request_count);
  _particle_grid_exchange_boundary_particles();

  free(global_requests);
  _particle_cell_point_list_free(&local_requests);

  particle_grid_update_cells();
  pg.dirty = false;
#else
  NOT_UNUSED(radius);
#endif
}

/*!
 * @brief Refresh selected ghost-particle fields from owner-rank values.
 *
 * This uses the current boundary exchange lists. Call
 * @ref particle_grid_update_boundary_broadcast_cells first to create ghosts and
 * populate @c pg.snd_boundary / @c pg.rcv_boundary. The selected fields are
 * copied from each owner particle to all ghost copies; values are overwritten,
 * not summed.
 */
trace void particle_boundary(Pscalar *fields) {
#if _MPI
  if (!pg.pool.pool || !fields)
    return;

  const int field_count = plist_len(fields);
  if (field_count <= 0)
    return;

  const int nproc = npe();
  const size_t record_size = _particle_boundary_record_size(field_count);
  assert(record_size <= (size_t)INT_MAX);

  int *send_counts = (int *)calloc((size_t)nproc, sizeof(*send_counts));
  int *recv_counts = (int *)calloc((size_t)nproc, sizeof(*recv_counts));
  int *send_bytes = (int *)calloc((size_t)nproc, sizeof(*send_bytes));
  int *recv_bytes = (int *)calloc((size_t)nproc, sizeof(*recv_bytes));
  int *send_displs = (int *)calloc((size_t)nproc, sizeof(*send_displs));
  int *recv_displs = (int *)calloc((size_t)nproc, sizeof(*recv_displs));
  assert(send_counts);
  assert(recv_counts);
  assert(send_bytes);
  assert(recv_bytes);
  assert(send_displs);
  assert(recv_displs);

  for (int peer = 0; peer < nproc; peer++) {
    assert(pg.snd_boundary[peer].particles.size <= (size_t)INT_MAX);
    send_counts[peer] = (int)pg.snd_boundary[peer].particles.size;
  }

  MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT,
               MPI_COMM_WORLD);

  size_t total_send_bytes = 0;
  size_t total_recv_bytes = 0;
  for (int peer = 0; peer < nproc; peer++) {
    const size_t peer_send_bytes = (size_t)send_counts[peer] * record_size;
    const size_t peer_recv_bytes = (size_t)recv_counts[peer] * record_size;
    assert(peer_send_bytes <= (size_t)INT_MAX);
    assert(peer_recv_bytes <= (size_t)INT_MAX);
    assert(total_send_bytes <= (size_t)INT_MAX);
    assert(total_recv_bytes <= (size_t)INT_MAX);

    send_bytes[peer] = (int)peer_send_bytes;
    recv_bytes[peer] = (int)peer_recv_bytes;
    send_displs[peer] = (int)total_send_bytes;
    recv_displs[peer] = (int)total_recv_bytes;

    total_send_bytes += peer_send_bytes;
    total_recv_bytes += peer_recv_bytes;
  }

  unsigned char *send_buffer =
      total_send_bytes > 0 ? (unsigned char *)malloc(total_send_bytes) : NULL;
  unsigned char *recv_buffer =
      total_recv_bytes > 0 ? (unsigned char *)malloc(total_recv_bytes) : NULL;
  assert(total_send_bytes == 0 || send_buffer);
  assert(total_recv_bytes == 0 || recv_buffer);

  for (int peer = 0; peer < nproc; peer++) {
    if (send_counts[peer] == 0)
      continue;

    unsigned char *dst = send_buffer + send_displs[peer];
    for (size_t idx = 0; idx < pg.snd_boundary[peer].particles.size; idx++) {
      Particle *particle = pg.snd_boundary[peer].particles.ptrs[idx];
      assert(particle);
      _particle_boundary_pack_particle_fields(dst + idx * record_size, particle,
                                              fields);
    }
  }

  MPI_Alltoallv(send_buffer, send_bytes, send_displs, MPI_BYTE, recv_buffer,
                recv_bytes, recv_displs, MPI_BYTE, MPI_COMM_WORLD);

  for (int peer = 0; peer < nproc; peer++) {
    if (recv_counts[peer] == 0)
      continue;

    const unsigned char *src = recv_buffer + recv_displs[peer];
    for (int idx = 0; idx < recv_counts[peer]; idx++) {
      const unsigned char *record = src + (size_t)idx * record_size;
      uint64_t gid = 0;
      memcpy(&gid, record, sizeof(gid));

      Particle *particle = _particle_grid_find_boundary_ghost(peer, gid);
      assert(particle);
      _particle_boundary_unpack_particle_fields(particle, record, fields);
    }
  }

  free(recv_buffer);
  free(send_buffer);
  free(recv_displs);
  free(send_displs);
  free(recv_bytes);
  free(send_bytes);
  free(recv_counts);
  free(send_counts);
#else
  NOT_UNUSED(fields);
#endif
}

/*!
 * @brief
 */
trace void particle_grid_update_pid() {
#if _MPI
  if (!pg.pool.pool)
    return;

  const int nproc = npe();
  const size_t record_size = _particle_grid_migration_record_size();
  assert(record_size <= (size_t)INT_MAX);

  int *particle_counts = (int *)calloc((size_t)nproc, sizeof(*particle_counts));
  int *send_counts = (int *)calloc((size_t)nproc, sizeof(*send_counts));
  int *recv_counts = (int *)calloc((size_t)nproc, sizeof(*recv_counts));
  int *send_bytes = (int *)calloc((size_t)nproc, sizeof(*send_bytes));
  int *recv_bytes = (int *)calloc((size_t)nproc, sizeof(*recv_bytes));
  int *send_displs = (int *)calloc((size_t)nproc, sizeof(*send_displs));
  int *recv_displs = (int *)calloc((size_t)nproc, sizeof(*recv_displs));
  assert(particle_counts);
  assert(send_counts);
  assert(recv_counts);
  assert(send_bytes);
  assert(recv_bytes);
  assert(send_displs);
  assert(recv_displs);

  // assert (pg.pool.active.size <= (size_t) INT_MAX);
  // const int local_particle_count = (int) pg.pool.active.size;
  // MPI_Allgather (&local_particle_count, 1, MPI_INT, particle_counts, 1,
  // MPI_INT, MPI_COMM_WORLD);

  for (int peer = 0; peer < nproc; peer++) {
    particle_exchange_list_clear(&pg.snd_migrate[peer]);
    particle_exchange_list_clear(&pg.rcv_migrate[peer]);
  }

  for (size_t idx = 0; idx < pg.pool.active.size; idx++) {
    Particle *particle = pg.pool.active.ptrs[idx];
    if (!particle)
      continue;

    const int pid_old = particle->pid;
    if (pid_old >= 0 && pid_old != pid())
      continue;

    const int pid_new = _particle_grid_get_pid(particle);
    particle->pid = pid_new;

    if (pid_new >= 0 && pid_new != pid()) {
      assert(pid_new < nproc);
      particle_exchange_list_push_unique(&pg.snd_migrate[pid_new], particle);
    }
  }

  for (int peer = 0; peer < nproc; peer++) {
    assert(pg.snd_migrate[peer].particles.size <= (size_t)INT_MAX);
    send_counts[peer] = (int)pg.snd_migrate[peer].particles.size;
  }

  MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT,
               MPI_COMM_WORLD);

  size_t total_send_bytes = 0;
  size_t total_recv_bytes = 0;
  for (int peer = 0; peer < nproc; peer++) {
    const size_t peer_send_bytes = (size_t)send_counts[peer] * record_size;
    const size_t peer_recv_bytes = (size_t)recv_counts[peer] * record_size;
    assert(peer_send_bytes <= (size_t)INT_MAX);
    assert(peer_recv_bytes <= (size_t)INT_MAX);
    assert(total_send_bytes <= (size_t)INT_MAX);
    assert(total_recv_bytes <= (size_t)INT_MAX);

    send_bytes[peer] = (int)peer_send_bytes;
    recv_bytes[peer] = (int)peer_recv_bytes;
    send_displs[peer] = (int)total_send_bytes;
    recv_displs[peer] = (int)total_recv_bytes;

    total_send_bytes += peer_send_bytes;
    total_recv_bytes += peer_recv_bytes;
  }

  unsigned char *send_buffer =
      total_send_bytes > 0 ? (unsigned char *)malloc(total_send_bytes) : NULL;
  unsigned char *recv_buffer =
      total_recv_bytes > 0 ? (unsigned char *)malloc(total_recv_bytes) : NULL;
  assert(total_send_bytes == 0 || send_buffer);
  assert(total_recv_bytes == 0 || recv_buffer);

  for (int peer = 0; peer < nproc; peer++) {
    if (send_counts[peer] == 0)
      continue;
    unsigned char *dst = send_buffer + send_displs[peer];
    for (size_t idx = 0; idx < pg.snd_migrate[peer].particles.size; idx++) {
      Particle *particle = pg.snd_migrate[peer].particles.ptrs[idx];
      _particle_grid_pack_particle(dst + idx * record_size, particle);
    }
  }

  MPI_Alltoallv(send_buffer, send_bytes, send_displs, MPI_BYTE, recv_buffer,
                recv_bytes, recv_displs, MPI_BYTE, MPI_COMM_WORLD);

  for (int peer = 0; peer < nproc; peer++) {
    if (recv_counts[peer] == 0)
      continue;
    const unsigned char *src = recv_buffer + recv_displs[peer];
    for (int idx = 0; idx < recv_counts[peer]; idx++) {
      Particle *particle = particle_mempool_alloc_particle(&pg.pool);
      _particle_grid_unpack_particle(particle, src + (size_t)idx * record_size);
      particle->pid = pid();
      particle_exchange_list_push(&pg.rcv_migrate[peer], particle);
    }
  }

  for (int peer = 0; peer < nproc; peer++) {
    for (size_t idx = 0; idx < pg.snd_migrate[peer].particles.size; idx++) {
      Particle *particle = pg.snd_migrate[peer].particles.ptrs[idx];
      particle_mempool_free_particle_ptr(&pg.pool, particle);
    }
  }

  free(recv_buffer);
  free(send_buffer);
  free(recv_displs);
  free(send_displs);
  free(recv_bytes);
  free(send_bytes);
  free(recv_counts);
  free(send_counts);
  free(particle_counts);

  particle_grid_update_cells();
  pg.dirty = false;
#else
  for (size_t idx = 0; idx < pg.pool.active.size; idx++) {
    Particle *particle = pg.pool.active.ptrs[idx];
    if (particle)
      particle->pid = 0;
  }
#endif
}

/*!
 * @brief
 */
trace inline int _particle_grid_get_pid(Particle *particle) {
  if (!pg.cells)
    return -1;
#if _MPI
#if TREE
  {
    PARTICLE_VARIABLES();
    Point point = domain_locate_nonlocal(px, py, pz);
    int ig = 0, jg = 0, kg = 0;
    NOT_UNUSED(ig);
    NOT_UNUSED(jg);
    NOT_UNUSED(kg);

    if (allocated(0))
      return cell.pid;
    else
      return -1;
  }
#else // !TREE
  PARTICLE_VARIABLES();
  Point point = {0};
  point.level = depth();
  SET_DIMENSIONS();

  int coords[dimension];
  const int local_nx = point.n.x;
  const int global_nx = local_nx * Dimensions.x;
  const double Lx = domain_extent_x();

  double ppx = px;
  if (Period.x)
    ppx = X0 + fmod(fmod(ppx - X0, Lx) + Lx, Lx);

  int gi = (int)floor((ppx - X0) / Lx * global_nx);
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

  double ppy = py;
  if (Period.y)
    ppy = Y0 + fmod(fmod(ppy - Y0, Ly) + Ly, Ly);

  int gj = (int)floor((ppy - Y0) / Ly * global_ny);
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

  double ppz = pz;
  if (Period.z)
    ppz = Z0 + fmod(fmod(ppz - Z0, Lz) + Lz, Lz);

  int gk = (int)floor((ppz - Z0) / Lz * global_nz);
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

#endif
#else // !_MPI
  return 0;
#endif
}

/*!
 * @brief
 */
trace inline int _particle_grid_get_cell_pid(ParticleCellPoint pcp) {
  if (!pg.cells)
    return -1;

#if _MPI
  if (pcp.level < 0 || pcp.level > depth())
    return -1;

#if TREE
  {
    Point point = {0};
    point.level = pcp.level;
#if dimension >= 1
    point.i = pcp.i;
#endif
#if dimension >= 2
    point.j = pcp.j;
#endif
#if dimension >= 3
    point.k = pcp.k;
#endif

    int ig = 0, jg = 0, kg = 0;
    NOT_UNUSED(ig);
    NOT_UNUSED(jg);
    NOT_UNUSED(kg);
    POINT_VARIABLES();

    return allocated(0) ? cell.pid : -1;
  }
#else
  {
    Point point = {0};
    point.level = depth();
    SET_DIMENSIONS();

    const int shift = depth() - pcp.level;

#if dimension >= 1
    point.i = ((pcp.i - GHOSTS) << shift) + GHOSTS - mpi_coords[0] * point.n.x;
#endif
#if dimension >= 2
    point.j = ((pcp.j - GHOSTS) << shift) + GHOSTS - mpi_coords[1] * point.n.y;
#endif
#if dimension >= 3
    point.k = ((pcp.k - GHOSTS) << shift) + GHOSTS - mpi_coords[2] * point.n.z;
#endif

    int coords[dimension];
    coords[0] = mpi_coords[0];

    if (point.i < GHOSTS)
      coords[0]--;
    else if (point.i >= point.n.x + GHOSTS)
      coords[0]++;

#if dimension > 1
    coords[1] = mpi_coords[1];
    if (point.j < GHOSTS)
      coords[1]--;
    else if (point.j >= point.n.y + GHOSTS)
      coords[1]++;
#endif

#if dimension > 2
    coords[2] = mpi_coords[2];
    if (point.k < GHOSTS)
      coords[2]--;
    else if (point.k >= point.n.z + GHOSTS)
      coords[2]++;
#endif

    if (!Period.x && (coords[0] < 0 || coords[0] >= Dimensions.x))
      return -1;
    if (Period.x) {
      if (coords[0] < 0)
        coords[0] += Dimensions.x;
      else if (coords[0] >= Dimensions.x)
        coords[0] -= Dimensions.x;
    }

#if dimension > 1
    if (!Period.y && (coords[1] < 0 || coords[1] >= Dimensions.y))
      return -1;
    if (Period.y) {
      if (coords[1] < 0)
        coords[1] += Dimensions.y;
      else if (coords[1] >= Dimensions.y)
        coords[1] -= Dimensions.y;
    }
#endif

#if dimension > 2
    if (!Period.z && (coords[2] < 0 || coords[2] >= Dimensions.z))
      return -1;
    if (Period.z) {
      if (coords[2] < 0)
        coords[2] += Dimensions.z;
      else if (coords[2] >= Dimensions.z)
        coords[2] -= Dimensions.z;
    }
#endif

    int owner = -1;
    MPI_Cart_rank(pg.cartcomm, coords, &owner);
    return owner;
  }
#endif
#else
  return 0;
#endif
}
