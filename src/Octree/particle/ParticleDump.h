#pragma once

/*!
 * @file ParticleDump.h
 *
 * Particle checkpoint sidecar for library/io/output-dump.h.
 *
 * The Basilisk dump() file stores Eulerian fields and grid state, but it does
 * not know about this library's particle mempool. This header registers a
 * checkpointer that writes a raw particle sidecar next to the dump file:
 *
 *   <checkpoint-directory>/dump.particles
 *
 * The file layout is the same in serial and MPI:
 *
 *   ParticleDumpHeader
 *   uint64_t rank_counts[header.nranks]
 *   rank 0 owned particle records
 *   rank 1 owned particle records
 *   ...
 *
 * MPI writes and reads the record section collectively with MPI-IO. Boundary
 * ghost particles are deliberately excluded; they are reconstructed after
 * restore when a boundary radius is configured. This is a checkpoint format for
 * the same executable, not a portable interchange format across particle
 * layouts.
 *
 * Typical use:
 *
 *   #include "particle/ParticleDump.h"
 *
 *   static ParticleDumpConfig pdump = {
 *     .boundary_radius = 1
 *   };
 *
 *   event defaults (i = 0) {
 *     init_psolver();
 *     particle_grid_init(PARTICLE_GRID_LEVEL);
 *     particle_dump_register(&pdump);
 *   }
 *
 *   event init (i = 0) {
 *     if (restore_handler("checkpoint", all))
 *       return 0;
 *     ...
 *   }
 */

#include <limits.h>
#include <stdint.h>
#include <string.h>

#include "io/output-dump.h"
#if _MPI
#include "particle/ParticleGridMPI.h"
#else
#include "particle/ParticleGrid.h"
#endif

// ============================================================================
// Macro definitions
// ============================================================================

#define PARTICLE_DUMP_FILENAME ".particles"
#define PARTICLE_DUMP_MAGIC 0x5041525444554d50ULL
#define PARTICLE_DUMP_VERSION 2
#define PARTICLE_DUMP_NO_BOUNDARY_REFRESH (-1)

// ============================================================================
// Type definitions
// ============================================================================

typedef int (*ParticleDumpBoundaryRadiusFunction)(void *ctx);

typedef struct {
  uint64_t magic;
  uint32_t version;
  uint32_t particle_size;
  uint64_t datasize;
  uint64_t record_size;
  uint64_t count;
  uint32_t nranks;
  uint32_t reserved;
} ParticleDumpHeader;

/*!
 * @brief Particle dump registration options.
 *
 * @param boundary_radius Fixed MPI boundary refresh radius. Set to
 * PARTICLE_DUMP_NO_BOUNDARY_REFRESH to skip boundary ghost reconstruction.
 * @param boundary_radius_fn Optional callback used instead of boundary_radius.
 * This is useful when the interaction radius depends on runtime state.
 * @param boundary_ctx User context passed to boundary_radius_fn.
 *
 * A config passed to particle_dump_register() must outlive the simulation,
 * since output-dump.h stores the pointer as checkpointer context.
 */
typedef struct {
  int boundary_radius;
  ParticleDumpBoundaryRadiusFunction boundary_radius_fn;
  void *boundary_ctx;
} ParticleDumpConfig;

// ============================================================================
// Globals
// ============================================================================

static ParticleDumpConfig particle_dump_default_config = {
    .boundary_radius = PARTICLE_DUMP_NO_BOUNDARY_REFRESH,
    .boundary_radius_fn = NULL,
    .boundary_ctx = NULL};

// ============================================================================
// Function declarations
// ============================================================================

static uint64_t particle_dump_local_count();
static uint64_t particle_dump_prefix_count(uint64_t *counts, int rank);
static int particle_dump_record_bytes(uint64_t count, size_t record_size,
                                      size_t *bytes);
static unsigned char *
particle_dump_pack_owned(uint64_t count, size_t record_size, size_t *bytes);
static int particle_dump_unpack_records(unsigned char *buffer, uint64_t count,
                                        size_t record_size);
static int particle_dump_validate_header(ParticleDumpHeader *header,
                                         const char *path);
static int particle_dump_boundary_radius(ParticleDumpConfig *config);
#if _MPI
static int particle_dump_mpi_any_error(int failed);
#endif
static int particle_dump_write(const char *path, void *ctx);
static int particle_dump_restore(const char *path, void *ctx);
static inline Checkpointer
particle_dump_checkpointer(ParticleDumpConfig *config);
static inline int particle_dump_register(ParticleDumpConfig *config);
static inline int particle_dump_register_default();

// ============================================================================
// Function definitions
// ============================================================================

/*!
 * @brief Count locally-owned particles in the active mempool list.
 */
trace static uint64_t particle_dump_local_count() {
  uint64_t count = 0;
  for (size_t idx = 0; idx < pg.pool.active.size; idx++) {
    Particle *particle = pg.pool.active.ptrs[idx];
    if (particle && particle_is_local(particle))
      count++;
  }
  return count;
}

/*!
 * @brief Return the number of records belonging to ranks before @p rank.
 */
trace static uint64_t particle_dump_prefix_count(uint64_t *counts, int rank) {
  uint64_t prefix = 0;
  for (int peer = 0; peer < rank; peer++)
    prefix += counts[peer];
  return prefix;
}

/*!
 * @brief Convert a record count to bytes with overflow checking.
 */
trace static int particle_dump_record_bytes(uint64_t count, size_t record_size,
                                            size_t *bytes) {
  assert(bytes);
  if (record_size == 0 || count > (uint64_t)SIZE_MAX / record_size)
    return -1;
  *bytes = (size_t)count * record_size;
  return 0;
}

/*!
 * @brief Pack locally-owned particle records into a contiguous byte buffer.
 */
trace static unsigned char *
particle_dump_pack_owned(uint64_t count, size_t record_size, size_t *bytes) {
  if (particle_dump_record_bytes(count, record_size, bytes) != 0)
    return NULL;

  unsigned char *buffer = *bytes ? (unsigned char *)malloc(*bytes) : NULL;
  if (*bytes && !buffer)
    return NULL;

  unsigned char *cursor = buffer;
  for (size_t idx = 0; idx < pg.pool.active.size; idx++) {
    Particle *particle = pg.pool.active.ptrs[idx];
    if (!particle || !particle_is_local(particle))
      continue;

    memcpy(cursor, particle, record_size);
    cursor += record_size;
  }

  return buffer;
}

/*!
 * @brief Allocate particles and unpack raw records into the mempool.
 */
trace static int particle_dump_unpack_records(unsigned char *buffer,
                                              uint64_t count,
                                              size_t record_size) {
  unsigned char *cursor = buffer;
  for (uint64_t idx = 0; idx < count; idx++) {
    Particle *particle = particle_mempool_alloc_particle(&pg.pool);
    if (!particle)
      return -1;
    memcpy(particle, cursor, record_size);
    cursor += record_size;
  }
  return 0;
}

/*!
 * @brief Validate that @p header matches this executable's particle layout.
 */
trace static int particle_dump_validate_header(ParticleDumpHeader *header,
                                               const char *path) {
  assert(header);
  const size_t record_size = sizeof(Particle) + pg.pool.datasize;
  if (header->magic != PARTICLE_DUMP_MAGIC ||
      header->version != PARTICLE_DUMP_VERSION ||
      header->particle_size != sizeof(Particle) ||
      header->datasize != pg.pool.datasize ||
      header->record_size != record_size || header->nranks < 1) {
    fprintf(ferr,
            "particle_dump_restore(): incompatible particle checkpoint %s\n",
            path);
    return -1;
  }
  return 0;
}

/*!
 * @brief Return the configured boundary refresh radius.
 */
trace static int particle_dump_boundary_radius(ParticleDumpConfig *config) {
  if (!config)
    config = &particle_dump_default_config;

  if (config->boundary_radius_fn)
    return config->boundary_radius_fn(config->boundary_ctx);

  return config->boundary_radius;
}

#if _MPI
/*!
 * @brief Return non-zero if any rank reports failure.
 */
trace static int particle_dump_mpi_any_error(int failed) {
  int local = failed ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return global;
}
#endif

/*!
 * @brief Write locally-owned particles to a checkpoint sidecar.
 */
trace static int particle_dump_write(const char *path, void *ctx) {
  NOT_UNUSED(ctx);

  const size_t record_size = sizeof(Particle) + pg.pool.datasize;
  const uint64_t local_count = particle_dump_local_count();
  size_t local_bytes = 0;
  unsigned char *buffer =
      particle_dump_pack_owned(local_count, record_size, &local_bytes);

#if _MPI
  if (particle_dump_mpi_any_error(local_bytes && !buffer)) {
    free(buffer);
    return -1;
  }

  uint64_t *counts = (uint64_t *)calloc((size_t)npe(), sizeof(*counts));
  if (particle_dump_mpi_any_error(!counts)) {
    free(buffer);
    free(counts);
    return -1;
  }

  MPI_Allgather(&local_count, 1, MPI_UINT64_T, counts, 1, MPI_UINT64_T,
                MPI_COMM_WORLD);

  uint64_t global_count = 0;
  for (int peer = 0; peer < npe(); peer++)
    global_count += counts[peer];

  ParticleDumpHeader header = {.magic = PARTICLE_DUMP_MAGIC,
                               .version = PARTICLE_DUMP_VERSION,
                               .particle_size = sizeof(Particle),
                               .datasize = pg.pool.datasize,
                               .record_size = record_size,
                               .count = global_count,
                               .nranks = (uint32_t)npe(),
                               .reserved = 0};

  const MPI_Offset counts_offset = (MPI_Offset)sizeof(header);
  const MPI_Offset records_offset =
      counts_offset + (MPI_Offset)npe() * (MPI_Offset)sizeof(uint64_t);
  const MPI_Offset record_offset =
      records_offset + (MPI_Offset)particle_dump_prefix_count(counts, pid()) *
                           (MPI_Offset)record_size;
  const MPI_Offset total_size =
      records_offset + (MPI_Offset)global_count * (MPI_Offset)record_size;

  MPI_File fh;
  int ok = MPI_File_open(MPI_COMM_WORLD, path,
                         MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  const bool opened = ok == MPI_SUCCESS;
  if (particle_dump_mpi_any_error(ok != MPI_SUCCESS)) {
    free(counts);
    free(buffer);
    return -1;
  }

  if (ok == MPI_SUCCESS)
    ok = MPI_File_set_size(fh, total_size);
  int failed = particle_dump_mpi_any_error(ok != MPI_SUCCESS);

  MPI_Status status;
  const int header_bytes = pid() == 0 ? (int)sizeof(header) : 0;
  unsigned char dummy = 0;
  if (!failed)
    ok = MPI_File_write_at_all(fh, 0,
                               pid() == 0 ? (void *)&header : (void *)&dummy,
                               header_bytes, MPI_BYTE, &status);
  failed = particle_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  assert(sizeof(uint64_t) <= (size_t)INT_MAX);
  if (!failed)
    ok = MPI_File_write_at_all(
        fh, counts_offset + (MPI_Offset)pid() * (MPI_Offset)sizeof(uint64_t),
        &local_count, (int)sizeof(uint64_t), MPI_BYTE, &status);
  failed = particle_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  assert(local_bytes <= (size_t)INT_MAX);
  if (!failed)
    ok = MPI_File_write_at_all(fh, record_offset, local_bytes ? buffer : &dummy,
                               (int)local_bytes, MPI_BYTE, &status);
  failed = particle_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  if (opened)
    ok = MPI_File_close(&fh);
  failed = particle_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  free(counts);
  free(buffer);
  return failed ? -1 : 0;
#else
  if (local_bytes && !buffer)
    return -1;

  FILE *fp = fopen(path, "wb");
  if (!fp) {
    free(buffer);
    return -1;
  }

  ParticleDumpHeader header = {.magic = PARTICLE_DUMP_MAGIC,
                               .version = PARTICLE_DUMP_VERSION,
                               .particle_size = sizeof(Particle),
                               .datasize = pg.pool.datasize,
                               .record_size = record_size,
                               .count = local_count,
                               .nranks = 1,
                               .reserved = 0};

  int ok = fwrite(&header, sizeof(header), 1, fp) == 1 ? 0 : -1;
  if (ok == 0)
    ok = fwrite(&local_count, sizeof(local_count), 1, fp) == 1 ? 0 : -1;
  if (ok == 0 && local_bytes)
    ok = fwrite(buffer, local_bytes, 1, fp) == 1 ? 0 : -1;

  fclose(fp);
  free(buffer);
  return ok;
#endif
}

/*!
 * @brief Restore particles from a checkpoint sidecar.
 */
trace static int particle_dump_restore(const char *path, void *ctx) {
  ParticleDumpConfig *config = (ParticleDumpConfig *)ctx;

  if (!file_exists(path))
    return 0;

#if _MPI
  MPI_File fh;
  int ok =
      MPI_File_open(MPI_COMM_WORLD, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  const bool opened = ok == MPI_SUCCESS;
  if (particle_dump_mpi_any_error(ok != MPI_SUCCESS))
    return -1;

  ParticleDumpHeader header = {0};
  MPI_Status status;
  if (pid() == 0)
    ok = MPI_File_read_at(fh, 0, &header, (int)sizeof(header), MPI_BYTE,
                          &status);

  MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&header, (int)sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);

  if (ok == MPI_SUCCESS && particle_dump_validate_header(&header, path) != 0)
    ok = MPI_ERR_OTHER;
  if (ok == MPI_SUCCESS && header.nranks != (uint32_t)npe()) {
    fprintf(ferr,
            "particle_dump_restore(): rank count mismatch in %s: %u != %d\n",
            path, header.nranks, npe());
    ok = MPI_ERR_OTHER;
  }
  int failed = particle_dump_mpi_any_error(ok != MPI_SUCCESS);

  uint64_t *counts = NULL;
  if (!failed) {
    counts = (uint64_t *)calloc((size_t)header.nranks, sizeof(*counts));
    if (!counts)
      ok = MPI_ERR_OTHER;
  }
  failed = particle_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  const MPI_Offset counts_offset = (MPI_Offset)sizeof(header);
  const MPI_Offset records_offset =
      counts_offset + (MPI_Offset)header.nranks * (MPI_Offset)sizeof(uint64_t);

  if (!failed) {
    assert((size_t)header.nranks * sizeof(uint64_t) <= (size_t)INT_MAX);
    ok = MPI_File_read_at_all(fh, counts_offset, counts,
                              (int)((size_t)header.nranks * sizeof(uint64_t)),
                              MPI_BYTE, &status);
  }
  failed = particle_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  size_t local_bytes = 0;
  unsigned char *buffer = NULL;
  if (!failed) {
    ok = particle_dump_record_bytes(counts[pid()], (size_t)header.record_size,
                                    &local_bytes) == 0
             ? MPI_SUCCESS
             : MPI_ERR_OTHER;
  }
  if (!failed && ok == MPI_SUCCESS) {
    buffer = local_bytes ? (unsigned char *)malloc(local_bytes) : NULL;
    if (local_bytes && !buffer)
      ok = MPI_ERR_OTHER;
  }
  failed = particle_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  unsigned char dummy = 0;
  if (!failed) {
    assert(local_bytes <= (size_t)INT_MAX);
    const MPI_Offset record_offset =
        records_offset + (MPI_Offset)particle_dump_prefix_count(counts, pid()) *
                             (MPI_Offset)header.record_size;
    ok = MPI_File_read_at_all(fh, record_offset, local_bytes ? buffer : &dummy,
                              (int)local_bytes, MPI_BYTE, &status);
  }
  failed = particle_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  if (opened)
    ok = MPI_File_close(&fh);
  failed = particle_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  if (!failed) {
    particle_grid_delete_all_particles();
    if (particle_dump_unpack_records(buffer, counts[pid()],
                                     (size_t)header.record_size) != 0)
      ok = MPI_ERR_OTHER;
  }
  failed = particle_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  free(counts);
  free(buffer);
  if (failed)
    return -1;
#else
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return -1;

  ParticleDumpHeader header = {0};
  int ok = fread(&header, sizeof(header), 1, fp) == 1 ? 0 : -1;
  if (ok == 0)
    ok = particle_dump_validate_header(&header, path);
  if (ok == 0 && header.nranks != 1) {
    fprintf(ferr,
            "particle_dump_restore(): rank count mismatch in %s: %u != 1\n",
            path, header.nranks);
    ok = -1;
  }

  uint64_t local_count = 0;
  if (ok == 0)
    ok = fread(&local_count, sizeof(local_count), 1, fp) == 1 ? 0 : -1;

  size_t local_bytes = 0;
  if (ok == 0)
    ok = particle_dump_record_bytes(local_count, (size_t)header.record_size,
                                    &local_bytes);

  unsigned char *buffer =
      local_bytes ? (unsigned char *)malloc(local_bytes) : NULL;
  if (ok == 0 && local_bytes && !buffer)
    ok = -1;

  if (ok == 0 && local_bytes)
    ok = fread(buffer, local_bytes, 1, fp) == 1 ? 0 : -1;

  fclose(fp);

  if (ok == 0) {
    particle_grid_delete_all_particles();
    ok = particle_dump_unpack_records(buffer, local_count,
                                      (size_t)header.record_size);
  }

  free(buffer);
  if (ok != 0)
    return -1;
#endif

  particle_grid_update_cells();
#if _MPI
  particle_grid_update_pid();
  const int boundary_radius = particle_dump_boundary_radius(config);
  if (boundary_radius >= 0)
    particle_grid_update_boundary_broadcast_cells(boundary_radius);
#endif

  return 0;
}

/*!
 * @brief Build the Checkpointer used by output-dump.h.
 */
static inline Checkpointer
particle_dump_checkpointer(ParticleDumpConfig *config) {
  if (!config)
    config = &particle_dump_default_config;

  return (Checkpointer){.filename = PARTICLE_DUMP_FILENAME,
                        .dump_phase = CKPT_PHASE_POST_DUMP,
                        .dump = particle_dump_write,
                        .restore_phase = CKPT_PHASE_POST_RESTORE,
                        .restore = particle_dump_restore,
                        .ctx = config};
}

/*!
 * @brief Register particle checkpointing with output-dump.h.
 */
static inline int particle_dump_register(ParticleDumpConfig *config) {
  return checkpointer_register(particle_dump_checkpointer(config));
}

/*!
 * @brief Register particle checkpointing without MPI boundary ghost refresh.
 */
static inline int particle_dump_register_default() {
  return particle_dump_register(&particle_dump_default_config);
}
