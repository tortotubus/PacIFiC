#pragma once

/**
# Capsule checkpoint sidecar for `output-dump.h`

This sidecar assumes full capsule and IBM initialization has already run before
`restore_handler()`. Capsule restore writes only mutable capsule solver state
into the init-built capsule objects. Sparse capsule-backed IBM nodes are
ephemeral local-support caches and are recreated from capsule state after the
Basilisk grid restore.
*/

#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "io/output-dump.h"

// ============================================================================
// Macro definitions
// ============================================================================

#define CAPSULE_DUMP_FILENAME ".caps"
#define CAPSULE_DUMP_MAGIC 0x43415053444d5031ULL
#define CAPSULE_DUMP_VERSION 3

// ============================================================================
// Type definitions
// ============================================================================

typedef struct {
  uint64_t magic;
  uint32_t version;
  uint32_t nranks;
  uint64_t total_capsules;
  uint64_t total_owned_capsules;
} CapsuleDumpHeader;

typedef struct {
  int cap_id;
  int cap_type;
  int owner_pid;
  int mesh_id;
  int nln;
  int nle;
  int nlt;
  coord centroid;
  coord ang_vel;
  double initial_volume;
  double volume;
  double circum_radius;
  double taylor_deform;
  int updated_stretches;
  int updated_normals;
  int updated_curvatures;
  int isactive;
} CapsuleDumpRecordHeader;

typedef struct {
  coord pos;
  coord lagVel;
  coord lagForce;
  coord normal;
  double curv;
  double gcurv;
  double ref_curv;
  int nb_fit_iterations;
} CapsuleDumpNodeState;

typedef struct {
  unsigned char *data;
  size_t size;
  size_t capacity;
  size_t offset;
} CapsuleDumpBuffer;

// ============================================================================
// Function declarations
// ============================================================================

static uint64_t capsule_dump_local_count(void);
static uint64_t capsule_dump_prefix_count(uint64_t *counts, int rank);
static int capsule_dump_record_bytes(uint64_t count, size_t record_size,
                                     size_t *bytes);
static int capsule_dump_buffer_reserve(CapsuleDumpBuffer *buffer, size_t extra);
static int capsule_dump_buffer_write(CapsuleDumpBuffer *buffer,
                                     const void *data, size_t len);
static int capsule_dump_buffer_read(CapsuleDumpBuffer *buffer, void *data,
                                    size_t len);
static void capsule_dump_buffer_reset(CapsuleDumpBuffer *buffer);
static size_t capsule_dump_mesh_record_size(CapsuleMesh *mesh);
static void capsule_dump_node_state_from_mesh(CapsuleDumpNodeState *state,
                                              CapNode *node);
static void capsule_dump_node_state_to_mesh(CapNode *node,
                                            CapsuleDumpNodeState *state);
static int capsule_dump_pack_mesh(CapsuleDumpBuffer *buffer, CapsuleMesh *mesh);
static int capsule_dump_unpack_mesh(CapsuleDumpBuffer *buffer);
static int capsule_dump_pack_local(CapsuleDumpBuffer *buffer,
                                   uint64_t *local_count);
static int capsule_dump_validate_header(CapsuleDumpHeader *header,
                                        const char *path);
static int capsule_dump_write_stream(FILE *fp);
#if _MPI
static int capsule_dump_mpi_any_error(int failed);
#endif
static int capsule_dump_write(const char *path, void *ctx);
static int capsule_dump_restore(const char *path, void *ctx);
static inline Checkpointer capsule_dump_checkpointer(void);
static inline int capsule_dump_register(void);
static inline int capsule_dump_register_default(void);

void dump_capsules(const char *fname, FILE *fp);
void restore_capsules(const char *filename);

// ============================================================================
// Events
// ============================================================================

event defaults(i = 0) { capsule_dump_register_default(); }

// ============================================================================
// Function definitions
// ============================================================================

static uint64_t capsule_dump_local_count(void) {
  uint64_t count = 0;
  for (int i = 0; i < allCaps.nbcaps; i++)
    if (CAPS(i).isactive && capsule_mesh_is_local_owner(&CAPS(i)))
      count++;
  return count;
}

static uint64_t capsule_dump_prefix_count(uint64_t *counts, int rank) {
  uint64_t prefix = 0;
  for (int peer = 0; peer < rank; peer++)
    prefix += counts[peer];
  return prefix;
}

static int capsule_dump_record_bytes(uint64_t count, size_t record_size,
                                     size_t *bytes) {
  assert(bytes);
  if (record_size == 0 || count > (uint64_t)SIZE_MAX / record_size)
    return -1;
  *bytes = (size_t)count * record_size;
  return 0;
}

static int capsule_dump_buffer_reserve(CapsuleDumpBuffer *buffer,
                                       size_t extra) {
  assert(buffer);
  if (extra == 0)
    return 0;
  if (buffer->size > SIZE_MAX - extra)
    return -1;

  size_t required = buffer->size + extra;
  if (required <= buffer->capacity)
    return 0;

  size_t next = buffer->capacity ? buffer->capacity : 256;
  while (next < required) {
    if (next > SIZE_MAX / 2)
      next = required;
    else
      next *= 2;
    if (next < required)
      return -1;
  }

  unsigned char *data =
      (unsigned char *)realloc(buffer->data, next * sizeof(unsigned char));
  if (!data)
    return -1;
  buffer->data = data;
  buffer->capacity = next;
  return 0;
}

static int capsule_dump_buffer_write(CapsuleDumpBuffer *buffer,
                                     const void *data, size_t len) {
  if (capsule_dump_buffer_reserve(buffer, len) != 0)
    return -1;
  if (len)
    memcpy(buffer->data + buffer->size, data, len);
  buffer->size += len;
  return 0;
}

static int capsule_dump_buffer_read(CapsuleDumpBuffer *buffer, void *data,
                                    size_t len) {
  assert(buffer);
  if (buffer->offset > buffer->size || len > buffer->size - buffer->offset)
    return -1;
  if (len)
    memcpy(data, buffer->data + buffer->offset, len);
  buffer->offset += len;
  return 0;
}

static void capsule_dump_buffer_reset(CapsuleDumpBuffer *buffer) {
  if (!buffer)
    return;
  free(buffer->data);
  buffer->data = NULL;
  buffer->size = 0;
  buffer->capacity = 0;
  buffer->offset = 0;
}

static size_t capsule_dump_mesh_record_size(CapsuleMesh *mesh) {
  if (!mesh || mesh->nln < 0 || mesh->nle < 0 || mesh->nlt < 0)
    return 0;

  size_t size = sizeof(CapsuleDumpRecordHeader);
  size_t bytes = 0;
  if (capsule_dump_record_bytes((uint64_t)mesh->nln,
                                sizeof(CapsuleDumpNodeState), &bytes) != 0)
    return 0;
  size += bytes;
  if (capsule_dump_record_bytes((uint64_t)mesh->nle, sizeof(EdgeState),
                                &bytes) != 0 ||
      size > SIZE_MAX - bytes)
    return 0;
  size += bytes;
  if (capsule_dump_record_bytes((uint64_t)mesh->nlt, sizeof(TriangleState),
                                &bytes) != 0 ||
      size > SIZE_MAX - bytes)
    return 0;
  size += bytes;
  return size;
}

static void capsule_dump_node_state_from_mesh(CapsuleDumpNodeState *state,
                                              CapNode *node) {
  assert(state);
  assert(node);
  state->pos = node->pos;
  state->lagVel = node->lagVel;
  state->lagForce = node->lagForce;
  state->normal = node->normal;
  state->curv = node->curv;
  state->gcurv = node->gcurv;
  state->ref_curv = node->ref_curv;
  state->nb_fit_iterations = node->nb_fit_iterations;
}

static void capsule_dump_node_state_to_mesh(CapNode *node,
                                            CapsuleDumpNodeState *state) {
  assert(node);
  assert(state);
  node->pos = state->pos;
  node->lagVel = state->lagVel;
  node->lagForce = state->lagForce;
  node->normal = state->normal;
  node->curv = state->curv;
  node->gcurv = state->gcurv;
  node->ref_curv = state->ref_curv;
  node->nb_fit_iterations = state->nb_fit_iterations;
}

static int capsule_dump_pack_mesh(CapsuleDumpBuffer *buffer,
                                  CapsuleMesh *mesh) {
  assert(buffer);
  assert(mesh);

  int mesh_id = capsule_ibm_mesh_id(mesh);
  if (mesh_id < 0 || mesh_id >= ibmm.nm)
    return -1;

  IBMesh *ibmesh = &ibmm.meshes[mesh_id];
  if (ibmesh->pid != mesh->owner_pid)
    return -1;

  CapsuleDumpRecordHeader header = {
      .cap_id = mesh->cap_id,
      .cap_type = mesh->cap_type,
      .owner_pid = mesh->owner_pid,
      .mesh_id = mesh_id,
      .nln = mesh->nln,
      .nle = mesh->nle,
      .nlt = mesh->nlt,
      .centroid = mesh->centroid,
      .ang_vel = mesh->ang_vel,
      .initial_volume = mesh->initial_volume,
      .volume = mesh->volume,
      .circum_radius = mesh->circum_radius,
      .taylor_deform = mesh->taylor_deform,
      .updated_stretches = mesh->updated_stretches ? 1 : 0,
      .updated_normals = mesh->updated_normals ? 1 : 0,
      .updated_curvatures = mesh->updated_curvatures ? 1 : 0,
      .isactive = mesh->isactive ? 1 : 0};

  if (capsule_dump_buffer_write(buffer, &header, sizeof(header)) != 0)
    return -1;

  for (int i = 0; i < mesh->nln; i++) {
    CapsuleDumpNodeState node_state = {0};
    capsule_dump_node_state_from_mesh(&node_state, &mesh->nodes[i]);
    if (capsule_dump_buffer_write(buffer, &node_state, sizeof(node_state)) != 0)
      return -1;
  }

  for (int i = 0; i < mesh->nle; i++) {
    EdgeState state = *LAG_EDGE_STATE(mesh, i);
    if (capsule_dump_buffer_write(buffer, &state, sizeof(state)) != 0)
      return -1;
  }

  for (int i = 0; i < mesh->nlt; i++) {
    TriangleState state = *LAG_TRI_STATE(mesh, i);
    if (capsule_dump_buffer_write(buffer, &state, sizeof(state)) != 0)
      return -1;
  }

  return 0;
}

static int capsule_dump_unpack_mesh(CapsuleDumpBuffer *buffer) {
  CapsuleDumpRecordHeader record = {0};
  if (capsule_dump_buffer_read(buffer, &record, sizeof(record)) != 0)
    return -1;

  if (record.mesh_id < 0 || record.mesh_id >= ibmm.nm) {
    fprintf(ferr,
            "capsule_dump_restore(): invalid mesh id %d in capsule sidecar\n",
            record.mesh_id);
    return -1;
  }

  IBMesh *ibmesh = &ibmm.meshes[record.mesh_id];
  if (ibmesh->pid != record.owner_pid) {
    fprintf(
        ferr,
        "capsule_dump_restore(): mesh owner mismatch for mesh %d: %d != %d\n",
        record.mesh_id, ibmesh->pid, record.owner_pid);
    return -1;
  }

#if _MPI
  if (record.owner_pid != pid()) {
    fprintf(ferr,
            "capsule_dump_restore(): rank %d received capsule owned by %d\n",
            pid(), record.owner_pid);
    return -1;
  }
#endif

  CapsuleMesh *mesh = capsule_ibm_lag_from_mesh_id(record.mesh_id);
  if (!mesh) {
    fprintf(ferr, "capsule_dump_restore(): missing capsule for mesh id %d\n",
            record.mesh_id);
    return -1;
  }

  if (mesh->cap_id != record.cap_id || mesh->cap_type != record.cap_type ||
      mesh->owner_pid != record.owner_pid || mesh->nln != record.nln ||
      mesh->nle != record.nle || mesh->nlt != record.nlt) {
    fprintf(
        ferr,
        "capsule_dump_restore(): capsule metadata mismatch for mesh id %d\n",
        record.mesh_id);
    return -1;
  }

  if (!mesh->nodes || (mesh->nle > 0 && !(mesh->edge_states || mesh->edges)) ||
      (mesh->nlt > 0 && !(mesh->triangle_states || mesh->triangles))) {
    fprintf(ferr,
            "capsule_dump_restore(): capsule storage missing for mesh id %d\n",
            record.mesh_id);
    return -1;
  }

  for (int i = 0; i < mesh->nln; i++) {
    CapsuleDumpNodeState node_state = {0};
    if (capsule_dump_buffer_read(buffer, &node_state, sizeof(node_state)) != 0)
      return -1;
    capsule_dump_node_state_to_mesh(&mesh->nodes[i], &node_state);
  }

  for (int i = 0; i < mesh->nle; i++) {
    EdgeState state = {0};
    if (capsule_dump_buffer_read(buffer, &state, sizeof(state)) != 0)
      return -1;
    *LAG_EDGE_STATE(mesh, i) = state;
  }

  for (int i = 0; i < mesh->nlt; i++) {
    TriangleState state = {0};
    if (capsule_dump_buffer_read(buffer, &state, sizeof(state)) != 0)
      return -1;
    *LAG_TRI_STATE(mesh, i) = state;
  }

  mesh->centroid = record.centroid;
  mesh->ang_vel = record.ang_vel;
  mesh->initial_volume = record.initial_volume;
  mesh->volume = record.volume;
  mesh->circum_radius = record.circum_radius;
  mesh->taylor_deform = record.taylor_deform;
  mesh->updated_stretches = record.updated_stretches != 0;
  mesh->updated_normals = record.updated_normals != 0;
  mesh->updated_curvatures = record.updated_curvatures != 0;
  mesh->isactive = record.isactive != 0;

  return 0;
}

static int capsule_dump_pack_local(CapsuleDumpBuffer *buffer,
                                   uint64_t *local_count) {
  assert(buffer);
  assert(local_count);

  *local_count = 0;
  for (int i = 0; i < allCaps.nbcaps; i++) {
    CapsuleMesh *mesh = &CAPS(i);
    if (!mesh->isactive || !capsule_mesh_is_local_owner(mesh))
      continue;

    if (capsule_dump_mesh_record_size(mesh) == 0)
      return -1;
    if (capsule_dump_pack_mesh(buffer, mesh) != 0)
      return -1;
    (*local_count)++;
  }

  return 0;
}

static int capsule_dump_validate_header(CapsuleDumpHeader *header,
                                        const char *path) {
  assert(header);
  if (header->magic != CAPSULE_DUMP_MAGIC ||
      header->version != CAPSULE_DUMP_VERSION || header->nranks < 1) {
    fprintf(ferr,
            "capsule_dump_restore(): incompatible capsule checkpoint %s\n",
            path);
    return -1;
  }
  return 0;
}

static int capsule_dump_write_stream(FILE *fp) {
  CapsuleDumpBuffer payload = {0};
  uint64_t local_count = 0;
  if (capsule_dump_pack_local(&payload, &local_count) != 0) {
    capsule_dump_buffer_reset(&payload);
    return -1;
  }

  CapsuleDumpHeader header = {.magic = CAPSULE_DUMP_MAGIC,
                              .version = CAPSULE_DUMP_VERSION,
                              .nranks = 1,
                              .total_capsules = (uint64_t)allCaps.nbcaps,
                              .total_owned_capsules = local_count};

  uint64_t rank_count = local_count;
  uint64_t rank_bytes = (uint64_t)payload.size;

  int ok = fwrite(&header, sizeof(header), 1, fp) == 1 ? 0 : -1;
  if (ok == 0)
    ok = fwrite(&rank_count, sizeof(rank_count), 1, fp) == 1 ? 0 : -1;
  if (ok == 0)
    ok = fwrite(&rank_bytes, sizeof(rank_bytes), 1, fp) == 1 ? 0 : -1;
  if (ok == 0 && payload.size)
    ok = fwrite(payload.data, payload.size, 1, fp) == 1 ? 0 : -1;

  capsule_dump_buffer_reset(&payload);
  return ok;
}

#if _MPI
static int capsule_dump_mpi_any_error(int failed) {
  int local = failed ? 1 : 0;
  int global = 0;
  MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return global;
}
#endif

static int capsule_dump_write(const char *path, void *ctx) {
  NOT_UNUSED(ctx);

  CapsuleDumpBuffer payload = {0};
  uint64_t local_count = 0;
  if (capsule_dump_pack_local(&payload, &local_count) != 0) {
    capsule_dump_buffer_reset(&payload);
    return -1;
  }

#if _MPI
  uint64_t *counts = (uint64_t *)calloc((size_t)npe(), sizeof(*counts));
  uint64_t *bytes = (uint64_t *)calloc((size_t)npe(), sizeof(*bytes));
  if (capsule_dump_mpi_any_error(!counts || !bytes)) {
    free(counts);
    free(bytes);
    capsule_dump_buffer_reset(&payload);
    return -1;
  }

  uint64_t local_bytes = (uint64_t)payload.size;
  MPI_Allgather(&local_count, 1, MPI_UINT64_T, counts, 1, MPI_UINT64_T,
                MPI_COMM_WORLD);
  MPI_Allgather(&local_bytes, 1, MPI_UINT64_T, bytes, 1, MPI_UINT64_T,
                MPI_COMM_WORLD);

  uint64_t total_owned = 0;
  uint64_t total_bytes = 0;
  for (int peer = 0; peer < npe(); peer++) {
    total_owned += counts[peer];
    total_bytes += bytes[peer];
  }

  CapsuleDumpHeader header = {.magic = CAPSULE_DUMP_MAGIC,
                              .version = CAPSULE_DUMP_VERSION,
                              .nranks = (uint32_t)npe(),
                              .total_capsules = (uint64_t)allCaps.nbcaps,
                              .total_owned_capsules = total_owned};

  const MPI_Offset counts_offset = (MPI_Offset)sizeof(header);
  const MPI_Offset bytes_offset =
      counts_offset + (MPI_Offset)npe() * (MPI_Offset)sizeof(uint64_t);
  const MPI_Offset payload_offset =
      bytes_offset + (MPI_Offset)npe() * (MPI_Offset)sizeof(uint64_t);
  const MPI_Offset local_payload_offset =
      payload_offset + (MPI_Offset)capsule_dump_prefix_count(bytes, pid());
  const MPI_Offset total_size = payload_offset + (MPI_Offset)total_bytes;

  MPI_File fh;
  int ok = MPI_File_open(MPI_COMM_WORLD, path,
                         MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  const bool opened = ok == MPI_SUCCESS;
  if (capsule_dump_mpi_any_error(ok != MPI_SUCCESS)) {
    free(counts);
    free(bytes);
    capsule_dump_buffer_reset(&payload);
    return -1;
  }

  if (ok == MPI_SUCCESS)
    ok = MPI_File_set_size(fh, total_size);
  int failed = capsule_dump_mpi_any_error(ok != MPI_SUCCESS);

  MPI_Status status;
  unsigned char dummy = 0;
  const int header_bytes = pid() == 0 ? (int)sizeof(header) : 0;
  if (!failed)
    ok = MPI_File_write_at_all(fh, 0,
                               pid() == 0 ? (void *)&header : (void *)&dummy,
                               header_bytes, MPI_BYTE, &status);
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  if (!failed)
    ok = MPI_File_write_at_all(
        fh, counts_offset + (MPI_Offset)pid() * (MPI_Offset)sizeof(uint64_t),
        &local_count, (int)sizeof(uint64_t), MPI_BYTE, &status);
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  if (!failed)
    ok = MPI_File_write_at_all(
        fh, bytes_offset + (MPI_Offset)pid() * (MPI_Offset)sizeof(uint64_t),
        &local_bytes, (int)sizeof(uint64_t), MPI_BYTE, &status);
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  assert(payload.size <= (size_t)INT_MAX);
  if (!failed)
    ok = MPI_File_write_at_all(fh, local_payload_offset,
                               payload.size ? payload.data : &dummy,
                               (int)payload.size, MPI_BYTE, &status);
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  if (opened)
    ok = MPI_File_close(&fh);
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  free(counts);
  free(bytes);
  capsule_dump_buffer_reset(&payload);
  if (!failed)
    ibmeshmanager_sync_velocity_coupled_model_outputs();
  return failed ? -1 : 0;
#else
  FILE *fp = fopen(path, "wb");
  if (!fp) {
    capsule_dump_buffer_reset(&payload);
    return -1;
  }

  CapsuleDumpHeader header = {.magic = CAPSULE_DUMP_MAGIC,
                              .version = CAPSULE_DUMP_VERSION,
                              .nranks = 1,
                              .total_capsules = (uint64_t)allCaps.nbcaps,
                              .total_owned_capsules = local_count};
  uint64_t rank_count = local_count;
  uint64_t rank_bytes = (uint64_t)payload.size;

  int ok = fwrite(&header, sizeof(header), 1, fp) == 1 ? 0 : -1;
  if (ok == 0)
    ok = fwrite(&rank_count, sizeof(rank_count), 1, fp) == 1 ? 0 : -1;
  if (ok == 0)
    ok = fwrite(&rank_bytes, sizeof(rank_bytes), 1, fp) == 1 ? 0 : -1;
  if (ok == 0 && payload.size)
    ok = fwrite(payload.data, payload.size, 1, fp) == 1 ? 0 : -1;

  fclose(fp);
  capsule_dump_buffer_reset(&payload);
  if (ok == 0)
    ibmeshmanager_sync_velocity_coupled_model_outputs();
  return ok;
#endif
}

static int capsule_dump_restore(const char *path, void *ctx) {
  NOT_UNUSED(ctx);

  if (!file_exists(path))
    return 0;

#if _MPI
  MPI_File fh;
  int ok =
      MPI_File_open(MPI_COMM_WORLD, path, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  const bool opened = ok == MPI_SUCCESS;
  if (capsule_dump_mpi_any_error(ok != MPI_SUCCESS))
    return -1;

  CapsuleDumpHeader header = {0};
  MPI_Status status;
  if (pid() == 0)
    ok = MPI_File_read_at(fh, 0, &header, (int)sizeof(header), MPI_BYTE,
                          &status);

  MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&header, (int)sizeof(header), MPI_BYTE, 0, MPI_COMM_WORLD);

  if (ok == MPI_SUCCESS && capsule_dump_validate_header(&header, path) != 0)
    ok = MPI_ERR_OTHER;
  if (ok == MPI_SUCCESS && header.nranks != (uint32_t)npe()) {
    fprintf(ferr,
            "capsule_dump_restore(): rank count mismatch in %s: %u != %d\n",
            path, header.nranks, npe());
    ok = MPI_ERR_OTHER;
  }
  if (ok == MPI_SUCCESS && header.total_capsules != (uint64_t)allCaps.nbcaps) {
    fprintf(ferr, "capsule_dump_restore(): capsule count mismatch");
    ok = MPI_ERR_OTHER;
  }
  int failed = capsule_dump_mpi_any_error(ok != MPI_SUCCESS);

  uint64_t *counts = NULL;
  uint64_t *bytes = NULL;
  if (!failed) {
    counts = (uint64_t *)calloc((size_t)header.nranks, sizeof(*counts));
    bytes = (uint64_t *)calloc((size_t)header.nranks, sizeof(*bytes));
    if (!counts || !bytes)
      ok = MPI_ERR_OTHER;
  }
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  const MPI_Offset counts_offset = (MPI_Offset)sizeof(header);
  const MPI_Offset bytes_offset =
      counts_offset + (MPI_Offset)header.nranks * (MPI_Offset)sizeof(uint64_t);
  const MPI_Offset payload_offset =
      bytes_offset + (MPI_Offset)header.nranks * (MPI_Offset)sizeof(uint64_t);

  if (!failed) {
    assert((size_t)header.nranks * sizeof(uint64_t) <= (size_t)INT_MAX);
    ok = MPI_File_read_at_all(fh, counts_offset, counts,
                              (int)((size_t)header.nranks * sizeof(uint64_t)),
                              MPI_BYTE, &status);
  }
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  if (!failed) {
    assert((size_t)header.nranks * sizeof(uint64_t) <= (size_t)INT_MAX);
    ok = MPI_File_read_at_all(fh, bytes_offset, bytes,
                              (int)((size_t)header.nranks * sizeof(uint64_t)),
                              MPI_BYTE, &status);
  }
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  CapsuleDumpBuffer payload = {0};
  if (!failed && bytes[pid()] > SIZE_MAX)
    ok = MPI_ERR_OTHER;
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  if (!failed && bytes[pid()] > 0) {
    payload.data = (unsigned char *)malloc((size_t)bytes[pid()]);
    if (!payload.data)
      ok = MPI_ERR_OTHER;
    else {
      payload.size = (size_t)bytes[pid()];
      payload.capacity = payload.size;
    }
  }
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  unsigned char dummy = 0;
  if (!failed) {
    assert(payload.size <= (size_t)INT_MAX);
    const MPI_Offset local_payload_offset =
        payload_offset + (MPI_Offset)capsule_dump_prefix_count(bytes, pid());
    ok = MPI_File_read_at_all(fh, local_payload_offset,
                              payload.size ? payload.data : &dummy,
                              (int)payload.size, MPI_BYTE, &status);
  }
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  if (opened)
    ok = MPI_File_close(&fh);
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  if (!failed) {
    uint64_t counted_caps = 0;
    for (int peer = 0; peer < npe(); peer++)
      counted_caps += counts[peer];
    if (counted_caps != header.total_owned_capsules) {
      fprintf(ferr,
              "capsule_dump_restore(): owned capsule count mismatch in %s\n",
              path);
      ok = MPI_ERR_OTHER;
    }
  }
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  if (!failed) {
    payload.offset = 0;
    for (uint64_t i = 0; i < counts[pid()]; i++) {
      if (capsule_dump_unpack_mesh(&payload) != 0) {
        ok = MPI_ERR_OTHER;
        break;
      }
    }
    if (ok == MPI_SUCCESS && payload.offset != payload.size)
      ok = MPI_ERR_OTHER;
  }
  failed = capsule_dump_mpi_any_error(failed || ok != MPI_SUCCESS);

  free(counts);
  free(bytes);
  capsule_dump_buffer_reset(&payload);
  if (!failed)
    ibmeshmanager_sync_velocity_coupled_model_outputs();
  return failed ? -1 : 0;
#else
  FILE *fp = fopen(path, "rb");
  if (!fp)
    return -1;

  CapsuleDumpHeader header = {0};
  int ok = fread(&header, sizeof(header), 1, fp) == 1 ? 0 : -1;
  if (ok == 0)
    ok = capsule_dump_validate_header(&header, path);
  if (ok == 0 && header.nranks != 1) {
    fprintf(ferr,
            "capsule_dump_restore(): rank count mismatch in %s: %u != 1\n",
            path, header.nranks);
    ok = -1;
  }
  if (ok == 0 && header.total_capsules != (uint64_t)allCaps.nbcaps) {
    fprintf(ferr, "capsule_dump_restore(): capsule count mismatch");
    ok = -1;
  }

  uint64_t local_count = 0;
  uint64_t local_bytes = 0;
  if (ok == 0)
    ok = fread(&local_count, sizeof(local_count), 1, fp) == 1 ? 0 : -1;
  if (ok == 0)
    ok = fread(&local_bytes, sizeof(local_bytes), 1, fp) == 1 ? 0 : -1;

  CapsuleDumpBuffer payload = {0};
  if (ok == 0 && local_bytes > 0) {
    if (local_bytes > SIZE_MAX)
      ok = -1;
    else {
      payload.data = (unsigned char *)malloc((size_t)local_bytes);
      if (!payload.data)
        ok = -1;
      else {
        payload.size = (size_t)local_bytes;
        payload.capacity = payload.size;
      }
    }
  }

  if (ok == 0 && payload.size)
    ok = fread(payload.data, payload.size, 1, fp) == 1 ? 0 : -1;

  fclose(fp);

  if (ok == 0) {
    payload.offset = 0;
    for (uint64_t i = 0; i < local_count; i++) {
      if (capsule_dump_unpack_mesh(&payload) != 0) {
        ok = -1;
        break;
      }
    }
    if (ok == 0 && payload.offset != payload.size)
      ok = -1;
  }

  capsule_dump_buffer_reset(&payload);
  if (ok == 0)
    ibmeshmanager_sync_velocity_coupled_model_outputs();
  return ok;
#endif
}

static inline Checkpointer capsule_dump_checkpointer(void) {
  return (Checkpointer){.filename = CAPSULE_DUMP_FILENAME,
                        .dump_phase = CKPT_PHASE_POST_DUMP,
                        .dump = capsule_dump_write,
                        .restore_phase = CKPT_PHASE_POST_RESTORE,
                        .restore = capsule_dump_restore,
                        .ctx = NULL};
}

static inline int capsule_dump_register(void) {
  return checkpointer_register(capsule_dump_checkpointer());
}

static inline int capsule_dump_register_default(void) {
  return capsule_dump_register();
}

void dump_capsules(const char *fname, FILE *fp) {
  const char default_name[] = "caps.dump";
  const char *name = fname ? fname : default_name;

  if (fp) {
#if _MPI
    assert(npe() == 1 &&
           "dump_capsules(FILE*) is only supported for serial execution");
#endif
    if (capsule_dump_write_stream(fp) != 0) {
      fprintf(ferr, "dump_capsules(): failed to write capsule dump\n");
      exit(1);
    }
    return;
  }

  if (capsule_dump_write(name, NULL) != 0) {
    fprintf(ferr, "dump_capsules(): failed to write capsule dump %s\n", name);
    exit(1);
  }
}

void restore_capsules(const char *filename) {
  if (capsule_dump_restore(filename, NULL) != 0) {
    fprintf(ferr, "restore_capsules(): failed to restore %s\n", filename);
    exit(1);
  }
}
