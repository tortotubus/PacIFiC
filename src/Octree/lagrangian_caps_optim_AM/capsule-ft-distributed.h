#ifndef LAGRANGIAN_CAPS_OPTIM_AM_CAPSULE_FT_DISTRIBUTED_H
#define LAGRANGIAN_CAPS_OPTIM_AM_CAPSULE_FT_DISTRIBUTED_H

/**
# Distributed capsule MPI helpers

This header owns the distributed-capsule routing, sparse lagVel exchange,
owner-to-ghost geometry synchronization, soft lifecycle, and output gather
helpers used by capsule-ft.h.
*/

coord proc_max = {-HUGE, -HUGE, -HUGE};
coord proc_min = {HUGE, HUGE, HUGE};
#if _MPI
coord* all_proc_max = NULL;
coord* all_proc_min = NULL;
int* ncaps_for_proc = NULL;
int* proc_cap_offsets = NULL;
int* proc_cap_ids = NULL;
int proc_cap_ids_nm = 0;
int* owner_to_ghost_send_caps = NULL;
int* owner_to_ghost_send_int_counts = NULL;
int* owner_to_ghost_send_double_counts = NULL;
int* owner_to_ghost_recv_caps = NULL;
int* owner_to_ghost_recv_int_counts = NULL;
int* owner_to_ghost_recv_double_counts = NULL;
int* owner_to_ghost_send_int_offsets = NULL;
int* owner_to_ghost_send_double_offsets = NULL;
int* owner_to_ghost_recv_int_offsets = NULL;
int* owner_to_ghost_recv_double_offsets = NULL;
int* owner_to_ghost_send_int_buffer = NULL;
double* owner_to_ghost_send_double_buffer = NULL;
int* owner_to_ghost_recv_int_buffer = NULL;
double* owner_to_ghost_recv_double_buffer = NULL;
int* ghost_to_owner_send_caps = NULL;
int* ghost_to_owner_send_int_counts = NULL;
int* ghost_to_owner_send_double_counts = NULL;
int* ghost_to_owner_recv_caps = NULL;
int* ghost_to_owner_recv_int_counts = NULL;
int* ghost_to_owner_recv_double_counts = NULL;
int* ghost_to_owner_send_int_offsets = NULL;
int* ghost_to_owner_send_double_offsets = NULL;
int* ghost_to_owner_recv_int_offsets = NULL;
int* ghost_to_owner_recv_double_offsets = NULL;
int* ghost_to_owner_send_int_buffer = NULL;
double* ghost_to_owner_send_double_buffer = NULL;
int* ghost_to_owner_recv_int_buffer = NULL;
double* ghost_to_owner_recv_double_buffer = NULL;
coord** debug_pre_reduce_lagVel = NULL;
int* debug_pre_reduce_lagVel_nln = NULL;
coord** debug_pre_advect_pos = NULL;
int* debug_pre_advect_pos_nln = NULL;
int* debug_capsule_advected_on_this_rank = NULL;
#endif

#ifndef DEBUG_AABB
  #define DEBUG_AABB 0
#endif
#ifndef DEBUG_AABB_FREQ
  #define DEBUG_AABB_FREQ 10000
#endif
#ifndef DEBUG_APPLY_SPARSE_OWNER_LAGVEL
  #define DEBUG_APPLY_SPARSE_OWNER_LAGVEL 0
#endif
#ifndef DEBUG_OWNER_ONLY_ADVECTION
  #define DEBUG_OWNER_ONLY_ADVECTION 0
#endif
#ifndef DEBUG_OWNER_GEOM_TO_ALL_RANKS
  #define DEBUG_OWNER_GEOM_TO_ALL_RANKS 0
#endif
#ifndef DEBUG_OUTPUT_OWNER_GEOM_TO_RANK0
  #define DEBUG_OUTPUT_OWNER_GEOM_TO_RANK0 0
#endif
#ifndef DEBUG_CAPSULE_LIFECYCLE_DRYRUN
  #define DEBUG_CAPSULE_LIFECYCLE_DRYRUN 0
#endif
#ifndef DEBUG_APPLY_SOFT_CAPSULE_LIFECYCLE
  #define DEBUG_APPLY_SOFT_CAPSULE_LIFECYCLE 0
#endif
#ifndef DEBUG_FREE_INACTIVE_CAPSULE_STORAGE
  #define DEBUG_FREE_INACTIVE_CAPSULE_STORAGE 0
#endif
#ifndef DEBUG_SOFT_CAPSULE_LIFECYCLE
  #define DEBUG_SOFT_CAPSULE_LIFECYCLE 0
#endif
#ifndef DEBUG_DISTRIBUTED_GEOMETRY_EXCHANGE
  #define DEBUG_DISTRIBUTED_GEOMETRY_EXCHANGE 0
#endif
#ifndef DEBUG_DRAW_SOFT_INACTIVE_CAPS
  #define DEBUG_DRAW_SOFT_INACTIVE_CAPS 1
#endif
#ifndef DEBUG_DISTRIBUTED_MEMORY_AUDIT
  #define DEBUG_DISTRIBUTED_MEMORY_AUDIT 0
#endif
#ifndef DEBUG_DISTRIBUTED_MEMORY_AUDIT_FREQ
  #define DEBUG_DISTRIBUTED_MEMORY_AUDIT_FREQ DEBUG_AABB_FREQ
#endif
#ifndef DEBUG_LUBRICATION_PAIR_ROUTING
  #define DEBUG_LUBRICATION_PAIR_ROUTING 0
#endif
#ifndef DEBUG_LUBRICATION_PAIR_ROUTING_FREQ
  #define DEBUG_LUBRICATION_PAIR_ROUTING_FREQ DEBUG_AABB_FREQ
#endif
#ifndef DEBUG_SPARSE_OWNER_LAGVEL_ACCUM
  #define DEBUG_SPARSE_OWNER_LAGVEL_ACCUM 0
#endif

#ifndef LAG_DISTRIBUTED_CAPSULES
  #define LAG_DISTRIBUTED_CAPSULES DEBUG_APPLY_SPARSE_OWNER_LAGVEL
#endif
#ifndef LAG_DISTRIBUTED_INITIAL_PLACEMENT
  #define LAG_DISTRIBUTED_INITIAL_PLACEMENT LAG_DISTRIBUTED_CAPSULES
#endif
#ifndef LAG_DISTRIBUTED_OWNER_ADVECTION
  #define LAG_DISTRIBUTED_OWNER_ADVECTION LAG_DISTRIBUTED_CAPSULES
#endif
#ifndef LAG_DISTRIBUTED_SOFT_LIFECYCLE
  #define LAG_DISTRIBUTED_SOFT_LIFECYCLE LAG_DISTRIBUTED_CAPSULES
#endif
#ifndef LAG_DISTRIBUTED_FREE_INACTIVE
  #define LAG_DISTRIBUTED_FREE_INACTIVE LAG_DISTRIBUTED_CAPSULES
#endif
#ifndef LAG_DISTRIBUTED_OUTPUT_GATHER_RANK0
  #define LAG_DISTRIBUTED_OUTPUT_GATHER_RANK0 LAG_DISTRIBUTED_CAPSULES
#endif
#ifndef LAG_DISTRIBUTED_ENSURE_RECEIVED_CAPSULE_TEMPLATE
  #define LAG_DISTRIBUTED_ENSURE_RECEIVED_CAPSULE_TEMPLATE(mesh, cap_id, \
    cap_type, nln, nle, nlt, cap_es, cap_radius) ((void) 0)
#endif

#if _MPI
void debug_ensure_mpi_routing_arrays()
{
  if (all_proc_min == NULL) {
    all_proc_min = (coord*)malloc(npe()*sizeof(coord));
    all_proc_max = (coord*)malloc(npe()*sizeof(coord));
    ncaps_for_proc = (int*)malloc(npe()*sizeof(int));
    proc_cap_offsets = (int*)malloc((npe() + 1)*sizeof(int));
    owner_to_ghost_send_caps = (int*)malloc(npe()*sizeof(int));
    owner_to_ghost_send_int_counts = (int*)malloc(npe()*sizeof(int));
    owner_to_ghost_send_double_counts = (int*)malloc(npe()*sizeof(int));
    owner_to_ghost_recv_caps = (int*)malloc(npe()*sizeof(int));
    owner_to_ghost_recv_int_counts = (int*)malloc(npe()*sizeof(int));
    owner_to_ghost_recv_double_counts = (int*)malloc(npe()*sizeof(int));
    owner_to_ghost_send_int_offsets = (int*)malloc(npe()*sizeof(int));
    owner_to_ghost_send_double_offsets = (int*)malloc(npe()*sizeof(int));
    owner_to_ghost_recv_int_offsets = (int*)malloc(npe()*sizeof(int));
    owner_to_ghost_recv_double_offsets = (int*)malloc(npe()*sizeof(int));
    ghost_to_owner_send_caps = (int*)malloc(npe()*sizeof(int));
    ghost_to_owner_send_int_counts = (int*)malloc(npe()*sizeof(int));
    ghost_to_owner_send_double_counts = (int*)malloc(npe()*sizeof(int));
    ghost_to_owner_recv_caps = (int*)malloc(npe()*sizeof(int));
    ghost_to_owner_recv_int_counts = (int*)malloc(npe()*sizeof(int));
    ghost_to_owner_recv_double_counts = (int*)malloc(npe()*sizeof(int));
    ghost_to_owner_send_int_offsets = (int*)malloc(npe()*sizeof(int));
    ghost_to_owner_send_double_offsets = (int*)malloc(npe()*sizeof(int));
    ghost_to_owner_recv_int_offsets = (int*)malloc(npe()*sizeof(int));
    ghost_to_owner_recv_double_offsets = (int*)malloc(npe()*sizeof(int));
    debug_capsule_advected_on_this_rank =
      (int*)calloc(NCAPS, sizeof(int));
  }
}

void distributed_store_pre_reduce_lagVel(int iter)
{
  if (iter % DEBUG_AABB_FREQ != 0 && !LAG_DISTRIBUTED_CAPSULES)
    return;

  if (debug_pre_reduce_lagVel == NULL) {
    debug_pre_reduce_lagVel = (coord**)calloc(NCAPS, sizeof(coord*));
    debug_pre_reduce_lagVel_nln = (int*)calloc(NCAPS, sizeof(int));
    assert(debug_pre_reduce_lagVel);
    assert(debug_pre_reduce_lagVel_nln);
  }
  for(int cap=0; cap<NCAPS; cap++) {
    if (CAPS(cap).isactive) {
      if (debug_pre_reduce_lagVel_nln[cap] != CAPS(cap).nln) {
        debug_pre_reduce_lagVel[cap] =
          (coord*)realloc(debug_pre_reduce_lagVel[cap],
            CAPS(cap).nln*sizeof(coord));
        assert(debug_pre_reduce_lagVel[cap]);
        debug_pre_reduce_lagVel_nln[cap] = CAPS(cap).nln;
      }
      for(int node_id=0; node_id<CAPS(cap).nln; node_id++)
        debug_pre_reduce_lagVel[cap][node_id] =
          CAPS(cap).nodes[node_id].lagVel;
    }
  }
}

void distributed_prepare_owner_advection()
{
  compute_proc_borders(&proc_max, &proc_min);
  debug_ensure_mpi_routing_arrays();
  gather_all_proc_borders(proc_min, proc_max, all_proc_min, all_proc_max);
  for(int cap=0; cap<NCAPS; cap++)
    debug_capsule_advected_on_this_rank[cap] = false;
}

int distributed_capsule_is_owned_here(lagMesh* mesh)
{
  int owner_proc = find_capsule_owner_proc(mesh, all_proc_min, all_proc_max);
  return owner_proc == pid();
}

int distributed_capsule_has_local_fluid_support(lagMesh* mesh)
{
  compute_proc_borders(&proc_max, &proc_min);
  return lagmesh_bounding_sphere_intersects_box(mesh, proc_min, proc_max);
}

void distributed_mark_capsule_advected(int cap)
{
  debug_capsule_advected_on_this_rank[cap] = true;
}

void debug_lubrication_pair_routing(int iter)
{
  #if DEBUG_LUBRICATION_PAIR_ROUTING && LUBR_VEL == 1
    if (iter % DEBUG_LUBRICATION_PAIR_ROUTING_FREQ != 0)
      return;

    compute_proc_borders(&proc_max, &proc_min);
    debug_ensure_mpi_routing_arrays();
    gather_all_proc_borders(proc_min, proc_max, all_proc_min, all_proc_max);

    int* local_owner_valid = (int*)calloc(NCAPS, sizeof(int));
    int* all_owner_valid = (int*)calloc(NCAPS, sizeof(int));
    int* local_owner_proc = (int*)malloc(NCAPS*sizeof(int));
    int* all_owner_proc = (int*)malloc(NCAPS*sizeof(int));
    double* local_meta = (double*)calloc(4*NCAPS, sizeof(double));
    double* all_meta = (double*)calloc(4*NCAPS, sizeof(double));
    int* local_pair_resident = (int*)calloc(NCAPS*NCAPS, sizeof(int));
    int* all_pair_resident = (int*)calloc(NCAPS*NCAPS, sizeof(int));
    assert(local_owner_valid);
    assert(all_owner_valid);
    assert(local_owner_proc);
    assert(all_owner_proc);
    assert(local_meta);
    assert(all_meta);
    assert(local_pair_resident);
    assert(all_pair_resident);

    for(int cap=0; cap<NCAPS; cap++)
      local_owner_proc[cap] = -1;

    int* local_lubr_available = (int*)calloc(NCAPS, sizeof(int));
    assert(local_lubr_available);

    for(int cap=0; cap<NCAPS; cap++) {
      if (!CAPS(cap).isactive || CAPS(cap).nodes == NULL)
        continue;

      int owner_proc = find_capsule_owner_proc(&CAPS(cap),
        all_proc_min, all_proc_max);
      if (owner_proc == pid()) {
        local_owner_valid[cap] = 1;
        local_owner_proc[cap] = owner_proc;
        local_meta[4*cap] = CAPS(cap).centroid.x;
        local_meta[4*cap + 1] = CAPS(cap).centroid.y;
        local_meta[4*cap + 2] = CAPS(cap).centroid.z;
        local_meta[4*cap + 3] = CAPS(cap).circum_radius;
      }

      local_lubr_available[cap] =
        lagmesh_bounding_sphere_intersects_box(&CAPS(cap),
          proc_min, proc_max);
    }

    for(int a=0; a<NCAPS; a++)
      if (local_lubr_available[a])
        for(int b=a + 1; b<NCAPS; b++)
          if (local_lubr_available[b])
            local_pair_resident[a*NCAPS + b] = 1;

    MPI_Reduce(local_owner_valid, all_owner_valid, NCAPS, MPI_INT, MPI_MAX,
      0, MPI_COMM_WORLD);
    MPI_Reduce(local_owner_proc, all_owner_proc, NCAPS, MPI_INT, MPI_MAX,
      0, MPI_COMM_WORLD);
    MPI_Reduce(local_meta, all_meta, 4*NCAPS, MPI_DOUBLE, MPI_SUM,
      0, MPI_COMM_WORLD);
    MPI_Reduce(local_pair_resident, all_pair_resident, NCAPS*NCAPS, MPI_INT,
      MPI_MAX, 0, MPI_COMM_WORLD);

    if (pid() == 0) {
      int n_sphere_pairs = 0;
      int n_missing_pairs = 0;
      fprintf(stderr, "DEBUG_LUBRICATION_PAIR_ROUTING iter %d missing_pairs=",
        iter);
      for(int a=0; a<NCAPS; a++) {
        if (!all_owner_valid[a])
          continue;
        coord ac = {
          all_meta[4*a],
          all_meta[4*a + 1],
          all_meta[4*a + 2]
        };
        double ar = all_meta[4*a + 3];
        for(int b=a + 1; b<NCAPS; b++) {
          if (!all_owner_valid[b])
            continue;
          coord bc = {
            all_meta[4*b],
            all_meta[4*b + 1],
            all_meta[4*b + 2]
          };
          double br = all_meta[4*b + 3];
          if (bounding_spheres_intersect(ac, ar, bc, br)) {
            n_sphere_pairs++;
            if (!all_pair_resident[a*NCAPS + b]) {
              fprintf(stderr, " (%d,%d owners=%d,%d)", a, b,
                all_owner_proc[a], all_owner_proc[b]);
              n_missing_pairs++;
            }
          }
        }
      }
      if (n_missing_pairs == 0)
        fprintf(stderr, " none");
      fprintf(stderr, " n_sphere_pairs=%d n_missing_pairs=%d\n",
        n_sphere_pairs, n_missing_pairs);
    }

    free(local_owner_valid);
    free(all_owner_valid);
    free(local_owner_proc);
    free(all_owner_proc);
    free(local_meta);
    free(all_meta);
    free(local_pair_resident);
    free(all_pair_resident);
    free(local_lubr_available);
  #else
    (void) iter;
  #endif
}

#if DEBUG_DISTRIBUTED_MEMORY_AUDIT
size_t debug_lagmesh_current_storage_bytes(lagMesh* mesh)
{
  if (!mesh->isactive || mesh->nodes == NULL)
    return 0;

  size_t bytes = sizeof(lagMesh);
  bytes += mesh->nln*sizeof(lagNode);
  bytes += mesh->nle*sizeof(Edge);
  #if dimension > 2
    bytes += mesh->nlt*sizeof(Triangle);
  #endif
  bytes += mesh->nln*STENCIL_SIZE*sizeof(Index);
  bytes += mesh->nln*sizeof(Index);
  return bytes;
}

size_t debug_lag_topology_storage_bytes(lagTopology* topology)
{
  if (!topology)
    return 0;

  size_t bytes = sizeof(lagTopology);
  bytes += topology->nle*sizeof(int[2]);
  #if dimension < 3
    bytes += topology->nln*sizeof(int[2]);
  #else
    bytes += topology->nln*sizeof(int);
    bytes += topology->nln*sizeof(int[6]);
    bytes += topology->nln*sizeof(int[6]);
    bytes += topology->nln*sizeof(int);
    bytes += topology->nln*sizeof(int[6]);
    bytes += topology->nle*sizeof(int[2]);
    bytes += topology->nlt*sizeof(int[3]);
    bytes += topology->nlt*sizeof(int[3]);
  #endif
  return bytes;
}

size_t debug_lag_ref_geometry_storage_bytes(lagReferenceGeometry* ref)
{
  if (!ref)
    return 0;

  size_t bytes = sizeof(lagReferenceGeometry);
  bytes += ref->nle*sizeof(double);
  #if dimension > 2
    bytes += ref->nlt*sizeof(coord[2]);
    bytes += ref->nlt*sizeof(double[3][2]);
  #endif
  return bytes;
}

size_t debug_distributed_mpi_buffer_bytes()
{
  size_t bytes = 0;
  if (all_proc_min) {
    bytes += 2*npe()*sizeof(coord);
    bytes += npe()*sizeof(int);
    bytes += (npe() + 1)*sizeof(int);
    bytes += 24*npe()*sizeof(int);
    bytes += NCAPS*sizeof(int);
  }
  if (proc_cap_ids)
    bytes += proc_cap_ids_nm*sizeof(int);

  int owner_to_ghost_ints = 0, owner_to_ghost_doubles = 0;
  int owner_to_ghost_recv_ints = 0, owner_to_ghost_recv_doubles = 0;
  int ghost_to_owner_ints = 0, ghost_to_owner_doubles = 0;
  int ghost_to_owner_recv_ints = 0, ghost_to_owner_recv_doubles = 0;
  if (owner_to_ghost_send_int_counts)
    for(int p=0; p<npe(); p++) {
      owner_to_ghost_ints += owner_to_ghost_send_int_counts[p];
      owner_to_ghost_doubles += owner_to_ghost_send_double_counts[p];
      owner_to_ghost_recv_ints += owner_to_ghost_recv_int_counts[p];
      owner_to_ghost_recv_doubles += owner_to_ghost_recv_double_counts[p];
      ghost_to_owner_ints += ghost_to_owner_send_int_counts[p];
      ghost_to_owner_doubles += ghost_to_owner_send_double_counts[p];
      ghost_to_owner_recv_ints += ghost_to_owner_recv_int_counts[p];
      ghost_to_owner_recv_doubles += ghost_to_owner_recv_double_counts[p];
    }

  bytes += (owner_to_ghost_ints + owner_to_ghost_recv_ints +
    ghost_to_owner_ints + ghost_to_owner_recv_ints)*sizeof(int);
  bytes += (owner_to_ghost_doubles + owner_to_ghost_recv_doubles +
    ghost_to_owner_doubles + ghost_to_owner_recv_doubles)*sizeof(double);
  return bytes;
}

void debug_print_distributed_memory_audit(int iter)
{
  if (iter % DEBUG_DISTRIBUTED_MEMORY_AUDIT_FREQ != 0)
    return;

  compute_proc_borders(&proc_max, &proc_min);
  debug_ensure_mpi_routing_arrays();
  gather_all_proc_borders(proc_min, proc_max, all_proc_min, all_proc_max);

  int nactive = 0, nowner = 0, nghost = 0, ninactive = 0;
  int nfreed_storage = 0;
  size_t current_bytes = 0;
  size_t snapshot_bytes = 0;

  for(int cap=0; cap<NCAPS; cap++) {
    if (!CAPS(cap).isactive) {
      ninactive++;
      if (CAPS(cap).nodes == NULL && CAPS(cap).edges == NULL)
        nfreed_storage++;
      continue;
    }
    nactive++;
    int owner_proc = find_capsule_owner_proc(&CAPS(cap),
      all_proc_min, all_proc_max);
    if (owner_proc == pid())
      nowner++;
    else
      nghost++;
    current_bytes += debug_lagmesh_current_storage_bytes(&CAPS(cap));
  }

  if (debug_pre_reduce_lagVel) {
    snapshot_bytes += NCAPS*(sizeof(coord*) + sizeof(int));
    for(int cap=0; cap<NCAPS; cap++)
      snapshot_bytes += debug_pre_reduce_lagVel_nln[cap]*sizeof(coord);
  }
  if (debug_pre_advect_pos) {
    snapshot_bytes += NCAPS*(sizeof(coord*) + sizeof(int));
    for(int cap=0; cap<NCAPS; cap++)
      snapshot_bytes += debug_pre_advect_pos_nln[cap]*sizeof(coord);
  }

  size_t topology_bytes = 0;
  for(int i=0; i<lag_topologies.n; i++)
    topology_bytes += debug_lag_topology_storage_bytes(lag_topologies.items[i]);

  size_t ref_bytes = 0;
  for(int i=0; i<lag_ref_geometries.n; i++)
    ref_bytes += debug_lag_ref_geometry_storage_bytes(lag_ref_geometries.items[i]);

  size_t mpi_bytes = debug_distributed_mpi_buffer_bytes();

  unsigned long long local[12] = {
    (unsigned long long)nactive,
    (unsigned long long)nowner,
    (unsigned long long)nghost,
    (unsigned long long)ninactive,
    (unsigned long long)nfreed_storage,
    (unsigned long long)current_bytes,
    (unsigned long long)topology_bytes,
    (unsigned long long)ref_bytes,
    (unsigned long long)snapshot_bytes,
    (unsigned long long)mpi_bytes,
    (unsigned long long)lag_topologies.n,
    (unsigned long long)lag_ref_geometries.n
  };

  unsigned long long* all = NULL;
  if (pid() == 0) {
    all = (unsigned long long*)malloc(npe()*12*sizeof(unsigned long long));
    assert(all);
  }
  MPI_Gather(local, 12, MPI_UNSIGNED_LONG_LONG, all, 12,
    MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);

  if (pid() == 0) {
    unsigned long long totals[12] = {0};
    fprintf(stderr, "DEBUG_DISTRIBUTED_MEMORY_AUDIT iter %d npe %d\n",
      iter, npe());
    for(int p=0; p<npe(); p++) {
      unsigned long long* row = &all[12*p];
      for(int j=0; j<12; j++)
        totals[j] += row[j];
      fprintf(stderr,
        "DEBUG_DISTRIBUTED_MEMORY_AUDIT iter %d rank %d active=%llu owner=%llu ghost=%llu inactive=%llu freed_storage=%llu current_bytes=%llu topology_bytes=%llu ref_bytes=%llu snapshot_bytes=%llu mpi_bytes=%llu ntopologies=%llu nrefs=%llu\n",
        iter, p, row[0], row[1], row[2], row[3], row[4], row[5],
        row[6], row[7], row[8], row[9], row[10], row[11]);
    }
    fprintf(stderr,
      "DEBUG_DISTRIBUTED_MEMORY_AUDIT_TOTAL iter %d active=%llu owner=%llu ghost=%llu inactive=%llu freed_storage=%llu current_bytes=%llu topology_bytes=%llu ref_bytes=%llu snapshot_bytes=%llu mpi_bytes=%llu ntopologies=%llu nrefs=%llu\n",
      iter, totals[0], totals[1], totals[2], totals[3], totals[4],
      totals[5], totals[6], totals[7], totals[8], totals[9],
      totals[10], totals[11]);
    free(all);
  }
}
#else
void debug_print_distributed_memory_audit(int iter)
{
  (void) iter;
}
#endif

void debug_store_pre_advect_pos(int iter)
{
  if (iter % DEBUG_AABB_FREQ != 0)
    return;

  if (debug_pre_advect_pos == NULL) {
    debug_pre_advect_pos = (coord**)calloc(NCAPS, sizeof(coord*));
    debug_pre_advect_pos_nln = (int*)calloc(NCAPS, sizeof(int));
    assert(debug_pre_advect_pos);
    assert(debug_pre_advect_pos_nln);
  }
  for(int cap=0; cap<NCAPS; cap++) {
    if (CAPS(cap).isactive) {
      if (debug_pre_advect_pos_nln[cap] != CAPS(cap).nln) {
        debug_pre_advect_pos[cap] =
          (coord*)realloc(debug_pre_advect_pos[cap],
            CAPS(cap).nln*sizeof(coord));
        assert(debug_pre_advect_pos[cap]);
        debug_pre_advect_pos_nln[cap] = CAPS(cap).nln;
      }
      for(int node_id=0; node_id<CAPS(cap).nln; node_id++)
        debug_pre_advect_pos[cap][node_id] =
          CAPS(cap).nodes[node_id].pos;
    }
  }
}

void debug_sparse_owner_lagVel_from_recv_payloads(int iter,
  int total_ghost_to_owner_recv_caps, int compare_euler_move,
  int print_debug)
{
  if (total_ghost_to_owner_recv_caps > 0) {
    if (print_debug)
      fprintf(stderr,
        "DEBUG_OWNER_LAGVEL_ACCUM_DRYRUN pid %d/%d iter %d caps=",
        pid(), npe(), iter);
    for(int cap=0; cap<NCAPS; cap++) {
      if (!CAPS(cap).isactive)
        continue;
      int owner_proc = find_capsule_owner_proc(&CAPS(cap),
        all_proc_min, all_proc_max);
      if (owner_proc != pid())
        continue;

      coord* lagvel_sum = (coord*)calloc(CAPS(cap).nln, sizeof(coord));
      assert(lagvel_sum);
      coord* base_lagVel =
        debug_pre_reduce_lagVel &&
        debug_pre_reduce_lagVel_nln[cap] == CAPS(cap).nln ?
        debug_pre_reduce_lagVel[cap] : NULL;
      double base_lagvel_abs_sum = 0.;
      for(int node_id=0; node_id<CAPS(cap).nln; node_id++)
        foreach_dimension() {
          lagvel_sum[node_id].x = base_lagVel ?
            base_lagVel[node_id].x : CAPS(cap).nodes[node_id].lagVel.x;
          base_lagvel_abs_sum += fabs(lagvel_sum[node_id].x);
        }

      int recv_caps_added = 0;
      int recv_sources_added = 0;
      for(int p=0; p<npe(); p++) {
        int int_pos = ghost_to_owner_recv_int_offsets[p];
        int double_pos = ghost_to_owner_recv_double_offsets[p];
        for(int q=0; q<ghost_to_owner_recv_caps[p]; q++) {
          int recv_cap_id = ghost_to_owner_recv_int_buffer[int_pos++];
          int recv_nln = ghost_to_owner_recv_int_buffer[int_pos++];
          int add_payload = recv_cap_id == CAPS(cap).cap_id &&
            recv_nln == CAPS(cap).nln;
          if (add_payload) {
            recv_caps_added++;
            recv_sources_added += p != pid();
          }
          for(int node_id=0; node_id<recv_nln; node_id++)
            foreach_dimension() {
              double vel_comp =
                ghost_to_owner_recv_double_buffer[double_pos++];
              if (add_payload)
                lagvel_sum[node_id].x += vel_comp;
            }
        }
      }

      coord accum_vel_min = {HUGE, HUGE, HUGE};
      coord accum_vel_max = {-HUGE, -HUGE, -HUGE};
      double accum_lagvel_abs_sum = 0.;
      double reduced_lagvel_abs_sum = 0.;
      double reduced_vs_accum_max_abs_diff = 0.;
      double owner_euler_move_max_abs_diff = 0.;
      coord* pre_advect_pos =
        compare_euler_move && debug_pre_advect_pos &&
        debug_pre_advect_pos_nln[cap] == CAPS(cap).nln ?
        debug_pre_advect_pos[cap] : NULL;
      for(int node_id=0; node_id<CAPS(cap).nln; node_id++)
        foreach_dimension() {
          accum_vel_min.x = min(accum_vel_min.x, lagvel_sum[node_id].x);
          accum_vel_max.x = max(accum_vel_max.x, lagvel_sum[node_id].x);
          accum_lagvel_abs_sum += fabs(lagvel_sum[node_id].x);
          reduced_lagvel_abs_sum += fabs(CAPS(cap).nodes[node_id].lagVel.x);
          reduced_vs_accum_max_abs_diff = max(
            reduced_vs_accum_max_abs_diff,
            fabs(CAPS(cap).nodes[node_id].lagVel.x - lagvel_sum[node_id].x));
          if (pre_advect_pos) {
            double predicted_pos =
              pre_advect_pos[node_id].x + dt*lagvel_sum[node_id].x;
            owner_euler_move_max_abs_diff = max(
              owner_euler_move_max_abs_diff,
              fabs(CAPS(cap).nodes[node_id].pos.x - predicted_pos));
          }
        }

      int applied_sparse_owner_lagVel = false;
      #if LAG_DISTRIBUTED_CAPSULES
        for(int node_id=0; node_id<CAPS(cap).nln; node_id++)
          CAPS(cap).nodes[node_id].lagVel = lagvel_sum[node_id];
        applied_sparse_owner_lagVel = true;
      #endif
      if (print_debug)
        fprintf(stderr,
          " cap=%d,owner=%d,nln=%d,base_abs_sum=%g,recv_caps_added=%d,recv_sources_added=%d,accum_vel_min=(%g %g %g),accum_vel_max=(%g %g %g),accum_abs_sum=%g,reduced_abs_sum=%g,reduced_vs_accum_max_abs_diff=%g,owner_euler_move_max_abs_diff=%g,applied_sparse_owner_lagVel=%d",
          CAPS(cap).cap_id, owner_proc, CAPS(cap).nln,
          base_lagvel_abs_sum, recv_caps_added, recv_sources_added,
          accum_vel_min.x, accum_vel_min.y, accum_vel_min.z,
          accum_vel_max.x, accum_vel_max.y, accum_vel_max.z,
          accum_lagvel_abs_sum, reduced_lagvel_abs_sum,
          reduced_vs_accum_max_abs_diff, owner_euler_move_max_abs_diff,
          applied_sparse_owner_lagVel);
      free(lagvel_sum);
    }
    if (print_debug)
      fprintf(stderr, "\n");
  }

  if (!print_debug)
    return;

  int local_owner_accum_ints[6] = {-1, -1, 0, 0, 0,
    LAG_DISTRIBUTED_CAPSULES};
  double local_owner_accum_doubles[11] =
    {0., HUGE, HUGE, HUGE, -HUGE, -HUGE, -HUGE, 0., 0., 0., 0.};
  if (total_ghost_to_owner_recv_caps > 0) {
    for(int cap=0; cap<NCAPS; cap++) {
      if (!CAPS(cap).isactive)
        continue;
      int owner_proc = find_capsule_owner_proc(&CAPS(cap),
        all_proc_min, all_proc_max);
      if (owner_proc != pid())
        continue;

      coord* lagvel_sum = (coord*)calloc(CAPS(cap).nln, sizeof(coord));
      assert(lagvel_sum);
      coord* base_lagVel =
        debug_pre_reduce_lagVel &&
        debug_pre_reduce_lagVel_nln[cap] == CAPS(cap).nln ?
        debug_pre_reduce_lagVel[cap] : NULL;
      double base_lagvel_abs_sum = 0.;
      for(int node_id=0; node_id<CAPS(cap).nln; node_id++)
        foreach_dimension() {
          lagvel_sum[node_id].x = base_lagVel ?
            base_lagVel[node_id].x : CAPS(cap).nodes[node_id].lagVel.x;
          base_lagvel_abs_sum += fabs(lagvel_sum[node_id].x);
        }

      int recv_caps_added = 0;
      int recv_sources_added = 0;
      for(int p=0; p<npe(); p++) {
        int int_pos = ghost_to_owner_recv_int_offsets[p];
        int double_pos = ghost_to_owner_recv_double_offsets[p];
        for(int q=0; q<ghost_to_owner_recv_caps[p]; q++) {
          int recv_cap_id = ghost_to_owner_recv_int_buffer[int_pos++];
          int recv_nln = ghost_to_owner_recv_int_buffer[int_pos++];
          int add_payload = recv_cap_id == CAPS(cap).cap_id &&
            recv_nln == CAPS(cap).nln;
          if (add_payload) {
            recv_caps_added++;
            recv_sources_added += p != pid();
          }
          for(int node_id=0; node_id<recv_nln; node_id++)
            foreach_dimension() {
              double vel_comp =
                ghost_to_owner_recv_double_buffer[double_pos++];
              if (add_payload)
                lagvel_sum[node_id].x += vel_comp;
            }
        }
      }

      coord accum_vel_min = {HUGE, HUGE, HUGE};
      coord accum_vel_max = {-HUGE, -HUGE, -HUGE};
      double accum_lagvel_abs_sum = 0.;
      double reduced_lagvel_abs_sum = 0.;
      double reduced_vs_accum_max_abs_diff = 0.;
      double owner_euler_move_max_abs_diff = 0.;
      coord* pre_advect_pos =
        compare_euler_move && debug_pre_advect_pos &&
        debug_pre_advect_pos_nln[cap] == CAPS(cap).nln ?
        debug_pre_advect_pos[cap] : NULL;
      for(int node_id=0; node_id<CAPS(cap).nln; node_id++)
        foreach_dimension() {
          accum_vel_min.x = min(accum_vel_min.x, lagvel_sum[node_id].x);
          accum_vel_max.x = max(accum_vel_max.x, lagvel_sum[node_id].x);
          accum_lagvel_abs_sum += fabs(lagvel_sum[node_id].x);
          reduced_lagvel_abs_sum += fabs(CAPS(cap).nodes[node_id].lagVel.x);
          reduced_vs_accum_max_abs_diff = max(
            reduced_vs_accum_max_abs_diff,
            fabs(CAPS(cap).nodes[node_id].lagVel.x - lagvel_sum[node_id].x));
          if (pre_advect_pos) {
            double predicted_pos =
              pre_advect_pos[node_id].x + dt*lagvel_sum[node_id].x;
            owner_euler_move_max_abs_diff = max(
              owner_euler_move_max_abs_diff,
              fabs(CAPS(cap).nodes[node_id].pos.x - predicted_pos));
          }
        }

      local_owner_accum_ints[0] = CAPS(cap).cap_id;
      local_owner_accum_ints[1] = owner_proc;
      local_owner_accum_ints[2] = CAPS(cap).nln;
      local_owner_accum_ints[3] = recv_caps_added;
      local_owner_accum_ints[4] = recv_sources_added;
      local_owner_accum_ints[5] = LAG_DISTRIBUTED_CAPSULES;
      local_owner_accum_doubles[0] = base_lagvel_abs_sum;
      local_owner_accum_doubles[1] = accum_vel_min.x;
      local_owner_accum_doubles[2] = accum_vel_min.y;
      local_owner_accum_doubles[3] = accum_vel_min.z;
      local_owner_accum_doubles[4] = accum_vel_max.x;
      local_owner_accum_doubles[5] = accum_vel_max.y;
      local_owner_accum_doubles[6] = accum_vel_max.z;
      local_owner_accum_doubles[7] = accum_lagvel_abs_sum;
      local_owner_accum_doubles[8] = reduced_lagvel_abs_sum;
      local_owner_accum_doubles[9] = reduced_vs_accum_max_abs_diff;
      local_owner_accum_doubles[10] = owner_euler_move_max_abs_diff;
      free(lagvel_sum);
      break;
    }
  }

  int* all_owner_accum_ints = NULL;
  double* all_owner_accum_doubles = NULL;
  if (pid() == 0) {
    all_owner_accum_ints = (int*)malloc(npe()*6*sizeof(int));
    all_owner_accum_doubles = (double*)malloc(npe()*11*sizeof(double));
    assert(all_owner_accum_ints);
    assert(all_owner_accum_doubles);
  }
  MPI_Gather(local_owner_accum_ints, 6, MPI_INT,
    all_owner_accum_ints, 6, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(local_owner_accum_doubles, 11, MPI_DOUBLE,
    all_owner_accum_doubles, 11, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (pid() == 0) {
    int printed_owner_accum = false;
    for(int rank=0; rank<npe(); rank++) {
      int* row_i = all_owner_accum_ints + rank*6;
      double* row_d = all_owner_accum_doubles + rank*11;
      if (row_i[0] < 0)
        continue;
      if (!printed_owner_accum) {
        fprintf(stderr,
          "DEBUG_ALL_OWNER_LAGVEL_ACCUM_DRYRUN iter %d", iter);
        printed_owner_accum = true;
      }
      fprintf(stderr,
        " rank=%d,cap=%d,owner=%d,nln=%d,base_abs_sum=%g,recv_caps_added=%d,recv_sources_added=%d,accum_vel_min=(%g %g %g),accum_vel_max=(%g %g %g),accum_abs_sum=%g,reduced_abs_sum=%g,reduced_vs_accum_max_abs_diff=%g,owner_euler_move_max_abs_diff=%g,applied_sparse_owner_lagVel=%d",
        rank, row_i[0], row_i[1], row_i[2], row_d[0],
        row_i[3], row_i[4], row_d[1], row_d[2], row_d[3],
        row_d[4], row_d[5], row_d[6], row_d[7], row_d[8],
        row_d[9], row_d[10], row_i[5]);
    }
    if (printed_owner_accum)
      fprintf(stderr, "\n");
  }
  free(all_owner_accum_ints);
  free(all_owner_accum_doubles);
}

void debug_sparse_owner_lagVel_exchange_before_advection(int iter)
{
  compute_proc_borders(&proc_max, &proc_min);
  debug_ensure_mpi_routing_arrays();
  gather_all_proc_borders(proc_min, proc_max, all_proc_min, all_proc_max);
  for(int p=0; p<npe(); p++) {
    ghost_to_owner_send_caps[p] = 0;
    ghost_to_owner_send_int_counts[p] = 0;
    ghost_to_owner_send_double_counts[p] = 0;
    ghost_to_owner_recv_caps[p] = 0;
    ghost_to_owner_recv_int_counts[p] = 0;
    ghost_to_owner_recv_double_counts[p] = 0;
  }

  for(int cap=0; cap<NCAPS; cap++) {
    if (CAPS(cap).isactive) {
      int owner_proc = find_capsule_owner_proc(&CAPS(cap),
        all_proc_min, all_proc_max);
      if (owner_proc >= 0 && owner_proc != pid()) {
        int intersects_proc = lagmesh_bounding_sphere_intersects_box(
          &CAPS(cap), proc_min, proc_max);
        if (intersects_proc) {
          ghost_to_owner_send_caps[owner_proc]++;
          ghost_to_owner_send_int_counts[owner_proc] +=
            estimate_ghost_to_owner_nints(&CAPS(cap));
          ghost_to_owner_send_double_counts[owner_proc] +=
            estimate_ghost_to_owner_ndoubles(&CAPS(cap));
        }
      }
    }
  }

  MPI_Alltoall(ghost_to_owner_send_caps, 1, MPI_INT,
    ghost_to_owner_recv_caps, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Alltoall(ghost_to_owner_send_int_counts, 1, MPI_INT,
    ghost_to_owner_recv_int_counts, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Alltoall(ghost_to_owner_send_double_counts, 1, MPI_INT,
    ghost_to_owner_recv_double_counts, 1, MPI_INT, MPI_COMM_WORLD);

  int total_ghost_to_owner_ints = 0;
  int total_ghost_to_owner_doubles = 0;
  int total_ghost_to_owner_recv_caps = 0;
  int total_ghost_to_owner_recv_ints = 0;
  int total_ghost_to_owner_recv_doubles = 0;
  for(int p=0; p<npe(); p++) {
    ghost_to_owner_send_int_offsets[p] = total_ghost_to_owner_ints;
    ghost_to_owner_send_double_offsets[p] = total_ghost_to_owner_doubles;
    ghost_to_owner_recv_int_offsets[p] = total_ghost_to_owner_recv_ints;
    ghost_to_owner_recv_double_offsets[p] = total_ghost_to_owner_recv_doubles;
    total_ghost_to_owner_ints += ghost_to_owner_send_int_counts[p];
    total_ghost_to_owner_doubles += ghost_to_owner_send_double_counts[p];
    total_ghost_to_owner_recv_caps += ghost_to_owner_recv_caps[p];
    total_ghost_to_owner_recv_ints += ghost_to_owner_recv_int_counts[p];
    total_ghost_to_owner_recv_doubles += ghost_to_owner_recv_double_counts[p];
  }

  if (total_ghost_to_owner_ints > 0)
    ghost_to_owner_send_int_buffer = (int*)realloc(
      ghost_to_owner_send_int_buffer,
      total_ghost_to_owner_ints*sizeof(int));
  else {
    free(ghost_to_owner_send_int_buffer);
    ghost_to_owner_send_int_buffer = NULL;
  }
  if (total_ghost_to_owner_doubles > 0)
    ghost_to_owner_send_double_buffer = (double*)realloc(
      ghost_to_owner_send_double_buffer,
      total_ghost_to_owner_doubles*sizeof(double));
  else {
    free(ghost_to_owner_send_double_buffer);
    ghost_to_owner_send_double_buffer = NULL;
  }
  if (total_ghost_to_owner_recv_ints > 0)
    ghost_to_owner_recv_int_buffer = (int*)realloc(
      ghost_to_owner_recv_int_buffer,
      total_ghost_to_owner_recv_ints*sizeof(int));
  else {
    free(ghost_to_owner_recv_int_buffer);
    ghost_to_owner_recv_int_buffer = NULL;
  }
  if (total_ghost_to_owner_recv_doubles > 0)
    ghost_to_owner_recv_double_buffer = (double*)realloc(
      ghost_to_owner_recv_double_buffer,
      total_ghost_to_owner_recv_doubles*sizeof(double));
  else {
    free(ghost_to_owner_recv_double_buffer);
    ghost_to_owner_recv_double_buffer = NULL;
  }

  int* ghost_to_owner_pack_int_pos = (int*)malloc(npe()*sizeof(int));
  int* ghost_to_owner_pack_double_pos = (int*)malloc(npe()*sizeof(int));
  assert(ghost_to_owner_pack_int_pos);
  assert(ghost_to_owner_pack_double_pos);
  for(int p=0; p<npe(); p++) {
    ghost_to_owner_pack_int_pos[p] = ghost_to_owner_send_int_offsets[p];
    ghost_to_owner_pack_double_pos[p] = ghost_to_owner_send_double_offsets[p];
  }

  for(int cap=0; cap<NCAPS; cap++) {
    if (CAPS(cap).isactive) {
      int owner_proc = find_capsule_owner_proc(&CAPS(cap),
        all_proc_min, all_proc_max);
      if (owner_proc >= 0 && owner_proc != pid()) {
        int intersects_proc = lagmesh_bounding_sphere_intersects_box(
          &CAPS(cap), proc_min, proc_max);
        if (intersects_proc) {
          coord* pre_reduce_lagVel =
            debug_pre_reduce_lagVel &&
            debug_pre_reduce_lagVel_nln[cap] == CAPS(cap).nln ?
            debug_pre_reduce_lagVel[cap] : NULL;
          pack_ghost_to_owner_capsule_lagVel(&CAPS(cap), pre_reduce_lagVel,
            ghost_to_owner_send_int_buffer,
            &ghost_to_owner_pack_int_pos[owner_proc],
            ghost_to_owner_send_double_buffer,
            &ghost_to_owner_pack_double_pos[owner_proc]);
        }
      }
    }
  }
  free(ghost_to_owner_pack_int_pos);
  free(ghost_to_owner_pack_double_pos);

  int dummy_int_buffer = 0;
  double dummy_double_buffer = 0.;
  int* ghost_to_owner_send_int_exchange =
    ghost_to_owner_send_int_buffer ?
    ghost_to_owner_send_int_buffer : &dummy_int_buffer;
  double* ghost_to_owner_send_double_exchange =
    ghost_to_owner_send_double_buffer ?
    ghost_to_owner_send_double_buffer : &dummy_double_buffer;
  int* ghost_to_owner_recv_int_exchange =
    ghost_to_owner_recv_int_buffer ?
    ghost_to_owner_recv_int_buffer : &dummy_int_buffer;
  double* ghost_to_owner_recv_double_exchange =
    ghost_to_owner_recv_double_buffer ?
    ghost_to_owner_recv_double_buffer : &dummy_double_buffer;

  MPI_Alltoallv(ghost_to_owner_send_int_exchange,
    ghost_to_owner_send_int_counts, ghost_to_owner_send_int_offsets,
    MPI_INT, ghost_to_owner_recv_int_exchange,
    ghost_to_owner_recv_int_counts, ghost_to_owner_recv_int_offsets,
    MPI_INT, MPI_COMM_WORLD);
  MPI_Alltoallv(ghost_to_owner_send_double_exchange,
    ghost_to_owner_send_double_counts, ghost_to_owner_send_double_offsets,
    MPI_DOUBLE, ghost_to_owner_recv_double_exchange,
    ghost_to_owner_recv_double_counts, ghost_to_owner_recv_double_offsets,
    MPI_DOUBLE, MPI_COMM_WORLD);

  debug_sparse_owner_lagVel_from_recv_payloads(iter,
    total_ghost_to_owner_recv_caps, false,
    DEBUG_SPARSE_OWNER_LAGVEL_ACCUM && iter % DEBUG_AABB_FREQ == 0);
}

void debug_owner_advection_dryrun(int iter)
{
  if (iter % DEBUG_AABB_FREQ != 0)
    return;

  compute_proc_borders(&proc_max, &proc_min);
  debug_ensure_mpi_routing_arrays();
  gather_all_proc_borders(proc_min, proc_max, all_proc_min, all_proc_max);

  if (pid() == 0) {
    fprintf(stderr, "DEBUG_OWNER_ADVECTION_DRYRUN iter %d cap_owners=", iter);
    for(int cap=0; cap<NCAPS; cap++) {
      if (CAPS(cap).isactive) {
        int owner_proc = find_capsule_owner_proc(&CAPS(cap),
          all_proc_min, all_proc_max);
        fprintf(stderr, " cap=%d:owner=%d", CAPS(cap).cap_id, owner_proc);
      }
    }
    fprintf(stderr, "\n");
  }

  int nowned = 0;
  fprintf(stderr,
    "DEBUG_LOCAL_OWNER_ADVECTION_DRYRUN pid %d/%d iter %d would_advect=",
    pid(), npe(), iter);
  for(int cap=0; cap<NCAPS; cap++) {
    if (CAPS(cap).isactive) {
      int owner_proc = find_capsule_owner_proc(&CAPS(cap),
        all_proc_min, all_proc_max);
      if (owner_proc == pid()) {
        fprintf(stderr, " cap=%d", CAPS(cap).cap_id);
        nowned++;
      }
    }
  }
  if (nowned == 0)
    fprintf(stderr, " none");
  fprintf(stderr, " nowned=%d\n", nowned);
}

void debug_capsule_lifecycle_dryrun(int iter)
{
  #if DEBUG_CAPSULE_LIFECYCLE_DRYRUN
    if (iter % DEBUG_AABB_FREQ != 0)
      return;

    compute_proc_borders(&proc_max, &proc_min);
    debug_ensure_mpi_routing_arrays();
    gather_all_proc_borders(proc_min, proc_max, all_proc_min, all_proc_max);

    int nkeep_owner = 0;
    int nkeep_ghost = 0;
    int ncreate_ghost = 0;
    int ndestroy = 0;

    fprintf(stderr,
      "DEBUG_CAPSULE_LIFECYCLE_DRYRUN pid %d/%d iter %d actions=",
      pid(), npe(), iter);
    for(int cap=0; cap<NCAPS; cap++) {
      int local_exists = CAPS(cap).isactive;
      int owner_proc = -1;
      int intersects_local = false;
      int should_exist = false;

      if (local_exists) {
        owner_proc = find_capsule_owner_proc(&CAPS(cap),
          all_proc_min, all_proc_max);
        intersects_local = lagmesh_bounding_sphere_intersects_box(
          &CAPS(cap), proc_min, proc_max);
        should_exist = owner_proc == pid() || intersects_local;
      }

      const char* action = "inactive";
      if (local_exists && owner_proc == pid()) {
        action = "keep_owner";
        nkeep_owner++;
      }
      else if (local_exists && should_exist) {
        action = "keep_ghost";
        nkeep_ghost++;
      }
      else if (local_exists) {
        action = "destroy_local_copy";
        ndestroy++;
      }
      else if (should_exist) {
        action = "create_ghost";
        ncreate_ghost++;
      }

      fprintf(stderr,
        " cap=%d,owner=%d,intersects=%d,local_exists=%d,action=%s",
        local_exists ? CAPS(cap).cap_id : cap, owner_proc,
        intersects_local, local_exists, action);
    }
    fprintf(stderr,
      " totals=(owner=%d ghost=%d create=%d destroy=%d)\n",
      nkeep_owner, nkeep_ghost, ncreate_ghost, ndestroy);
  #endif
}

void debug_apply_soft_capsule_lifecycle(int iter)
{
  #if LAG_DISTRIBUTED_SOFT_LIFECYCLE
    compute_proc_borders(&proc_max, &proc_min);
    debug_ensure_mpi_routing_arrays();
    gather_all_proc_borders(proc_min, proc_max, all_proc_min, all_proc_max);

    int ndeactivated = 0;
    int print_debug = DEBUG_SOFT_CAPSULE_LIFECYCLE &&
      iter % DEBUG_AABB_FREQ == 0;
    if (print_debug)
      fprintf(stderr,
        "DEBUG_SOFT_CAPSULE_LIFECYCLE_APPLY pid %d/%d iter %d deactivate=",
        pid(), npe(), iter);
    for(int cap=0; cap<NCAPS; cap++) {
      if (!CAPS(cap).isactive)
        continue;
      int owner_proc = find_capsule_owner_proc(&CAPS(cap),
        all_proc_min, all_proc_max);
      int intersects_local = lagmesh_bounding_sphere_intersects_box(
        &CAPS(cap), proc_min, proc_max);
      if (owner_proc >= 0 && owner_proc != pid() && !intersects_local) {
        if (print_debug)
          fprintf(stderr, " cap=%d,owner=%d,local_index=%d",
            CAPS(cap).cap_id, owner_proc, cap);
        #if LAG_DISTRIBUTED_FREE_INACTIVE
          free_lagMesh_current_storage_keep_slot(&CAPS(cap));
        #else
          CAPS(cap).isactive = false;
        #endif
        ndeactivated++;
      }
    }
    if (print_debug) {
      if (ndeactivated == 0)
        fprintf(stderr, " none");
      fprintf(stderr, " ndeactivated=%d\n", ndeactivated);
    }
  #endif
}

void debug_print_local_capsule_lifecycle_counts(int iter)
{
  #if LAG_DISTRIBUTED_SOFT_LIFECYCLE
    if (!DEBUG_SOFT_CAPSULE_LIFECYCLE)
      return;
    if (iter % DEBUG_AABB_FREQ != 0)
      return;

    compute_proc_borders(&proc_max, &proc_min);
    debug_ensure_mpi_routing_arrays();
    gather_all_proc_borders(proc_min, proc_max, all_proc_min, all_proc_max);

    int nactive = 0;
    int nowner = 0;
    int nghost = 0;
    int ninactive = 0;
    int nfreed_storage = 0;
    for(int cap=0; cap<NCAPS; cap++) {
      if (!CAPS(cap).isactive) {
        ninactive++;
        if (CAPS(cap).nodes == NULL && CAPS(cap).edges == NULL)
          nfreed_storage++;
        continue;
      }
      nactive++;
      int owner_proc = find_capsule_owner_proc(&CAPS(cap),
        all_proc_min, all_proc_max);
      if (owner_proc == pid())
        nowner++;
      else
        nghost++;
    }

    fprintf(stderr,
      "DEBUG_SOFT_CAPSULE_COUNTS pid %d/%d iter %d active=%d owner=%d ghost=%d inactive=%d freed_storage=%d all_caps=%d\n",
      pid(), npe(), iter, nactive, nowner, nghost, ninactive,
      nfreed_storage, NCAPS);

    int local_counts[5] = {
      nactive, nowner, nghost, ninactive, nfreed_storage
    };
    int* all_counts = NULL;
    if (pid() == 0) {
      all_counts = (int*)malloc(npe()*5*sizeof(int));
      assert(all_counts);
    }
    MPI_Gather(local_counts, 5, MPI_INT, all_counts, 5, MPI_INT,
      0, MPI_COMM_WORLD);
    if (pid() == 0) {
      fprintf(stderr, "DEBUG_ALL_SOFT_CAPSULE_COUNTS iter %d npe %d\n",
        iter, npe());
      for(int p=0; p<npe(); p++) {
        int off = 5*p;
        fprintf(stderr,
          "DEBUG_ALL_SOFT_CAPSULE_COUNTS iter %d rank %d active=%d owner=%d ghost=%d inactive=%d freed_storage=%d all_caps=%d\n",
          iter, p, all_counts[off], all_counts[off + 1],
          all_counts[off + 2], all_counts[off + 3],
          all_counts[off + 4], NCAPS);
      }
      free(all_counts);
    }
  #endif
}

void debug_update_local_capsule_from_owner_payload(int* int_data, int* int_pos,
  double* double_data, int* double_pos, int iter)
{
  int cap_id, cap_type, nln, nle, nlt;
  double cap_es, cap_radius, circum_radius;
  unpack_owner_to_ghost_header(int_data, int_pos, double_data, double_pos,
    &cap_id, &cap_type, &nln, &nle, &nlt, &cap_es, &cap_radius,
    &circum_radius);

  int local_cap = -1;
  for(int cap=0; cap<NCAPS; cap++)
    if (CAPS(cap).cap_id == cap_id) {
      local_cap = cap;
      break;
    }

  coord recv_centroid = {0};
  foreach_dimension()
    recv_centroid.x = double_data[(*double_pos)++];

  int can_update = local_cap >= 0;
  if (can_update) {
    LAG_DISTRIBUTED_ENSURE_RECEIVED_CAPSULE_TEMPLATE(&CAPS(local_cap),
      cap_id, cap_type, nln, nle, nlt, cap_es, cap_radius);
    ensure_lagMesh_current_storage(&CAPS(local_cap), nln, nle, nlt);
  }

  can_update = can_update &&
    CAPS(local_cap).nln == nln && CAPS(local_cap).nle == nle &&
    CAPS(local_cap).nodes != NULL && CAPS(local_cap).edges != NULL;
  #if dimension > 2
    can_update = can_update && CAPS(local_cap).nlt == nlt &&
      CAPS(local_cap).triangles != NULL;
  #endif

  for(int node_id=0; node_id<nln; node_id++)
    foreach_dimension() {
      double pos_comp = double_data[(*double_pos)++];
      if (can_update)
        CAPS(local_cap).nodes[node_id].pos.x = pos_comp;
    }
  for(int node_id=0; node_id<nln; node_id++)
    foreach_dimension() {
      double force_comp = double_data[(*double_pos)++];
      if (can_update)
        CAPS(local_cap).nodes[node_id].lagForce.x = force_comp;
    }

  if (can_update) {
    int was_active = CAPS(local_cap).isactive;
    CAPS(local_cap).isactive = true;
    CAPS(local_cap).cap_type = cap_type;
    CAPS(local_cap).cap_es = cap_es;
    CAPS(local_cap).cap_radius = cap_radius;
    CAPS(local_cap).centroid = recv_centroid;
    CAPS(local_cap).circum_radius = circum_radius;
    CAPS(local_cap).updated_stretches = false;
    CAPS(local_cap).updated_normals = false;
    CAPS(local_cap).updated_curvatures = false;
    comp_capsule_geodynamics(&CAPS(local_cap));
    #if LAG_DISTRIBUTED_SOFT_LIFECYCLE
      if (DEBUG_SOFT_CAPSULE_LIFECYCLE && !was_active &&
        iter % DEBUG_AABB_FREQ == 0)
        fprintf(stderr,
          "DEBUG_SOFT_CAPSULE_LIFECYCLE_REACTIVATE pid %d/%d cap=%d local_index=%d\n",
          pid(), npe(), cap_id, local_cap);
    #endif
  }
}

void debug_owner_to_ghost_geometry_exchange_after_advection(int iter)
{
  compute_proc_borders(&proc_max, &proc_min);
  debug_ensure_mpi_routing_arrays();
  gather_all_proc_borders(proc_min, proc_max, all_proc_min, all_proc_max);

  for(int p=0; p<npe(); p++) {
    owner_to_ghost_send_caps[p] = 0;
    owner_to_ghost_send_int_counts[p] = 0;
    owner_to_ghost_send_double_counts[p] = 0;
    owner_to_ghost_recv_caps[p] = 0;
    owner_to_ghost_recv_int_counts[p] = 0;
    owner_to_ghost_recv_double_counts[p] = 0;
  }

  for(int cap=0; cap<NCAPS; cap++) {
    if (CAPS(cap).isactive) {
      int owner_proc = find_capsule_owner_proc(&CAPS(cap),
        all_proc_min, all_proc_max);
      int sender_proc = debug_capsule_advected_on_this_rank != NULL &&
        debug_capsule_advected_on_this_rank[cap];
      if (sender_proc) {
        int nints = estimate_owner_to_ghost_nints(&CAPS(cap));
        int ndoubles = estimate_owner_to_ghost_ndoubles(&CAPS(cap));
        for(int p=0; p<npe(); p++) {
          int intersects_proc = lagmesh_bounding_sphere_intersects_box(
            &CAPS(cap), all_proc_min[p], all_proc_max[p]);
          int send_geometry = DEBUG_OWNER_GEOM_TO_ALL_RANKS ||
            intersects_proc || p == owner_proc;
          if (send_geometry && p != pid()) {
            owner_to_ghost_send_caps[p]++;
            owner_to_ghost_send_int_counts[p] += nints;
            owner_to_ghost_send_double_counts[p] += ndoubles;
          }
        }
      }
    }
  }

  MPI_Alltoall(owner_to_ghost_send_caps, 1, MPI_INT,
    owner_to_ghost_recv_caps, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Alltoall(owner_to_ghost_send_int_counts, 1, MPI_INT,
    owner_to_ghost_recv_int_counts, 1, MPI_INT, MPI_COMM_WORLD);
  MPI_Alltoall(owner_to_ghost_send_double_counts, 1, MPI_INT,
    owner_to_ghost_recv_double_counts, 1, MPI_INT, MPI_COMM_WORLD);

  int total_owner_to_ghost_ints = 0;
  int total_owner_to_ghost_doubles = 0;
  int total_owner_to_ghost_recv_caps = 0;
  int total_owner_to_ghost_recv_ints = 0;
  int total_owner_to_ghost_recv_doubles = 0;
  for(int p=0; p<npe(); p++) {
    owner_to_ghost_send_int_offsets[p] = total_owner_to_ghost_ints;
    owner_to_ghost_send_double_offsets[p] = total_owner_to_ghost_doubles;
    owner_to_ghost_recv_int_offsets[p] = total_owner_to_ghost_recv_ints;
    owner_to_ghost_recv_double_offsets[p] = total_owner_to_ghost_recv_doubles;
    total_owner_to_ghost_ints += owner_to_ghost_send_int_counts[p];
    total_owner_to_ghost_doubles += owner_to_ghost_send_double_counts[p];
    total_owner_to_ghost_recv_caps += owner_to_ghost_recv_caps[p];
    total_owner_to_ghost_recv_ints += owner_to_ghost_recv_int_counts[p];
    total_owner_to_ghost_recv_doubles += owner_to_ghost_recv_double_counts[p];
  }

  if (total_owner_to_ghost_ints > 0)
    owner_to_ghost_send_int_buffer = (int*)realloc(
      owner_to_ghost_send_int_buffer,
      total_owner_to_ghost_ints*sizeof(int));
  else {
    free(owner_to_ghost_send_int_buffer);
    owner_to_ghost_send_int_buffer = NULL;
  }
  if (total_owner_to_ghost_doubles > 0)
    owner_to_ghost_send_double_buffer = (double*)realloc(
      owner_to_ghost_send_double_buffer,
      total_owner_to_ghost_doubles*sizeof(double));
  else {
    free(owner_to_ghost_send_double_buffer);
    owner_to_ghost_send_double_buffer = NULL;
  }
  if (total_owner_to_ghost_recv_ints > 0)
    owner_to_ghost_recv_int_buffer = (int*)realloc(
      owner_to_ghost_recv_int_buffer,
      total_owner_to_ghost_recv_ints*sizeof(int));
  else {
    free(owner_to_ghost_recv_int_buffer);
    owner_to_ghost_recv_int_buffer = NULL;
  }
  if (total_owner_to_ghost_recv_doubles > 0)
    owner_to_ghost_recv_double_buffer = (double*)realloc(
      owner_to_ghost_recv_double_buffer,
      total_owner_to_ghost_recv_doubles*sizeof(double));
  else {
    free(owner_to_ghost_recv_double_buffer);
    owner_to_ghost_recv_double_buffer = NULL;
  }

  int* pack_int_pos = (int*)malloc(npe()*sizeof(int));
  int* pack_double_pos = (int*)malloc(npe()*sizeof(int));
  assert(pack_int_pos);
  assert(pack_double_pos);
  for(int p=0; p<npe(); p++) {
    pack_int_pos[p] = owner_to_ghost_send_int_offsets[p];
    pack_double_pos[p] = owner_to_ghost_send_double_offsets[p];
  }
  for(int cap=0; cap<NCAPS; cap++) {
    if (CAPS(cap).isactive) {
      int owner_proc = find_capsule_owner_proc(&CAPS(cap),
        all_proc_min, all_proc_max);
      int sender_proc = debug_capsule_advected_on_this_rank != NULL &&
        debug_capsule_advected_on_this_rank[cap];
      if (sender_proc) {
        for(int p=0; p<npe(); p++) {
          int intersects_proc = lagmesh_bounding_sphere_intersects_box(
            &CAPS(cap), all_proc_min[p], all_proc_max[p]);
          int send_geometry = DEBUG_OWNER_GEOM_TO_ALL_RANKS ||
            intersects_proc || p == owner_proc;
          if (send_geometry && p != pid())
            pack_owner_to_ghost_capsule(&CAPS(cap),
              owner_to_ghost_send_int_buffer, &pack_int_pos[p],
              owner_to_ghost_send_double_buffer, &pack_double_pos[p]);
        }
      }
    }
  }
  free(pack_int_pos);
  free(pack_double_pos);

  int dummy_int_buffer = 0;
  double dummy_double_buffer = 0.;
  int* send_int_exchange = owner_to_ghost_send_int_buffer ?
    owner_to_ghost_send_int_buffer : &dummy_int_buffer;
  double* send_double_exchange = owner_to_ghost_send_double_buffer ?
    owner_to_ghost_send_double_buffer : &dummy_double_buffer;
  int* recv_int_exchange = owner_to_ghost_recv_int_buffer ?
    owner_to_ghost_recv_int_buffer : &dummy_int_buffer;
  double* recv_double_exchange = owner_to_ghost_recv_double_buffer ?
    owner_to_ghost_recv_double_buffer : &dummy_double_buffer;

  MPI_Alltoallv(send_int_exchange,
    owner_to_ghost_send_int_counts, owner_to_ghost_send_int_offsets,
    MPI_INT, recv_int_exchange,
    owner_to_ghost_recv_int_counts, owner_to_ghost_recv_int_offsets,
    MPI_INT, MPI_COMM_WORLD);
  MPI_Alltoallv(send_double_exchange,
    owner_to_ghost_send_double_counts, owner_to_ghost_send_double_offsets,
    MPI_DOUBLE, recv_double_exchange,
    owner_to_ghost_recv_double_counts, owner_to_ghost_recv_double_offsets,
    MPI_DOUBLE, MPI_COMM_WORLD);

  for(int p=0; p<npe(); p++) {
    int int_pos = owner_to_ghost_recv_int_offsets[p];
    int double_pos = owner_to_ghost_recv_double_offsets[p];
    for(int q=0; q<owner_to_ghost_recv_caps[p]; q++)
      debug_update_local_capsule_from_owner_payload(
        owner_to_ghost_recv_int_buffer, &int_pos,
        owner_to_ghost_recv_double_buffer, &double_pos, iter);
  }

  if (DEBUG_DISTRIBUTED_GEOMETRY_EXCHANGE &&
    iter % DEBUG_AABB_FREQ == 0) {
    fprintf(stderr,
      "DEBUG_OWNER_TO_GHOST_GEOM_APPLY pid %d/%d iter %d recv_caps=%d\n",
      pid(), npe(), iter, total_owner_to_ghost_recv_caps);
  }
}

void debug_owner_geometry_exchange_to_rank0_for_output(int iter)
{
  #if LAG_DISTRIBUTED_OUTPUT_GATHER_RANK0
    compute_proc_borders(&proc_max, &proc_min);
    debug_ensure_mpi_routing_arrays();
    gather_all_proc_borders(proc_min, proc_max, all_proc_min, all_proc_max);

    int send_caps = 0;
    int send_int_count = 0;
    int send_double_count = 0;

    for(int cap=0; cap<NCAPS; cap++) {
      if (CAPS(cap).isactive) {
        int owner_proc = find_capsule_owner_proc(&CAPS(cap),
          all_proc_min, all_proc_max);
        if (owner_proc == pid()) {
          send_caps++;
          send_int_count += estimate_owner_to_ghost_nints(&CAPS(cap));
          send_double_count += estimate_owner_to_ghost_ndoubles(&CAPS(cap));
        }
      }
    }

    int* recv_caps = (int*)calloc(npe(), sizeof(int));
    int* recv_int_counts = (int*)calloc(npe(), sizeof(int));
    int* recv_double_counts = (int*)calloc(npe(), sizeof(int));
    assert(recv_caps);
    assert(recv_int_counts);
    assert(recv_double_counts);

    MPI_Gather(&send_caps, 1, MPI_INT, recv_caps, 1, MPI_INT,
      0, MPI_COMM_WORLD);
    MPI_Gather(&send_int_count, 1, MPI_INT, recv_int_counts, 1, MPI_INT,
      0, MPI_COMM_WORLD);
    MPI_Gather(&send_double_count, 1, MPI_INT, recv_double_counts, 1,
      MPI_INT, 0, MPI_COMM_WORLD);

    int* recv_int_offsets = (int*)calloc(npe(), sizeof(int));
    int* recv_double_offsets = (int*)calloc(npe(), sizeof(int));
    int total_recv_caps = 0;
    int total_recv_ints = 0;
    int total_recv_doubles = 0;
    assert(recv_int_offsets);
    assert(recv_double_offsets);
    if (pid() == 0) {
      for(int p=0; p<npe(); p++) {
        recv_int_offsets[p] = total_recv_ints;
        recv_double_offsets[p] = total_recv_doubles;
        total_recv_caps += recv_caps[p];
        total_recv_ints += recv_int_counts[p];
        total_recv_doubles += recv_double_counts[p];
      }
    }

    int* send_int_buffer = NULL;
    double* send_double_buffer = NULL;
    if (send_int_count > 0) {
      send_int_buffer = (int*)malloc(send_int_count*sizeof(int));
      assert(send_int_buffer);
    }
    if (send_double_count > 0) {
      send_double_buffer = (double*)malloc(send_double_count*sizeof(double));
      assert(send_double_buffer);
    }

    int send_int_pos = 0;
    int send_double_pos = 0;
    for(int cap=0; cap<NCAPS; cap++) {
      if (CAPS(cap).isactive) {
        int owner_proc = find_capsule_owner_proc(&CAPS(cap),
          all_proc_min, all_proc_max);
        if (owner_proc == pid())
          pack_owner_to_ghost_capsule(&CAPS(cap), send_int_buffer,
            &send_int_pos, send_double_buffer, &send_double_pos);
      }
    }

    int* recv_int_buffer = NULL;
    double* recv_double_buffer = NULL;
    if (pid() == 0 && total_recv_ints > 0) {
      recv_int_buffer = (int*)malloc(total_recv_ints*sizeof(int));
      assert(recv_int_buffer);
    }
    if (pid() == 0 && total_recv_doubles > 0) {
      recv_double_buffer = (double*)malloc(total_recv_doubles*sizeof(double));
      assert(recv_double_buffer);
    }

    int dummy_int = 0;
    double dummy_double = 0.;
    int* send_int_exchange = send_int_buffer ? send_int_buffer : &dummy_int;
    double* send_double_exchange = send_double_buffer ?
      send_double_buffer : &dummy_double;
    int* recv_int_exchange = recv_int_buffer ? recv_int_buffer : &dummy_int;
    double* recv_double_exchange = recv_double_buffer ?
      recv_double_buffer : &dummy_double;

    MPI_Gatherv(send_int_exchange, send_int_count, MPI_INT,
      recv_int_exchange, recv_int_counts, recv_int_offsets, MPI_INT,
      0, MPI_COMM_WORLD);
    MPI_Gatherv(send_double_exchange, send_double_count, MPI_DOUBLE,
      recv_double_exchange, recv_double_counts, recv_double_offsets,
      MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (pid() == 0) {
      for(int p=0; p<npe(); p++) {
        int int_pos = recv_int_offsets[p];
        int double_pos = recv_double_offsets[p];
        for(int q=0; q<recv_caps[p]; q++)
          debug_update_local_capsule_from_owner_payload(recv_int_buffer,
            &int_pos, recv_double_buffer, &double_pos, iter);
      }
      if (DEBUG_DISTRIBUTED_GEOMETRY_EXCHANGE &&
        iter % DEBUG_AABB_FREQ == 0) {
        int total_bytes = total_recv_ints*((int) sizeof(int)) +
          total_recv_doubles*((int) sizeof(double));
        fprintf(stderr,
          "DEBUG_OUTPUT_OWNER_GEOM_TO_RANK0 iter %d recv_caps=%d recv_ints=%d recv_doubles=%d bytes=%d\n",
          iter, total_recv_caps, total_recv_ints, total_recv_doubles,
          total_bytes);
      }
    }

    free(send_int_buffer);
    free(send_double_buffer);
    free(recv_int_buffer);
    free(recv_double_buffer);
    free(recv_caps);
    free(recv_int_counts);
    free(recv_double_counts);
    free(recv_int_offsets);
    free(recv_double_offsets);
  #endif
}

void free_distributed_capsule_mpi_storage()
{
  free(all_proc_min);
  free(all_proc_max);
  free(ncaps_for_proc);
  free(proc_cap_offsets);
  free(proc_cap_ids);
  free(owner_to_ghost_send_caps);
  free(owner_to_ghost_send_int_counts);
  free(owner_to_ghost_send_double_counts);
  free(owner_to_ghost_recv_caps);
  free(owner_to_ghost_recv_int_counts);
  free(owner_to_ghost_recv_double_counts);
  free(owner_to_ghost_send_int_offsets);
  free(owner_to_ghost_send_double_offsets);
  free(owner_to_ghost_recv_int_offsets);
  free(owner_to_ghost_recv_double_offsets);
  free(owner_to_ghost_send_int_buffer);
  free(owner_to_ghost_send_double_buffer);
  free(owner_to_ghost_recv_int_buffer);
  free(owner_to_ghost_recv_double_buffer);
  free(ghost_to_owner_send_caps);
  free(ghost_to_owner_send_int_counts);
  free(ghost_to_owner_send_double_counts);
  free(ghost_to_owner_recv_caps);
  free(ghost_to_owner_recv_int_counts);
  free(ghost_to_owner_recv_double_counts);
  free(ghost_to_owner_send_int_offsets);
  free(ghost_to_owner_send_double_offsets);
  free(ghost_to_owner_recv_int_offsets);
  free(ghost_to_owner_recv_double_offsets);
  free(ghost_to_owner_send_int_buffer);
  free(ghost_to_owner_send_double_buffer);
  free(ghost_to_owner_recv_int_buffer);
  free(ghost_to_owner_recv_double_buffer);
  if (debug_pre_reduce_lagVel) {
    for(int cap=0; cap<NCAPS; cap++)
      free(debug_pre_reduce_lagVel[cap]);
  }
  free(debug_pre_reduce_lagVel);
  free(debug_pre_reduce_lagVel_nln);
  if (debug_pre_advect_pos) {
    for(int cap=0; cap<NCAPS; cap++)
      free(debug_pre_advect_pos[cap]);
  }
  free(debug_pre_advect_pos);
  free(debug_pre_advect_pos_nln);
  free(debug_capsule_advected_on_this_rank);
}

void distributed_lagVel_exchange_before_advection(int iter)
{
  debug_sparse_owner_lagVel_exchange_before_advection(iter);
}

void distributed_geometry_exchange_after_advection(int iter)
{
  debug_owner_to_ghost_geometry_exchange_after_advection(iter);
}

void distributed_lifecycle_after_advection(int iter)
{
  debug_capsule_lifecycle_dryrun(iter);
  debug_apply_soft_capsule_lifecycle(iter);
  debug_print_local_capsule_lifecycle_counts(iter);
  debug_print_distributed_memory_audit(iter);
}
#endif


void debug_print_distributed_aabb_state(int i)
{
  /* Compute borders of the curren proc */
  compute_proc_borders(&proc_max, &proc_min);
  #if DEBUG_AABB
    if (i % DEBUG_AABB_FREQ == 0) {
      #if _MPI
        debug_ensure_mpi_routing_arrays();
        gather_all_proc_borders(proc_min, proc_max, all_proc_min, all_proc_max);
        for(int p=0; p<npe(); p++)
          ncaps_for_proc[p] = 0;
        for(int cap=0; cap<NCAPS; cap++) {
          if (CAPS(cap).isactive) {
            for(int p=0; p<npe(); p++) {
              bool intersects_proc = lagmesh_bounding_sphere_intersects_box(
                &CAPS(cap), all_proc_min[p], all_proc_max[p]);
              if (intersects_proc)
                ncaps_for_proc[p]++;
            }
          }
        }
        proc_cap_offsets[0] = 0;
        for(int p=0; p<npe(); p++)
          proc_cap_offsets[p + 1] = proc_cap_offsets[p] + ncaps_for_proc[p];
        int total_proc_cap_routes = proc_cap_offsets[npe()];
        if (total_proc_cap_routes > proc_cap_ids_nm) {
          proc_cap_ids_nm = total_proc_cap_routes;
          proc_cap_ids = (int*)realloc(proc_cap_ids,
            proc_cap_ids_nm*sizeof(int));
        }
        for(int p=0; p<npe(); p++)
          ncaps_for_proc[p] = 0;
        for(int cap=0; cap<NCAPS; cap++) {
          if (CAPS(cap).isactive) {
            for(int p=0; p<npe(); p++) {
              bool intersects_proc = lagmesh_bounding_sphere_intersects_box(
                &CAPS(cap), all_proc_min[p], all_proc_max[p]);
              if (intersects_proc) {
                int slot = proc_cap_offsets[p] + ncaps_for_proc[p]++;
                proc_cap_ids[slot] = cap;
              }
            }
          }
        }
        for(int p=0; p<npe(); p++) {
          owner_to_ghost_send_caps[p] = 0;
          owner_to_ghost_send_int_counts[p] = 0;
          owner_to_ghost_send_double_counts[p] = 0;
          owner_to_ghost_recv_caps[p] = 0;
          owner_to_ghost_recv_int_counts[p] = 0;
          owner_to_ghost_recv_double_counts[p] = 0;
          ghost_to_owner_send_caps[p] = 0;
          ghost_to_owner_send_int_counts[p] = 0;
          ghost_to_owner_send_double_counts[p] = 0;
          ghost_to_owner_recv_caps[p] = 0;
          ghost_to_owner_recv_int_counts[p] = 0;
          ghost_to_owner_recv_double_counts[p] = 0;
        }
        for(int cap=0; cap<NCAPS; cap++) {
          if (CAPS(cap).isactive) {
            int owner_proc = find_capsule_owner_proc(&CAPS(cap),
              all_proc_min, all_proc_max);
            if (owner_proc >= 0) {
              int owner_to_ghost_nints = estimate_owner_to_ghost_nints(&CAPS(cap));
              int owner_to_ghost_ndoubles = estimate_owner_to_ghost_ndoubles(&CAPS(cap));
              int ghost_to_owner_nints = estimate_ghost_to_owner_nints(&CAPS(cap));
              int ghost_to_owner_ndoubles = estimate_ghost_to_owner_ndoubles(&CAPS(cap));
              for(int p=0; p<npe(); p++) {
                bool intersects_proc = lagmesh_bounding_sphere_intersects_box(
                  &CAPS(cap), all_proc_min[p], all_proc_max[p]);
                if (intersects_proc && p != owner_proc) {
                  if (owner_proc == pid()) {
                    owner_to_ghost_send_caps[p]++;
                    owner_to_ghost_send_int_counts[p] += owner_to_ghost_nints;
                    owner_to_ghost_send_double_counts[p] += owner_to_ghost_ndoubles;
                  }
                  if (p == pid()) {
                    ghost_to_owner_send_caps[owner_proc]++;
                    ghost_to_owner_send_int_counts[owner_proc] += ghost_to_owner_nints;
                    ghost_to_owner_send_double_counts[owner_proc] += ghost_to_owner_ndoubles;
                  }
                }
              }
            }
          }
        }
        MPI_Alltoall(owner_to_ghost_send_caps, 1, MPI_INT,
          owner_to_ghost_recv_caps, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Alltoall(owner_to_ghost_send_int_counts, 1, MPI_INT,
          owner_to_ghost_recv_int_counts, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Alltoall(owner_to_ghost_send_double_counts, 1, MPI_INT,
          owner_to_ghost_recv_double_counts, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Alltoall(ghost_to_owner_send_caps, 1, MPI_INT,
          ghost_to_owner_recv_caps, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Alltoall(ghost_to_owner_send_int_counts, 1, MPI_INT,
          ghost_to_owner_recv_int_counts, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Alltoall(ghost_to_owner_send_double_counts, 1, MPI_INT,
          ghost_to_owner_recv_double_counts, 1, MPI_INT, MPI_COMM_WORLD);
      #endif
      fprintf(stderr,
        "DEBUG_AABB pid %d/%d iter %d proc_min=(%g %g %g) proc_max=(%g %g %g)\n",
        pid(), npe(), i,
        proc_min.x, proc_min.y, proc_min.z,
        proc_max.x, proc_max.y, proc_max.z);
      #if _MPI
        int total_owner_to_ghost_caps = 0;
        int total_owner_to_ghost_ints = 0;
        int total_owner_to_ghost_doubles = 0;
        int total_ghost_to_owner_caps = 0;
        int total_ghost_to_owner_ints = 0;
        int total_ghost_to_owner_doubles = 0;
        int total_owner_to_ghost_recv_caps = 0;
        int total_owner_to_ghost_recv_ints = 0;
        int total_owner_to_ghost_recv_doubles = 0;
        int total_ghost_to_owner_recv_caps = 0;
        int total_ghost_to_owner_recv_ints = 0;
        int total_ghost_to_owner_recv_doubles = 0;
        for(int p=0; p<npe(); p++) {
          owner_to_ghost_send_int_offsets[p] = total_owner_to_ghost_ints;
          owner_to_ghost_send_double_offsets[p] = total_owner_to_ghost_doubles;
          ghost_to_owner_send_int_offsets[p] = total_ghost_to_owner_ints;
          ghost_to_owner_send_double_offsets[p] = total_ghost_to_owner_doubles;
          owner_to_ghost_recv_int_offsets[p] = total_owner_to_ghost_recv_ints;
          owner_to_ghost_recv_double_offsets[p] = total_owner_to_ghost_recv_doubles;
          ghost_to_owner_recv_int_offsets[p] = total_ghost_to_owner_recv_ints;
          ghost_to_owner_recv_double_offsets[p] = total_ghost_to_owner_recv_doubles;
          total_owner_to_ghost_caps += owner_to_ghost_send_caps[p];
          total_owner_to_ghost_ints += owner_to_ghost_send_int_counts[p];
          total_owner_to_ghost_doubles += owner_to_ghost_send_double_counts[p];
          total_ghost_to_owner_caps += ghost_to_owner_send_caps[p];
          total_ghost_to_owner_ints += ghost_to_owner_send_int_counts[p];
          total_ghost_to_owner_doubles += ghost_to_owner_send_double_counts[p];
          total_owner_to_ghost_recv_caps += owner_to_ghost_recv_caps[p];
          total_owner_to_ghost_recv_ints += owner_to_ghost_recv_int_counts[p];
          total_owner_to_ghost_recv_doubles += owner_to_ghost_recv_double_counts[p];
          total_ghost_to_owner_recv_caps += ghost_to_owner_recv_caps[p];
          total_ghost_to_owner_recv_ints += ghost_to_owner_recv_int_counts[p];
          total_ghost_to_owner_recv_doubles += ghost_to_owner_recv_double_counts[p];
        }
        if (total_owner_to_ghost_ints > 0)
          owner_to_ghost_send_int_buffer = (int*)realloc(
            owner_to_ghost_send_int_buffer,
            total_owner_to_ghost_ints*sizeof(int));
        else {
          free(owner_to_ghost_send_int_buffer);
          owner_to_ghost_send_int_buffer = NULL;
        }
        if (total_owner_to_ghost_doubles > 0)
          owner_to_ghost_send_double_buffer = (double*)realloc(
            owner_to_ghost_send_double_buffer,
            total_owner_to_ghost_doubles*sizeof(double));
        else {
          free(owner_to_ghost_send_double_buffer);
          owner_to_ghost_send_double_buffer = NULL;
        }
        if (total_owner_to_ghost_recv_ints > 0)
          owner_to_ghost_recv_int_buffer = (int*)realloc(
            owner_to_ghost_recv_int_buffer,
            total_owner_to_ghost_recv_ints*sizeof(int));
        else {
          free(owner_to_ghost_recv_int_buffer);
          owner_to_ghost_recv_int_buffer = NULL;
        }
        if (total_owner_to_ghost_recv_doubles > 0)
          owner_to_ghost_recv_double_buffer = (double*)realloc(
            owner_to_ghost_recv_double_buffer,
            total_owner_to_ghost_recv_doubles*sizeof(double));
        else {
          free(owner_to_ghost_recv_double_buffer);
          owner_to_ghost_recv_double_buffer = NULL;
        }
        if (total_ghost_to_owner_ints > 0)
          ghost_to_owner_send_int_buffer = (int*)realloc(
            ghost_to_owner_send_int_buffer,
            total_ghost_to_owner_ints*sizeof(int));
        else {
          free(ghost_to_owner_send_int_buffer);
          ghost_to_owner_send_int_buffer = NULL;
        }
        if (total_ghost_to_owner_doubles > 0)
          ghost_to_owner_send_double_buffer = (double*)realloc(
            ghost_to_owner_send_double_buffer,
            total_ghost_to_owner_doubles*sizeof(double));
        else {
          free(ghost_to_owner_send_double_buffer);
          ghost_to_owner_send_double_buffer = NULL;
        }
        if (total_ghost_to_owner_recv_ints > 0)
          ghost_to_owner_recv_int_buffer = (int*)realloc(
            ghost_to_owner_recv_int_buffer,
            total_ghost_to_owner_recv_ints*sizeof(int));
        else {
          free(ghost_to_owner_recv_int_buffer);
          ghost_to_owner_recv_int_buffer = NULL;
        }
        if (total_ghost_to_owner_recv_doubles > 0)
          ghost_to_owner_recv_double_buffer = (double*)realloc(
            ghost_to_owner_recv_double_buffer,
            total_ghost_to_owner_recv_doubles*sizeof(double));
        else {
          free(ghost_to_owner_recv_double_buffer);
          ghost_to_owner_recv_double_buffer = NULL;
        }
        int* owner_to_ghost_pack_int_pos = (int*)malloc(npe()*sizeof(int));
        int* owner_to_ghost_pack_double_pos = (int*)malloc(npe()*sizeof(int));
        for(int p=0; p<npe(); p++) {
          owner_to_ghost_pack_int_pos[p] = owner_to_ghost_send_int_offsets[p];
          owner_to_ghost_pack_double_pos[p] = owner_to_ghost_send_double_offsets[p];
        }
        if (total_owner_to_ghost_caps > 0) {
          for(int cap=0; cap<NCAPS; cap++) {
            if (CAPS(cap).isactive) {
              int owner_proc = find_capsule_owner_proc(&CAPS(cap),
                all_proc_min, all_proc_max);
              if (owner_proc == pid()) {
                for(int p=0; p<npe(); p++) {
                  bool intersects_proc = lagmesh_bounding_sphere_intersects_box(
                    &CAPS(cap), all_proc_min[p], all_proc_max[p]);
                  if (intersects_proc && p != owner_proc)
                    pack_owner_to_ghost_capsule(&CAPS(cap),
                      owner_to_ghost_send_int_buffer,
                      &owner_to_ghost_pack_int_pos[p],
                      owner_to_ghost_send_double_buffer,
                      &owner_to_ghost_pack_double_pos[p]);
                }
              }
            }
          }
          fprintf(stderr,
            "DEBUG_OWNER_TO_GHOST_PACKED pid %d/%d iter %d dests=",
            pid(), npe(), i);
          for(int p=0; p<npe(); p++) {
            int int_used = owner_to_ghost_pack_int_pos[p]
              - owner_to_ghost_send_int_offsets[p];
            int double_used = owner_to_ghost_pack_double_pos[p]
              - owner_to_ghost_send_double_offsets[p];
            if (owner_to_ghost_send_caps[p] > 0)
              fprintf(stderr,
                " %d:caps=%d,int_used=%d/%d,double_used=%d/%d",
                p, owner_to_ghost_send_caps[p],
                int_used, owner_to_ghost_send_int_counts[p],
                double_used, owner_to_ghost_send_double_counts[p]);
          }
          fprintf(stderr, "\n");
        }
        int local_pack_row_len = 2*npe();
        int* local_pack_row = (int*)calloc(local_pack_row_len, sizeof(int));
        int* all_pack_rows = NULL;
        for(int p=0; p<npe(); p++) {
          local_pack_row[p] = owner_to_ghost_pack_int_pos[p]
            - owner_to_ghost_send_int_offsets[p];
          local_pack_row[npe() + p] = owner_to_ghost_pack_double_pos[p]
            - owner_to_ghost_send_double_offsets[p];
        }
        if (pid() == 0)
          all_pack_rows = (int*)malloc(npe()*local_pack_row_len*sizeof(int));
        MPI_Gather(local_pack_row, local_pack_row_len, MPI_INT,
          all_pack_rows, local_pack_row_len, MPI_INT, 0, MPI_COMM_WORLD);
        if (pid() == 0) {
          for(int rank=0; rank<npe(); rank++) {
            int* row = all_pack_rows + rank*local_pack_row_len;
            int packed_ints = 0, packed_doubles = 0;
            for(int p=0; p<npe(); p++) {
              packed_ints += row[p];
              packed_doubles += row[npe() + p];
            }
            if (packed_ints > 0 || packed_doubles > 0) {
              fprintf(stderr,
                "DEBUG_ALL_OWNER_TO_GHOST_PACKED iter %d rank %d dests=",
                i, rank);
              for(int p=0; p<npe(); p++)
                if (row[p] > 0 || row[npe() + p] > 0)
                  fprintf(stderr, " %d:int_used=%d,double_used=%d",
                    p, row[p], row[npe() + p]);
              fprintf(stderr, " total_int_used=%d total_double_used=%d\n",
                packed_ints, packed_doubles);
            }
          }
        }
        free(local_pack_row);
        free(all_pack_rows);
        free(owner_to_ghost_pack_int_pos);
        free(owner_to_ghost_pack_double_pos);
        int* ghost_to_owner_pack_int_pos = (int*)malloc(npe()*sizeof(int));
        int* ghost_to_owner_pack_double_pos = (int*)malloc(npe()*sizeof(int));
        for(int p=0; p<npe(); p++) {
          ghost_to_owner_pack_int_pos[p] = ghost_to_owner_send_int_offsets[p];
          ghost_to_owner_pack_double_pos[p] = ghost_to_owner_send_double_offsets[p];
        }
        if (total_ghost_to_owner_caps > 0) {
          for(int cap=0; cap<NCAPS; cap++) {
            if (CAPS(cap).isactive) {
              int owner_proc = find_capsule_owner_proc(&CAPS(cap),
                all_proc_min, all_proc_max);
              if (owner_proc >= 0 && owner_proc != pid()) {
                bool intersects_proc = lagmesh_bounding_sphere_intersects_box(
                  &CAPS(cap), proc_min, proc_max);
                if (intersects_proc) {
                  coord* pre_reduce_lagVel =
                    debug_pre_reduce_lagVel &&
                    debug_pre_reduce_lagVel_nln[cap] == CAPS(cap).nln ?
                    debug_pre_reduce_lagVel[cap] : NULL;
                  pack_ghost_to_owner_capsule_lagVel(&CAPS(cap),
                    pre_reduce_lagVel,
                    ghost_to_owner_send_int_buffer,
                    &ghost_to_owner_pack_int_pos[owner_proc],
                    ghost_to_owner_send_double_buffer,
                    &ghost_to_owner_pack_double_pos[owner_proc]);
                }
              }
            }
          }
          fprintf(stderr,
            "DEBUG_GHOST_TO_OWNER_PACKED pid %d/%d iter %d owners=",
            pid(), npe(), i);
          for(int p=0; p<npe(); p++) {
            int int_used = ghost_to_owner_pack_int_pos[p]
              - ghost_to_owner_send_int_offsets[p];
            int double_used = ghost_to_owner_pack_double_pos[p]
              - ghost_to_owner_send_double_offsets[p];
            if (ghost_to_owner_send_caps[p] > 0)
              fprintf(stderr,
                " %d:caps=%d,int_used=%d/%d,double_used=%d/%d",
                p, ghost_to_owner_send_caps[p],
                int_used, ghost_to_owner_send_int_counts[p],
                double_used, ghost_to_owner_send_double_counts[p]);
          }
          fprintf(stderr, "\n");
        }
        int local_ghost_pack_row_len = 2*npe();
        int* local_ghost_pack_row = (int*)calloc(local_ghost_pack_row_len, sizeof(int));
        int* all_ghost_pack_rows = NULL;
        for(int p=0; p<npe(); p++) {
          local_ghost_pack_row[p] = ghost_to_owner_pack_int_pos[p]
            - ghost_to_owner_send_int_offsets[p];
          local_ghost_pack_row[npe() + p] = ghost_to_owner_pack_double_pos[p]
            - ghost_to_owner_send_double_offsets[p];
        }
        if (pid() == 0)
          all_ghost_pack_rows = (int*)malloc(
            npe()*local_ghost_pack_row_len*sizeof(int));
        MPI_Gather(local_ghost_pack_row, local_ghost_pack_row_len, MPI_INT,
          all_ghost_pack_rows, local_ghost_pack_row_len, MPI_INT, 0,
          MPI_COMM_WORLD);
        if (pid() == 0) {
          for(int rank=0; rank<npe(); rank++) {
            int* row = all_ghost_pack_rows + rank*local_ghost_pack_row_len;
            int packed_ints = 0, packed_doubles = 0;
            for(int p=0; p<npe(); p++) {
              packed_ints += row[p];
              packed_doubles += row[npe() + p];
            }
            if (packed_ints > 0 || packed_doubles > 0) {
              fprintf(stderr,
                "DEBUG_ALL_GHOST_TO_OWNER_PACKED iter %d rank %d owners=",
                i, rank);
              for(int p=0; p<npe(); p++)
                if (row[p] > 0 || row[npe() + p] > 0)
                  fprintf(stderr, " %d:int_used=%d,double_used=%d",
                    p, row[p], row[npe() + p]);
              fprintf(stderr, " total_int_used=%d total_double_used=%d\n",
                packed_ints, packed_doubles);
            }
          }
        }
        free(local_ghost_pack_row);
        free(all_ghost_pack_rows);
        free(ghost_to_owner_pack_int_pos);
        free(ghost_to_owner_pack_double_pos);
        int dummy_int_buffer = 0;
        double dummy_double_buffer = 0.;
        int* owner_to_ghost_send_int_exchange =
          owner_to_ghost_send_int_buffer ?
          owner_to_ghost_send_int_buffer : &dummy_int_buffer;
        double* owner_to_ghost_send_double_exchange =
          owner_to_ghost_send_double_buffer ?
          owner_to_ghost_send_double_buffer : &dummy_double_buffer;
        int* owner_to_ghost_recv_int_exchange =
          owner_to_ghost_recv_int_buffer ?
          owner_to_ghost_recv_int_buffer : &dummy_int_buffer;
        double* owner_to_ghost_recv_double_exchange =
          owner_to_ghost_recv_double_buffer ?
          owner_to_ghost_recv_double_buffer : &dummy_double_buffer;
        int* ghost_to_owner_send_int_exchange =
          ghost_to_owner_send_int_buffer ?
          ghost_to_owner_send_int_buffer : &dummy_int_buffer;
        double* ghost_to_owner_send_double_exchange =
          ghost_to_owner_send_double_buffer ?
          ghost_to_owner_send_double_buffer : &dummy_double_buffer;
        int* ghost_to_owner_recv_int_exchange =
          ghost_to_owner_recv_int_buffer ?
          ghost_to_owner_recv_int_buffer : &dummy_int_buffer;
        double* ghost_to_owner_recv_double_exchange =
          ghost_to_owner_recv_double_buffer ?
          ghost_to_owner_recv_double_buffer : &dummy_double_buffer;

        MPI_Alltoallv(owner_to_ghost_send_int_exchange,
          owner_to_ghost_send_int_counts, owner_to_ghost_send_int_offsets,
          MPI_INT, owner_to_ghost_recv_int_exchange,
          owner_to_ghost_recv_int_counts, owner_to_ghost_recv_int_offsets,
          MPI_INT, MPI_COMM_WORLD);
        MPI_Alltoallv(owner_to_ghost_send_double_exchange,
          owner_to_ghost_send_double_counts, owner_to_ghost_send_double_offsets,
          MPI_DOUBLE, owner_to_ghost_recv_double_exchange,
          owner_to_ghost_recv_double_counts, owner_to_ghost_recv_double_offsets,
          MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Alltoallv(ghost_to_owner_send_int_exchange,
          ghost_to_owner_send_int_counts, ghost_to_owner_send_int_offsets,
          MPI_INT, ghost_to_owner_recv_int_exchange,
          ghost_to_owner_recv_int_counts, ghost_to_owner_recv_int_offsets,
          MPI_INT, MPI_COMM_WORLD);
        MPI_Alltoallv(ghost_to_owner_send_double_exchange,
          ghost_to_owner_send_double_counts, ghost_to_owner_send_double_offsets,
          MPI_DOUBLE, ghost_to_owner_recv_double_exchange,
          ghost_to_owner_recv_double_counts, ghost_to_owner_recv_double_offsets,
          MPI_DOUBLE, MPI_COMM_WORLD);

        int received_ghost_caps_nm = total_owner_to_ghost_recv_caps > 0 ?
          total_owner_to_ghost_recv_caps : 1;
        int* received_ghost_cap_ids =
          malloc(received_ghost_caps_nm*sizeof(int));
        int* received_ghost_owner_procs =
          malloc(received_ghost_caps_nm*sizeof(int));
        assert(received_ghost_cap_ids);
        assert(received_ghost_owner_procs);
        int received_ghost_caps = 0;

        if (total_owner_to_ghost_recv_caps > 0) {
          fprintf(stderr,
            "DEBUG_RECV_GHOST_CAP_HEADERS pid %d/%d iter %d caps=",
            pid(), npe(), i);
          for(int p=0; p<npe(); p++) {
            int int_pos = owner_to_ghost_recv_int_offsets[p];
            int double_pos = owner_to_ghost_recv_double_offsets[p];
            for(int q=0; q<owner_to_ghost_recv_caps[p]; q++) {
              int cap_id, cap_type, nln, nle, nlt;
              double cap_es, cap_radius, circum_radius;
              unpack_owner_to_ghost_header(owner_to_ghost_recv_int_buffer,
                &int_pos, owner_to_ghost_recv_double_buffer, &double_pos,
                &cap_id, &cap_type, &nln, &nle, &nlt, &cap_es, &cap_radius,
                &circum_radius);
              received_ghost_cap_ids[received_ghost_caps] = cap_id;
              received_ghost_owner_procs[received_ghost_caps] = p;
              received_ghost_caps++;
              fprintf(stderr,
                " from=%d,cap=%d,type=%d,nln=%d,nle=%d",
                p, cap_id, cap_type, nln, nle);
              #if dimension > 2
                fprintf(stderr, ",nlt=%d", nlt);
              #endif
              fprintf(stderr,
                ",cap_es=%g,cap_radius=%g,circum_radius=%g",
                cap_es, cap_radius, circum_radius);
              coord recv_centroid;
              foreach_dimension()
                recv_centroid.x =
                  owner_to_ghost_recv_double_buffer[double_pos++];
              coord recv_pos_min = {HUGE, HUGE, HUGE};
              coord recv_pos_max = {-HUGE, -HUGE, -HUGE};
              for(int node_id=0; node_id<nln; node_id++) {
                foreach_dimension() {
                  double pos_comp =
                    owner_to_ghost_recv_double_buffer[double_pos++];
                  recv_pos_min.x = min(recv_pos_min.x, pos_comp);
                  recv_pos_max.x = max(recv_pos_max.x, pos_comp);
                }
              }
              double recv_lagforce_abs_sum = 0.;
              for(int node_id=0; node_id<nln; node_id++)
                foreach_dimension()
                  recv_lagforce_abs_sum +=
                    fabs(owner_to_ghost_recv_double_buffer[double_pos++]);
              fprintf(stderr,
                ",recv_centroid=(%g %g %g),recv_pos_min=(%g %g %g),recv_pos_max=(%g %g %g),recv_force_abs_sum=%g",
                recv_centroid.x, recv_centroid.y, recv_centroid.z,
                recv_pos_min.x, recv_pos_min.y, recv_pos_min.z,
                recv_pos_max.x, recv_pos_max.y, recv_pos_max.z,
                recv_lagforce_abs_sum);
            }
          }
          fprintf(stderr, "\n");
        }

        fprintf(stderr,
          "DEBUG_GHOST_CAP_LIFECYCLE_PLAN pid %d/%d iter %d receive_actions=",
          pid(), npe(), i);
        if (received_ghost_caps == 0)
          fprintf(stderr, " none");
        for(int r=0; r<received_ghost_caps; r++) {
          int local_index = -1;
          for(int cap=0; cap<NCAPS; cap++)
            if (CAPS(cap).isactive &&
              CAPS(cap).cap_id == received_ghost_cap_ids[r]) {
              local_index = cap;
              break;
            }
          fprintf(stderr, " cap=%d,owner=%d,local_index=%d,action=%s",
            received_ghost_cap_ids[r], received_ghost_owner_procs[r],
            local_index, local_index >= 0 ? "update_ghost" : "create_ghost");
        }
        fprintf(stderr, " destroy_candidates=");
        int ndestroy_candidates = 0;
        for(int cap=0; cap<NCAPS; cap++) {
          if (!CAPS(cap).isactive)
            continue;
          int owner_proc =
            find_capsule_owner_proc(&CAPS(cap), all_proc_min, all_proc_max);
          if (owner_proc == pid())
            continue;
          int received_here = false;
          for(int r=0; r<received_ghost_caps; r++)
            if (received_ghost_cap_ids[r] == CAPS(cap).cap_id) {
              received_here = true;
              break;
            }
          if (!received_here) {
            fprintf(stderr, " cap=%d,owner=%d,local_index=%d",
              CAPS(cap).cap_id, owner_proc, cap);
            ndestroy_candidates++;
          }
        }
        if (ndestroy_candidates == 0)
          fprintf(stderr, " none");
        fprintf(stderr, "\n");

        if (total_ghost_to_owner_recv_caps > 0) {
          fprintf(stderr,
            "DEBUG_RECV_OWNER_LAGVEL_PAYLOAD pid %d/%d iter %d caps=",
            pid(), npe(), i);
          for(int p=0; p<npe(); p++) {
            int int_pos = ghost_to_owner_recv_int_offsets[p];
            int double_pos = ghost_to_owner_recv_double_offsets[p];
            for(int q=0; q<ghost_to_owner_recv_caps[p]; q++) {
              int cap_id = ghost_to_owner_recv_int_buffer[int_pos++];
              int nln = ghost_to_owner_recv_int_buffer[int_pos++];
              coord recv_vel_min = {HUGE, HUGE, HUGE};
              coord recv_vel_max = {-HUGE, -HUGE, -HUGE};
              double recv_lagvel_abs_sum = 0.;
              for(int node_id=0; node_id<nln; node_id++) {
                foreach_dimension() {
                  double vel_comp =
                    ghost_to_owner_recv_double_buffer[double_pos++];
                  recv_vel_min.x = min(recv_vel_min.x, vel_comp);
                  recv_vel_max.x = max(recv_vel_max.x, vel_comp);
                  recv_lagvel_abs_sum += fabs(vel_comp);
                }
              }
              fprintf(stderr,
                " from=%d,cap=%d,nln=%d,recv_vel_min=(%g %g %g),recv_vel_max=(%g %g %g),recv_lagvel_abs_sum=%g",
                p, cap_id, nln,
                recv_vel_min.x, recv_vel_min.y, recv_vel_min.z,
                recv_vel_max.x, recv_vel_max.y, recv_vel_max.z,
                recv_lagvel_abs_sum);
            }
          }
          fprintf(stderr, "\n");
        }

        int local_lagvel_payload_ints[3] = {-1, -1, 0};
        double local_lagvel_payload_doubles[7] =
          {HUGE, HUGE, HUGE, -HUGE, -HUGE, -HUGE, 0.};
        if (total_ghost_to_owner_recv_caps > 0) {
          for(int p=0; p<npe(); p++) {
            if (ghost_to_owner_recv_caps[p] > 0) {
              int int_pos = ghost_to_owner_recv_int_offsets[p];
              int double_pos = ghost_to_owner_recv_double_offsets[p];
              int cap_id = ghost_to_owner_recv_int_buffer[int_pos++];
              int nln = ghost_to_owner_recv_int_buffer[int_pos++];
              coord recv_vel_min = {HUGE, HUGE, HUGE};
              coord recv_vel_max = {-HUGE, -HUGE, -HUGE};
              double recv_lagvel_abs_sum = 0.;
              for(int node_id=0; node_id<nln; node_id++) {
                foreach_dimension() {
                  double vel_comp =
                    ghost_to_owner_recv_double_buffer[double_pos++];
                  recv_vel_min.x = min(recv_vel_min.x, vel_comp);
                  recv_vel_max.x = max(recv_vel_max.x, vel_comp);
                  recv_lagvel_abs_sum += fabs(vel_comp);
                }
              }
              local_lagvel_payload_ints[0] = p;
              local_lagvel_payload_ints[1] = cap_id;
              local_lagvel_payload_ints[2] = nln;
              local_lagvel_payload_doubles[0] = recv_vel_min.x;
              local_lagvel_payload_doubles[1] = recv_vel_min.y;
              local_lagvel_payload_doubles[2] = recv_vel_min.z;
              local_lagvel_payload_doubles[3] = recv_vel_max.x;
              local_lagvel_payload_doubles[4] = recv_vel_max.y;
              local_lagvel_payload_doubles[5] = recv_vel_max.z;
              local_lagvel_payload_doubles[6] = recv_lagvel_abs_sum;
              break;
            }
          }
        }
        int* all_lagvel_payload_ints = NULL;
        double* all_lagvel_payload_doubles = NULL;
        if (pid() == 0) {
          all_lagvel_payload_ints = (int*)malloc(npe()*3*sizeof(int));
          all_lagvel_payload_doubles =
            (double*)malloc(npe()*7*sizeof(double));
          assert(all_lagvel_payload_ints);
          assert(all_lagvel_payload_doubles);
        }
        MPI_Gather(local_lagvel_payload_ints, 3, MPI_INT,
          all_lagvel_payload_ints, 3, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(local_lagvel_payload_doubles, 7, MPI_DOUBLE,
          all_lagvel_payload_doubles, 7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (pid() == 0) {
          int printed_lagvel_payload = false;
          for(int rank=0; rank<npe(); rank++) {
            int* row_i = all_lagvel_payload_ints + rank*3;
            double* row_d = all_lagvel_payload_doubles + rank*7;
            if (row_i[0] < 0)
              continue;
            if (!printed_lagvel_payload) {
              fprintf(stderr,
                "DEBUG_ALL_RECV_OWNER_LAGVEL_PAYLOAD iter %d", i);
              printed_lagvel_payload = true;
            }
            fprintf(stderr,
              " rank=%d,from=%d,cap=%d,nln=%d,recv_vel_min=(%g %g %g),recv_vel_max=(%g %g %g),recv_lagvel_abs_sum=%g",
              rank, row_i[0], row_i[1], row_i[2],
              row_d[0], row_d[1], row_d[2],
              row_d[3], row_d[4], row_d[5], row_d[6]);
          }
          if (printed_lagvel_payload)
            fprintf(stderr, "\n");
        }
        free(all_lagvel_payload_ints);
        free(all_lagvel_payload_doubles);

        #if !LAG_DISTRIBUTED_CAPSULES
          debug_sparse_owner_lagVel_from_recv_payloads(i,
            total_ghost_to_owner_recv_caps, true, true);
        #endif

        free(received_ghost_cap_ids);
        free(received_ghost_owner_procs);

        if (total_owner_to_ghost_recv_caps > 0 ||
          total_ghost_to_owner_recv_caps > 0) {
          fprintf(stderr,
            "DEBUG_MPI_PAYLOAD_EXCHANGE pid %d/%d iter %d owner_to_ghost_recv=",
            pid(), npe(), i);
          for(int p=0; p<npe(); p++) {
            if (owner_to_ghost_recv_caps[p] > 0) {
              int int_offset = owner_to_ghost_recv_int_offsets[p];
              int double_offset = owner_to_ghost_recv_double_offsets[p];
              fprintf(stderr,
                " %d:caps=%d,first_cap=%d,first_type=%d,nln=%d,nle=%d",
                p, owner_to_ghost_recv_caps[p],
                owner_to_ghost_recv_int_buffer[int_offset],
                owner_to_ghost_recv_int_buffer[int_offset + 1],
                owner_to_ghost_recv_int_buffer[int_offset + 2],
                owner_to_ghost_recv_int_buffer[int_offset + 3]);
              #if dimension > 2
                fprintf(stderr, ",nlt=%d",
                  owner_to_ghost_recv_int_buffer[int_offset + 4]);
              #endif
              fprintf(stderr, ",cap_es=%g,cap_radius=%g,circum_radius=%g",
                owner_to_ghost_recv_double_buffer[double_offset],
                owner_to_ghost_recv_double_buffer[double_offset + 1],
                owner_to_ghost_recv_double_buffer[double_offset + 2]);
            }
          }
          fprintf(stderr, " ghost_to_owner_recv=");
          for(int p=0; p<npe(); p++) {
            if (ghost_to_owner_recv_caps[p] > 0) {
              int int_offset = ghost_to_owner_recv_int_offsets[p];
              int double_offset = ghost_to_owner_recv_double_offsets[p];
              fprintf(stderr,
                " %d:caps=%d,first_cap=%d,nln=%d,first_lagVel=%g",
                p, ghost_to_owner_recv_caps[p],
                ghost_to_owner_recv_int_buffer[int_offset],
                ghost_to_owner_recv_int_buffer[int_offset + 1],
                ghost_to_owner_recv_double_buffer[double_offset]);
            }
          }
          fprintf(stderr, "\n");
        }
        int local_exchange_row_len = 4*npe();
        int* local_exchange_row = (int*)malloc(local_exchange_row_len*sizeof(int));
        int* all_exchange_rows = NULL;
        for(int p=0; p<npe(); p++) {
          local_exchange_row[p] = owner_to_ghost_recv_caps[p];
          local_exchange_row[npe() + p] = ghost_to_owner_recv_caps[p];
          local_exchange_row[2*npe() + p] = -1;
          local_exchange_row[3*npe() + p] = -1;
          if (owner_to_ghost_recv_caps[p] > 0)
            local_exchange_row[2*npe() + p] =
              owner_to_ghost_recv_int_buffer[owner_to_ghost_recv_int_offsets[p]];
          if (ghost_to_owner_recv_caps[p] > 0)
            local_exchange_row[3*npe() + p] =
              ghost_to_owner_recv_int_buffer[ghost_to_owner_recv_int_offsets[p]];
        }
        if (pid() == 0)
          all_exchange_rows = (int*)malloc(
            npe()*local_exchange_row_len*sizeof(int));
        MPI_Gather(local_exchange_row, local_exchange_row_len, MPI_INT,
          all_exchange_rows, local_exchange_row_len, MPI_INT, 0,
          MPI_COMM_WORLD);
        if (pid() == 0) {
          for(int rank=0; rank<npe(); rank++) {
            int* row = all_exchange_rows + rank*local_exchange_row_len;
            int owner_to_ghost_recv_total = 0;
            int ghost_to_owner_recv_total = 0;
            for(int p=0; p<npe(); p++) {
              owner_to_ghost_recv_total += row[p];
              ghost_to_owner_recv_total += row[npe() + p];
            }
            if (owner_to_ghost_recv_total > 0 ||
              ghost_to_owner_recv_total > 0) {
              fprintf(stderr,
                "DEBUG_ALL_MPI_PAYLOAD_EXCHANGE iter %d rank %d owner_to_ghost_recv=",
                i, rank);
              for(int p=0; p<npe(); p++)
                if (row[p] > 0)
                  fprintf(stderr, " %d:caps=%d,first_cap=%d",
                    p, row[p], row[2*npe() + p]);
              fprintf(stderr, " ghost_to_owner_recv=");
              for(int p=0; p<npe(); p++)
                if (row[npe() + p] > 0)
                  fprintf(stderr, " %d:caps=%d,first_cap=%d",
                    p, row[npe() + p], row[3*npe() + p]);
              fprintf(stderr,
                " total_owner_to_ghost_recv=%d total_ghost_to_owner_recv=%d\n",
                owner_to_ghost_recv_total, ghost_to_owner_recv_total);
            }
          }
        }
        free(local_exchange_row);
        free(all_exchange_rows);
        if (total_owner_to_ghost_caps > 0 || total_ghost_to_owner_caps > 0 ||
          total_owner_to_ghost_recv_caps > 0 || total_ghost_to_owner_recv_caps > 0) {
          size_t owner_to_ghost_send_bytes =
            total_owner_to_ghost_ints*sizeof(int)
            + total_owner_to_ghost_doubles*sizeof(double);
          size_t owner_to_ghost_recv_bytes =
            total_owner_to_ghost_recv_ints*sizeof(int)
            + total_owner_to_ghost_recv_doubles*sizeof(double);
          size_t ghost_to_owner_send_bytes =
            total_ghost_to_owner_ints*sizeof(int)
            + total_ghost_to_owner_doubles*sizeof(double);
          size_t ghost_to_owner_recv_bytes =
            total_ghost_to_owner_recv_ints*sizeof(int)
            + total_ghost_to_owner_recv_doubles*sizeof(double);
          fprintf(stderr,
            "DEBUG_LOCAL_PAYLOAD_BUFFERS pid %d/%d iter %d owner_to_ghost_send=(caps=%d,int=%d,double=%d,bytes=%zu) owner_to_ghost_recv=(caps=%d,int=%d,double=%d,bytes=%zu) ghost_to_owner_send=(caps=%d,int=%d,double=%d,bytes=%zu) ghost_to_owner_recv=(caps=%d,int=%d,double=%d,bytes=%zu)\n",
            pid(), npe(), i,
            total_owner_to_ghost_caps, total_owner_to_ghost_ints,
            total_owner_to_ghost_doubles, owner_to_ghost_send_bytes,
            total_owner_to_ghost_recv_caps, total_owner_to_ghost_recv_ints,
            total_owner_to_ghost_recv_doubles, owner_to_ghost_recv_bytes,
            total_ghost_to_owner_caps, total_ghost_to_owner_ints,
            total_ghost_to_owner_doubles, ghost_to_owner_send_bytes,
            total_ghost_to_owner_recv_caps, total_ghost_to_owner_recv_ints,
            total_ghost_to_owner_recv_doubles, ghost_to_owner_recv_bytes);
        }
        if (total_owner_to_ghost_caps > 0) {
          fprintf(stderr,
            "DEBUG_LOCAL_OWNER_TO_GHOST_SEND_COUNTS pid %d/%d iter %d sends=",
            pid(), npe(), i);
          for(int p=0; p<npe(); p++)
            if (owner_to_ghost_send_caps[p] > 0)
              fprintf(stderr, " %d:caps=%d,int=%d,double=%d",
                p, owner_to_ghost_send_caps[p],
                owner_to_ghost_send_int_counts[p],
                owner_to_ghost_send_double_counts[p]);
          fprintf(stderr, " total_caps=%d total_int=%d total_double=%d\n",
            total_owner_to_ghost_caps, total_owner_to_ghost_ints,
            total_owner_to_ghost_doubles);
        }
        if (total_ghost_to_owner_caps > 0) {
          fprintf(stderr,
            "DEBUG_LOCAL_GHOST_TO_OWNER_SEND_COUNTS pid %d/%d iter %d sends=",
            pid(), npe(), i);
          for(int p=0; p<npe(); p++)
            if (ghost_to_owner_send_caps[p] > 0)
              fprintf(stderr, " %d:caps=%d,int=%d,double=%d",
                p, ghost_to_owner_send_caps[p],
                ghost_to_owner_send_int_counts[p],
                ghost_to_owner_send_double_counts[p]);
          fprintf(stderr, " total_caps=%d total_int=%d total_double=%d\n",
            total_ghost_to_owner_caps, total_ghost_to_owner_ints,
            total_ghost_to_owner_doubles);
        }
        if (total_owner_to_ghost_recv_caps > 0) {
          fprintf(stderr,
            "DEBUG_LOCAL_OWNER_TO_GHOST_RECV_COUNTS pid %d/%d iter %d recvs=",
            pid(), npe(), i);
          for(int p=0; p<npe(); p++)
            if (owner_to_ghost_recv_caps[p] > 0)
              fprintf(stderr, " %d:caps=%d,int=%d,double=%d",
                p, owner_to_ghost_recv_caps[p],
                owner_to_ghost_recv_int_counts[p],
                owner_to_ghost_recv_double_counts[p]);
          fprintf(stderr, " total_caps=%d total_int=%d total_double=%d\n",
            total_owner_to_ghost_recv_caps, total_owner_to_ghost_recv_ints,
            total_owner_to_ghost_recv_doubles);
        }
        if (total_ghost_to_owner_recv_caps > 0) {
          fprintf(stderr,
            "DEBUG_LOCAL_GHOST_TO_OWNER_RECV_COUNTS pid %d/%d iter %d recvs=",
            pid(), npe(), i);
          for(int p=0; p<npe(); p++)
            if (ghost_to_owner_recv_caps[p] > 0)
              fprintf(stderr, " %d:caps=%d,int=%d,double=%d",
                p, ghost_to_owner_recv_caps[p],
                ghost_to_owner_recv_int_counts[p],
                ghost_to_owner_recv_double_counts[p]);
          fprintf(stderr, " total_caps=%d total_int=%d total_double=%d\n",
            total_ghost_to_owner_recv_caps, total_ghost_to_owner_recv_ints,
            total_ghost_to_owner_recv_doubles);
        }
        int local_count_row_len = 6*npe();
        int* local_count_row = (int*)calloc(local_count_row_len, sizeof(int));
        int* all_local_count_rows = NULL;
        for(int p=0; p<npe(); p++) {
          local_count_row[p] = owner_to_ghost_send_caps[p];
          local_count_row[npe() + p] = owner_to_ghost_send_int_counts[p];
          local_count_row[2*npe() + p] = owner_to_ghost_send_double_counts[p];
          local_count_row[3*npe() + p] = ghost_to_owner_send_caps[p];
          local_count_row[4*npe() + p] = ghost_to_owner_send_int_counts[p];
          local_count_row[5*npe() + p] = ghost_to_owner_send_double_counts[p];
        }
        if (pid() == 0)
          all_local_count_rows = (int*)malloc(npe()*local_count_row_len*sizeof(int));
        MPI_Gather(local_count_row, local_count_row_len, MPI_INT,
          all_local_count_rows, local_count_row_len, MPI_INT, 0, MPI_COMM_WORLD);
        if (pid() == 0) {
          fprintf(stderr, "DEBUG_ALL_LOCAL_SEND_COUNTS iter %d npe %d\n",
            i, npe());
          for(int rank=0; rank<npe(); rank++) {
            int* row = all_local_count_rows + rank*local_count_row_len;
            int owner_caps = 0, owner_ints = 0, owner_doubles = 0;
            int ghost_caps = 0, ghost_ints = 0, ghost_doubles = 0;
            for(int p=0; p<npe(); p++) {
              owner_caps += row[p];
              owner_ints += row[npe() + p];
              owner_doubles += row[2*npe() + p];
              ghost_caps += row[3*npe() + p];
              ghost_ints += row[4*npe() + p];
              ghost_doubles += row[5*npe() + p];
            }
            if (owner_caps > 0) {
              fprintf(stderr,
                "DEBUG_ALL_LOCAL_OWNER_TO_GHOST_SEND_COUNTS iter %d rank %d sends=",
                i, rank);
              for(int p=0; p<npe(); p++)
                if (row[p] > 0)
                  fprintf(stderr, " %d:caps=%d,int=%d,double=%d",
                    p, row[p], row[npe() + p], row[2*npe() + p]);
              fprintf(stderr, " total_caps=%d total_int=%d total_double=%d\n",
                owner_caps, owner_ints, owner_doubles);
            }
            if (ghost_caps > 0) {
              fprintf(stderr,
                "DEBUG_ALL_LOCAL_GHOST_TO_OWNER_SEND_COUNTS iter %d rank %d sends=",
                i, rank);
              for(int p=0; p<npe(); p++)
                if (row[3*npe() + p] > 0)
                  fprintf(stderr, " %d:caps=%d,int=%d,double=%d",
                    p, row[3*npe() + p],
                    row[4*npe() + p], row[5*npe() + p]);
              fprintf(stderr, " total_caps=%d total_int=%d total_double=%d\n",
                ghost_caps, ghost_ints, ghost_doubles);
            }
          }
        }
        free(local_count_row);
        free(all_local_count_rows);
        if (pid() == 0) {
          fprintf(stderr, "DEBUG_AABB_TABLE iter %d npe %d\n", i, npe());
          for(int p=0; p<npe(); p++)
            fprintf(stderr,
              "DEBUG_AABB_TABLE iter %d proc %d proc_min=(%g %g %g) proc_max=(%g %g %g)\n",
              i, p,
              all_proc_min[p].x, all_proc_min[p].y, all_proc_min[p].z,
              all_proc_max[p].x, all_proc_max[p].y, all_proc_max[p].z);
          for(int cap=0; cap<NCAPS; cap++) {
            if (CAPS(cap).isactive) {
              int nintersections = 0;
              int owner_proc = find_capsule_owner_proc(&CAPS(cap),
                all_proc_min, all_proc_max);
              fprintf(stderr,
                "DEBUG_CAP_OWNER iter %d cap %d centroid=(%g %g %g) owner=%d\n",
                i, cap,
                CAPS(cap).centroid.x, CAPS(cap).centroid.y, CAPS(cap).centroid.z,
                owner_proc);
              fprintf(stderr,
                "DEBUG_CAP_PROC_TABLE iter %d cap %d centroid=(%g %g %g) circum_radius=%g intersecting_procs=",
                i, cap,
                CAPS(cap).centroid.x, CAPS(cap).centroid.y, CAPS(cap).centroid.z,
                CAPS(cap).circum_radius);
              for(int p=0; p<npe(); p++) {
                bool intersects_proc = lagmesh_bounding_sphere_intersects_box(
                  &CAPS(cap), all_proc_min[p], all_proc_max[p]);
                if (intersects_proc) {
                  fprintf(stderr, " %d", p);
                  nintersections++;
                }
              }
              fprintf(stderr, " nintersections=%d\n", nintersections);
              int nsends = 0;
              fprintf(stderr,
                "DEBUG_CAP_SEND_PLAN iter %d cap %d owner=%d send_to=",
                i, cap, owner_proc);
              if (owner_proc >= 0) {
                for(int p=0; p<npe(); p++) {
                  bool intersects_proc = lagmesh_bounding_sphere_intersects_box(
                    &CAPS(cap), all_proc_min[p], all_proc_max[p]);
                  if (intersects_proc && p != owner_proc) {
                    fprintf(stderr, " %d", p);
                    nsends++;
                  }
                }
              }
              fprintf(stderr, " nsends=%d\n", nsends);
            }
          }
          int* owner_send_counts = (int*)calloc(npe()*npe(), sizeof(int));
          int* owner_to_ghost_int_counts = (int*)calloc(npe()*npe(), sizeof(int));
          int* owner_to_ghost_double_counts = (int*)calloc(npe()*npe(), sizeof(int));
          int* ghost_to_owner_int_counts = (int*)calloc(npe()*npe(), sizeof(int));
          int* ghost_to_owner_double_counts = (int*)calloc(npe()*npe(), sizeof(int));
          for(int cap=0; cap<NCAPS; cap++) {
            if (CAPS(cap).isactive) {
              int owner_proc = find_capsule_owner_proc(&CAPS(cap),
                all_proc_min, all_proc_max);
              if (owner_proc >= 0) {
                int owner_to_ghost_nints = estimate_owner_to_ghost_nints(&CAPS(cap));
                int owner_to_ghost_ndoubles = estimate_owner_to_ghost_ndoubles(&CAPS(cap));
                int ghost_to_owner_nints = estimate_ghost_to_owner_nints(&CAPS(cap));
                int ghost_to_owner_ndoubles = estimate_ghost_to_owner_ndoubles(&CAPS(cap));
                for(int p=0; p<npe(); p++) {
                  bool intersects_proc = lagmesh_bounding_sphere_intersects_box(
                    &CAPS(cap), all_proc_min[p], all_proc_max[p]);
                  if (intersects_proc && p != owner_proc) {
                    owner_send_counts[owner_proc*npe() + p]++;
                    owner_to_ghost_int_counts[owner_proc*npe() + p] += owner_to_ghost_nints;
                    owner_to_ghost_double_counts[owner_proc*npe() + p] += owner_to_ghost_ndoubles;
                    ghost_to_owner_int_counts[p*npe() + owner_proc] += ghost_to_owner_nints;
                    ghost_to_owner_double_counts[p*npe() + owner_proc] += ghost_to_owner_ndoubles;
                  }
                }
              }
            }
          }
          fprintf(stderr, "DEBUG_OWNER_SEND_COUNTS iter %d npe %d\n", i, npe());
          for(int owner=0; owner<npe(); owner++) {
            int total_owner_sends = 0;
            for(int dest=0; dest<npe(); dest++)
              total_owner_sends += owner_send_counts[owner*npe() + dest];
            if (total_owner_sends > 0) {
              fprintf(stderr,
                "DEBUG_OWNER_SEND_COUNTS iter %d owner %d send_counts=",
                i, owner);
              for(int dest=0; dest<npe(); dest++)
                if (owner_send_counts[owner*npe() + dest] > 0)
                  fprintf(stderr, " %d:%d", dest,
                    owner_send_counts[owner*npe() + dest]);
              fprintf(stderr, " total=%d\n", total_owner_sends);
            }
          }
          fprintf(stderr, "DEBUG_OWNER_TO_GHOST_BUFFER_EST iter %d npe %d\n",
            i, npe());
          for(int owner=0; owner<npe(); owner++) {
            int total_ints = 0;
            int total_doubles = 0;
            for(int dest=0; dest<npe(); dest++) {
              total_ints += owner_to_ghost_int_counts[owner*npe() + dest];
              total_doubles += owner_to_ghost_double_counts[owner*npe() + dest];
            }
            if (total_ints > 0 || total_doubles > 0) {
              size_t total_bytes = total_ints*sizeof(int)
                + total_doubles*sizeof(double);
              fprintf(stderr,
                "DEBUG_OWNER_TO_GHOST_BUFFER_EST iter %d owner %d payloads=",
                i, owner);
              for(int dest=0; dest<npe(); dest++) {
                int nints = owner_to_ghost_int_counts[owner*npe() + dest];
                int ndoubles = owner_to_ghost_double_counts[owner*npe() + dest];
                if (nints > 0 || ndoubles > 0) {
                  size_t bytes = nints*sizeof(int) + ndoubles*sizeof(double);
                  fprintf(stderr, " %d:int=%d,double=%d,bytes=%zu",
                    dest, nints, ndoubles, bytes);
                }
              }
              fprintf(stderr, " total_int=%d total_double=%d total_bytes=%zu\n",
                total_ints, total_doubles, total_bytes);
            }
          }
          fprintf(stderr, "DEBUG_GHOST_TO_OWNER_BUFFER_EST iter %d npe %d\n",
            i, npe());
          for(int ghost=0; ghost<npe(); ghost++) {
            int total_ints = 0;
            int total_doubles = 0;
            for(int owner=0; owner<npe(); owner++) {
              total_ints += ghost_to_owner_int_counts[ghost*npe() + owner];
              total_doubles += ghost_to_owner_double_counts[ghost*npe() + owner];
            }
            if (total_ints > 0 || total_doubles > 0) {
              size_t total_bytes = total_ints*sizeof(int)
                + total_doubles*sizeof(double);
              fprintf(stderr,
                "DEBUG_GHOST_TO_OWNER_BUFFER_EST iter %d ghost %d payloads=",
                i, ghost);
              for(int owner=0; owner<npe(); owner++) {
                int nints = ghost_to_owner_int_counts[ghost*npe() + owner];
                int ndoubles = ghost_to_owner_double_counts[ghost*npe() + owner];
                if (nints > 0 || ndoubles > 0) {
                  size_t bytes = nints*sizeof(int) + ndoubles*sizeof(double);
                  fprintf(stderr, " %d:int=%d,double=%d,bytes=%zu",
                    owner, nints, ndoubles, bytes);
                }
              }
              fprintf(stderr, " total_int=%d total_double=%d total_bytes=%zu\n",
                total_ints, total_doubles, total_bytes);
            }
          }
          free(owner_send_counts);
          free(owner_to_ghost_int_counts);
          free(owner_to_ghost_double_counts);
          free(ghost_to_owner_int_counts);
          free(ghost_to_owner_double_counts);
          fprintf(stderr, "DEBUG_PROC_CAP_LIST iter %d npe %d\n", i, npe());
          for(int p=0; p<npe(); p++) {
            fprintf(stderr, "DEBUG_PROC_CAP_LIST iter %d proc %d ncaps=%d offset=%d cap_ids=",
              i, p, ncaps_for_proc[p], proc_cap_offsets[p]);
            for(int q=0; q<ncaps_for_proc[p]; q++)
              fprintf(stderr, " %d", proc_cap_ids[proc_cap_offsets[p] + q]);
            fprintf(stderr, "\n");
          }
          fprintf(stderr, "DEBUG_PROC_CAP_LIST iter %d total_routes=%d\n",
            i, proc_cap_offsets[npe()]);
        }
      #endif
    }
  #endif

}

#endif
