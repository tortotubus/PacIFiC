#pragma once

#include "ibm/IBMeshManager.h"
#include "ibm/IBKernels.h"

#define CAPSULE_IBM_MESH_ID_INVALID (-1)

typedef struct CapsuleIBMProxy {
  CapsuleMesh* lag;
  int mesh_id;
  int node_count;
  int owner_pid;
  bool midpoint_output;
  double midpoint_dt;
  bool active;
} CapsuleIBMProxy;

typedef void (*CapsuleIBMForceAssembler)(CapsuleMesh* mesh);

static CapsuleIBMProxy* capsule_ibm_proxies = NULL;
static int capsule_ibm_proxy_capacity = 0;
static bool capsule_ibm_proxy_table_initialized = false;
static CapsuleIBMForceAssembler capsule_ibm_force_assembler = NULL;
static IBvector capsule_viscosity_area_normal;
static bool capsule_ibm_viscosity_fields_initialized = false;

static void capsule_ibm_init_viscosity_fields(void);

static void capsule_ibm_set_force_assembler(CapsuleIBMForceAssembler assembler) {
  init_ibsolver();
  capsule_ibm_init_viscosity_fields();
  capsule_ibm_force_assembler = assembler;
}

static void capsule_ibm_init_viscosity_fields(void) {
  if (capsule_ibm_viscosity_fields_initialized)
    return;

  new_ibvector(capsule_viscosity_area_normal);
  foreach_dimension()
    ibnodump(capsule_viscosity_area_normal.x) = true;
  capsule_ibm_viscosity_fields_initialized = true;
}

static void capsule_ibm_zero_lag_forces(CapsuleMesh* lag) {
  if (!lag)
    return;

  for (int i = 0; i < lag->nln; i++)
    foreach_dimension()
      lag->nodes[i].lagForce.x = 0.;
}

static void capsule_ibm_assemble_lag_forces(CapsuleMesh* lag) {
  capsule_ibm_zero_lag_forces(lag);
  if (capsule_ibm_force_assembler)
    capsule_ibm_force_assembler(lag);
}

static void capsule_ibm_clear_proxy_slot(int i) {
  capsule_ibm_proxies[i].lag = NULL;
  capsule_ibm_proxies[i].mesh_id = CAPSULE_IBM_MESH_ID_INVALID;
  capsule_ibm_proxies[i].node_count = 0;
  capsule_ibm_proxies[i].owner_pid = 0;
  capsule_ibm_proxies[i].midpoint_output = false;
  capsule_ibm_proxies[i].midpoint_dt = 0.;
  capsule_ibm_proxies[i].active = false;
}

static void capsule_ibm_proxy_table_reserve(int capacity) {
  if (capacity <= capsule_ibm_proxy_capacity)
    return;

  capsule_ibm_proxies =
    (CapsuleIBMProxy*) realloc(capsule_ibm_proxies,
                               (size_t) capacity*sizeof(CapsuleIBMProxy));
  assert(capsule_ibm_proxies);
  for (int i = capsule_ibm_proxy_capacity; i < capacity; i++)
    capsule_ibm_clear_proxy_slot(i);
  capsule_ibm_proxy_capacity = capacity;
}

static void capsule_ibm_proxy_table_init(void) {
  if (capsule_ibm_proxy_table_initialized)
    return;

  capsule_ibm_proxy_table_reserve((allCaps.nbcaps > NCAPS ? allCaps.nbcaps : NCAPS) + 4);
  capsule_ibm_proxy_table_initialized = true;
}

static void capsule_ibm_ensure_manager(void) {
  init_ibsolver();
  capsule_ibm_init_viscosity_fields();
  ibmeshmanager_init(0);
  capsule_ibm_proxy_table_init();
}

static IBscalar* capsule_ibm_model_output_fields(void* ctx) {
  (void) ctx;

  IBscalar* fields = NULL;
  foreach_dimension() {
    fields = iblist_add(fields, npos.x);
    fields = iblist_add(fields, nvel.x);
    fields = iblist_add(fields, nforce.x);
    fields = iblist_add(fields, capsule_viscosity_area_normal.x);
  }
  fields = iblist_add(fields, nweight);
  return fields;
}

static IBscalar* capsule_ibm_model_input_fields(void* ctx) {
  (void) ctx;

  IBscalar* fields = NULL;
  foreach_dimension()
    fields = iblist_add(fields, nvel.x);
  return fields;
}

static int capsule_ibm_find_proxy(CapsuleMesh* lag) {
  capsule_ibm_proxy_table_init();
  for (int i = 0; i < capsule_ibm_proxy_capacity; i++)
    if (capsule_ibm_proxies[i].active && capsule_ibm_proxies[i].lag == lag)
      return i;
  return -1;
}

static int capsule_ibm_first_free_proxy(void) {
  capsule_ibm_proxy_table_init();
  for (int i = 0; i < capsule_ibm_proxy_capacity; i++)
    if (!capsule_ibm_proxies[i].active)
      return i;
  int first_new = capsule_ibm_proxy_capacity;
  capsule_ibm_proxy_table_reserve(2*capsule_ibm_proxy_capacity);
  return first_new;
}

static int capsule_ibm_requested_node_depth(IBMesh* ibmesh) {
#if TREE
  return ibmesh && ibmesh->depth > 0 ? ibmesh->depth : grid->maxdepth;
#else
  return depth();
#endif
}

static void capsule_ibm_copy_lag_to_node(CapsuleMesh* lag, IBMesh* ibmesh,
                                         int node_id, IBNode* node) {
  double node_area = compute_node_area(lag, node_id);
  if (node_area <= SEPS)
    node_area = 1.;

  foreach_dimension() {
    ibval(npos.x) = lag->nodes[node_id].pos.x;
    ibval(nvel.x) = lag->nodes[node_id].lagVel.x;
    ibval(nforce.x) = lag->nodes[node_id].lagForce.x/node_area;
  }
  ibval(nweight) = node_area;
  node->depth = capsule_ibm_requested_node_depth(ibmesh);
}

static coord capsule_ibm_node_viscosity_area_normal(CapsuleMesh* lag,
                                                    int node_id) {
  coord area_normal = {0};

#if dimension == 3
  if (!lag || node_id < 0 || node_id >= lag->nln)
    return area_normal;

  if (!lag->updated_normals)
    comp_normals(lag);

  CapNode* cap_node = &lag->nodes[node_id];
  for (int j = 0; j < cap_node->nb_triangles; j++) {
    int tid = cap_node->triangle_ids[j];
    if (tid < 0 || tid >= lag->nlt)
      continue;

    foreach_dimension()
      area_normal.x +=
        LAG_TRI_STATE(lag, tid)->normal.x*LAG_TRI_STATE(lag, tid)->area/3.;
  }
#else
  (void) lag;
  (void) node_id;
#endif

  return area_normal;
}

static void capsule_ibm_copy_lag_to_node_keyed(CapsuleIBMProxy* proxy,
                                               IBMesh* ibmesh,
                                               int node_id,
                                               IBNode* node) {
  if (!proxy || !proxy->lag || !ibmesh || !node ||
      node_id < 0 || node_id >= proxy->lag->nln)
    return;

  CapsuleMesh* lag = proxy->lag;
  double node_area = compute_node_area(lag, node_id);
  if (node_area <= SEPS)
    node_area = 1.;

  node->mesh_gid = ibmesh->gid;
  node->node_lid = node_id;
  node->depth = capsule_ibm_requested_node_depth(ibmesh);

  foreach_dimension() {
    ibval(npos.x) = lag->nodes[node_id].pos.x;
    if (proxy->midpoint_output)
      ibval(npos.x) += 0.5*proxy->midpoint_dt*lag->nodes[node_id].lagVel.x;
    ibval(nvel.x) = lag->nodes[node_id].lagVel.x;
    ibval(nforce.x) = lag->nodes[node_id].lagForce.x/node_area;
  }
  ibval(nweight) = node_area;

  coord visc_area_normal =
    capsule_ibm_node_viscosity_area_normal(lag, node_id);
  foreach_dimension()
    ibval(capsule_viscosity_area_normal.x) = visc_area_normal.x;
}

static void capsule_ibm_copy_node_velocity_to_lag(CapsuleMesh* lag, int node_id, IBNode* node) {
  if (!lag || !node || node_id < 0 || node_id >= lag->nln)
    return;

  foreach_dimension()
    lag->nodes[node_id].lagVel.x = ibval(nvel.x);
}

static void capsule_ibm_copy_lag_positions_to_nodes(CapsuleIBMProxy* proxy,
                                                    IBMesh* ibmesh) {
  if (!proxy || !proxy->lag || !ibmesh)
    return;

  int nn = proxy->lag->nln;
  if ((int) ibmesh->nodes.size < nn)
    nn = (int) ibmesh->nodes.size;

  for (int i = 0; i < nn; i++) {
    IBNode* node = ibmesh->nodes.ptrs[i];
    foreach_dimension()
      ibval(npos.x) = proxy->lag->nodes[i].pos.x;
    node->depth = capsule_ibm_requested_node_depth(ibmesh);
  }
}

static coord capsule_ibm_triangle_area_normal(CapsuleMesh* mesh, int tid) {
  coord area_normal = {0};

#if dimension == 3
  if (!mesh || tid < 0 || tid >= mesh->nlt)
    return area_normal;

  int* node_ids = LAG_TRI_TOPO(mesh, tid)->node_ids;
  coord e10 = {0}, e20 = {0};
  foreach_dimension() {
    e10.x = GENERAL_1DIST(mesh->nodes[node_ids[1]].pos.x,
                          mesh->nodes[node_ids[0]].pos.x);
    e20.x = GENERAL_1DIST(mesh->nodes[node_ids[2]].pos.x,
                          mesh->nodes[node_ids[0]].pos.x);
  }

  area_normal.x = 0.5*(e10.y*e20.z - e10.z*e20.y);
  area_normal.y = 0.5*(e10.z*e20.x - e10.x*e20.z);
  area_normal.z = 0.5*(e10.x*e20.y - e10.y*e20.x);
#endif

  return area_normal;
}

static void capsule_ibm_zero_viscosity_area_normals(IBMesh* ibmesh) {
  if (!ibmesh)
    return;

  for (size_t i = 0; i < ibmesh->nodes.size; i++) {
    IBNode* node = ibmesh->nodes.ptrs[i];
    foreach_dimension()
      ibval(capsule_viscosity_area_normal.x) = 0.;
  }
}

static void capsule_ibm_compute_viscosity_area_normals(CapsuleIBMProxy* proxy,
                                                       IBMesh* ibmesh) {
  if (!proxy || !proxy->lag || !ibmesh)
    return;

  capsule_ibm_zero_viscosity_area_normals(ibmesh);

#if dimension == 3
  CapsuleMesh* lag = proxy->lag;
  comp_normals(lag);
  for (int tid = 0; tid < lag->nlt; tid++) {
    int* node_ids = LAG_TRI_TOPO(lag, tid)->node_ids;
    coord area_normal = {0};
    foreach_dimension()
      area_normal.x = LAG_TRI_STATE(lag, tid)->normal.x*
                      LAG_TRI_STATE(lag, tid)->area;

    for (int j = 0; j < 3; j++) {
      int node_id = node_ids[j];
      if (node_id < 0 || (size_t) node_id >= ibmesh->nodes.size)
        continue;

      IBNode* node = ibmesh->nodes.ptrs[(size_t) node_id];
      foreach_dimension()
        ibval(capsule_viscosity_area_normal.x) += area_normal.x/3.;
    }
  }
#endif
}

static void capsule_ibm_sync_viscosity_area_normals(CapsuleIBMProxy* proxy,
                                                    IBMesh* ibmesh) {
  if (!ibmesh)
    return;

  int owner = 0;
#if _MPI
  owner = ibmesh->pid >= 0 ? ibmesh->pid : 0;
#endif

  if (owner == pid() && proxy && proxy->lag)
    capsule_ibm_compute_viscosity_area_normals(proxy, ibmesh);

#if _MPI
  {
    IBscalar* fields = NULL;
    foreach_dimension()
      fields = iblist_add(fields, capsule_viscosity_area_normal.x);
    ibmeshmanager_broadcast_fields_from_mesh_owner(ibmesh, fields);
    free(fields);
  }
#endif
}

static CapsuleIBMProxy* capsule_ibm_proxy_for_lag(CapsuleMesh* lag) {
  if (!lag)
    return NULL;

  capsule_ibm_proxy_table_init();

  int proxy_id = capsule_ibm_find_proxy(lag);
  if (proxy_id < 0) {
    proxy_id = capsule_ibm_first_free_proxy();
    if (proxy_id < 0) {
      fprintf(stderr, "Error: no free capsule IBM proxy slots.\n");
      assert(false);
    }

    CapsuleIBMProxy* proxy = &capsule_ibm_proxies[proxy_id];
    proxy->lag = lag;
    proxy->mesh_id = CAPSULE_IBM_MESH_ID_INVALID;
    proxy->node_count = lag->nln;
    proxy->owner_pid = lag->owner_pid;
    proxy->midpoint_output = false;
    proxy->midpoint_dt = 0.;
    proxy->active = true;
  }

  return &capsule_ibm_proxies[proxy_id];
}

static int capsule_ibm_model_node_count(void* ctx) {
  CapsuleIBMProxy* proxy = (CapsuleIBMProxy*) ctx;
  if (!proxy)
    return 0;
  return proxy->lag ? proxy->lag->nln : proxy->node_count;
}

static int capsule_ibm_model_sync(void* ctx, void* ibmesh_ptr) {
  CapsuleIBMProxy* proxy = (CapsuleIBMProxy*) ctx;
  IBMesh* ibmesh = (IBMesh*) ibmesh_ptr;
  if (!proxy || !ibmesh)
    return -1;

  ibmesh->pid = proxy->owner_pid;
  proxy->midpoint_output = false;
  proxy->midpoint_dt = 0.;

  if (!proxy->lag) {
    for (size_t i = 0; i < ibmesh->nodes.size; i++)
      ibmesh->nodes.ptrs[i]->depth = capsule_ibm_requested_node_depth(ibmesh);
    proxy->node_count = (int) ibmesh->nodes.size;
    return 0;
  }

  int nn = proxy->lag->nln;
  if ((int) ibmesh->nodes.size < nn)
    nn = (int) ibmesh->nodes.size;

  for (int i = 0; i < nn; i++) {
    IBNode* node = ibmesh->nodes.ptrs[i];
    int node_id = node->node_lid >= 0 ? node->node_lid : i;
    capsule_ibm_copy_lag_to_node_keyed(proxy, ibmesh, node_id, node);
  }

  if (!ibmesh->sparse_managed)
    capsule_ibm_compute_viscosity_area_normals(proxy, ibmesh);
  proxy->node_count = proxy->lag->nln;
  return 0;
}

static int capsule_ibm_model_sync_node(void* ctx,
                                       void* ibmesh_ptr,
                                       int node_lid,
                                       void* node_ptr) {
  CapsuleIBMProxy* proxy = (CapsuleIBMProxy*) ctx;
  IBMesh* ibmesh = (IBMesh*) ibmesh_ptr;
  IBNode* node = (IBNode*) node_ptr;
  if (!proxy || !proxy->lag || !ibmesh || !node)
    return -1;

  capsule_ibm_copy_lag_to_node_keyed(proxy, ibmesh, node_lid, node);
  return 0;
}

static int capsule_ibm_model_input_node(void* ctx,
                                        void* ibmesh_ptr,
                                        int node_lid,
                                        void* node_ptr) {
  (void) ibmesh_ptr;
  CapsuleIBMProxy* proxy = (CapsuleIBMProxy*) ctx;
  IBNode* node = (IBNode*) node_ptr;
  if (!proxy || !proxy->lag || !node ||
      node_lid < 0 || node_lid >= proxy->lag->nln)
    return -1;

  capsule_ibm_copy_node_velocity_to_lag(proxy->lag, node_lid, node);
  return 0;
}

static int capsule_ibm_model_midpoint(void* ctx, void* ibmesh_ptr, double local_dt) {
  CapsuleIBMProxy* proxy = (CapsuleIBMProxy*) ctx;
  IBMesh* ibmesh = (IBMesh*) ibmesh_ptr;
  if (!proxy || !proxy->lag || !ibmesh)
    return -1;

  if (!ibmesh->sparse_managed) {
    for (size_t i = 0; i < ibmesh->nodes.size; i++) {
      IBNode* node = ibmesh->nodes.ptrs[i];
      int node_id = node->node_lid >= 0 ? node->node_lid : (int) i;
      capsule_ibm_copy_node_velocity_to_lag(proxy->lag, node_id, node);
    }
  }

  proxy->midpoint_output = true;
  proxy->midpoint_dt = local_dt;
  capsule_ibm_assemble_lag_forces(proxy->lag);

  for (size_t i = 0; i < ibmesh->nodes.size; i++) {
    IBNode* node = ibmesh->nodes.ptrs[i];
    int node_id = node->node_lid >= 0 ? node->node_lid : (int) i;
    capsule_ibm_copy_lag_to_node_keyed(proxy, ibmesh, node_id, node);
  }

  if (!ibmesh->sparse_managed)
    capsule_ibm_compute_viscosity_area_normals(proxy, ibmesh);
  return 0;
}

static int capsule_ibm_model_advance(void* ctx, void* ibmesh_ptr, double local_dt) {
  CapsuleIBMProxy* proxy = (CapsuleIBMProxy*) ctx;
  IBMesh* ibmesh = (IBMesh*) ibmesh_ptr;
  if (!proxy || !proxy->lag || !ibmesh)
    return -1;

  if (!ibmesh->sparse_managed) {
    for (size_t i = 0; i < ibmesh->nodes.size; i++) {
      IBNode* node = ibmesh->nodes.ptrs[i];
      int node_id = node->node_lid >= 0 ? node->node_lid : (int) i;
      capsule_ibm_copy_node_velocity_to_lag(proxy->lag, node_id, node);
    }
  }

  proxy->midpoint_output = false;
  proxy->midpoint_dt = 0.;
  for (int i = 0; i < proxy->lag->nln; i++) {
    foreach_dimension()
      proxy->lag->nodes[i].pos.x += local_dt*proxy->lag->nodes[i].lagVel.x;
  }

  correct_lag_pos(proxy->lag);
#if CONSERVE_VOLUME
  enforce_optimal_volume_conservation(proxy->lag);
#endif
  comp_centroid(proxy->lag);
  comp_volume(proxy->lag);
  comp_capsule_geodynamics(proxy->lag);
  capsule_ibm_model_sync(ctx, ibmesh_ptr);
  return 0;
}

static void capsule_ibm_model_destroy(void* ctx) {
  CapsuleIBMProxy* proxy = (CapsuleIBMProxy*) ctx;
  if (!proxy)
    return;

  proxy->lag = NULL;
  proxy->mesh_id = CAPSULE_IBM_MESH_ID_INVALID;
  proxy->node_count = 0;
  proxy->owner_pid = 0;
  proxy->active = false;
}

static IBMeshModel capsule_ibm_model_from_proxy(CapsuleIBMProxy* proxy) {
  IBMeshModel model = ibmeshmodel_velocity_coupled_init();
  model.ctx = proxy;
  model.velocity_ops->node_count = capsule_ibm_model_node_count;
  model.velocity_ops->sync = capsule_ibm_model_sync;
  model.velocity_ops->sync_node = capsule_ibm_model_sync_node;
  model.velocity_ops->input_node = capsule_ibm_model_input_node;
  model.velocity_ops->midpoint = capsule_ibm_model_midpoint;
  model.velocity_ops->advance = capsule_ibm_model_advance;
  model.velocity_ops->input_fields = capsule_ibm_model_input_fields;
  model.velocity_ops->output_fields = capsule_ibm_model_output_fields;
  /* The MPI mesh manager destroys temporary active models when replacing them
     with passive replicas on non-owner ranks. Proxy lifetime is therefore
     managed explicitly by capsule_ibm_unregister_mesh()/cleanup(). */
  model.velocity_ops->destroy = NULL;
  return model;
}

/** Public adapter interface: returns a model suitable for ibmeshmanager_set_model(). */
static IBMeshModel capsule_ibm_mesh_model(CapsuleMesh* lag) {
  CapsuleIBMProxy* proxy = capsule_ibm_proxy_for_lag(lag);
  if (!proxy)
    return ibmeshmodel_init();
  return capsule_ibm_model_from_proxy(proxy);
}

static IBMeshModel capsule_ibm_model(CapsuleMesh* lag) {
  return capsule_ibm_mesh_model(lag);
}

static int capsule_ibm_find_mesh_id_for_proxy(CapsuleIBMProxy* proxy) {
  if (!proxy)
    return CAPSULE_IBM_MESH_ID_INVALID;

  if (proxy->mesh_id >= 0 && proxy->mesh_id < ibmm.nm &&
      ibmm.meshes[proxy->mesh_id].pid == proxy->owner_pid)
    return proxy->mesh_id;

  /* Non-owner ranks replace capsule models with passive models, so model.ctx
     is only a reliable reverse lookup on the owning rank. */
  for (int i = 0; i < ibmm.nm; i++) {
    if (ibmm.meshes[i].model.ctx == proxy) {
      proxy->mesh_id = i;
      return i;
    }
  }

  return CAPSULE_IBM_MESH_ID_INVALID;
}

static int capsule_ibm_mesh_id(CapsuleMesh* lag) {
  if (!lag)
    return CAPSULE_IBM_MESH_ID_INVALID;

  int proxy_id = capsule_ibm_find_proxy(lag);
  if (proxy_id < 0)
    return CAPSULE_IBM_MESH_ID_INVALID;

  return capsule_ibm_find_mesh_id_for_proxy(&capsule_ibm_proxies[proxy_id]);
}

static CapsuleIBMProxy* capsule_ibm_proxy_from_mesh_id(int mesh_id) {
  if (mesh_id < 0 || mesh_id >= ibmm.nm)
    return NULL;

  capsule_ibm_proxy_table_init();
  for (int i = 0; i < capsule_ibm_proxy_capacity; i++) {
    CapsuleIBMProxy* proxy = &capsule_ibm_proxies[i];
    if (!proxy->active)
      continue;
    if (proxy->mesh_id == mesh_id)
      return proxy;
  }
  return NULL;
}

static CapsuleMesh* capsule_ibm_lag_from_mesh_id(int mesh_id) {
  CapsuleIBMProxy* proxy = capsule_ibm_proxy_from_mesh_id(mesh_id);
  return proxy ? proxy->lag : NULL;
}

static void capsule_ibm_install_velocity_model(CapsuleIBMProxy* proxy, int mesh_id) {
  proxy->mesh_id = mesh_id;
  IBMesh* ibmesh = &ibmm.meshes[mesh_id];
  if (ibmesh->model.type != IB_MODEL_NONE) {
    IBMeshModel old_model = ibmesh->model;
    if (old_model.type == IB_MODEL_VELOCITY_COUPLED && old_model.velocity_ops)
      old_model.velocity_ops->destroy = NULL;
    if (old_model.type == IB_MODEL_FORCE_COUPLED && old_model.force_ops)
      old_model.force_ops->destroy = NULL;
    ibmeshmodel_destroy(&old_model);
    ibmesh->model = ibmeshmodel_init();
  }

  IBMeshModel model = capsule_ibm_model_from_proxy(proxy);
  int mesh_gid = proxy->lag ? proxy->lag->cap_id : mesh_id;
  ibmeshmanager_set_model_sparse(mesh_id, mesh_gid, model, proxy->owner_pid);
}

static CapsuleIBMProxy* capsule_ibm_register_mesh_internal(CapsuleMesh* lag,
                                                          bool sync_now) {
  if (!lag || !lag->isactive || lag->nln <= 0)
    return NULL;

  capsule_ibm_ensure_manager();

  CapsuleIBMProxy* proxy = capsule_ibm_proxy_for_lag(lag);
  int mesh_id = capsule_ibm_find_mesh_id_for_proxy(proxy);
  if (mesh_id < 0) {
    mesh_id = ibmeshmanager_add_mesh();
    capsule_ibm_install_velocity_model(proxy, mesh_id);
  } else if (proxy->node_count != lag->nln) {
    capsule_ibm_install_velocity_model(proxy, mesh_id);
  }

  if (sync_now && capsule_mesh_is_local_owner(lag))
    capsule_ibm_model_sync(proxy, &ibmm.meshes[mesh_id]);
  return proxy;
}

static CapsuleIBMProxy* capsule_ibm_register_mesh(CapsuleMesh* lag) {
  return capsule_ibm_register_mesh_internal(lag, true);
}

static CapsuleIBMProxy* capsule_ibm_register_mesh_no_sync(CapsuleMesh* lag) {
  return capsule_ibm_register_mesh_internal(lag, false);
}

static void capsule_ibm_unregister_mesh(CapsuleMesh* lag) {
  int proxy_id = capsule_ibm_find_proxy(lag);
  if (proxy_id < 0)
    return;

  CapsuleIBMProxy* proxy = &capsule_ibm_proxies[proxy_id];
  int mesh_id = capsule_ibm_find_mesh_id_for_proxy(proxy);
  if (mesh_id >= 0 && mesh_id < ibmm.nm) {
    int deleted_mesh_id = mesh_id;
    int last_mesh_id = ibmm.nm - 1;
    ibmeshmanager_delete_mesh(deleted_mesh_id);
    if (deleted_mesh_id != last_mesh_id) {
      for (int i = 0; i < capsule_ibm_proxy_capacity; i++)
        if (i != proxy_id && capsule_ibm_proxies[i].active &&
            capsule_ibm_proxies[i].mesh_id == last_mesh_id)
          capsule_ibm_proxies[i].mesh_id = deleted_mesh_id;
    }
  }

  capsule_ibm_clear_proxy_slot(proxy_id);
}

static void capsule_ibm_register_active_capsules(void) {
  for (int i = 0; i < allCaps.nbcaps; i++)
    if (CAPS(i).isactive)
      capsule_ibm_register_mesh(&CAPS(i));
#if _MPI
  ibmeshmanager_update_pid();
#endif
}

static void capsule_ibm_gather_velocity(CapsuleMesh* lag) {
  CapsuleIBMProxy* proxy = capsule_ibm_register_mesh(lag);
  if (!proxy)
    return;

  IBMesh* ibmesh = &ibmm.meshes[proxy->mesh_id];
  for (size_t i = 0; i < ibmesh->nodes.size; i++) {
    IBNode* node = ibmesh->nodes.ptrs[i];
    foreach_dimension()
      ibval(nvel.x) = 0.;
    peskin_cosine_kernel_gather_dimensionless(node) {
      foreach_dimension()
        ibval(nvel.x) += weight*u.x[];
    }
  }

#if _MPI
  {
    IBscalar* list = NULL;
    foreach_dimension()
      list = iblist_add(list, nvel.x);
    ibmeshmanager_boundary(list);
    free(list);
  }
#endif

  for (size_t i = 0; i < ibmesh->nodes.size; i++) {
    IBNode *node = ibmesh->nodes.ptrs[i];
    if (!node)
      continue;
    int node_id = node->node_lid >= 0 ? node->node_lid : (int) i;
    if (node_id >= 0 && node_id < lag->nln)
      capsule_ibm_copy_node_velocity_to_lag(lag, node_id, node);
  }
}

static void capsule_ibm_spread_force(vector forcing, CapsuleMesh* lag) {
  CapsuleIBMProxy* proxy = capsule_ibm_register_mesh(lag);
  if (!proxy)
    return;

  IBMesh* ibmesh = &ibmm.meshes[proxy->mesh_id];
  capsule_ibm_model_sync(proxy, ibmesh);

  for (size_t i = 0; i < ibmesh->nodes.size; i++) {
    IBNode* node = ibmesh->nodes.ptrs[i];
    peskin_cosine_kernel_spread_dimensionless(node) {
      foreach_dimension()
        forcing.x[] += weight*ibval(nforce.x)*ibval(nweight)/dv();
    }
  }
}

static void capsule_ibm_cleanup(void) {
  if (!capsule_ibm_proxy_table_initialized)
    return;

  for (int i = 0; i < capsule_ibm_proxy_capacity; i++)
    capsule_ibm_clear_proxy_slot(i);
  free(capsule_ibm_proxies);
  capsule_ibm_proxies = NULL;
  capsule_ibm_proxy_capacity = 0;
  capsule_ibm_proxy_table_initialized = false;
}
