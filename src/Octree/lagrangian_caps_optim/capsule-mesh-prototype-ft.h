#pragma once

#include "capsule-mesh-ft.h"

// ============================================================================
// Type Definitions
// ============================================================================

struct CapsuleMeshPrototype {
  int refcount;
  int nln;
  CapNode *nodes;
  int nle;
  EdgeTopology *edges;
  EdgeState *edge_states;
  int nlt;
  TriangleTopology *triangles;
  TriangleState *triangle_states;
};

// ============================================================================
// Macros
// ============================================================================

#define LAG_EDGE_TOPO(mesh, i)                                                 \
  ((mesh)->prototype ? &((mesh)->prototype->edges[(i)])                        \
                     : &((mesh)->edges[(i)].topo))
#define LAG_EDGE_STATE(mesh, i)                                                \
  ((mesh)->edge_states ? &((mesh)->edge_states[(i)])                           \
                       : &((mesh)->edges[(i)].state))
#define LAG_TRI_TOPO(mesh, i)                                                  \
  ((mesh)->prototype ? &((mesh)->prototype->triangles[(i)])                    \
                     : &((mesh)->triangles[(i)].topo))
#define LAG_TRI_STATE(mesh, i)                                                 \
  ((mesh)->triangle_states ? &((mesh)->triangle_states[(i)])                   \
                           : &((mesh)->triangles[(i)].state))

// ============================================================================
// Function Declarations
// ============================================================================

CapsuleMeshPrototype *capsule_mesh_prototype_new(int nln, int nle, int nlt);
CapsuleMeshPrototype *capsule_mesh_prototype_retain(CapsuleMeshPrototype *prototype);
void capsule_mesh_prototype_release(CapsuleMeshPrototype *prototype);
void capsule_mesh_clear_storage(CapsuleMesh *mesh);
void capsule_mesh_extract_prototype(CapsuleMesh *mesh);
void capsule_mesh_update_prototype_template(CapsuleMesh *mesh);
void capsule_mesh_attach_prototype(CapsuleMesh *mesh, CapsuleMeshPrototype *prototype);

// ============================================================================
// Function Definitions
// ============================================================================

/*!
 * @brief
 */
CapsuleMeshPrototype *capsule_mesh_prototype_new(int nln, int nle, int nlt) {
  CapsuleMeshPrototype *prototype =
      (CapsuleMeshPrototype *)calloc(1, sizeof(CapsuleMeshPrototype));
  assert(prototype);
  prototype->refcount = 1;
  prototype->nln = nln;
  prototype->nodes =
      nln > 0 ? (CapNode *)calloc((size_t)nln, sizeof(CapNode)) : NULL;
  assert(nln == 0 || prototype->nodes);
  prototype->nle = nle;
  prototype->edges =
      nle > 0 ? (EdgeTopology *)calloc((size_t)nle, sizeof(EdgeTopology))
              : NULL;
  assert(nle == 0 || prototype->edges);
  prototype->edge_states =
      nle > 0 ? (EdgeState *)calloc((size_t)nle, sizeof(EdgeState)) : NULL;
  assert(nle == 0 || prototype->edge_states);
  prototype->nlt = nlt;
  prototype->triangles =
      nlt > 0
          ? (TriangleTopology *)calloc((size_t)nlt, sizeof(TriangleTopology))
          : NULL;
  assert(nlt == 0 || prototype->triangles);
  prototype->triangle_states =
      nlt > 0 ? (TriangleState *)calloc((size_t)nlt, sizeof(TriangleState))
              : NULL;
  assert(nlt == 0 || prototype->triangle_states);
  return prototype;
}

/*!
 * @brief
 */
CapsuleMeshPrototype *capsule_mesh_prototype_retain(CapsuleMeshPrototype *prototype) {
  if (prototype)
    prototype->refcount++;
  return prototype;
}

/*!
 * @brief
 */
void capsule_mesh_prototype_release(CapsuleMeshPrototype *prototype) {
  if (!prototype)
    return;
  prototype->refcount--;
  if (prototype->refcount > 0)
    return;
  free(prototype->nodes);
  free(prototype->edges);
  free(prototype->edge_states);
#if dimension > 2
  free(prototype->triangles);
  free(prototype->triangle_states);
#endif
  free(prototype);
}

/*!
 * @brief
 */
void capsule_mesh_clear_storage(CapsuleMesh *mesh) {
  if (!mesh)
    return;
  capsule_mesh_prototype_release(mesh->prototype);
  free(mesh->nodes);
  free(mesh->edges);
  free(mesh->edge_states);
#if dimension > 2
  free(mesh->triangles);
  free(mesh->triangle_states);
#endif
  mesh->prototype = NULL;
  mesh->nodes = NULL;
  mesh->edges = NULL;
  mesh->edge_states = NULL;
#if dimension > 2
  mesh->triangles = NULL;
  mesh->triangle_states = NULL;
#endif
}

/*!
 * @brief
 */
void capsule_mesh_extract_prototype(CapsuleMesh *mesh) {
  assert(mesh);
  assert(!mesh->prototype);
  assert(mesh->edges || mesh->nle == 0);
  assert(mesh->triangles || mesh->nlt == 0);

  CapsuleMeshPrototype *prototype =
      capsule_mesh_prototype_new(mesh->nln, mesh->nle, mesh->nlt);

  if (mesh->nln > 0)
    memcpy(prototype->nodes, mesh->nodes, (size_t)mesh->nln * sizeof(CapNode));

  if (mesh->nle > 0) {
    mesh->edge_states =
        (EdgeState *)calloc((size_t)mesh->nle, sizeof(EdgeState));
    assert(mesh->edge_states);
    for (int i = 0; i < mesh->nle; i++) {
      prototype->edges[i] = mesh->edges[i].topo;
      mesh->edge_states[i] = mesh->edges[i].state;
      prototype->edge_states[i] = mesh->edges[i].state;
    }
  }

  if (mesh->nlt > 0) {
    mesh->triangle_states =
        (TriangleState *)calloc((size_t)mesh->nlt, sizeof(TriangleState));
    assert(mesh->triangle_states);
    for (int i = 0; i < mesh->nlt; i++) {
      prototype->triangles[i] = mesh->triangles[i].topo;
      mesh->triangle_states[i] = mesh->triangles[i].state;
      prototype->triangle_states[i] = mesh->triangles[i].state;
    }
  }

  free(mesh->edges);
  mesh->edges = NULL;
  free(mesh->triangles);
  mesh->triangles = NULL;
  mesh->prototype = prototype;
}

/*!
 * @brief
 */
void capsule_mesh_update_prototype_template(CapsuleMesh *mesh) {
  if (!mesh || !mesh->prototype)
    return;
  CapsuleMeshPrototype *prototype = mesh->prototype;
  assert(prototype->nln == mesh->nln);
  assert(prototype->nle == mesh->nle);
  if (mesh->nln > 0)
    memcpy(prototype->nodes, mesh->nodes, (size_t)mesh->nln * sizeof(CapNode));
  if (mesh->nle > 0)
    memcpy(prototype->edge_states, mesh->edge_states,
           (size_t)mesh->nle * sizeof(EdgeState));
  assert(prototype->nlt == mesh->nlt);
  if (mesh->nlt > 0)
    memcpy(prototype->triangle_states, mesh->triangle_states,
           (size_t)mesh->nlt * sizeof(TriangleState));
}

/*!
 * @brief
 */
void capsule_mesh_attach_prototype(CapsuleMesh *mesh, CapsuleMeshPrototype *prototype) {
  assert(mesh);
  assert(prototype);
  capsule_mesh_clear_storage(mesh);
  mesh->prototype = capsule_mesh_prototype_retain(prototype);
  mesh->nln = prototype->nln;
  mesh->nodes =
      prototype->nln > 0
          ? (CapNode *)malloc((size_t)prototype->nln * sizeof(CapNode))
          : NULL;
  assert(prototype->nln == 0 || mesh->nodes);
  if (prototype->nln > 0) {
    memcpy(mesh->nodes, prototype->nodes,
           (size_t)prototype->nln * sizeof(CapNode));
    for (int i = 0; i < mesh->nln; i++) {
      mesh->nodes[i].lagVel = (coord){0};
      mesh->nodes[i].lagForce = (coord){0};
    }
  }
  mesh->nle = prototype->nle;
  mesh->edge_states =
      prototype->nle > 0
          ? (EdgeState *)malloc((size_t)prototype->nle * sizeof(EdgeState))
          : NULL;
  assert(prototype->nle == 0 || mesh->edge_states);
  if (prototype->nle > 0)
    memcpy(mesh->edge_states, prototype->edge_states,
           (size_t)prototype->nle * sizeof(EdgeState));
  mesh->nlt = prototype->nlt;
  mesh->triangle_states = prototype->nlt > 0
                              ? (TriangleState *)malloc((size_t)prototype->nlt *
                                                        sizeof(TriangleState))
                              : NULL;
  assert(prototype->nlt == 0 || mesh->triangle_states);
  if (prototype->nlt > 0)
    memcpy(mesh->triangle_states, prototype->triangle_states,
           (size_t)prototype->nlt * sizeof(TriangleState));
}
