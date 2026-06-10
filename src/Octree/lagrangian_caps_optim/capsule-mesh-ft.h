#pragma once

// ============================================================================
// Type Definitions
// ============================================================================

/*!
 * @brief
 */
typedef struct CapsuleMeshPrototype CapsuleMeshPrototype;

/*!
 * @brief
 */
typedef struct CapNode {
  coord pos;
  coord lagVel;
  coord lagForce;
  coord normal;
  int nb_neighbors;
  int nb_triangles;
  int nb_fit_iterations;
  int neighbor_ids[6];
  int edge_ids[6];
  int triangle_ids[6];
  double curv;
  double ref_curv;
  double gcurv;
} CapNode;

/*!
 * @brief
 */
typedef struct EdgeTopology {
  int node_ids[2];
  int triangle_ids[2];
  double l0;
} EdgeTopology;

/*!
 * @brief
 */
typedef struct EdgeState {
  double length;
  coord normal;
} EdgeState;

/*!
 * @brief
 */
typedef struct Edge {
  union {
    EdgeTopology topo;
    struct {
      int node_ids[2];
      int triangle_ids[2];
      double l0;
    };
  };
  union {
    EdgeState state;
    struct {
      double length;
      coord normal;
    };
  };
} Edge;

/*!
 * @brief
 */
typedef struct TriangleTopology {
  int node_ids[3];
  int edge_ids[3];
  coord refShape[2];
  double sfc[3][2];
} TriangleTopology;

/*!
 * @brief
 */
typedef struct TriangleState {
  double area;
  coord normal;
  coord centroid;
  double stretch[2];
  double tension[2];
} TriangleState;

/*!
 * @brief
 */
typedef struct Triangle {
  union {
    TriangleTopology topo;
    struct {
      int node_ids[3];
      int edge_ids[3];
      coord refShape[2];
      double sfc[3][2];
    };
  };
  union {
    TriangleState state;
    struct {
      double area;
      coord normal;
      coord centroid;
      double stretch[2];
      double tension[2];
    };
  };
} Triangle;

typedef struct CapsuleMesh {
  int cap_id;
  int cap_type;
  int owner_pid;
  double cap_es;
  double cap_radius;
  CapsuleMeshPrototype *prototype;
  int nln;
  CapNode *nodes;
  int nle;
  Edge *edges;
  EdgeState *edge_states;
  int nlt;
  Triangle *triangles;
  TriangleState *triangle_states;
  coord centroid;
  coord ang_vel;
  double initial_volume;
  double volume;
  double circum_radius;
  double taylor_deform;
  bool updated_stretches;
  bool updated_normals;
  bool updated_curvatures;
  bool isactive;
} CapsuleMesh;

// ============================================================================
// Function Declarations
// ============================================================================

void initialize_empty_capsule(CapsuleMesh *mesh);
static bool capsule_mesh_is_local_owner(CapsuleMesh *mesh);

// ============================================================================
// Function Definitions
// ============================================================================

/*!
 * @brief
 */
void initialize_empty_capsule(CapsuleMesh *mesh) {
  mesh->cap_es = 1.;
  mesh->cap_radius = 1.;
  mesh->cap_id = -1;
  mesh->cap_type = -1;
  mesh->owner_pid = 0;
  mesh->prototype = NULL;
  mesh->nln = 0;
  mesh->nle = 0;
  mesh->nodes = NULL;
  mesh->edges = NULL;
  mesh->edge_states = NULL;
  mesh->nlt = 0;
  mesh->triangles = NULL;
  mesh->triangle_states = NULL;
  mesh->centroid = (coord){0};
  mesh->ang_vel = (coord){0};
  mesh->initial_volume = 0.;
  mesh->volume = 0.;
  mesh->circum_radius = 0.;
  mesh->taylor_deform = 0.;
  mesh->updated_stretches = false;
  mesh->updated_normals = false;
  mesh->updated_curvatures = false;
  mesh->isactive = false;
}

/*!
 * @brief
 */
static bool capsule_mesh_is_local_owner(CapsuleMesh *mesh) {
  return mesh && mesh->owner_pid == pid();
}
