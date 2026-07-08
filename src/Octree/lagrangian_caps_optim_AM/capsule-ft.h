/**
# Triangulated mesh of a capsule

In this file, we implement a Lagrangian mesh for the front-tracking method,
meant to track the position and compute the stresses of an elasitc membrane.
*/

#ifndef dimension
  #define dimension 3
#endif


#ifndef MULT_GRID
  #define MULT_GRID 0
#endif

/*The stencil type is chosen as 3 if not 5 */
#ifndef STENCIL_TYPE
  #define STENCIL_TYPE 3
#endif

/* Repulsive lubrication force to avoid overlapping */
#ifndef LUBR_FORCE
  #define LUBR_FORCE 0
#endif

#ifndef LUBR_VEL
  #define LUBR_VEL 0
#endif

#ifndef LAG_SHARED_TOPOLOGY
  #define LAG_SHARED_TOPOLOGY 1
#endif

#ifndef LAG_REF_GEOMETRY
  #define LAG_REF_GEOMETRY 1
#endif

/*Create the Index_lag*/
scalar Index_lagnode[];
vector Index_lag_id[];

/*Create the dimension ratios for non-cubic domain*/
#ifndef L0_X
  #define L0_X L0
#endif
#ifndef L0_Y
  #define L0_Y L0
#endif
#ifndef L0_Z
  #define L0_Z L0
#endif

coord L0_ratio = {L0_X/L0, L0_Y/L0, L0_Z/L0};

/**
## Structure of the mesh

In the Lagrangian mesh, each node is assigned several attributes:

* ```pos```, the node coordinates
* ```lagVel```, the node velocity vector
* ```normal```, the vector normal to the membrane at the node coordinates
* ```curv```, the (mean) curvature of the membrane at the node coordinates
* ```ref_curv```, the reference (mean) curvature of the membrane at the node coordinates (default is zero)
* ```lagForce```, an elastic and/or bending force exerted by the membrane on its surrounding fluid
* ```stencil```, a [Cache](http://basilisk.fr/src/grid/tree.h#82) structure used for averaging the neighboring velocities and spreading the Lagrangian force as a body force on the neighboring Eulerian nodes
* in case of MPI simulations, ```pid```, the rank of the processor owning the Eulerian cell which contains the lagNode
* ```edge_ids```, the IDs of its connecting edges: 2 in 2D, up to 6 in 3D (because every considered shape is derived by subdividing in icosahedron, leading to 5 or 6 neighbors only).

In case of 3D simulations, other attributes are introduced:

* ```nb_neighbors```, the number of neighbors of the node, which should only be 5 or 6
* ```neighbor_ids```, the IDs of the node neighbors
* ```nb_triangles```, the number of triangles having the node as a vertex
* ```triangle_ids```, the ids of the above triangles
* ```gcurv```, the Gaussian curvature of the membrane at the node coordinates
* ```nb_fit_iterations```, the number of iterations needed to compute the membrane curvature and normal vector to the desired convergence threshold.

*/
typedef struct lagNode {
  coord pos;
  coord lagVel;
  coord normal;
  double curv;
  double ref_curv;
  coord lagForce;
  Cache stencil;
  Cache eulcell;
  #if _MPI
    int pid;
  #endif
  #if dimension < 3
    #if !LAG_SHARED_TOPOLOGY
      int edge_ids[2];
    #endif
  #else
    #if !LAG_SHARED_TOPOLOGY
      int nb_neighbors;       // 1 ints
      int neighbor_ids[6];    // 6 ints
      int edge_ids[6];        // 6 ints
      int nb_triangles;       // 1 ints
      int triangle_ids[6];    // 6 ints
    #endif
    double gcurv;
    int nb_fit_iterations;
  #endif
} lagNode;

/** We specify the size of the 3x3(x3) or 5x5(x5) stencil in 2D(3D). */
#if dimension < 3
  #if STENCIL_TYPE == 3
    #define STENCIL_SIZE 9
  #else 
    #define STENCIL_SIZE 25
  #endif
#else
  #if STENCIL_TYPE == 3
    #define STENCIL_SIZE 27
  #else 
    #define STENCIL_SIZE 125
  #endif
#endif

/** Similarly, the edges of the mesh are assigned:

* ```node_ids```, the IDs of the two nodes they connect
* In case of 3D simulations, ```triangle_ids```, the IDs of the two triangles they separate
* ```l0```, the length of the edge in the initial stress-free configuration
* ```length```, the current length of the edge
* ```normal```, a vector normal to the membrane at the edge midpoint.

*/
typedef struct Edge {
  #if !LAG_SHARED_TOPOLOGY
    int node_ids[2];    // 2 ints
    #if dimension > 2
      int triangle_ids[2];  //2 ints
    #endif
  #endif
  #if !LAG_REF_GEOMETRY
    double l0;
  #endif
  double length;
  coord normal;
} Edge;

/** In case of 3D simulations, we define triangle faces. Each ```Triangle``` structure has the following attributes:

* ```node_ids```, the IDs of the triangle vertices
* ```edge_ids```, the IDs of the triangle edges
* ```area```, the current area of the triangle
* ```normal```, the normal vector to the triangle pointing outside the membrane
* ```centroid```, the coordinates of the triangle centroid
* ```refShape```, the coordinates of the vertices in the Common Plane in the stress-free configuration. By convention the first vertex is always placed at $(0, 0, 0)$, so only the coordinates of the second and third vertex are stored in ```refShape```
* ```sfc```, the shape function coefficients used in the finite elements method and computed in [store_initial_configuration](elasticity-ft.h#store_initial_configuration). There are two coefficients per vertex, hence the $3 \times 2$ array structure of ```sfc```.

 */
#if dimension > 2
  typedef struct Triangle {
    #if !LAG_SHARED_TOPOLOGY
      int node_ids[3];   // 3 ints
      int edge_ids[3];   // 3 ints
    #endif
    double area;
    coord normal;
    coord centroid;
    #if !LAG_REF_GEOMETRY
       coord refShape[2];
       double sfc[3][2]; // sfc for "shape function coefficients"
    #endif
    double stretch[2];
    double tension[2];
  } Triangle;
#endif

/** Connectivity shared by capsules with the same membrane topology.

The arrays below mirror the integer connectivity currently stored inside
lagNode, Edge and Triangle, but are grouped in one topology object so that
multiple capsules can reference the same mesh connectivity.
*/
typedef struct lagTopology {
  int nln;
  int nle;
  #if dimension > 2
    int nlt;
  #endif
  #if dimension < 3
    int (*node_edge_ids)[2];
  #else
    int *node_nb_neighbors;
    int (*node_neighbor_ids)[6];
    int (*node_edge_ids)[6];
    int *node_nb_triangles;
    int (*node_triangle_ids)[6];
    int (*edge_triangle_ids)[2];
    int (*triangle_node_ids)[3];
    int (*triangle_edge_ids)[3];
  #endif
  int (*edge_node_ids)[2];
} lagTopology;

typedef struct lagReferenceGeometry {
  int nle;
  double *edge_l0;
  #if dimension > 2
    int nlt;
    coord (*triangle_refShape)[2];
    double (*triangle_sfc)[3][2];
  #endif
  double initial_volume;
} lagReferenceGeometry;

/** The ```lagMesh``` structure contains arrays of the previously introduced nodes, edges and triangles. It defines an unstructured mesh, the membrane of our capsule. Its attributes are:

*  ```nln``` the number of Lagrangian nodes
* ```nodes```, the array of Lagrangian nodes
* ```nle```, the number of Lagrangian edges
* ```edges```, the array of Lagrangian edges
* In case of 3D simulations:
    * ```nlt```, the number of Lagrangian triangles
    * ```triangles```, the array of Lagrangian triangles
* ```updated_stretches```, a boolean used to check if the current length of the edges has been updated since the last advection of the Lagrangian nodes
* ```updated_normals```, a similar boolean telling if the nodal normal vectors should be recomputed
* ```updated_curvatures```, a last boolean telling if the nodal curvatures should be recomputed.
* ```isactive```, a boolean indicating if the capsule exists in the flow (useful
when capsules are introduced during a simulation)
*/

typedef struct lagMesh {
  int cap_id;
  int cap_type;
  double cap_es;
  double cap_radius;
  lagTopology* topology;
  lagReferenceGeometry* ref_geometry;
  int nln;
  lagNode* nodes;
  int nle;
  Edge* edges;
  #if dimension > 2
    int nlt;
    Triangle* triangles;
  #endif
  coord centroid;
  coord ang_vel;
  #if !LAG_REF_GEOMETRY
    double initial_volume;
  #endif
  double volume;
  double circum_radius;
  double taylor_deform;
  bool updated_stretches;
  bool updated_normals;
  bool updated_curvatures;
  bool isactive;
} lagMesh;

/** We denote by ```NCAPS``` the number of Lagrangian meshes, or capsules, in
the simulation. It is one by default. */
#ifndef NCAPS
  #define NCAPS 1
#endif
#ifndef RESTART_CASE
  #define RESTART_CASE 0
#endif
#ifndef LAG_TOPOLOGY_DEBUG
  #define LAG_TOPOLOGY_DEBUG 0
#endif

/** The Lagrangian mesh is accessible in the code thanks to the structure
below, which is simply an array of Lagrangian meshes (useful when several of
them are considered). The macro $CAPS(k)$ can be used as a shortcut to access the
$k^{th}$ membrane. */
typedef struct Capsules {
  lagMesh caps[NCAPS];
  int nbcaps;
} Capsules;
Capsules allCaps;
#define CAPS(i) (allCaps.caps[i])

#if LAG_SHARED_TOPOLOGY
  #if dimension < 3
    #define LAG_NODE_EDGE_ID(M, I, J) ((M)->topology->node_edge_ids[I][J])
    #define SET_LAG_NODE_EDGE_ID(M, I, J, V) do { \
      (M)->topology->node_edge_ids[I][J] = (V); \
    } while (0)
  #else
    #define LAG_NODE_NB_NEIGHBORS(M, I) ((M)->topology->node_nb_neighbors[I])
    #define LAG_NODE_NEIGHBOR_ID(M, I, J) ((M)->topology->node_neighbor_ids[I][J])
    #define LAG_NODE_EDGE_ID(M, I, J) ((M)->topology->node_edge_ids[I][J])
    #define LAG_NODE_NB_TRIANGLES(M, I) ((M)->topology->node_nb_triangles[I])
    #define LAG_NODE_TRIANGLE_ID(M, I, J) ((M)->topology->node_triangle_ids[I][J])
    #define LAG_EDGE_TRIANGLE_ID(M, I, J) ((M)->topology->edge_triangle_ids[I][J])
    #define LAG_TRIANGLE_NODE_ID(M, I, J) ((M)->topology->triangle_node_ids[I][J])
    #define LAG_TRIANGLE_EDGE_ID(M, I, J) ((M)->topology->triangle_edge_ids[I][J])
    #define SET_LAG_NODE_NB_NEIGHBORS(M, I, V) do { \
      (M)->topology->node_nb_neighbors[I] = (V); \
    } while (0)
    #define SET_LAG_NODE_NEIGHBOR_ID(M, I, J, V) do { \
      (M)->topology->node_neighbor_ids[I][J] = (V); \
    } while (0)
    #define SET_LAG_NODE_EDGE_ID(M, I, J, V) do { \
      (M)->topology->node_edge_ids[I][J] = (V); \
    } while (0)
    #define SET_LAG_NODE_NB_TRIANGLES(M, I, V) do { \
      (M)->topology->node_nb_triangles[I] = (V); \
    } while (0)
    #define SET_LAG_NODE_TRIANGLE_ID(M, I, J, V) do { \
      (M)->topology->node_triangle_ids[I][J] = (V); \
    } while (0)
    #define SET_LAG_EDGE_TRIANGLE_ID(M, I, J, V) do { \
      (M)->topology->edge_triangle_ids[I][J] = (V); \
    } while (0)
    #define SET_LAG_TRIANGLE_NODE_ID(M, I, J, V) do { \
      (M)->topology->triangle_node_ids[I][J] = (V); \
    } while (0)
    #define SET_LAG_TRIANGLE_EDGE_ID(M, I, J, V) do { \
      (M)->topology->triangle_edge_ids[I][J] = (V); \
    } while (0)
  #endif
  #define LAG_EDGE_NODE_ID(M, I, J) ((M)->topology->edge_node_ids[I][J])
  #define SET_LAG_EDGE_NODE_ID(M, I, J, V) do { \
    (M)->topology->edge_node_ids[I][J] = (V); \
  } while (0)
#else
  #if dimension < 3
    #define LAG_NODE_EDGE_ID(M, I, J) ((M)->topology ? \
      (M)->topology->node_edge_ids[I][J] : (M)->nodes[I].edge_ids[J])
    #define SET_LAG_NODE_EDGE_ID(M, I, J, V) do { \
      if ((M)->topology) (M)->topology->node_edge_ids[I][J] = (V); \
      else (M)->nodes[I].edge_ids[J] = (V); \
    } while (0)
  #else
    #define LAG_NODE_NB_NEIGHBORS(M, I) ((M)->topology ? \
      (M)->topology->node_nb_neighbors[I] : (M)->nodes[I].nb_neighbors)
    #define LAG_NODE_NEIGHBOR_ID(M, I, J) ((M)->topology ? \
      (M)->topology->node_neighbor_ids[I][J] : (M)->nodes[I].neighbor_ids[J])
    #define LAG_NODE_EDGE_ID(M, I, J) ((M)->topology ? \
      (M)->topology->node_edge_ids[I][J] : (M)->nodes[I].edge_ids[J])
    #define LAG_NODE_NB_TRIANGLES(M, I) ((M)->topology ? \
      (M)->topology->node_nb_triangles[I] : (M)->nodes[I].nb_triangles)
    #define LAG_NODE_TRIANGLE_ID(M, I, J) ((M)->topology ? \
      (M)->topology->node_triangle_ids[I][J] : (M)->nodes[I].triangle_ids[J])
    #define LAG_EDGE_TRIANGLE_ID(M, I, J) ((M)->topology ? \
      (M)->topology->edge_triangle_ids[I][J] : (M)->edges[I].triangle_ids[J])
    #define LAG_TRIANGLE_NODE_ID(M, I, J) ((M)->topology ? \
      (M)->topology->triangle_node_ids[I][J] : (M)->triangles[I].node_ids[J])
    #define LAG_TRIANGLE_EDGE_ID(M, I, J) ((M)->topology ? \
      (M)->topology->triangle_edge_ids[I][J] : (M)->triangles[I].edge_ids[J])
    #define SET_LAG_NODE_NB_NEIGHBORS(M, I, V) do { \
      if ((M)->topology) (M)->topology->node_nb_neighbors[I] = (V); \
      else (M)->nodes[I].nb_neighbors = (V); \
    } while (0)
    #define SET_LAG_NODE_NEIGHBOR_ID(M, I, J, V) do { \
      if ((M)->topology) (M)->topology->node_neighbor_ids[I][J] = (V); \
      else (M)->nodes[I].neighbor_ids[J] = (V); \
    } while (0)
    #define SET_LAG_NODE_EDGE_ID(M, I, J, V) do { \
      if ((M)->topology) (M)->topology->node_edge_ids[I][J] = (V); \
      else (M)->nodes[I].edge_ids[J] = (V); \
    } while (0)
    #define SET_LAG_NODE_NB_TRIANGLES(M, I, V) do { \
      if ((M)->topology) (M)->topology->node_nb_triangles[I] = (V); \
      else (M)->nodes[I].nb_triangles = (V); \
    } while (0)
    #define SET_LAG_NODE_TRIANGLE_ID(M, I, J, V) do { \
      if ((M)->topology) (M)->topology->node_triangle_ids[I][J] = (V); \
      else (M)->nodes[I].triangle_ids[J] = (V); \
    } while (0)
    #define SET_LAG_EDGE_TRIANGLE_ID(M, I, J, V) do { \
      if ((M)->topology) (M)->topology->edge_triangle_ids[I][J] = (V); \
      else (M)->edges[I].triangle_ids[J] = (V); \
    } while (0)
    #define SET_LAG_TRIANGLE_NODE_ID(M, I, J, V) do { \
      if ((M)->topology) (M)->topology->triangle_node_ids[I][J] = (V); \
      else (M)->triangles[I].node_ids[J] = (V); \
    } while (0)
    #define SET_LAG_TRIANGLE_EDGE_ID(M, I, J, V) do { \
      if ((M)->topology) (M)->topology->triangle_edge_ids[I][J] = (V); \
      else (M)->triangles[I].edge_ids[J] = (V); \
    } while (0)
  #endif

  #define LAG_EDGE_NODE_ID(M, I, J) ((M)->topology ? \
    (M)->topology->edge_node_ids[I][J] : (M)->edges[I].node_ids[J])
  #define SET_LAG_EDGE_NODE_ID(M, I, J, V) do { \
    if ((M)->topology) (M)->topology->edge_node_ids[I][J] = (V); \
    else (M)->edges[I].node_ids[J] = (V); \
  } while (0)
#endif

#if LAG_REF_GEOMETRY
  #define LAG_EDGE_L0(M, I) ((M)->ref_geometry->edge_l0[I])
  #define SET_LAG_EDGE_L0(M, I, V) do { \
    (M)->ref_geometry->edge_l0[I] = (V); \
  } while (0)
  #define LAG_INITIAL_VOLUME(M) ((M)->ref_geometry->initial_volume)
  #define SET_LAG_INITIAL_VOLUME(M, V) do { \
    (M)->ref_geometry->initial_volume = (V); \
  } while (0)
  #if dimension > 2
    #define LAG_TRIANGLE_REFSHAPE(M, I) ((M)->ref_geometry->triangle_refShape[I])
    #define LAG_TRIANGLE_REFSHAPE_COMPONENT(M, I, J, C) \
      ((M)->ref_geometry->triangle_refShape[I][J].C)
    #define SET_LAG_TRIANGLE_REFSHAPE_COMPONENT(M, I, J, C, V) do { \
      (M)->ref_geometry->triangle_refShape[I][J].C = (V); \
    } while (0)
    #define LAG_TRIANGLE_SFC(M, I, J, K) ((M)->ref_geometry->triangle_sfc[I][J][K])
    #define SET_LAG_TRIANGLE_SFC(M, I, J, K, V) do { \
      (M)->ref_geometry->triangle_sfc[I][J][K] = (V); \
    } while (0)
  #endif
#else
  #define LAG_EDGE_L0(M, I) ((M)->edges[I].l0)
  #define SET_LAG_EDGE_L0(M, I, V) do { \
    (M)->edges[I].l0 = (V); \
  } while (0)
  #define LAG_INITIAL_VOLUME(M) ((M)->initial_volume)
  #define SET_LAG_INITIAL_VOLUME(M, V) do { \
    (M)->initial_volume = (V); \
  } while (0)
  #if dimension > 2
    #define LAG_TRIANGLE_REFSHAPE(M, I) ((M)->triangles[I].refShape)
    #define LAG_TRIANGLE_REFSHAPE_COMPONENT(M, I, J, C) \
      ((M)->triangles[I].refShape[J].C)
    #define SET_LAG_TRIANGLE_REFSHAPE_COMPONENT(M, I, J, C, V) do { \
      (M)->triangles[I].refShape[J].C = (V); \
    } while (0)
    #define LAG_TRIANGLE_SFC(M, I, J, K) ((M)->triangles[I].sfc[J][K])
    #define SET_LAG_TRIANGLE_SFC(M, I, J, K, V) do { \
      (M)->triangles[I].sfc[J][K] = (V); \
    } while (0)
  #endif
#endif


/**
## Initialization, memory management and useful macros.
*/
void initialize_empty_capsule(lagMesh* mesh) {
  mesh->cap_es = 1.;
  mesh->cap_radius = 1.;
  mesh->cap_id = -1;
  mesh->cap_type = -1;
  mesh->topology = NULL;
  mesh->ref_geometry = NULL;
  mesh->nln = 0;
  mesh->nle = 0;
  mesh->nodes = NULL;
  mesh->edges = NULL;
  #if dimension > 2
    mesh->nlt = 0;
    mesh->triangles = NULL;
  #endif
  mesh->updated_stretches = false;
  mesh->updated_normals = false;
  mesh->updated_curvatures = false;
  mesh->isactive = false;
}

lagTopology* allocate_lag_topology(int nln, int nle, int nlt) {
  lagTopology* topology = calloc(1, sizeof(lagTopology));
  assert(topology);
  topology->nln = nln;
  topology->nle = nle;
  topology->edge_node_ids = malloc(nle*sizeof(int[2]));
  assert(topology->edge_node_ids || nle == 0);

  #if dimension < 3
    topology->node_edge_ids = malloc(nln*sizeof(int[2]));
    assert(topology->node_edge_ids || nln == 0);
  #else
    topology->nlt = nlt;
    topology->node_nb_neighbors = malloc(nln*sizeof(int));
    topology->node_neighbor_ids = malloc(nln*sizeof(int[6]));
    topology->node_edge_ids = malloc(nln*sizeof(int[6]));
    topology->node_nb_triangles = malloc(nln*sizeof(int));
    topology->node_triangle_ids = malloc(nln*sizeof(int[6]));
    topology->edge_triangle_ids = malloc(nle*sizeof(int[2]));
    topology->triangle_node_ids = malloc(nlt*sizeof(int[3]));
    topology->triangle_edge_ids = malloc(nlt*sizeof(int[3]));
    assert(topology->node_nb_neighbors || nln == 0);
    assert(topology->node_neighbor_ids || nln == 0);
    assert(topology->node_edge_ids || nln == 0);
    assert(topology->node_nb_triangles || nln == 0);
    assert(topology->node_triangle_ids || nln == 0);
    assert(topology->edge_triangle_ids || nle == 0);
    assert(topology->triangle_node_ids || nlt == 0);
    assert(topology->triangle_edge_ids || nlt == 0);
  #endif

  for(int i=0; i<nln; i++) {
    #if dimension < 3
      for(int j=0; j<2; j++) topology->node_edge_ids[i][j] = -1;
    #else
      topology->node_nb_neighbors[i] = 0;
      topology->node_nb_triangles[i] = 0;
      for(int j=0; j<6; j++) {
        topology->node_neighbor_ids[i][j] = -1;
        topology->node_edge_ids[i][j] = -1;
        topology->node_triangle_ids[i][j] = -1;
      }
    #endif
  }
  for(int i=0; i<nle; i++) {
    for(int j=0; j<2; j++) topology->edge_node_ids[i][j] = -1;
    #if dimension > 2
      for(int j=0; j<2; j++) topology->edge_triangle_ids[i][j] = -1;
    #endif
  }
  #if dimension > 2
    for(int i=0; i<nlt; i++)
      for(int j=0; j<3; j++) {
        topology->triangle_node_ids[i][j] = -1;
        topology->triangle_edge_ids[i][j] = -1;
      }
  #endif

  return topology;
}

void resize_lag_topology(lagTopology* topology, int nln, int nle, int nlt) {
  assert(topology);
  int old_nln = topology->nln;
  int old_nle = topology->nle;
  #if dimension > 2
    int old_nlt = topology->nlt;
  #endif

  topology->nln = nln;
  topology->nle = nle;
  topology->edge_node_ids = realloc(topology->edge_node_ids, nle*sizeof(int[2]));
  assert(topology->edge_node_ids || nle == 0);

  #if dimension < 3
    topology->node_edge_ids = realloc(topology->node_edge_ids, nln*sizeof(int[2]));
    assert(topology->node_edge_ids || nln == 0);
  #else
    topology->nlt = nlt;
    topology->node_nb_neighbors = realloc(topology->node_nb_neighbors, nln*sizeof(int));
    topology->node_neighbor_ids = realloc(topology->node_neighbor_ids, nln*sizeof(int[6]));
    topology->node_edge_ids = realloc(topology->node_edge_ids, nln*sizeof(int[6]));
    topology->node_nb_triangles = realloc(topology->node_nb_triangles, nln*sizeof(int));
    topology->node_triangle_ids = realloc(topology->node_triangle_ids, nln*sizeof(int[6]));
    topology->edge_triangle_ids = realloc(topology->edge_triangle_ids, nle*sizeof(int[2]));
    topology->triangle_node_ids = realloc(topology->triangle_node_ids, nlt*sizeof(int[3]));
    topology->triangle_edge_ids = realloc(topology->triangle_edge_ids, nlt*sizeof(int[3]));
    assert(topology->node_nb_neighbors || nln == 0);
    assert(topology->node_neighbor_ids || nln == 0);
    assert(topology->node_edge_ids || nln == 0);
    assert(topology->node_nb_triangles || nln == 0);
    assert(topology->node_triangle_ids || nln == 0);
    assert(topology->edge_triangle_ids || nle == 0);
    assert(topology->triangle_node_ids || nlt == 0);
    assert(topology->triangle_edge_ids || nlt == 0);
  #endif

  for(int i=old_nln; i<nln; i++) {
    #if dimension < 3
      for(int j=0; j<2; j++) topology->node_edge_ids[i][j] = -1;
    #else
      topology->node_nb_neighbors[i] = 0;
      topology->node_nb_triangles[i] = 0;
      for(int j=0; j<6; j++) {
        topology->node_neighbor_ids[i][j] = -1;
        topology->node_edge_ids[i][j] = -1;
        topology->node_triangle_ids[i][j] = -1;
      }
    #endif
  }
  for(int i=old_nle; i<nle; i++) {
    for(int j=0; j<2; j++) topology->edge_node_ids[i][j] = -1;
    #if dimension > 2
      for(int j=0; j<2; j++) topology->edge_triangle_ids[i][j] = -1;
    #endif
  }
  #if dimension > 2
    for(int i=old_nlt; i<nlt; i++)
      for(int j=0; j<3; j++) {
        topology->triangle_node_ids[i][j] = -1;
        topology->triangle_edge_ids[i][j] = -1;
      }
  #else
    (void) nlt;
  #endif
}

lagTopology* copy_lag_topology_from_mesh(lagMesh* mesh) {
  #if dimension > 2
    int nlt = mesh->nlt;
  #else
    int nlt = 0;
  #endif
  lagTopology* topology = allocate_lag_topology(mesh->nln, mesh->nle, nlt);

  for(int i=0; i<mesh->nln; i++) {
    #if dimension < 3
      for(int j=0; j<2; j++)
        topology->node_edge_ids[i][j] = LAG_NODE_EDGE_ID(mesh, i, j);
    #else
      topology->node_nb_neighbors[i] = LAG_NODE_NB_NEIGHBORS(mesh, i);
      topology->node_nb_triangles[i] = LAG_NODE_NB_TRIANGLES(mesh, i);
      for(int j=0; j<6; j++) {
        topology->node_neighbor_ids[i][j] = LAG_NODE_NEIGHBOR_ID(mesh, i, j);
        topology->node_edge_ids[i][j] = LAG_NODE_EDGE_ID(mesh, i, j);
        topology->node_triangle_ids[i][j] = LAG_NODE_TRIANGLE_ID(mesh, i, j);
      }
    #endif
  }
  for(int i=0; i<mesh->nle; i++) {
    for(int j=0; j<2; j++)
      topology->edge_node_ids[i][j] = LAG_EDGE_NODE_ID(mesh, i, j);
    #if dimension > 2
      for(int j=0; j<2; j++)
        topology->edge_triangle_ids[i][j] = LAG_EDGE_TRIANGLE_ID(mesh, i, j);
    #endif
  }
  #if dimension > 2
    for(int i=0; i<mesh->nlt; i++) {
      for(int j=0; j<3; j++) {
        topology->triangle_node_ids[i][j] = LAG_TRIANGLE_NODE_ID(mesh, i, j);
        topology->triangle_edge_ids[i][j] = LAG_TRIANGLE_EDGE_ID(mesh, i, j);
      }
    }
  #endif
  return topology;
}

void free_lag_topology(lagTopology* topology) {
  if (!topology) return;
  free(topology->edge_node_ids);
  #if dimension < 3
    free(topology->node_edge_ids);
  #else
    free(topology->node_nb_neighbors);
    free(topology->node_neighbor_ids);
    free(topology->node_edge_ids);
    free(topology->node_nb_triangles);
    free(topology->node_triangle_ids);
    free(topology->edge_triangle_ids);
    free(topology->triangle_node_ids);
    free(topology->triangle_edge_ids);
  #endif
  free(topology);
}

lagReferenceGeometry* allocate_lag_ref_geometry(int nle, int nlt) {
  lagReferenceGeometry* ref = calloc(1, sizeof(lagReferenceGeometry));
  assert(ref);
  ref->nle = nle;
  ref->edge_l0 = malloc(nle*sizeof(double));
  assert(ref->edge_l0 || nle == 0);
  for(int i=0; i<nle; i++)
    ref->edge_l0[i] = 0.;

  #if dimension > 2
    ref->nlt = nlt;
    ref->triangle_refShape = malloc(nlt*sizeof(coord[2]));
    ref->triangle_sfc = malloc(nlt*sizeof(double[3][2]));
    assert(ref->triangle_refShape || nlt == 0);
    assert(ref->triangle_sfc || nlt == 0);
    for(int i=0; i<nlt; i++) {
      for(int j=0; j<2; j++)
        foreach_dimension()
          ref->triangle_refShape[i][j].x = 0.;
      for(int j=0; j<3; j++)
        for(int k=0; k<2; k++)
          ref->triangle_sfc[i][j][k] = 0.;
    }
  #else
    (void) nlt;
  #endif

  ref->initial_volume = 0.;
  return ref;
}

void resize_lag_ref_geometry(lagReferenceGeometry* ref, int nle, int nlt) {
  assert(ref);
  int old_nle = ref->nle;
  #if dimension > 2
    int old_nlt = ref->nlt;
  #endif

  ref->nle = nle;
  ref->edge_l0 = realloc(ref->edge_l0, nle*sizeof(double));
  assert(ref->edge_l0 || nle == 0);
  for(int i=old_nle; i<nle; i++)
    ref->edge_l0[i] = 0.;

  #if dimension > 2
    ref->nlt = nlt;
    ref->triangle_refShape = realloc(ref->triangle_refShape, nlt*sizeof(coord[2]));
    ref->triangle_sfc = realloc(ref->triangle_sfc, nlt*sizeof(double[3][2]));
    assert(ref->triangle_refShape || nlt == 0);
    assert(ref->triangle_sfc || nlt == 0);
    for(int i=old_nlt; i<nlt; i++) {
      for(int j=0; j<2; j++)
        foreach_dimension()
          ref->triangle_refShape[i][j].x = 0.;
      for(int j=0; j<3; j++)
        for(int k=0; k<2; k++)
          ref->triangle_sfc[i][j][k] = 0.;
    }
  #else
    (void) nlt;
  #endif
}

void free_lag_ref_geometry(lagReferenceGeometry* ref) {
  if (!ref) return;
  free(ref->edge_l0);
  #if dimension > 2
    free(ref->triangle_refShape);
    free(ref->triangle_sfc);
  #endif
  free(ref);
}

bool lag_ref_geometry_matches_mesh(lagReferenceGeometry* ref, lagMesh* mesh) {
  if (!ref || ref->nle != mesh->nle)
    return false;
  #if dimension > 2
    if (ref->nlt != mesh->nlt)
      return false;
  #endif

  for(int i=0; i<mesh->nle; i++)
    if (fabs(ref->edge_l0[i] - LAG_EDGE_L0(mesh, i)) > 1.e-12)
      return false;
  if (fabs(ref->initial_volume - LAG_INITIAL_VOLUME(mesh)) > 1.e-12)
    return false;

  #if dimension > 2
    for(int i=0; i<mesh->nlt; i++) {
      for(int j=0; j<2; j++)
        foreach_dimension()
          if (fabs(ref->triangle_refShape[i][j].x -
            LAG_TRIANGLE_REFSHAPE_COMPONENT(mesh, i, j, x)) > 1.e-12)
            return false;
      for(int j=0; j<3; j++)
        for(int k=0; k<2; k++)
          if (fabs(ref->triangle_sfc[i][j][k] -
            LAG_TRIANGLE_SFC(mesh, i, j, k)) > 1.e-12)
            return false;
    }
  #endif
  return true;
}

typedef struct lagTopologyRegistry {
  int n;
  int nm;
  lagTopology** items;
} lagTopologyRegistry;

static lagTopologyRegistry lag_topologies = {0, 0, NULL};

typedef struct lagRefGeometryRegistry {
  int n;
  int nm;
  lagReferenceGeometry** items;
} lagRefGeometryRegistry;

static lagRefGeometryRegistry lag_ref_geometries = {0, 0, NULL};

#if LAG_TOPOLOGY_DEBUG
  void lag_topology_debug_error(int* errors, const char* what, int a, int b, int c) {
    if (*errors < 20)
      fprintf(stderr, "lagTopology debug error: %s (%d, %d, %d)\n", what, a, b, c);
    (*errors)++;
  }

  int lag_topology_validate(lagMesh* mesh) {
    int errors = 0;
    if (!mesh->topology)
      lag_topology_debug_error(&errors, "missing topology", mesh->cap_id, -1, -1);
    if (mesh->topology && mesh->topology->nln != mesh->nln)
      lag_topology_debug_error(&errors, "nln mismatch", mesh->topology->nln, mesh->nln, -1);
    if (mesh->topology && mesh->topology->nle != mesh->nle)
      lag_topology_debug_error(&errors, "nle mismatch", mesh->topology->nle, mesh->nle, -1);
    #if dimension > 2
      if (mesh->topology && mesh->topology->nlt != mesh->nlt)
        lag_topology_debug_error(&errors, "nlt mismatch", mesh->topology->nlt, mesh->nlt, -1);
    #endif

    for(int i=0; i<mesh->nle; i++) {
      int n0 = LAG_EDGE_NODE_ID(mesh, i, 0);
      int n1 = LAG_EDGE_NODE_ID(mesh, i, 1);
      if (n0 < 0 || n0 >= mesh->nln)
        lag_topology_debug_error(&errors, "edge node 0 out of range", i, n0, mesh->nln);
      if (n1 < 0 || n1 >= mesh->nln)
        lag_topology_debug_error(&errors, "edge node 1 out of range", i, n1, mesh->nln);
    }

    #if dimension < 3
      for(int i=0; i<mesh->nln; i++)
        for(int j=0; j<2; j++) {
          int eid = LAG_NODE_EDGE_ID(mesh, i, j);
          if (eid < 0 || eid >= mesh->nle)
            lag_topology_debug_error(&errors, "node edge out of range", i, j, eid);
        }
    #else
      for(int i=0; i<mesh->nln; i++) {
        int nbn = LAG_NODE_NB_NEIGHBORS(mesh, i);
        int nbt = LAG_NODE_NB_TRIANGLES(mesh, i);
        if (nbn < 0 || nbn > 6)
          lag_topology_debug_error(&errors, "node neighbor count out of range", i, nbn, -1);
        if (nbt < 0 || nbt > 6)
          lag_topology_debug_error(&errors, "node triangle count out of range", i, nbt, -1);
        for(int j=0; j<nbn; j++) {
          int ngb = LAG_NODE_NEIGHBOR_ID(mesh, i, j);
          int eid = LAG_NODE_EDGE_ID(mesh, i, j);
          if (ngb < 0 || ngb >= mesh->nln) {
            lag_topology_debug_error(&errors, "node neighbor out of range", i, j, ngb);
            continue;
          }
          if (eid < 0 || eid >= mesh->nle)
            lag_topology_debug_error(&errors, "node edge out of range", i, j, eid);
          bool reciprocal = false;
          for(int k=0; k<LAG_NODE_NB_NEIGHBORS(mesh, ngb); k++)
            if (LAG_NODE_NEIGHBOR_ID(mesh, ngb, k) == i)
              reciprocal = true;
          if (!reciprocal)
            lag_topology_debug_error(&errors, "node neighbor not reciprocal", i, ngb, -1);
        }
        for(int j=0; j<nbt; j++) {
          int tid = LAG_NODE_TRIANGLE_ID(mesh, i, j);
          if (tid < 0 || tid >= mesh->nlt) {
            lag_topology_debug_error(&errors, "node triangle out of range", i, j, tid);
            continue;
          }
          bool found = false;
          for(int k=0; k<3; k++)
            if (LAG_TRIANGLE_NODE_ID(mesh, tid, k) == i)
              found = true;
          if (!found)
            lag_topology_debug_error(&errors, "node triangle does not contain node", i, tid, -1);
        }
      }

      for(int i=0; i<mesh->nlt; i++) {
        int tn[3];
        for(int j=0; j<3; j++) {
          tn[j] = LAG_TRIANGLE_NODE_ID(mesh, i, j);
          if (tn[j] < 0 || tn[j] >= mesh->nln)
            lag_topology_debug_error(&errors, "triangle node out of range", i, j, tn[j]);
        }
        for(int j=0; j<3; j++) {
          int eid = LAG_TRIANGLE_EDGE_ID(mesh, i, j);
          if (eid < 0 || eid >= mesh->nle) {
            lag_topology_debug_error(&errors, "triangle edge out of range", i, j, eid);
            continue;
          }
          int en0 = LAG_EDGE_NODE_ID(mesh, eid, 0);
          int en1 = LAG_EDGE_NODE_ID(mesh, eid, 1);
          bool en0_in_tri = false, en1_in_tri = false;
          for(int k=0; k<3; k++) {
            if (tn[k] == en0) en0_in_tri = true;
            if (tn[k] == en1) en1_in_tri = true;
          }
          if (!en0_in_tri || !en1_in_tri)
            lag_topology_debug_error(&errors, "triangle edge does not match triangle nodes", i, j, eid);
        }
      }

      for(int i=0; i<mesh->nle; i++)
        for(int j=0; j<2; j++) {
          int tid = LAG_EDGE_TRIANGLE_ID(mesh, i, j);
          if (tid < 0 || tid >= mesh->nlt)
            lag_topology_debug_error(&errors, "edge triangle out of range", i, j, tid);
        }
    #endif

    return errors;
  }

  void debug_lag_topology(lagMesh* mesh, const char* context) {
    static bool printed_sizes = false;
    if (!printed_sizes) {
      fprintf(stderr,
        "lagStruct sizes: rank %d LAG_SHARED_TOPOLOGY %d LAG_REF_GEOMETRY %d sizeof(coord) %zu sizeof(lagNode) %zu sizeof(Edge) %zu"
        #if dimension > 2
          " sizeof(Triangle) %zu"
        #endif
        "\n",
        pid(), LAG_SHARED_TOPOLOGY, LAG_REF_GEOMETRY, sizeof(coord),
        sizeof(lagNode), sizeof(Edge)
        #if dimension > 2
          , sizeof(Triangle)
        #endif
      );
      printed_sizes = true;
    }
    int errors = lag_topology_validate(mesh);
    #if dimension > 2
      fprintf(stderr,
        "lagTopology debug: rank %d context %s cap %d topology %p registry_n %d ref_geometry %p ref_registry_n %d nln %d nle %d nlt %d validation_errors %d\n",
        pid(), context ? context : "none", mesh->cap_id, (void*) mesh->topology,
        lag_topologies.n, (void*) mesh->ref_geometry, lag_ref_geometries.n,
        mesh->nln, mesh->nle, mesh->nlt, errors);
    #else
      fprintf(stderr,
        "lagTopology debug: rank %d context %s cap %d topology %p registry_n %d ref_geometry %p ref_registry_n %d nln %d nle %d validation_errors %d\n",
        pid(), context ? context : "none", mesh->cap_id, (void*) mesh->topology,
        lag_topologies.n, (void*) mesh->ref_geometry, lag_ref_geometries.n,
        mesh->nln, mesh->nle, errors);
    #endif
    if (errors > 20)
      fprintf(stderr, "lagTopology debug: %d additional errors suppressed\n", errors - 20);
  }
#else
  #define debug_lag_topology(mesh, context) ((void) 0)
#endif

bool lag_topology_matches_mesh(lagTopology* topology, lagMesh* mesh) {
  if (!topology || topology->nln != mesh->nln || topology->nle != mesh->nle)
    return false;
  #if dimension > 2
    if (topology->nlt != mesh->nlt)
      return false;
  #endif

  for(int i=0; i<mesh->nln; i++) {
    #if dimension < 3
      for(int j=0; j<2; j++)
        if (topology->node_edge_ids[i][j] != LAG_NODE_EDGE_ID(mesh, i, j))
          return false;
    #else
      if (topology->node_nb_neighbors[i] != LAG_NODE_NB_NEIGHBORS(mesh, i) ||
        topology->node_nb_triangles[i] != LAG_NODE_NB_TRIANGLES(mesh, i))
        return false;
      for(int j=0; j<6; j++) {
        if (topology->node_neighbor_ids[i][j] != LAG_NODE_NEIGHBOR_ID(mesh, i, j) ||
          topology->node_edge_ids[i][j] != LAG_NODE_EDGE_ID(mesh, i, j) ||
          topology->node_triangle_ids[i][j] != LAG_NODE_TRIANGLE_ID(mesh, i, j))
          return false;
      }
    #endif
  }
  for(int i=0; i<mesh->nle; i++) {
    for(int j=0; j<2; j++)
      if (topology->edge_node_ids[i][j] != LAG_EDGE_NODE_ID(mesh, i, j))
        return false;
    #if dimension > 2
      for(int j=0; j<2; j++)
        if (topology->edge_triangle_ids[i][j] != LAG_EDGE_TRIANGLE_ID(mesh, i, j))
          return false;
    #endif
  }
  #if dimension > 2
    for(int i=0; i<mesh->nlt; i++) {
      for(int j=0; j<3; j++) {
        if (topology->triangle_node_ids[i][j] != LAG_TRIANGLE_NODE_ID(mesh, i, j) ||
          topology->triangle_edge_ids[i][j] != LAG_TRIANGLE_EDGE_ID(mesh, i, j))
          return false;
      }
    }
  #endif
  return true;
}

void attach_shared_lag_topology(lagMesh* mesh) {
  #if LAG_SHARED_TOPOLOGY
    assert(mesh->topology);
  #endif
  for(int i=0; i<lag_topologies.n; i++) {
    if (lag_topology_matches_mesh(lag_topologies.items[i], mesh)) {
      if (mesh->topology && mesh->topology != lag_topologies.items[i])
        free_lag_topology(mesh->topology);
      mesh->topology = lag_topologies.items[i];
      return;
    }
  }
  if (lag_topologies.n >= lag_topologies.nm) {
    lag_topologies.nm = lag_topologies.nm ? 2*lag_topologies.nm : 4;
    lag_topologies.items =
      realloc(lag_topologies.items, lag_topologies.nm*sizeof(lagTopology*));
    assert(lag_topologies.items);
  }
  #if !LAG_SHARED_TOPOLOGY
    if (!mesh->topology)
      mesh->topology = copy_lag_topology_from_mesh(mesh);
  #endif
  lag_topologies.items[lag_topologies.n++] = mesh->topology;
}

void attach_shared_lag_ref_geometry(lagMesh* mesh) {
  #if LAG_REF_GEOMETRY
    assert(mesh->ref_geometry);
    for(int i=0; i<lag_ref_geometries.n; i++) {
      if (lag_ref_geometry_matches_mesh(lag_ref_geometries.items[i], mesh)) {
        if (mesh->ref_geometry != lag_ref_geometries.items[i])
          free_lag_ref_geometry(mesh->ref_geometry);
        mesh->ref_geometry = lag_ref_geometries.items[i];
        return;
      }
    }
    if (lag_ref_geometries.n >= lag_ref_geometries.nm) {
      lag_ref_geometries.nm = lag_ref_geometries.nm ? 2*lag_ref_geometries.nm : 4;
      lag_ref_geometries.items = realloc(lag_ref_geometries.items,
        lag_ref_geometries.nm*sizeof(lagReferenceGeometry*));
      assert(lag_ref_geometries.items);
    }
    lag_ref_geometries.items[lag_ref_geometries.n++] = mesh->ref_geometry;
  #else
    (void) mesh;
  #endif
}

void free_shared_lag_topologies() {
  for(int i=0; i<lag_topologies.n; i++)
    free_lag_topology(lag_topologies.items[i]);
  free(lag_topologies.items);
  lag_topologies.n = 0;
  lag_topologies.nm = 0;
  lag_topologies.items = NULL;
}

void free_shared_lag_ref_geometries() {
  for(int i=0; i<lag_ref_geometries.n; i++)
    free_lag_ref_geometry(lag_ref_geometries.items[i]);
  free(lag_ref_geometries.items);
  lag_ref_geometries.n = 0;
  lag_ref_geometries.nm = 0;
  lag_ref_geometries.items = NULL;
}

void free_one_caps(lagMesh* mesh) {
  for(int i=0; i<mesh->nln; i++) free(mesh->nodes[i].stencil.p);
  for(int i=0; i<mesh->nln; i++) free(mesh->nodes[i].eulcell.p);
  free(mesh->nodes);
  free(mesh->edges);
  #if dimension > 2
    free(mesh->triangles);
  #endif
  #if LAG_REF_GEOMETRY
    (void) mesh;
  #else
    free_lag_ref_geometry(mesh->ref_geometry);
  #endif
}

void free_all_caps(Capsules* caps) {
  for(int i=0; i<caps->nbcaps; i++)
    if (CAPS(i).isactive)
      free_one_caps(&(caps->caps[i]));
  free_shared_lag_topologies();
  free_shared_lag_ref_geometries();
}

/**
## Adding capsules after restart

In case we add membranes to a simulation after restart, we need to call
the function below. */
void initialize_capsules() {
  allCaps.nbcaps = NCAPS;
  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive) initialize_empty_capsule(&CAPS(i));
  }
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() a.x[] = 0.;
  }
}

void initialize_capsule_stencils(lagMesh* mesh) {
  for(int j=0; j<mesh->nln; j++) {
    mesh->nodes[j].stencil.n = STENCIL_SIZE;
    mesh->nodes[j].stencil.nm = STENCIL_SIZE;
    mesh->nodes[j].stencil.p = (Index*) malloc(STENCIL_SIZE*sizeof(Index));
    mesh->nodes[j].eulcell.n = 1;
    mesh->nodes[j].eulcell.nm = 1;
    mesh->nodes[j].eulcell.p = (Index*) malloc(sizeof(Index));
  }
  
}

void initialize_all_capsules_stencils() {
  for(int i=0; i<NCAPS; i++)
    if (CAPS(i).isactive)
      initialize_capsule_stencils(&CAPS(i));
}

void initialize_active_capsule(lagMesh* mesh, int cap_id, int cap_type) {
  initialize_empty_capsule(mesh);
  mesh->cap_type = cap_type;
  mesh->cap_id = cap_id;
  mesh->isactive = true;
  initialize_capsule_stencils(mesh);
}

/** The next few macros are useful to compute signed distances and averages
across periodic boundaries. We assume for this purpose that the length of
the edges are less that half the domain size, which in practice should always
be the case. */
#define ACROSS_PERIODIC(a,b,L) (fabs(a - b) > L/2.)
#define PERIODIC_1DIST(a,b,L) (fabs(a - L - b) > L/2. ? a + L - b : a - L - b)
#define GENERAL_1DIST(a,b,L) (ACROSS_PERIODIC(a,b,L) ? PERIODIC_1DIST(a,b,L) : a - b)
#define PERIODIC_1DAVG(a,b,L) (fabs(a - L - b) > L/2. ? a + L + b : a - L + b)
#define GENERAL_1DAVG(a,b,L) (ACROSS_PERIODIC(a,b,L) ? PERIODIC_1DAVG(a,b,L) : a + b)
#define GENERAL_SQNORM(a,b) (sq(GENERAL_1DIST(a.x, b.x, L0)) + \
  sq(GENERAL_1DIST(a.y, b.y, L0*Dimensions.y/Dimensions.x)) + sq(GENERAL_1DIST(a.z, b.z, L0*Dimensions.z/Dimensions.x)))

#if dimension < 3
  #define cnorm(a) (sqrt(sq(a.x) + sq(a.y)))
  #define cdot(a,b) (a.x*b.x + a.y*b.y)
#else
  #define cnorm(a) (sqrt(sq(a.x) + sq(a.y) + sq(a.z)))
  #define cdot(a,b) (a.x*b.x + a.y*b.y + a.z*b.z)
#endif

/**
## Advection of the mesh

The advection have to be compatible with MPI, so we start by including the MPI
communication functions. We also include [reg-dirac.h](reg-dirac.h), which
implements the interpolation of the Eulerian velocities onto the Lagrangian
nodes.
*/

/** By default, the mesh is advected using a second-order two-step Runge Kutta
scheme. If the following macro is set to 0, a first-order forward Euler schme
is used instead. */
#ifndef ADVECT_LAG_RK2
  #define ADVECT_LAG_RK2 1
#endif

#if _MPI
  #include "capsule-ft-mpi.h"
#endif
#include "ibm-ft.h"
#ifndef CONSERVE_VOLUME
  #define CONSERVE_VOLUME 1
#endif

/**
A lot of triangulated mesh-related functions are defined in a separate header
file.
*/
#include "geometry-ft.h"

#if CONSERVE_VOLUME
  #include "volume-conservation-ft.h"
#endif

/* Utilities for different numerical simulations */
#include "plugins-ft.h"



//----------------------------------------------------------------------------
trace void synchronize (scalar * list)
//----------------------------------------------------------------------------
{
  for (scalar s in list)
    s.dirty = true;
  boundary(list);
}


/**
The function below advects each Lagrangian node by
interpolating the velocities around the node of interest. By default, a
second-order Runge Kutta scheme is used. By setting the macro
```ADVECT_LAG_RK2``` to 0, a simple forward Euler scheme is used.
*/
trace
void advect_lagMesh(lagMesh* mesh) {

  #if !(ADVECT_LAG_RK2)
    for(int i=0; i < mesh->nln; i++) {
      foreach_dimension() {
        mesh->nodes[i].pos.x += dt*mesh->nodes[i].lagVel.x;
      }
    }
  #else
    lagMesh buffer_mesh;
    buffer_mesh.isactive = true;
    buffer_mesh.nln = mesh->nln;
    buffer_mesh.nodes = malloc(mesh->nln*sizeof(lagNode));
    for(int i=0; i<mesh->nln; i++) {
      // Step 1 of RK2
      foreach_dimension()
        buffer_mesh.nodes[i].pos.x = mesh->nodes[i].pos.x +
          .5*dt*mesh->nodes[i].lagVel.x;
    }
    correct_lag_pos(&buffer_mesh);
    for(int j=0; j<buffer_mesh.nln; j++) {
      buffer_mesh.nodes[j].stencil.n = STENCIL_SIZE;
      buffer_mesh.nodes[j].stencil.nm = STENCIL_SIZE;
      buffer_mesh.nodes[j].stencil.p = malloc(STENCIL_SIZE*sizeof(Index));
      buffer_mesh.nodes[j].eulcell.n = 1;
      buffer_mesh.nodes[j].eulcell.nm = 1;
      buffer_mesh.nodes[j].eulcell.p = malloc(sizeof(Index));
    }
    
    generate_lag_stencils_one_caps(&buffer_mesh, true);
    eul2lag(&buffer_mesh);
    for(int i=0; i<mesh->nln; i++) {
      // Step 2 of RK2
      foreach_dimension()
        mesh->nodes[i].pos.x += dt*buffer_mesh.nodes[i].lagVel.x;
    }
    for(int i=0; i<buffer_mesh.nln; i++) free(buffer_mesh.nodes[i].stencil.p);
    for(int i=0; i<buffer_mesh.nln; i++) free(buffer_mesh.nodes[i].eulcell.p);
    free(buffer_mesh.nodes);
  #endif

  correct_lag_pos(mesh);
  #if CONSERVE_VOLUME
    enforce_optimal_volume_conservation(mesh);
  #endif
  comp_centroid(mesh);
  comp_volume(mesh);
  comp_capsule_geodynamics(mesh);
}





/**
## Putting the pieces together

Below, we call the above functions at the appropriate time using the Basilisk
event syntax.

We start by creating empty Lagrangian meshes, and allocating an acceleration
field in case it isn't done yet by another Basilisk solver.
*/
event defaults (i = 0) {
  allCaps.nbcaps = NCAPS;
  // #if (RESTART_CASE == 0)
    for(int i=0; i<NCAPS; i++) {
      initialize_empty_capsule(&CAPS(i));
    }
  // #endif
  if (is_constant(a.x)) {
    a = new face vector;
    foreach_face() a.x[] = 0.;
  }
}



/*Repulsive lubrication nodal force*/
void repulsive_vel() 
{
  /*Compute the cell size in the grid*/
  #if MULT_GRID == 1   
    double delta = (L0/(1 << grid->maxdepth)/Dimensions.x);
  #else
    double delta = (L0/(1 << grid->maxdepth));
  #endif

  
  /*The value of K_lub is up to the */
  // double K_lub = 0.001/(E_S);

  for(int i = 0; i < NCAPS; i++) {
    if (CAPS(i).isactive) 
    {
      lagMesh* mesh = &(CAPS(i));

      for(int j=0; j<mesh->nln; j++) 
      { 
        foreach_cache(mesh->nodes[j].eulcell)
        {
          // int lagnode_id = (int)Index_lag_id.x[];
          int lagnode_id = j;
          coord lub_vel = {0};  
          double K_lub = 0.;
          foreach_dimension()
            K_lub += sq(mesh->nodes[lagnode_id].lagVel.x);
          K_lub = sqrt(K_lub)*0.25; //half the velocity as coefficient

          if(point.level>-1)
          {        
              coord lagpt = {0};
              lagpt.x = mesh->nodes[lagnode_id].pos.x;
              lagpt.y = mesh->nodes[lagnode_id].pos.y;
              lagpt.z = mesh->nodes[lagnode_id].pos.z;
                 
              foreach_neighbor()
              {
                if(point.level >-1)
                {
                  if(((int)Index_lagnode[] > -1) && ((mesh->cap_id) != (int)Index_lagnode[])) 
                  {        
                    coord checkpt = {0};
                    checkpt.x = CAPS((int)Index_lagnode[]).nodes[(int)Index_lag_id.x[]].pos.x;
                    checkpt.y = CAPS((int)Index_lagnode[]).nodes[(int)Index_lag_id.x[]].pos.y;
                    checkpt.z = CAPS((int)Index_lagnode[]).nodes[(int)Index_lag_id.x[]].pos.z;

                    coord lub_dir = {0};
                    
                    double lub_norm = sqrt(GENERAL_SQNORM(lagpt, checkpt));
            
                    foreach_dimension() lub_dir.x = GENERAL_1DIST(lagpt.x, checkpt.x, L0*L0_ratio.x)/lub_norm;
                    if(lub_norm < 2.*delta)
                    {
                      foreach_dimension() lub_vel.x += lub_dir.x * K_lub * (sq(2.*delta/lub_norm) - 1.);
                    }
                  }
                }
              }
            /** The velocity of the node is adjusted as the repulsive effect applies. */
            foreach_dimension() mesh->nodes[lagnode_id].lagVel.x += 0.5*lub_vel.x; 
      
            
            /*A special case where the two nodes from different caps lie in the same cell, we push the node of the other capsule as well */
            if(((int)Index_lagnode[] > -1) && ((int)Index_lag_id.y[] > -1) && ((mesh->cap_id) != (int)Index_lag_id.y[]))
            {
              coord checkpt = {0};
              checkpt.x = CAPS((int)Index_lag_id.y[]).nodes[(int)Index_lag_id.z[]].pos.x;
              checkpt.y = CAPS((int)Index_lag_id.y[]).nodes[(int)Index_lag_id.z[]].pos.y;
              checkpt.z = CAPS((int)Index_lag_id.y[]).nodes[(int)Index_lag_id.z[]].pos.z;

              coord lub_dir = {0};
              double lub_norm = sqrt(GENERAL_SQNORM(lagpt, checkpt));
              foreach_dimension() lub_dir.x = GENERAL_1DIST(lagpt.x, checkpt.x, L0*L0_ratio.x)/lub_norm;
              if(lub_norm < 2.*delta)
              {                 
                CAPS((int)Index_lag_id.y[]).nodes[(int)Index_lag_id.z[]].lagVel.x -= 0.5*lub_dir.x * K_lub * (sq(2.*delta/lub_norm) - 1.);
                CAPS((int)Index_lag_id.y[]).nodes[(int)Index_lag_id.z[]].lagVel.y -= 0.5*lub_dir.y * K_lub * (sq(2.*delta/lub_norm) - 1.);
                CAPS((int)Index_lag_id.y[]).nodes[(int)Index_lag_id.z[]].lagVel.z -= 0.5*lub_dir.z * K_lub * (sq(2.*delta/lub_norm) - 1.);
              }
            }
            
          }
        }
      }
    }
  }
}


/** Below, we advect each Lagrangian node using the interpolated Eulerian
velocities. We also use this loop as an opportunity to
re-initialize the Lagrangian forces to zero. */

coord proc_max = {-HUGE, -HUGE, -HUGE};
coord proc_min = {HUGE, HUGE, HUGE};
#if _MPI
coord* all_proc_max = NULL;
coord* all_proc_min = NULL;
int* ncaps_for_proc = NULL;
int* proc_cap_offsets = NULL;
int* proc_cap_ids = NULL;
int proc_cap_ids_nm = 0;
#endif

#ifndef DEBUG_AABB
  #define DEBUG_AABB 0
#endif
#ifndef DEBUG_AABB_FREQ
  #define DEBUG_AABB_FREQ 1
#endif

event tracer_advection(i++) {  

  /* Distribute velocity to the lagNodes */
  for(int i=0; i<NCAPS; i++) 
  {
      if (CAPS(i).isactive) 
        eul2lag(&CAPS(i));
  }   

 
  /*We synchronize the eul field and make sure that it is updated before applying repulsive velocity */
  
#if (LUBR_VEL==1)
  synchronize({Index_lagnode, Index_lag_id}); 
  repulsive_vel();
#endif 
  /**
  In case of parallel simulations, we communicate the Lagrangian velocity
  so that all processes have the same Lagrangian velocities.
  */
  reduce_alllagVel();


  /* Advection of the lagNode */
  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive) {
      advect_lagMesh(&CAPS(i));
      for(int j=0; j<CAPS(i).nln; j++)
        foreach_dimension() CAPS(i).nodes[j].lagForce.x = 0.;
    }
  }

  /* Compute borders of the curren proc */
  compute_proc_borders(&proc_max, &proc_min);
  #if DEBUG_AABB
    if (i % DEBUG_AABB_FREQ == 0) {
      #if _MPI
        if (all_proc_min == NULL) {
          all_proc_min = (coord*)malloc(npe()*sizeof(coord));
          all_proc_max = (coord*)malloc(npe()*sizeof(coord));
          ncaps_for_proc = (int*)malloc(npe()*sizeof(int));
          proc_cap_offsets = (int*)malloc((npe() + 1)*sizeof(int));
        }
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
      #endif
      fprintf(stderr,
        "DEBUG_AABB pid %d/%d iter %d proc_min=(%g %g %g) proc_max=(%g %g %g)\n",
        pid(), npe(), i,
        proc_min.x, proc_min.y, proc_min.z,
        proc_max.x, proc_max.y, proc_max.z);
      #if _MPI
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

  /*Clean the index field before generating the stencils*/
  foreach()
  {
    if (cm[] > 1.e-20) 
    { Index_lagnode[] = -1;
      foreach_dimension() Index_lag_id.x[] = -1;
    }
  }

  /* Generate new stencils in corresponding procs */
  for(int cap=0; cap<NCAPS; cap++) {
    if (CAPS(cap).isactive) {
      bool intersects_proc = is_capsule_in_boundingbox(proc_max, proc_min, &CAPS(cap));
      #if DEBUG_AABB
        if (i % DEBUG_AABB_FREQ == 0)
          fprintf(stderr,
            "DEBUG_AABB pid %d/%d iter %d cap %d centroid=(%g %g %g) circum_radius=%g intersects_proc=%d\n",
            pid(), npe(), i, cap,
            CAPS(cap).centroid.x, CAPS(cap).centroid.y, CAPS(cap).centroid.z,
            CAPS(cap).circum_radius, intersects_proc);
      #endif
      if(intersects_proc) 
        generate_lag_stencils_one_caps(&CAPS(cap), true);
    }
  }

}


/*Repulsive lubrication nodal force*/
void lubrication_force() 
{
  /*Compute the cell size in the grid*/
  #if MULT_GRID == 1   
    double delta = (L0/(1 << grid->maxdepth)/Dimensions.x);
    
  #else
    double delta = (L0/(1 << grid->maxdepth));
  #endif

  /*The value of K_lub is up to the */
  // double K_lub = 0.001/(E_S);

  for(int i = 0; i < NCAPS; i++) {
    if (CAPS(i).isactive) 
    {
      lagMesh* mesh = &(CAPS(i));

      for(int j=0; j<mesh->nln; j++) 
      { 
        foreach_cache(mesh->nodes[j].eulcell)
        {
          // int lagnode_id = (int)Index_lag_id.x[];
          int lagnode_id = j;
          coord lub_force = {0};  
          double K_lub = 0.;
          foreach_dimension()
            K_lub += sq(mesh->nodes[lagnode_id].lagForce.x);
          K_lub = 2.*sqrt(K_lub);

          if(point.level>-1)
          {        
              coord lagpt = {0};
              lagpt.x = mesh->nodes[lagnode_id].pos.x;
              lagpt.y = mesh->nodes[lagnode_id].pos.y;
              lagpt.z = mesh->nodes[lagnode_id].pos.z;
                 
              foreach_neighbor()
              {
                if(point.level >-1)
                {
                  if(((int)Index_lagnode[] > -1) && ((mesh->cap_id) != (int)Index_lagnode[])) 
                  {        
                    coord checkpt = {0};
                    checkpt.x = CAPS((int)Index_lagnode[]).nodes[(int)Index_lag_id.x[]].pos.x;
                    checkpt.y = CAPS((int)Index_lagnode[]).nodes[(int)Index_lag_id.x[]].pos.y;
                    checkpt.z = CAPS((int)Index_lagnode[]).nodes[(int)Index_lag_id.x[]].pos.z;

                    coord lub_dir = {0};
                    
                    double lub_norm = sqrt(GENERAL_SQNORM(lagpt, checkpt));
                    foreach_dimension() lub_dir.x = GENERAL_1DIST(lagpt.x, checkpt.x, L0*L0_ratio.x)/lub_norm;
                    if(lub_norm < 2*delta)
                    {
                      foreach_dimension() lub_force.x += lub_dir.x * K_lub * (sq(2*delta/lub_norm) - 1.);
                    }
                  }  
                }
              }
            /** The lubrication force is ready to be added to the Lagrangian force of the considered node. */
            foreach_dimension() mesh->nodes[lagnode_id].lagForce.x += lub_force.x;   
          }
        }
      }
    }
  }
}


/** In the acceleration event, we transfer the Lagrangian forces to the fluid
using a regularized Dirac function. The acceleration is stored on the cell
faces, and will be fed as a source term to the Navier-Stokes solver. */
vector forcing[];
event acceleration (i++) {

  /*We synchronize the eul field and make sure that it is updated before applying repulsive force */
  synchronize({Index_lagnode, Index_lag_id});

  /*We add the repulsive lubrication force for a better numerical stability*/
  # if LUBR_FORCE == 1  
  // lubrication_force(); 
  # endif

  face vector ae = a;
  foreach()
    if (cm[] > 1.e-20) foreach_dimension() forcing.x[] = 0.;
  for(int i=0; i<NCAPS; i++) {
    if (CAPS(i).isactive) lag2eul(forcing, &CAPS(i));
  }
  foreach_face()
    if (fm.x[] > 1.e-20) ae.x[] += .5*alpha.x[]*(forcing.x[] + forcing.x[-1]);
}

/** At the end of the simulation, we free the allocated memory.*/
event cleanup (t = end) {
  #if _MPI
    free(all_proc_min);
    free(all_proc_max);
    free(ncaps_for_proc);
    free(proc_cap_offsets);
    free(proc_cap_ids);
  #endif
  free_all_caps(&allCaps);
}


/**
## Additional functionalities
*/
#if dimension > 2
  #include "dump-ft.h"
  #include "post-processing-ft.h"
#endif


/**
## Tests
[advect_caps.c](../../tests/lagrangian_caps/advect_caps.c): Tests the
convergence of the advection scheme.


[curvature.c](../../tests/lagrangian_caps/curvature.c): Tests the computation of
the curvature at the Lagrangian nodes. Since the curvature depends on the
normals, this case also validates the computation of the normal vectors.
*/
