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
  #if _MPI && DEBUG_AABB
    if (i % DEBUG_AABB_FREQ == 0) {
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
  #endif
  /**
  In case of parallel simulations, we communicate the Lagrangian velocity
  so that all processes have the same Lagrangian velocities.
  */
  reduce_alllagVel();

  #if _MPI && DEBUG_AABB
    if (i % DEBUG_AABB_FREQ == 0) {
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
  #endif

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

        if (total_ghost_to_owner_recv_caps > 0) {
          fprintf(stderr,
            "DEBUG_OWNER_LAGVEL_ACCUM_DRYRUN pid %d/%d iter %d caps=",
            pid(), npe(), i);
          for(int cap=0; cap<NCAPS; cap++) {
            if (!CAPS(cap).isactive)
              continue;
            int owner_proc = find_capsule_owner_proc(&CAPS(cap),
              all_proc_min, all_proc_max);
            if (owner_proc != pid())
              continue;

            coord* lagvel_sum = (coord*)calloc(CAPS(cap).nln, sizeof(coord));
            assert(lagvel_sum);
            coord accum_vel_min = {HUGE, HUGE, HUGE};
            coord accum_vel_max = {-HUGE, -HUGE, -HUGE};
            double base_lagvel_abs_sum = 0.;
            coord* base_lagVel =
              debug_pre_reduce_lagVel &&
              debug_pre_reduce_lagVel_nln[cap] == CAPS(cap).nln ?
              debug_pre_reduce_lagVel[cap] : NULL;
            for(int node_id=0; node_id<CAPS(cap).nln; node_id++) {
              foreach_dimension() {
                lagvel_sum[node_id].x = base_lagVel ?
                  base_lagVel[node_id].x : CAPS(cap).nodes[node_id].lagVel.x;
                base_lagvel_abs_sum += fabs(lagvel_sum[node_id].x);
              }
            }

            int recv_caps_added = 0;
            int recv_sources_added = 0;
            for(int p=0; p<npe(); p++) {
              int int_pos = ghost_to_owner_recv_int_offsets[p];
              int double_pos = ghost_to_owner_recv_double_offsets[p];
              for(int q=0; q<ghost_to_owner_recv_caps[p]; q++) {
                int recv_cap_id = ghost_to_owner_recv_int_buffer[int_pos++];
                int recv_nln = ghost_to_owner_recv_int_buffer[int_pos++];
                bool add_payload = recv_cap_id == CAPS(cap).cap_id &&
                  recv_nln == CAPS(cap).nln;
                if (add_payload) {
                  recv_caps_added++;
                  recv_sources_added += p != pid();
                }
                for(int node_id=0; node_id<recv_nln; node_id++) {
                  foreach_dimension() {
                    double vel_comp =
                      ghost_to_owner_recv_double_buffer[double_pos++];
                    if (add_payload)
                      lagvel_sum[node_id].x += vel_comp;
                  }
                }
              }
            }

            double accum_lagvel_abs_sum = 0.;
            double reduced_lagvel_abs_sum = 0.;
            double reduced_vs_accum_max_abs_diff = 0.;
            double owner_euler_move_max_abs_diff = 0.;
            coord* pre_advect_pos =
              debug_pre_advect_pos &&
              debug_pre_advect_pos_nln[cap] == CAPS(cap).nln ?
              debug_pre_advect_pos[cap] : NULL;
            for(int node_id=0; node_id<CAPS(cap).nln; node_id++) {
              foreach_dimension() {
                accum_vel_min.x = min(accum_vel_min.x, lagvel_sum[node_id].x);
                accum_vel_max.x = max(accum_vel_max.x, lagvel_sum[node_id].x);
                accum_lagvel_abs_sum += fabs(lagvel_sum[node_id].x);
                reduced_lagvel_abs_sum +=
                  fabs(CAPS(cap).nodes[node_id].lagVel.x);
                reduced_vs_accum_max_abs_diff = max(
                  reduced_vs_accum_max_abs_diff,
                  fabs(CAPS(cap).nodes[node_id].lagVel.x -
                    lagvel_sum[node_id].x));
                if (pre_advect_pos) {
                  double predicted_pos =
                    pre_advect_pos[node_id].x + dt*lagvel_sum[node_id].x;
                  owner_euler_move_max_abs_diff = max(
                    owner_euler_move_max_abs_diff,
                    fabs(CAPS(cap).nodes[node_id].pos.x - predicted_pos));
                }
              }
            }
            fprintf(stderr,
              " cap=%d,owner=%d,nln=%d,base_abs_sum=%g,recv_caps_added=%d,recv_sources_added=%d,accum_vel_min=(%g %g %g),accum_vel_max=(%g %g %g),accum_abs_sum=%g,reduced_abs_sum=%g,reduced_vs_accum_max_abs_diff=%g,owner_euler_move_max_abs_diff=%g",
              CAPS(cap).cap_id, owner_proc, CAPS(cap).nln,
              base_lagvel_abs_sum, recv_caps_added, recv_sources_added,
              accum_vel_min.x, accum_vel_min.y, accum_vel_min.z,
              accum_vel_max.x, accum_vel_max.y, accum_vel_max.z,
              accum_lagvel_abs_sum, reduced_lagvel_abs_sum,
              reduced_vs_accum_max_abs_diff,
              owner_euler_move_max_abs_diff);
            free(lagvel_sum);
          }
          fprintf(stderr, "\n");
        }

        int local_owner_accum_ints[5] = {-1, -1, 0, 0, 0};
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
            double base_lagvel_abs_sum = 0.;
            coord* base_lagVel =
              debug_pre_reduce_lagVel &&
              debug_pre_reduce_lagVel_nln[cap] == CAPS(cap).nln ?
              debug_pre_reduce_lagVel[cap] : NULL;
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
                bool add_payload = recv_cap_id == CAPS(cap).cap_id &&
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
              debug_pre_advect_pos &&
              debug_pre_advect_pos_nln[cap] == CAPS(cap).nln ?
              debug_pre_advect_pos[cap] : NULL;
            for(int node_id=0; node_id<CAPS(cap).nln; node_id++)
              foreach_dimension() {
                accum_vel_min.x = min(accum_vel_min.x, lagvel_sum[node_id].x);
                accum_vel_max.x = max(accum_vel_max.x, lagvel_sum[node_id].x);
                accum_lagvel_abs_sum += fabs(lagvel_sum[node_id].x);
                reduced_lagvel_abs_sum +=
                  fabs(CAPS(cap).nodes[node_id].lagVel.x);
                reduced_vs_accum_max_abs_diff = max(
                  reduced_vs_accum_max_abs_diff,
                  fabs(CAPS(cap).nodes[node_id].lagVel.x -
                    lagvel_sum[node_id].x));
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
          all_owner_accum_ints = (int*)malloc(npe()*5*sizeof(int));
          all_owner_accum_doubles = (double*)malloc(npe()*11*sizeof(double));
          assert(all_owner_accum_ints);
          assert(all_owner_accum_doubles);
        }
        MPI_Gather(local_owner_accum_ints, 5, MPI_INT,
          all_owner_accum_ints, 5, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(local_owner_accum_doubles, 11, MPI_DOUBLE,
          all_owner_accum_doubles, 11, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (pid() == 0) {
          int printed_owner_accum = false;
          for(int rank=0; rank<npe(); rank++) {
            int* row_i = all_owner_accum_ints + rank*5;
            double* row_d = all_owner_accum_doubles + rank*11;
            if (row_i[0] < 0)
              continue;
            if (!printed_owner_accum) {
              fprintf(stderr,
                "DEBUG_ALL_OWNER_LAGVEL_ACCUM_DRYRUN iter %d", i);
              printed_owner_accum = true;
            }
            fprintf(stderr,
              " rank=%d,cap=%d,owner=%d,nln=%d,base_abs_sum=%g,recv_caps_added=%d,recv_sources_added=%d,accum_vel_min=(%g %g %g),accum_vel_max=(%g %g %g),accum_abs_sum=%g,reduced_abs_sum=%g,reduced_vs_accum_max_abs_diff=%g,owner_euler_move_max_abs_diff=%g",
              rank, row_i[0], row_i[1], row_i[2], row_d[0],
              row_i[3], row_i[4], row_d[1], row_d[2], row_d[3],
              row_d[4], row_d[5], row_d[6], row_d[7], row_d[8],
              row_d[9], row_d[10]);
          }
          if (printed_owner_accum)
            fprintf(stderr, "\n");
        }
        free(all_owner_accum_ints);
        free(all_owner_accum_doubles);

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
