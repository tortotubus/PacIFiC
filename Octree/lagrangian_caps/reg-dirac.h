/**
# Definition of the IBM regularized Dirac

The Lagrangian and Eulerian meshes communicate thanks to regularized Dirac-delta
functions: this is the core of the Immersed Boundary Method (IBM) as introduced
by [Peskin](Peskin1997).

This implementation of the IBM in Basilisk relies on the Cache structure to
efficiently loop through the Eulerian cells in the vicinity of the Lagrangian
nodes. Since the Cache structure is only defined for trees, the current
implementation is \textbf{only compatible with quadtree and octree grids, and
not with Cartesian nor multigrids}.
*/

#define BGHOSTS 2 // Having two layers of ghost cells should not be mandatory since the implementation relies on Caches

#define POS_PBC_X(X) ((u.x.boundary[left] != periodic_bc) ? (X) : (((X) > L0/2.) ? (X) - L0 : (X)))
#define POS_PBC_Y(Y) ((u.x.boundary[top] != periodic_bc) ? (Y) : (((Y) > L0/2.) ? (Y) - L0 : (Y)))
#define POS_PBC_Z(Z) ((u.x.boundary[top] != periodic_bc) ? (Z) : (((Z) > L0/2.) ? (Z) - L0 : (Z)))

struct _locate_lvl {int lvl; double x, y, z;};

Point locate_lvl(struct _locate_lvl p) {
  for (int l = p.lvl; l >= 0; l--) {
    Point point = {0};
    point.level = l;
    int n = 1 << point.level;
    point.i = (p.x - X0)/L0*n + GHOSTS;
    #if dimension >= 2
      point.j = (p.y - Y0)/L0*n + GHOSTS;
    #endif
    #if dimension >= 3
      point.k = (p.z - Z0)/L0*n + GHOSTS;
    #endif
      if (point.i >= 0 && point.i < n + 2*GHOSTS
    #if dimension >= 2
    && point.j >= 0 && point.j < n + 2*GHOSTS
    #endif
    #if dimension >= 3
    && point.k >= 0 && point.k < n + 2*GHOSTS
    #endif
    ) {
        if (allocated(0))
    return point;
      }
      else
        break;
  }
  Point point = {0};
  point.level = -1;
  return point;
}

/**
*/
#if CHANGING_STENCIL_LEVEL
  int get_level_IBM_stencil_one_step(lagNode* node) {
    int lvl = node->stenstat.slvl;
    node->stenstat.complete_stencil = true;
    double delta = (L0/(1 << lvl));
    for(int ni=-2; ni<=2; ni++) {
      for(int nj=-2; nj<=2; nj++) {
        #if dimension < 3
        Point point = locate(POS_PBC_X(node->pos.x + ni*delta),
          POS_PBC_Y(node->pos.y + nj*delta));
        #else
        /** In 3 dimensions, we implement a small acceleration: if the stencil
        is incomplete, there is no need to span the rest of it. */
        for(int nk=-2; nk<=2 && node->stenstat.complete_stencil; nk++) {
          Point point = locate(POS_PBC_X(node->pos.x + ni*delta),
            POS_PBC_Y(node->pos.y + nj*delta),
            POS_PBC_Z(node->pos.z + nk*delta));
        #endif
        if (point.level < lvl) {
          node->stenstat.slvl = point.level;
          node->stenstat.complete_stencil = false;
        }
      #if dimension > 2
        }
      #endif
      }
    }
    return node->stenstat.slvl;
  }
#endif

#if (_MPI && CHANGING_STENCIL_LEVEL)
  MPI_Op mpi_min_slvl;
  MPI_Datatype mpi_stenctil_status;
  stencil_status* stenstats_local;
  stencil_status* stenstats_global;

  /** The function below defines a custom MPI reduction operation for the
  stencil status of a node:
    * ```slvl``` is reduced to its minimum value
    * ```complete_stencil``` is reduced to the opposite of the logical OR */
  void min_slvl(void *invec, void *inoutvec, int *len, MPI_Datatype *dtype) {
    /** rl is the already reduced stencil level, nl is the new stencil level */
    int nl = ((stencil_status*)(invec))->slvl;
    int rl = ((stencil_status*)(inoutvec))->slvl;
    int ml;
    if (nl < rl) ml = nl;
    else ml = rl;
    ((stencil_status*)(inoutvec))->slvl = ml;

    /** rs is the already reduced stencil completeness status, ns is the new stencil completeness status */
    int ns = ((stencil_status*)(invec))->complete_stencil;
    int rs = ((stencil_status*)(inoutvec))->complete_stencil;
    int ms;
    if (!ns) ms = ns;
    else ms = rs;
    ((stencil_status*)(inoutvec))->complete_stencil = ms;
  }

  event defaults (i = 0) {
    int sizes[2] = {sizeof(int), sizeof(bool)};
    MPI_Aint offsets[2] = {0, sizeof(int)};
    MPI_Datatype types[2] = {MPI_INT, MPI_C_BOOL};
    MPI_Type_create_struct(1, sizes, offsets, types, &mpi_stenctil_status);
    MPI_Type_commit(&mpi_stenctil_status);
    MPI_Op_create(min_slvl, 1, &mpi_min_slvl);
    stenstats_local = malloc(total_nb_nodes(&mbs)*
      sizeof(stencil_status));
    stenstats_global = malloc(total_nb_nodes(&mbs)*
      sizeof(stencil_status));
  }

  event cleanup (t = end) {
    MPI_Type_free(&mpi_stenctil_status);
    MPI_Op_free(&mpi_min_slvl);
    free(stenstats_local);
    free(stenstats_global);
  }
#endif

/** In case of changing stencil levels, the levels of the stencils are
determined by the function below, which requires MPI communications for
parallel simulations (as many as the level difference between the minimum and
maximum levels in the membrane) */
#if CHANGING_STENCIL_LEVEL
  void get_level_IBM_stencils(Capsules* caps) {
    for(int k=0; k<caps->nbmb; k++) {
      for(int i=0; i<caps->mb[k].nlp; i++) {
        caps->mb[k].nodes[i].stenstat.slvl = MAX_STENCIL_LEVEL;
        caps->mb[k].nodes[i].stenstat.complete_stencil = false;
      }
    }
    bool complete_stencils = false;
    while (!complete_stencils) {
      complete_stencils = true;
      for(int k=0; k<caps->nbmb; k++) {
        for(int i=0; i<caps->mb[k].nlp; i++) {
          // Get local stencil level
          if (!(caps->mb[k].nodes[i].stenstat.complete_stencil))
            get_level_IBM_stencil_one_step(&(caps->mb[k].nodes[i]));
        }
      }
      /** In case of parallel simulations, we communicate the stencils' status,
        i.e. their level and their completeness. */
      #if _MPI
        int cni = 0;  // cni for "current node index"
        for(int k=0; k<caps->nbmb; k++) {
          for(int i=0; i<caps->mb[k].nlp; i++) {
            stenstats_local[cni].slvl = caps->mb[k].nodes[i].stenstat.slvl;
            stenstats_local[cni].complete_stencil =
              caps->mb[k].nodes[i].stenstat.complete_stencil;
            cni++;
          }
        }
        /** We now use our custom reduce datatype and function to find the
        lowest level and completeness of each stencil. */
        MPI_Allreduce(&stenstats_local, &stenstats_global, cni,
          mpi_stenctil_status, mpi_min_slvl, MPI_COMM_WORLD);

        /** The stencil status of each node is populated with the reduced
        information. */
        cni = 0;
        for(int k=0; k<caps->nbmb; k++) {
          for(int i=0; i<caps->mb[k].nlp; i++) {
            caps->mb[k].nodes[i].stenstat.slvl = stenstats_global[cni].slvl;
            caps->mb[k].nodes[i].stenstat.complete_stencil =
              stenstats_global[cni].complete_stencil;
            if (!stenstats_global[cni].complete_stencil)
              complete_stencils = false;
            cni++;
          }
        }
      #endif
    }
  }
#endif


/**
The function below loops through the Lagrangian nodes and "caches" the Eulerian
cells in a 5x5(x5) stencil around each node-> In case of parallel simulations,
the cached cells are tagged with the process id.
*/
scalar stencils[];
trace
void generate_lag_stencils_one_caps(lagMesh* mesh) {
  for(int i=0; i<mesh->nlp; i++) {
    mesh->nodes[i].stencil.n = 0;
    #if CONSTANT_STENCIL_LEVEL
      int lvl = MB_LEVEL;
    #else
      int lvl = mesh->nodes[i].stenstat.slvl;
    #endif
    double delta = L0/(1 << lvl);
    for(int ni=-2; ni<=2; ni++) {
      for(int nj=-2; nj<=2; nj++) {
        #if dimension < 3
          Point point = locate_lvl(lvl,
            POS_PBC_X(mesh->nodes[i].pos.x + ni*delta),
            POS_PBC_Y(mesh->nodes[i].pos.y + nj*delta));
        #else
          for(int nk=-2; nk<=2; nk++) {
            Point point = locate_lvl(lvl,
              POS_PBC_X(mesh->nodes[i].pos.x + ni*delta),
              POS_PBC_Y(mesh->nodes[i].pos.y + nj*delta),
              POS_PBC_Z(mesh->nodes[i].pos.z + nk*delta));
        #endif
        #if CONSTANT_STENCIL_LEVEL
          if (point.level >= 0 && point.level < MB_LEVEL)
        #else
          if (point.level >= 0 && point.level < MIN_STENCIL_LEVEL)
        #endif
          fprintf(stderr,
            "Warning: Lagrangian stencil not fully resolved.\n");
        cache_append(&(mesh->nodes[i].stencil), point, 0);
        #if _MPI
          #if dimension < 3
            if (ni == 0 && nj == 0) {
          #else
            if (ni == 0 && nj == 0 && nk == 0) {
          #endif
            if (point.level >= 0) mesh->nodes[i].pid = cell.pid;
            else mesh->nodes[i].pid = -1;
          }
        #endif
        #if dimension == 3
          }
        #endif
      }
    }
  }
}

trace
void generate_lag_stencils() {
  #if CHANGING_STENCIL_LEVEL
    get_level_IBM_stencils(&mbs);
  #endif
  for(int k=0; k<NCAPS; k++) generate_lag_stencils_one_caps(&MB(k));
}


/**
The two functions below are useful in case the cells of a given stencil are not
at the same level. The first function averages the values of a vector v into a
coord a, while the second function stores the content of a coord a into a
vector v.
*/
void read_vector_leaves(Point point, vector v, coord* a, int depth) {
  if (is_local(cell)) {
    if (is_leaf(cell))
      foreach_dimension() a->x += v.x[]/(1 << dimension*depth);
    else
      foreach_child() read_vector_leaves(point, v, a, depth + 1);
  }
}

/** Unlike the function above, we don't need to divide the values associated
to children cells since this correction is cancelled by the smaller delta in
those cells. */
void populate_vector_leaves(Point point, vector v, coord* a) {
  if (is_local(cell)) {
    if (is_leaf(cell))
      foreach_dimension() v.x[] += a->x;
    else
      foreach_child() populate_vector_leaves(point, v, a);
  }
}

/**
The function below fills the scalar $forcing$ with the regularized source term
of magnitude $force$ located at $x0$. This assumes that $forcing$ is "clean",
i.e. that it is zero everywhere (or that if there is already some non-zero terms
the intention is to include them in the forcing).
*/
trace
void lag2eul(vector forcing, lagMesh* mesh) {
  for(int i=0; i<mesh->nlp; i++) {
    double sdelta; // sdelta for "stencil delta"
    #if CONSTANT_STENCIL_LEVEL
      sdelta = L0/(1 << MB_LEVEL);
    #else
      sdelta = L0/(1 << mesh->nodes[i].stenstat.slvl);
    #endif
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[i].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[i].pos.y);
        #if dimension > 2
          dist.z = GENERAL_1DIST(z, mesh->nodes[i].pos.z);
        #endif
        #if dimension < 3
          if (fabs(dist.x) <= 2*sdelta && fabs(dist.y) <= 2*sdelta) {
            double weight =
              (1 + cos(.5*pi*dist.x/sdelta))*(1 + cos(.5*pi*dist.y/sdelta))
              /sq(4.*sdelta);
        #else
          if (fabs(dist.x) <= 2*sdelta && fabs(dist.y) <= 2*sdelta &&
            fabs(dist.z) <= 2*sdelta) {
            double weight =
              (1 + cos(.5*pi*dist.x/sdelta))*(1 + cos(.5*pi*dist.y/sdelta))
              *(1 + cos(.5*pi*dist.z/sdelta))/cube(4.*sdelta);
        #endif
          coord weighted_forcing;
          foreach_dimension()
            weighted_forcing.x = weight*mesh->nodes[i].lagForce.x;
          populate_vector_leaves(point, forcing, &weighted_forcing);
        }
      }
    }
  }
}

/**
The function below interpolates the eulerian velocities onto the nodes of
the Lagrangian mesh.
*/
trace
void eul2lag(lagMesh* mesh) {
  for(int i=0; i<mesh->nlp; i++) {
    double sdelta; // sdelta for "stencil delta"
    #if CONSTANT_STENCIL_LEVEL
      sdelta = L0/(1 << MB_LEVEL);
    #else
      sdelta = L0/(1 << mesh->nodes[i].stenstat.slvl);
    #endif
    foreach_dimension() mesh->nodes[i].lagVel.x = 0.;
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[i].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[i].pos.y);
        #if dimension > 2
        dist.z = GENERAL_1DIST(z, mesh->nodes[i].pos.z);
        #endif
        #if dimension < 3
        if (fabs(dist.x) <= 2*sdelta && fabs(dist.y) <= 2*sdelta) {
          double weight = (1 + cos(.5*pi*dist.x/sdelta))*
            (1 + cos(.5*pi*dist.y/sdelta))/16.;
        #else
        if (fabs(dist.x) <= 2*sdelta && fabs(dist.y) <= 2*sdelta
          && fabs(dist.z) <= 2*sdelta) {
          double weight = (1 + cos(.5*pi*dist.x/sdelta))*
            (1 + cos(.5*pi*dist.y/sdelta))*(1 + cos(.5*pi*dist.z/sdelta))/64.;
        #endif
          coord iu; // iu for "interpolated u"
          foreach_dimension() iu.x = 0.;
          read_vector_leaves(point, u, &iu, 0);
          foreach_dimension() mesh->nodes[i].lagVel.x += weight*iu.x;
        }
      }
    }
  }
  /**
  In case of parallel simulations, we communicate the Lagrangian velocity
  so that all processes have the same Lagrangian velocities.
  */
  #if _MPI
    if (mpi_npe > 1) reduce_lagVel(mesh);
  #endif
}

/**
The functions below fills a scalar field "stencils" with noise in all cached
cells. Passing this scalar to the \textit{adapt_wavelet} function ensure all
the 5x5(x5) stencils around the Lagrangian nodes are at the same level.
*/
// used for debugging
// #define STENCIL_TAG (point.level - 3)
#define STENCIL_TAG (noise())

void tag_stencil_leaves(Point point, coord dist, double sdelta) {
  if (is_local(cell)) {
    if (is_leaf(cell))
      foreach_dimension() stencils[] = STENCIL_TAG;
    else
    foreach_child() tag_stencil_leaves(point, dist, sdelta);
  }
}

trace
void tag_ibm_stencils_one_caps(lagMesh* mesh) {
  for(int i=0; i<mesh->nlp; i++) {
    foreach_dimension() mesh->nodes[i].lagVel.x = 0.;
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        double sdelta; // sdelta for "stencil delta"
        #if CONSTANT_STENCIL_LEVEL
          sdelta = L0/(1 << MB_LEVEL);
        #else
          sdelta = L0/(1 << mesh->nodes[i].stenstat.slvl);
        #endif
        coord dist;
        dist.x = GENERAL_1DIST(x, mesh->nodes[i].pos.x);
        dist.y = GENERAL_1DIST(y, mesh->nodes[i].pos.y);
        #if dimension > 2
        dist.z = GENERAL_1DIST(z, mesh->nodes[i].pos.z);
        #endif
        #if dimension < 3
          if (fabs(dist.x) <= 2*sdelta && fabs(dist.y) <= 2*sdelta) {
        #else
          if (fabs(dist.x) <= 2*sdelta && fabs(dist.y) <= 2*sdelta
            && fabs(dist.z) <= 2*sdelta) {
        #endif
          tag_stencil_leaves(point, dist, sdelta);
        }
      }
    }
  }
}

trace
void tag_ibm_stencils() {
  foreach() stencils[] = 0.;
  for(int k=0; k<NCAPS; k++) tag_ibm_stencils_one_caps(&MB(k));
}

/**
## References

~~~bib
@Article{Peskin1977,
  author    = {Peskin, C.S.},
  title     = {{Numerical analysis of blood flow in the heart}},
  journal   = {Journal of Computational Physics},
  year      = {1977},
  volume    = {25},
  number    = {3},
  pages     = {220--252},
  file      = {:files/Peskin1977 - Numerical Analysis of Blood Flow in the Heart.pdf:PDF},
  groups    = {FluidSolid flows},
  timestamp = {2013.07.18},
}
~~~
*/
