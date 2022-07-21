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

#ifndef CONSTANT_MB_LEVEL
  #define CONSTANT_MB_LEVEL 1
#endif

#if CONSTANT_MB_LEVEL
  #define MB_DELTA (L0/(1 << grid->maxdepth))
#endif

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

int get_level_IBM_stencil(lagNode* node) {
  #if CONSTANT_MB_LVL
    return grid->maxdepth;
  #else
    #if dimension < 3
      Point point = locate(node->pos.x, node->pos.y);
    #else
      Point point = locate(node->pos.x, node->pos.y, node->pos.z);
    #endif
    int lvl = point.level;
    bool complete_stencil = false;
    while (!complete_stencil) {
      complete_stencil = true;
      bool changed_lvl = false;
      double delta = (L0/(1 << lvl));
      for(int ni=-2; ni<=2; ni++) {
        for(int nj=-2; nj<=2 && !changed_lvl ; nj++) {
          #if dimension < 3
          point = locate(POS_PBC_X(node->pos.x + ni*delta),
            POS_PBC_Y(node->pos.y + nj*delta));
          #else
          for(int nk=-2; nk<=2; nk++) {
            point = locate(POS_PBC_X(node->pos.x + ni*delta),
              POS_PBC_Y(node->pos.y + nj*delta),
              POS_PBC_Z(node->pos.z + nk*delta));
          #endif
          if (point.level < lvl) {
            lvl = point.level;
            changed_lvl = true;
            complete_stencil = false;
          }
        #if dimension > 2
          }
        #endif
        }
      }
    }
    return lvl;
  #endif
}


/**
The function below loops through the Lagrangian nodes and "caches" the Eulerian
cells in a 5x5(x5) stencil around each node. In case of parallel simulations,
the cached cells are tagged with the process id.
*/
scalar stencils[];
trace
void generate_lag_stencils_one_caps(lagMesh* mesh) {
  for(int i=0; i<mesh->nlp; i++) {
    mesh->nodes[i].stencil.n = 0;
    /**
    The current implementation assumes that the Eulerian cells around Lagrangian
    node are all at the maximum level.
    */
    int lvl = get_level_IBM_stencil(&mesh->nodes[i]);
    double delta = L0/(1 << lvl);
    mesh->nodes[i].slvl = lvl;
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
        #if CONSTANT_MB_LEVEL
          if (point.level >= 0 && point.level != grid->maxdepth)
            fprintf(stderr,
              "Warning: Lagrangian stencil not fully resolved.\n");
        #endif
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
  for(int k=0; k<NCAPS; k++) generate_lag_stencils_one_caps(&MB(k));
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
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        double sdelta; // sdelta for "stencil delta"
        #if CONSTANT_MB_LEVEL
          sdelta = Delta;
        #else
          sdelta = L0/(1 << mesh->nodes[i].slvl);
        #endif
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
            /(sq(4*sdelta));
        #else
        if (fabs(dist.x) <= 2*sdelta && fabs(dist.y) <= 2*sdelta &&
          fabs(dist.z) <= 2*sdelta) {
          double weight =
            (1 + cos(.5*pi*dist.x/sdelta))*(1 + cos(.5*pi*dist.y/sdelta))
            *(1 + cos(.5*pi*dist.z/sdelta))/(cube(4*sdelta));
        #endif
        #if CONSTANT_MB_LEVEL
          foreach_dimension() forcing.x[] += weight*mesh->nodes[i].lagForce.x;
        #else
          if (is_leaf(cell)) {
            foreach_dimension() forcing.x[] += weight*mesh->nodes[i].lagForce.x;
          }
          else {
            foreach_child() {
              if (is_local(cell)) {
                if (is_leaf(cell)) {
                  foreach_dimension() forcing.x[] +=
                    weight*mesh->nodes[i].lagForce.x/(1 << dimension);
                }
                else {
                  foreach_child() {
                    if (is_local(cell)) {
                      if (is_leaf(cell)) {
                        foreach_dimension() forcing.x[] +=
                          weight*mesh->nodes[i].lagForce.x/sq(1 << dimension);
                      }
                      else fprintf(stderr, "Error: too many different levels in"
                        "one IBM stencil\n");
                    }
                  }
                }
              }
            }
          }
        #endif
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
    foreach_dimension() mesh->nodes[i].lagVel.x = 0.;
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        double sdelta; // sdelta for "stencil delta"
        #if CONSTANT_MB_LEVEL
          sdelta = Delta;
        #else
          sdelta = L0/(1 << mesh->nodes[i].slvl);
        #endif
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
        #if CONSTANT_MB_LEVEL
          foreach_dimension() iu.x = u.x[];
        #else
          foreach_dimension() iu.x = 0.;
          if (is_leaf(cell)) {
            foreach_dimension() iu.x = u.x[];
          }
          else {
            foreach_child() {
              if (is_local(cell)) {
                if (is_leaf(cell)) {
                  foreach_dimension() iu.x += u.x[]/(1 << dimension);
                }
                else {
                  foreach_child() {
                    if (is_local(cell)) {
                      if (is_leaf(cell)) {
                        foreach_dimension() iu.x += u.x[]/sq(1 << dimension);
                      }
                      else fprintf(stderr, "Error: too many different levels in"
                        "one IBM stencil\n");
                    }
                  }
                }
              }
            }
          }
        #endif
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
The function below fills a scalar field "stencils" with noise in all "cached"
cells. Passing this scalar to the \textit{adapt_wavelet} function ensure all
the 5x5(x5) stencils around the Lagrangian nodes are at the same level.
*/
trace
void tag_ibm_stencils_one_caps(lagMesh* mesh) {
  for(int i=0; i<mesh->nlp; i++) {
    foreach_dimension() mesh->nodes[i].lagVel.x = 0.;
    foreach_cache(mesh->nodes[i].stencil) {
      if (point.level >= 0) {
        double sdelta; // sdelta for "stencil delta"
        #if CONSTANT_MB_LEVEL
          sdelta = Delta;
        #else
          sdelta = L0/(1 << mesh->nodes[i].slvl);
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
        #if CONSTANT_MB_LEVEL
          #if dimension < 3
          stencils[] = sq(dist.x + dist.y )/
            sq(2.*sdelta)*(2.+noise());
          #else
            stencils[] = sq(dist.x + dist.y + dist.z)/
              cube(2.*sdelta)*(2.+noise());
          #endif
        #else
          if (is_leaf(cell)) {
            #if dimension < 3
              stencils[] = sq(dist.x + dist.y )/
                sq(2.*sdelta)*(2.+noise());
            #else
              stencils[] = sq(dist.x + dist.y + dist.z)/
                cube(2.*sdelta)*(2.+noise());
            #endif
          }
          else {
            foreach_child() {
              if (is_local(cell)) {
                if (is_leaf(cell)) {
                  #if dimension < 3
                    stencils[] = sq(dist.x + dist.y )/
                      sq(2.*sdelta)*(2.+noise());
                  #else
                    stencils[] = sq(dist.x + dist.y + dist.z)/
                      cube(2.*sdelta)*(2.+noise());
                  #endif
                }
                else {
                  foreach_child() {
                    if (is_local(cell)) {
                      if (is_leaf(cell)) {
                        #if dimension < 3
                          stencils[] = sq(dist.x + dist.y )/
                            sq(2.*sdelta)*(2.+noise());
                        #else
                          stencils[] = sq(dist.x + dist.y + dist.z)/
                            cube(2.*sdelta)*(2.+noise());
                        #endif
                      }
                      else fprintf(stderr, "Error: too many different levels in"
                        "one IBM stencil\n");
                    }
                  }
                }
              }
            }
          }
        #endif
        }
      }
    }
  }

  #if OLD_QCC
  boundary({stencils});
  #endif
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
