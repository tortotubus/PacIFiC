/**
 * @brief This macro performs a [breadth-first
 * search](https://en.wikipedia.org/wiki/Breadth-first_search) of the basilisk
 * tree, and is largely rewritten and based on the builtin basilisk macro
 * [foreach_cell_root](https://basilisk.fr/src/grid/foreach_cell.h), which
 * performs a [depth-first
 * search](https://en.wikipedia.org/wiki/Depth-first_search) of the tree.
 *
 * @param root The Point where the root of the tree is
 *
 * Aside from the basilisk-specific semantics, it is worthwhile to understand
 * what breadth-first search is and how queues inherently relate to the data
 * structure
 */
macro2 foreach_cell_root_BFS(Point root) {
  {
#if dimension == 1
    Point point = {0};
    int ig = 0;
    NOT_UNUSED(ig);
    typedef struct {
      int l, i;
    } queue_t;
#elif dimension == 2
    Point point = {0};
    int ig = 0, jg = 0;
    NOT_UNUSED(ig);
    NOT_UNUSED(jg);
    typedef struct {
      int l, i, j;
    } queue_t;
#else // dimension == 3
    Point point = {0};
    int ig = 0, jg = 0;
    int kg = 0;
    NOT_UNUSED(ig);
    NOT_UNUSED(jg);
    NOT_UNUSED(kg);
    typedef struct {
      int l, i, j, k;
    } queue_t;
#endif

    size_t QMAX = 0;
    foreach_cell() { QMAX++; }

    // queue_t queue[QMAX];
    queue_t *queue = calloc(QMAX, sizeof(queue_t));
    size_t front = 0, rear = 0;

#if dimension == 1
    queue[rear++] = (queue_t){root.level, root.i};
#elif dimension == 2
    queue[rear++] = (queue_t){root.level, root.i, root.j};
#else // dimension == 3
    queue[rear++] = (queue_t){root.level, root.i, root.j, root.k};
#endif

    while (front < rear) {
      queue_t current = queue[front++];
      point.level = current.l;
#if dimension == 1
      point.i = current.i;
#elif dimension == 2
      point.i = current.i;
      point.j = current.j;
#else // dimension == 3
      point.i = current.i;
      point.j = current.j;
      point.k = current.k;
#endif

      /* Check if the cell is allocated (macro 'allocated' expects 'point') */
      if (!allocated(0, 0, 0))
        continue;

      // clang-format off
      {...} // User defined code block
      // clang-format on

#if dimension == 1
      if (point.level < grid->depth && !is_leaf(cell)) {
        if (rear + 2 > QMAX)
          abort();
        queue[rear++] = (queue_t){point.level + 1, _LEFT};
        queue[rear++] = (queue_t){point.level + 1, _RIGHT};
      }
#elif dimension == 2
      if (point.level < grid->depth && !is_leaf(cell)) {
        if (rear + 4 > QMAX)
          abort();
        queue[rear++] = (queue_t){point.level + 1, _LEFT, _BOTTOM};
        queue[rear++] = (queue_t){point.level + 1, _RIGHT, _BOTTOM};
        queue[rear++] = (queue_t){point.level + 1, _LEFT, _TOP};
        queue[rear++] = (queue_t){point.level + 1, _RIGHT, _TOP};
      }
#else // dimension == 3
      if (point.level < grid->depth && !is_leaf(cell)) {
        if (rear + 8 > QMAX)
          abort();
        queue[rear++] = (queue_t){point.level + 1, _LEFT, _BOTTOM, _BACK};
        queue[rear++] = (queue_t){point.level + 1, _RIGHT, _BOTTOM, _BACK};
        queue[rear++] = (queue_t){point.level + 1, _LEFT, _TOP, _BACK};
        queue[rear++] = (queue_t){point.level + 1, _RIGHT, _TOP, _BACK};
        queue[rear++] = (queue_t){point.level + 1, _LEFT, _BOTTOM, _FRONT};
        queue[rear++] = (queue_t){point.level + 1, _RIGHT, _BOTTOM, _FRONT};
        queue[rear++] = (queue_t){point.level + 1, _LEFT, _TOP, _FRONT};
        queue[rear++] = (queue_t){point.level + 1, _RIGHT, _TOP, _FRONT};
      }
#endif
    } // End of loop over the current level

    free(queue);
  }
}

macro2 foreach_cell_BFS(ivec d = Dimensions) {
  {
    /*
     * Root trees must be emitted in vtk_hdf_tree_id_from_root() order:
     * x fastest, then y, then z. Descriptors, masks, and cell data are all
     * packed in this stream and decoded through TreeIds.
     */
#if dimension == 1
    for (int ox = 0; ox < d.x; ox++) {
      Point root = {GHOSTS + ox, 0};
      foreach_cell_root_BFS(root)
      // clang-format off
      {...}
      // clang-format on
    }
#elif dimension == 2
    for (int oy = 0; oy < d.y; oy++) {
      for (int ox = 0; ox < d.x; ox++) {
        Point root = {GHOSTS + ox, GHOSTS + oy, 0};
        foreach_cell_root_BFS(root)
        // clang-format off
        {...}
        // clang-format on
      }
    }
#else // dimension == 3
    for (int oz = 0; oz < d.z; oz++) {
      for (int oy = 0; oy < d.y; oy++) {
        for (int ox = 0; ox < d.x; ox++) {
          Point root = {GHOSTS + ox, GHOSTS + oy, GHOSTS + oz};
          foreach_cell_root_BFS(root)
          // clang-format off
          {...}
          // clang-format on
        }
      }
    }
#endif
  }
}
