#include "foreach_cell_bfs.h"
#include "vtkType.h"

/**
 * @brief This struct and the associated members collect descriptive information from basilisk to describe the tree as a
 * HyperTreeGrid.
 */
typedef struct {
    size_t max_vertices;    /**< Number of cells in the tree level with the greatest number of cells */
    int64_t number_of_depths; /**< Common number of depth levels in this Basilisk forest view. */
    int64_t *depth_per_tree;  /**< Per-root tree depth array expected by VTKHDF. */

    Bit_t *descriptors;      /**< Breadth-first search of our tree. See @ref vtk_hdf_get_descriptors */
    size_t descriptors_size; /**< Number of bytes in @ref descriptors */
    int64_t n_descriptors;   /**< Number of bits in @ref descriptors */
    int64_t *descriptors_size_per_tree; /**< Internal descriptor-bit count for each root tree. */

    int64_t number_of_cells;                 /**< Total number of cells in the local forest view. */
    int64_t *number_of_cells_per_tree;       /**< Number of cells in each root tree. */
    int64_t *number_of_cells_per_tree_depth; /**< Flattened per-root, per-depth cell counts. */
    int64_t number_of_trees;                 /**< Number of root trees in the local forest view. */
    int64_t *tree_ids;                       /**< Ids of each tree in VTK root-grid order. */

    bool has_mask;    /**< If we are using a mask or not */
    Bit_t *mask;      /**< Binary mask to hide cells in the tree */
    size_t mask_size; /**< Number of bytes in @ref mask */
    size_t n_masks;   /**< Number of bits in @ref mask */

    double *x;  /**< Bounds of the tree on the x-axis. */
    double *y;  /**< Bounds of the tree on the y-axis. */
    double *z;  /**< Bounds of the tree on the z-axis. */
    size_t n_x; /**< Number of bounds of the tree on the x-axis. See @ref x */
    size_t n_y; /**< Number of bounds of the tree on the y-axis. See @ref y */
    size_t n_z; /**< Number of bounds of the tree on the z-axis. See @ref z */

} vtkHDFHyperTreeGridData;

static size_t vtk_hdf_number_of_tree_depth_values(const vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {
    return (size_t)vtk_hdf_hypertreegrid_data->number_of_trees *
           (size_t)vtk_hdf_hypertreegrid_data->number_of_depths;
}

static int64_t vtk_hdf_tree_id_from_root(Point point) {
#if dimension == 1
    return (int64_t)(point.i - GHOSTS);
#elif dimension == 2
    int64_t ix = (int64_t)(point.i - GHOSTS);
    int64_t iy = (int64_t)(point.j - GHOSTS);
    return ix + (int64_t)Dimensions.x * iy;
#else
    int64_t ix = (int64_t)(point.i - GHOSTS);
    int64_t iy = (int64_t)(point.j - GHOSTS);
    int64_t iz = (int64_t)(point.k - GHOSTS);
    return ix + (int64_t)Dimensions.x * (iy + (int64_t)Dimensions.y * iz);
#endif
}

static double vtk_hdf_root_delta(void) {
    return L0 / (double)Dimensions.x;
}

static int64_t vtk_hdf_number_of_roots(void) {
#if dimension == 1
    return (int64_t)Dimensions.x;
#elif dimension == 2
    return (int64_t)Dimensions.x * (int64_t)Dimensions.y;
#else
    return (int64_t)Dimensions.x * (int64_t)Dimensions.y * (int64_t)Dimensions.z;
#endif
}

static size_t vtk_hdf_bitstream_nbytes(size_t nbits) {
    return (nbits + 7) / 8;
}

/**
 * @brief Get the X/Y/Z corners of the tree.
 * 
 * @memberof vtkHDFHyperTreeGridData
 */
void vtk_hdf_get_coordinates(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {

    double root_delta = vtk_hdf_root_delta();

#if dimension >= 1
    /*
     * VTKHDF HyperTreeGrid stores a regular grid of root trees using
     * coordinate arrays plus the Dimensions attribute. We assume Dimensions is
     * the number of coordinates, not the number of root cells; consequently a
     * forest with Dimensions.x root cells needs Dimensions.x + 1 x-coordinates.
     *
     * Basilisk forest roots are square/cubic. The root spacing in every active
     * direction is therefore L0/Dimensions.x, not L0/Dimensions.y or
     * L0/Dimensions.z.
     */
    vtk_hdf_hypertreegrid_data->n_x = (size_t)Dimensions.x + 1;
    vtk_hdf_hypertreegrid_data->x = malloc(vtk_hdf_hypertreegrid_data->n_x * sizeof(double));
    for (size_t i = 0; i < vtk_hdf_hypertreegrid_data->n_x; i++)
        vtk_hdf_hypertreegrid_data->x[i] = X0 + (double)i * root_delta;
#else
    vtk_hdf_hypertreegrid_data->n_x = 1;
    vtk_hdf_hypertreegrid_data->x = malloc(1 * sizeof(double));
    vtk_hdf_hypertreegrid_data->x[0] = 0;
#endif

#if dimension >= 2
    vtk_hdf_hypertreegrid_data->n_y = (size_t)Dimensions.y + 1;
    vtk_hdf_hypertreegrid_data->y = malloc(vtk_hdf_hypertreegrid_data->n_y * sizeof(double));
    for (size_t j = 0; j < vtk_hdf_hypertreegrid_data->n_y; j++)
        vtk_hdf_hypertreegrid_data->y[j] = Y0 + (double)j * root_delta;
#else
    vtk_hdf_hypertreegrid_data->n_y = 1;
    vtk_hdf_hypertreegrid_data->y = malloc(1 * sizeof(double));
    vtk_hdf_hypertreegrid_data->y[0] = 0;
#endif

#if dimension >= 3
    vtk_hdf_hypertreegrid_data->n_z = (size_t)Dimensions.z + 1;
    vtk_hdf_hypertreegrid_data->z = malloc(vtk_hdf_hypertreegrid_data->n_z * sizeof(double));
    for (size_t k = 0; k < vtk_hdf_hypertreegrid_data->n_z; k++)
        vtk_hdf_hypertreegrid_data->z[k] = Z0 + (double)k * root_delta;
#else
    vtk_hdf_hypertreegrid_data->n_z = 1;
    vtk_hdf_hypertreegrid_data->z = malloc(1 * sizeof(double));
    vtk_hdf_hypertreegrid_data->z[0] = 0;
#endif
}

/**
 * @brief Count the number of cells on each level of the (binary/quad/oct)tree. Also compute the maximum
 *
 * @memberof vtkHDFHyperTreeGridData
 */
void vtk_hdf_get_number_of_vertices_per_depth(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {

    int depth = grid->maxdepth;
    size_t number_of_trees = (size_t)vtk_hdf_hypertreegrid_data->number_of_trees;

    vtk_hdf_hypertreegrid_data->number_of_depths = (int64_t)(depth + 1);
    size_t depth_values = vtk_hdf_number_of_tree_depth_values(vtk_hdf_hypertreegrid_data);

    free(vtk_hdf_hypertreegrid_data->depth_per_tree);
    free(vtk_hdf_hypertreegrid_data->number_of_cells_per_tree);
    free(vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth);
    free(vtk_hdf_hypertreegrid_data->descriptors_size_per_tree);

    vtk_hdf_hypertreegrid_data->depth_per_tree = malloc(number_of_trees * sizeof(int64_t));
    vtk_hdf_hypertreegrid_data->number_of_cells_per_tree = calloc(number_of_trees, sizeof(int64_t));
    vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth = calloc(depth_values, sizeof(int64_t));
    vtk_hdf_hypertreegrid_data->descriptors_size_per_tree = calloc(number_of_trees, sizeof(int64_t));

    for (size_t root_id = 0; root_id < number_of_trees; root_id++)
        vtk_hdf_hypertreegrid_data->depth_per_tree[root_id] = vtk_hdf_hypertreegrid_data->number_of_depths;

    /*
     * Assumption: VTKHDF splits forest metadata by TreeIds order. We therefore
     * count every root tree independently and flatten NumberOfCellsPerTreeDepth
     * as tree-major rows:
     *
     *   tree0[level0..levelN], tree1[level0..levelN], ...
     *
     * This also assumes foreach_cell_BFS() behaves like a forest iterator: it
     * visits a level-0 root before that root's descendants, and each root's BFS
     * stream is contiguous.
     */
    int64_t current_tree = 0;
    foreach_cell_BFS() {
        if (level == 0)
            current_tree = vtk_hdf_tree_id_from_root(point);
        if (current_tree < 0 || current_tree >= vtk_hdf_hypertreegrid_data->number_of_trees)
            continue;

        size_t offset = (size_t)current_tree * (size_t)vtk_hdf_hypertreegrid_data->number_of_depths + (size_t)level;
        vtk_hdf_hypertreegrid_data->number_of_cells_per_tree[current_tree]++;
        vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth[offset]++;
    }

    // Compute the largest cell count on any level and total number of cells
    vtk_hdf_hypertreegrid_data->max_vertices = 0;
    vtk_hdf_hypertreegrid_data->number_of_cells = 0;

    for (size_t i = 0; i < depth_values; i++) {
        int64_t nv = vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth[i];
        if (nv > vtk_hdf_hypertreegrid_data->max_vertices)
            vtk_hdf_hypertreegrid_data->max_vertices = nv;
    }

    for (size_t root_id = 0; root_id < number_of_trees; root_id++) {
        size_t deepest_offset =
            root_id * (size_t)vtk_hdf_hypertreegrid_data->number_of_depths +
            (size_t)vtk_hdf_hypertreegrid_data->number_of_depths - 1;
        int64_t deepest = vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth[deepest_offset];

        vtk_hdf_hypertreegrid_data->number_of_cells += vtk_hdf_hypertreegrid_data->number_of_cells_per_tree[root_id];
        vtk_hdf_hypertreegrid_data->descriptors_size_per_tree[root_id] =
            vtk_hdf_hypertreegrid_data->number_of_cells_per_tree[root_id] - deepest;
    }
}

/**
 * @brief Create a mask for the grid to hide remote leaf cells.
 *
 * @memberof vtkHDFHyperTreeGridData
 */
void vtk_hdf_get_local_mask(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {

    size_t number_of_cells = vtk_hdf_hypertreegrid_data->number_of_cells;
    // size_t max_depth = vtk_hdf_hypertreegrid_data->depth_per_tree;

    size_t n_masks = number_of_cells;
    vtk_hdf_hypertreegrid_data->n_masks = n_masks;
    /*
     * VTKHDF stores Mask as a packed bitstream. The bit count is not stored in
     * a separate MaskSize dataset; readers derive it from NumberOfCells for
     * each piece. In MPI this byte count is therefore rounded up independently
     * on each rank/piece before the byte streams are concatenated.
     */
    size_t mask_size = vtk_hdf_bitstream_nbytes(n_masks);
    vtk_hdf_hypertreegrid_data->mask_size = mask_size;

    vtk_hdf_hypertreegrid_data->mask = malloc(mask_size * sizeof(Bit_t));
    memset(vtk_hdf_hypertreegrid_data->mask, 0, mask_size);

    uint8_t bit_count = 0;
    size_t byte_count = 0;

    foreach_cell_BFS() {
        if (!(is_local(cell)) && is_leaf(cell)) {
            bool val = true;
            foreach_neighbor(1) {
              if (is_local(cell)) {
                val = false;
              }
            }
            assert(byte_count < vtk_hdf_hypertreegrid_data->mask_size);
            vtk_hdf_hypertreegrid_data->mask[byte_count] |= (uint8_t)(val << (7 - bit_count));
        }

        bit_count++;

        if (bit_count == 8) {
            bit_count = 0;
            byte_count++;
        }
    }

    vtk_hdf_hypertreegrid_data->has_mask = true;
}

/**
 * @brief Traverse the (binary-/quad-/oct-)tree using breadth-first search and record if the node is refined or not.
 *
 * This function uses the custom @ref foreach_cell_BFS() iterator macro to perform breadth first search. If the nth
 * visited node is not a leaf (i.e. it is refined), then we set the corresponding bit in our byte array to 1 (true).
 * Otherwise we set it to 0 (false).
 *
 * @memberof vtkHDFHyperTreeGridData
 */
void vtk_hdf_get_descriptors(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {

    // Get the maximum depth of the tree
    size_t max_depth = (size_t)vtk_hdf_hypertreegrid_data->number_of_depths;

    /*
     * VTKHDF stores descriptor bits for all trees in the current piece as one
     * concatenated packed stream. DescriptorsSize is the total bit count for
     * the current piece, not a per-tree array. DepthPerTree/TreeIds split the
     * stream back into trees.
     *
     * As in the original single-tree writer, deepest-level cells do not need
     * descriptor bits because they cannot be refined further.
     */
    size_t n_descriptors = 0;
    for (size_t root_id = 0; root_id < (size_t)vtk_hdf_hypertreegrid_data->number_of_trees; root_id++)
        n_descriptors += (size_t)vtk_hdf_hypertreegrid_data->descriptors_size_per_tree[root_id];

    // Store the result in our object; we need this later
    vtk_hdf_hypertreegrid_data->n_descriptors = (int64_t)n_descriptors;

    // Calculate the number of bytes needed to pack n number of descriptors
    size_t descriptors_size = vtk_hdf_bitstream_nbytes(n_descriptors);

    // Store this result
    vtk_hdf_hypertreegrid_data->descriptors_size = descriptors_size;

    // Heap allocate the byte array
    vtk_hdf_hypertreegrid_data->descriptors = malloc(descriptors_size * sizeof(Bit_t));
    memset(vtk_hdf_hypertreegrid_data->descriptors, 0, descriptors_size);

    // Keep a count of the bit + byte where we are in the byte array
    uint8_t bit_count = 0;
    size_t byte_count = 0;

    int64_t current_tree = 0;
    foreach_cell_BFS() {
        if (level == 0)
            current_tree = vtk_hdf_tree_id_from_root(point);
        if (current_tree < 0 || current_tree >= vtk_hdf_hypertreegrid_data->number_of_trees)
            continue;
        if ((size_t)level == max_depth - 1)
            continue;

        // If the node is refined (not a leaf), set the correct bit of the byte to 1
        assert(byte_count < vtk_hdf_hypertreegrid_data->descriptors_size);
        if (!is_leaf(cell))
            vtk_hdf_hypertreegrid_data->descriptors[byte_count] |= (uint8_t)(1 << (7 - bit_count));

        // Increment bit
        bit_count++;

        // Increment byte
        if (bit_count == 8) {
            bit_count = 0;
            byte_count++;
        }
    }
}

/**
 * @brief Desctructor for the @ref vtkHDFHyperTreeGridData struct.
 *
 * @memberof vtkHDFHyperTreeGridData
 */
void vtk_hdf_hypertreegrid_data_free(vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data) {
    if (vtk_hdf_hypertreegrid_data) {
        free(vtk_hdf_hypertreegrid_data->x);
        free(vtk_hdf_hypertreegrid_data->y);
        free(vtk_hdf_hypertreegrid_data->z);
        free(vtk_hdf_hypertreegrid_data->depth_per_tree);
        free(vtk_hdf_hypertreegrid_data->number_of_cells_per_tree);
        free(vtk_hdf_hypertreegrid_data->number_of_cells_per_tree_depth);
        free(vtk_hdf_hypertreegrid_data->descriptors_size_per_tree);
        free(vtk_hdf_hypertreegrid_data->descriptors);
        free(vtk_hdf_hypertreegrid_data->mask);
        free(vtk_hdf_hypertreegrid_data->tree_ids);
    }
}

/**
 * @brief Constructor for the @ref vtkHDFHyperTreeGridData struct. 
 *
 * @memberof vtkHDFHyperTreeGridData
 */
vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data_init(void) {
    vtkHDFHyperTreeGridData *vtk_hdf_hypertreegrid_data = calloc(1, sizeof(vtkHDFHyperTreeGridData));

    if (!vtk_hdf_hypertreegrid_data) {
        // handle malloc error
    };

    if (!grid) {
        fprintf(stderr, "Grid has not yet been initialized: Please initialize the "
                        "grid first.\n");
        exit(1);
    }

    vtk_hdf_hypertreegrid_data->number_of_trees = vtk_hdf_number_of_roots();
    vtk_hdf_hypertreegrid_data->tree_ids = qcalloc(vtk_hdf_hypertreegrid_data->number_of_trees, int64_t);

    /*
     * Assumption: TreeIds use row-major root-grid order matching
     * vtk_hdf_tree_id_from_root() and the order in which the forest BFS emits
     * root trees.
     */
    for (int64_t i = 0; i < vtk_hdf_hypertreegrid_data->number_of_trees; i++) {
        vtk_hdf_hypertreegrid_data->tree_ids[i] = i;
    }

    vtk_hdf_get_coordinates(vtk_hdf_hypertreegrid_data);
    vtk_hdf_get_number_of_vertices_per_depth(vtk_hdf_hypertreegrid_data);
    vtk_hdf_get_descriptors(vtk_hdf_hypertreegrid_data);

#if MPI_SINGLE_FILE
    vtk_hdf_get_local_mask(vtk_hdf_hypertreegrid_data);
#endif

    return vtk_hdf_hypertreegrid_data;
}
