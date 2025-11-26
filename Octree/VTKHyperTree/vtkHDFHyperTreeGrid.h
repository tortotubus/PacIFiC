#include "vtkHDF.h"
#include "vtkHDFHyperTreeGridData.h"

#include <float.h>
#include <math.h>

#define COMPRESSION 1
#define COMPRESSION_LEVEL 7

// #define CHUNK_SIZE (1 << (4))

/**
 * @brief This struct holds various IDs needed by the HDF5 library to read and write HDF5 files for our particular
 * HyperTreeGrid/PHyperTreeGrid schema.
 */
typedef struct {
    /* Parent */
    vtkHDF vtk_hdf;

    hid_t grp_celldata_id;
    // hid_t grp_steps_id;

    /* Dataspace, datatype, and property-list identifiers */
    hid_t attr_space_id;
    hid_t attr_dtype_id;
    hid_t dset_space_id;
    hid_t dcpl_id;
    hid_t dset_dtype_id;

    /* Dataset and attribute identifiers */
    hid_t dset_id;
    hid_t attr_id;

    /* MPIIO */
    hid_t file_space;
    hid_t mem_space;
    hid_t xfer_plist;

} vtkHDFHyperTreeGrid;

/**
 * @brief This function
 */
void vtk_HDF_hypertreegrid_close(vtkHDFHyperTreeGrid *vtk_hdf_htg) {
    //
    if (vtk_hdf_htg->grp_celldata_id >= 0)
        H5Gclose(vtk_hdf_htg->grp_celldata_id);

    // if (vtk_hdf_htg->grp_steps_id >= 0)
    //   H5Gclose(vtk_hdf_htg->grp_steps_id);

    if (vtk_hdf_htg->xfer_plist >= 0)
        H5Pclose(vtk_hdf_htg->xfer_plist);

    if (vtk_hdf_htg->mem_space >= 0)
        H5Sclose(vtk_hdf_htg->mem_space);

    if (vtk_hdf_htg->file_space >= 0)
        H5Sclose(vtk_hdf_htg->file_space);

    vtk_HDF_close(&vtk_hdf_htg->vtk_hdf);
}

/**
 * @brief
 */
void vtk_HDF_hypertreegrid_error(vtkHDFHyperTreeGrid *vtk_hdf_htg) {
    vtk_HDF_hypertreegrid_close(vtk_hdf_htg);
    assert(1 == 2);
}

/**
 * @brief This function does the actual writing of the HyperTreeGrid/PHyperTreeGrid in the HDF5 format.
 *
 * @param scalar_list
 * @param vector_list
 * @param fname
 *
 * The HyperTreeGrid in [VTK](https://vtk.org/) is a
 * [class](https://vtk.org/doc/nightly/html/classvtkHyperTreeGrid.html#details) that is intended exactly for
 * representation of tree-based grids in VTK. It is similar to Basilisk, although the major difference is that while the
 * foreach_cell() iterator in Basilisk uses [depth-first search](https://en.wikipedia.org/wiki/Depth-first_search), the
 * HyperTreeGrid format requires us to describe the tree using [breadth-first
 * search](https://en.wikipedia.org/wiki/Breadth-first_search).
 *
 * Some other differences worth noting is that the HyperTreeGrid also supports multiple trees in a regular grid, even
 * though basilisk does not support this. Thus, there are some seemingly extraneous fields included that we do not make
 * use of, though are neccesary if we wanted to describe such a "forest" of trees.
 *
 * It is also worth noting that the HDF5 format itself is somewhat generic (like XML), even if VTK imposes a specific
 * schema for HyperTreeGrid formats. In HDF5, we have two basic building blocks: Dataset and Group objects. These are
 * self-descriptive. In addition, each may have "attributes" attached to them, which themselves are small pieces of data
 * used to help interpret and understand the Group and Dataset objects. In addition, it is worthwhile to understand HDF5
 * chunking, since this has some implications for writing appendable datasets (H5S_UNLIMITED) or when using compression
 * features of HDF5.
 *
 * In regards to what collection of Groups, Datasets, and Attributes in HDF5 define a HyperTreeGrid, the only official
 * documentation on the schema from Kitware is found
 * [here](https://web.archive.org/web/20250804033704/https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#hypertreegrid).
 * Note that even though the format supports the creation of temporal datasets, we write one separate file per time
 * step. This is to avoid the risk of corruption. In addition, this writer was made for version 2.4 (of the file
 * format). The format
 *
 * In summary, our schema in the case of a non-temporal (single-timestep) HyperTreeGrid is as follows
 *
 * - **Group**: "/VTKHDF"
 *
 *  - **Attribute**: "Version"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{2}`
 *    - **Description**: This is the VTKHDF file 'version' we use: The first entry is the major version and the second
 * entry represents the minor version.
 *
 *  - **Attribute**: "Type"
 *    - **Datatype**: `char`
 *    - **Dimension**: `{13}`
 *    - **Description**: This tells VTK what format we are using (HyperTreeGrid) as opposed to ImageGrid,
 * UnstructuredGrid, etc.
 *
 *  - **Attribute**: "BranchFactor"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{1}`
 *    - **Description**: When a tree cell is refined, it should have n^d children, where d is the dimension and n is the
 * branch factor. For basilisk, this is always 2.
 *
 *  - **Attribute**: "Dimensions"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{3}`
 *    - **Description**: The number of coordinates describing the trees in each dimension; this should match/describe
 * the number of entries in the XCoordinate, YCoordinate, and ZCoordinate datasets, when the tree is not a temporal tree
 * or PHyperTree.
 *
 *  - **Attribute**: "TransposedRootIndexing"
 *
 *  - **Dataset**: "XCoordinates"
 *    - **Datatype**: `double`
 *    - **Dimension**: `{2}`
 *    - **Description**: The x coordinate of two corners of our tree.
 *
 *  - **Dataset**: "YCoordinates"
 *    - **Datatype**: `double`
 *    - **Dimension**: `{2} or {1}`
 *    - **Description**: The y coordinate of two corners of our tree. If the tree is really a binary tree, just write 0
 * as the single entry.
 *
 *  - **Dataset**: ZCoordinates
 *    - **Datatype**: `double`
 *    - **Dimension**: `{2} or {1}`
 *    - **Description**: The z coordinate of two corners of our tree. If the tree is really a quad tree, just write 0 as
 * the single entry. *
 *
 *  - **Dataset**: "Descriptors"
 *    - **Datatype**: `char`
 *    - **Description**: This is the description of the actual structure of the tree by breadth-first search. \sa @ref
 * hdf_get_descriptors
 *
 *  - **Dataset**: "DescriptorsSize"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{1}`
 *    - **Description**: The size of the descriptors dataset. This only becomes useful when there are a forest of trees
 * (or parallel copies of the same tree), so that the VTK reader knows where the descriptors of each tree starts and
 * ends in the concatenation of all descriptors.
 *
 *  - **Dataset**: "TreeIds"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{1}`
 *    - **Value**: `0`
 *    - **Description**: TODO
 *
 *  - **Dataset**: "DepthPerTree"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{1}`
 *    - **Description**: The maximum depth of the tree. In the case of basilisk we only have one.
 *
 *  - **Dataset**: "NumberOfTrees"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{1}`
 *    - **Value**: `1`
 *    - **Description**:
 *
 *  - **Dataset**: "NumberOfDepths"
 *    - **Datatype**: `int64`
 *    - **Dimension**: `{1}`
 *    - **Value**: `1`
 *    - **Description**:
 *
 *  - **Dataset**: "Mask"
 *    - **Datatype**: `char`
 *    - **Dimension**: `{}`
 *    - **Description**: This field is entirely optional. If it is included, VTK will hide or show certain cells. The
 * format is nearly identical to that of the descriptors dataset.
 *
 *  - **Group**: "/VTKHDF/CellData"
 *
 *    - **Dataset**: "MyCellCenteredScalarField"
 *      - **DataType**: `float` or `double`
 *      - **Dimension**: `{}`
 *    - **Dataset**: "MyCellCenteredVectorField"
 *      - **DataType**: `float` or `double`
 *      - **Dimension**: `{}`
 *
 * If done correctly, we should be able to run `h5dump -H` on our file and get something like
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *  HDF5 "myfile.vtkhdf" {
 *  GROUP "/" {
 *     GROUP "VTKHDF" {
 *        ATTRIBUTE "BranchFactor" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
 *        }
 *        ATTRIBUTE "Dimensions" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 3 ) / ( 3 ) }
 *        }
 *        ATTRIBUTE "TransposedRootIndexing" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( 1 ) }
 *        }
 *        ATTRIBUTE "Type" {
 *           DATATYPE  H5T_STRING {
 *              STRSIZE 13;
 *              STRPAD H5T_STR_NULLPAD;
 *              CSET H5T_CSET_ASCII;
 *              CTYPE H5T_C_S1;
 *           }
 *           DATASPACE  SCALAR
 *        }
 *        ATTRIBUTE "Version" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 2 ) / ( 2 ) }
 *        }
 *        GROUP "CellData" {
 *           DATASET "p" {
 *              DATATYPE  H5T_IEEE_F32LE
 *              DATASPACE  SIMPLE { ( 87381 ) / ( H5S_UNLIMITED ) }
 *           }
 *           DATASET "u" {
 *              DATATYPE  H5T_IEEE_F32LE
 *              DATASPACE  SIMPLE { ( 87381, 2 ) / ( H5S_UNLIMITED, 2 ) }
 *           }
 *        }
 *        DATASET "DepthPerTree" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( H5S_UNLIMITED ) }
 *        }
 *        DATASET "Descriptors" {
 *           DATATYPE  H5T_STD_U8LE
 *           DATASPACE  SIMPLE { ( 2731 ) / ( H5S_UNLIMITED ) }
 *        }
 *        DATASET "DescriptorsSize" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( H5S_UNLIMITED ) }
 *        }
 *        DATASET "NumberOfCells" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( H5S_UNLIMITED ) }
 *        }
 *        DATASET "NumberOfCellsPerTreeDepth" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 9 ) / ( H5S_UNLIMITED ) }
 *        }
 *        DATASET "NumberOfDepths" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( H5S_UNLIMITED ) }
 *        }
 *        DATASET "NumberOfTrees" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( H5S_UNLIMITED ) }
 *        }
 *        DATASET "TreeIds" {
 *           DATATYPE  H5T_STD_I64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( H5S_UNLIMITED ) }
 *        }
 *        DATASET "XCoordinates" {
 *           DATATYPE  H5T_IEEE_F64LE
 *           DATASPACE  SIMPLE { ( 2 ) / ( H5S_UNLIMITED ) }
 *        }
 *        DATASET "YCoordinates" {
 *           DATATYPE  H5T_IEEE_F64LE
 *           DATASPACE  SIMPLE { ( 2 ) / ( H5S_UNLIMITED ) }
 *        }
 *        DATASET "ZCoordinates" {
 *           DATATYPE  H5T_IEEE_F64LE
 *           DATASPACE  SIMPLE { ( 1 ) / ( H5S_UNLIMITED ) }
 *        }
 *     }
 *  }
 *  }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * For the non-temporal parallel "PHyperTreeGrid", the schema is exactly identical, although we now write multiple trees
 * over the same X/Y/Z coordinates and (for the most part) VTK will stitch the trees from each process together in a
 * coherent way. This requires some MPI coordination on our part to know the offsets into each dataset--though
 * otherwise, the HDF5 library handles all other details of MPIIO coordination.
 *
 *
 *
 */
vtkHDFHyperTreeGrid vtk_HDF_hypertreegrid_init(scalar *scalar_list, vector *vector_list, const char *fname) {

    // Create the vtkHDF struct
    vtkHDF vtk_hdf = vtk_HDF_init(fname);

    // Initialize our object ids with -1 or H5I_INVALID_HID and the vtkHDF struct
    vtkHDFHyperTreeGrid vtk_hdf_htg = {
        .vtk_hdf = vtk_hdf,
        .grp_celldata_id = H5I_INVALID_HID,
        .attr_space_id = H5I_INVALID_HID,
        .attr_dtype_id = H5I_INVALID_HID,
        .attr_id = H5I_INVALID_HID,
        .dset_space_id = H5I_INVALID_HID,
        .dcpl_id = H5I_INVALID_HID,
        .dset_dtype_id = H5I_INVALID_HID,
        .dset_id = H5I_INVALID_HID,
        .file_space = H5I_INVALID_HID,
        .mem_space = H5I_INVALID_HID,
        .xfer_plist = H5I_INVALID_HID,
    };

    /*
     * Attribute: /VTKHDF/BranchFactor
     * Datatype: H5T_NATIVE_INT64 / int64_t
     * Dimension: {1}
     *
     * This attribute describes the number of children of a refined cell. This is \f(n^d\f) where \f(n\f) is the branch
     * factor and \f(d\f) is the dimension of the tree. In basilisk \f(n=2\f) always.
     */
    {
        // Value for the BranchFactor attrbiute
        int64_t bf_value = 2;
        hsize_t dims_attr[1] = {1};

        // Create the attribute space
        vtk_hdf_htg.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
        if (vtk_hdf_htg.attr_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Set the datatype
        vtk_hdf_htg.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
        if (vtk_hdf_htg.attr_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Create the attribute
        vtk_hdf_htg.attr_id = H5Acreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "BranchFactor", vtk_hdf_htg.attr_dtype_id,
                                         vtk_hdf_htg.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
        if (vtk_hdf_htg.attr_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Write the attribute
        if (H5Awrite(vtk_hdf_htg.attr_id, vtk_hdf_htg.attr_dtype_id, &bf_value) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Close
        H5Aclose(vtk_hdf_htg.attr_id);
        H5Tclose(vtk_hdf_htg.attr_dtype_id);
        H5Sclose(vtk_hdf_htg.attr_space_id);
    }

    /*
     * Attribute: /VTKHDF/Dimensions
     * Datatype: H5T_NATIVE_INT64 / int64_t
     * Dimension: {3}
     *
     * This attribute describes the number of coordinates in each cartesian X,Y,Z direction. For basilisk, each will
     * never be greater than 2, since we can only have 1 tree.
     */
    {
#if dimension == 1
        int64_t dims_value[3] = {2, 1, 1};
#elif dimension == 2
        int64_t dims_value[3] = {2, 2, 1};
#else
        int64_t dims_value[3] = {2, 2, 2};
#endif

        hsize_t dims_attr[1] = {3};

        /* Create a 1D dataspace of length 3 */
        vtk_hdf_htg.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
        if (vtk_hdf_htg.attr_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        /* Use a 64-bit‐int little‐endian type */
        vtk_hdf_htg.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
        if (vtk_hdf_htg.attr_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.attr_id = H5Acreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "Dimensions", vtk_hdf_htg.attr_dtype_id,
                                         vtk_hdf_htg.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
        if (vtk_hdf_htg.attr_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        /* Now actually write dims_value[3] into the attribute */
        if (H5Awrite(vtk_hdf_htg.attr_id, vtk_hdf_htg.attr_dtype_id, dims_value) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Aclose(vtk_hdf_htg.attr_id);
        H5Tclose(vtk_hdf_htg.attr_dtype_id);
        H5Sclose(vtk_hdf_htg.attr_space_id);
    }

    /*
     * Attribute: /VTKHDF/TransposedRootIndexing
     * Datatype: H5T_NATIVE_INT64 / int64_t
     * Dimension: {1}
     */
    {
        int64_t tri_value = 0;
        hsize_t dims_attr[1] = {1};
        vtk_hdf_htg.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
        if (vtk_hdf_htg.attr_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
        if (vtk_hdf_htg.attr_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.attr_id =
            H5Acreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "TransposedRootIndexing", vtk_hdf_htg.attr_dtype_id,
                       vtk_hdf_htg.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
        if (vtk_hdf_htg.attr_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Awrite(vtk_hdf_htg.attr_id, vtk_hdf_htg.attr_dtype_id, &tri_value) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Aclose(vtk_hdf_htg.attr_id);
        H5Tclose(vtk_hdf_htg.attr_dtype_id);
        H5Sclose(vtk_hdf_htg.attr_space_id);
    }

    /*
     * Attribute: /VTKHDF/Type
     * Datatype: H5T_C_S1 / char
     * Dimension: {13}
     *
     * This should be "HyperTreeGrid" to distinguish from other supported types (e.g. "ImageData")
     */
    {
        // Value of the attribute
        const char *type_str = "HyperTreeGrid";

        // Create the attribute space
        vtk_hdf_htg.attr_space_id = H5Screate(H5S_SCALAR);
        if (vtk_hdf_htg.attr_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Create a fixed-length string datatype of length 13, null-padded, ASCII
        vtk_hdf_htg.attr_dtype_id = H5Tcopy(H5T_C_S1);
        if (vtk_hdf_htg.attr_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Set the size
        if (H5Tset_size(vtk_hdf_htg.attr_dtype_id, (size_t)13) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Pad the string
        if (H5Tset_strpad(vtk_hdf_htg.attr_dtype_id, H5T_STR_NULLPAD) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Set the string value
        if (H5Tset_cset(vtk_hdf_htg.attr_dtype_id, H5T_CSET_ASCII) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Create the attribute
        vtk_hdf_htg.attr_id = H5Acreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "Type", vtk_hdf_htg.attr_dtype_id,
                                         vtk_hdf_htg.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
        if (vtk_hdf_htg.attr_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        /* Write the string (automatically null‐padded up to length 13) */
        if (H5Awrite(vtk_hdf_htg.attr_id, vtk_hdf_htg.attr_dtype_id, type_str) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Aclose(vtk_hdf_htg.attr_id);
        H5Tclose(vtk_hdf_htg.attr_dtype_id);
        H5Sclose(vtk_hdf_htg.attr_space_id);
    }

    /*
     * Attribute: /VTKHDF/Version
     * Datatype: H5T_NATIVE_INT64 / int64_t
     * Dimension: {2}
     */
    {
        int64_t vers_value[2] = {2, 4};
        hsize_t dims_attr[1] = {2};
        vtk_hdf_htg.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
        if (vtk_hdf_htg.attr_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
        if (vtk_hdf_htg.attr_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.attr_id = H5Acreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "Version", vtk_hdf_htg.attr_dtype_id,
                                         vtk_hdf_htg.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
        if (vtk_hdf_htg.attr_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Awrite(vtk_hdf_htg.attr_id, vtk_hdf_htg.attr_dtype_id, vers_value) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Aclose(vtk_hdf_htg.attr_id);
        H5Tclose(vtk_hdf_htg.attr_dtype_id);
        H5Sclose(vtk_hdf_htg.attr_space_id);
    }

    /*
     * Create the vtkHDFHyperTreeGridData object: This object contains a lot of local-view data that we will write after
     * this. See @ref vtk_hdf_hypertreegrid_data_init in @ref vtkHDFHyperTreeGridData.h
     */
    vtkHDFHyperTreeGridData *vtk_hdf_htg_data = vtk_hdf_hypertreegrid_data_init();

    /*
     * Group: /VTKHDF/CellData
     *
     * Here we place all of our grid, tree, and field data. This group does not require any attributes.
     */
    {
        vtk_hdf_htg.grp_celldata_id =
            H5Gcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "CellData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (vtk_hdf_htg.grp_celldata_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
    }

#if MPI_SINGLE_FILE
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    /*
     * Dataset: /VTKHDF/DepthPerTree
     * Datatype: H5T_NATIVE_INT64 / int64_t
     * Dimension: {1*npe()}
     *
     * The maximum depth of our single tree. In the MPI_SINGLE_FILE case, we must write this once per proc according to
     * the maximum (actual) depth of the tree on that proc.
     */
    {
#if MPI_SINGLE_FILE
        hsize_t local_size = (hsize_t)1;
        hsize_t global_size = (hsize_t)npe();
        hsize_t offset = pid() * local_size;
#else
        hsize_t local_size = (hsize_t)1;
        hsize_t global_size = (hsize_t)local_size;
#endif

        int64_t *depth_per_tree = &vtk_hdf_htg_data->depth_per_tree;
        hsize_t dims_d[1] = {global_size};
        hsize_t maxdims_d[1] = {global_size};

        vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
        if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t chunk_dims[1] = {local_size};
        if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
        if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_id = H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "DepthPerTree", vtk_hdf_htg.dset_dtype_id,
                                         vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
        if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE

        vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
        if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t start[1] = {offset};
        hsize_t count[1] = {local_size};
        if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        hsize_t m_dims[1] = {local_size};
        vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
        if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                     vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
                     vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
                     vtk_hdf_htg.file_space,    /* file dataspace with hyperslab selected */
                     vtk_hdf_htg.xfer_plist,    /* collective MPI‐IO transfer property */
                     depth_per_tree             /* pointer to local data */
                     ) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        H5Pclose(vtk_hdf_htg.xfer_plist);
        vtk_hdf_htg.xfer_plist = -1;
        H5Sclose(vtk_hdf_htg.mem_space);
        vtk_hdf_htg.mem_space = -1;
        H5Sclose(vtk_hdf_htg.file_space);
        vtk_hdf_htg.file_space = -1;
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#else
        if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, depth_per_tree) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Dclose(vtk_hdf_htg.dset_id);
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
    }

    /*
     * Dataset: /VTKHDF/Descriptors
     * Datatype: H5T_STD_U8LE / UInt8_t
     * Dimension: {global_size}
     *
     * The [breadth-first search](https://en.wikipedia.org/wiki/Breadth-first_search) description of the tree. Each bit
     * represents a node of the tree, starting with the root. If a given node is refined, then a '1' is written. If a
     * given node is not refined, then a '0' is written instead. For example, the binary tree
     *
     * [Level 0]                      o
     *                               / \
     *                              /   \
     * [Level 1]                   o     o
     *                            / \
     *                           /   \
     * [Level 2]                o     o
     *
     * would be written as "1 10 00" or "11000000" formatted as a byte. The actual values of this are computed in @ref
     * hdf_get_descriptors() in @ref vtkHDFHyperTreeGridData.h and the actual breadth-first traversal (in the correct
     * order) is performed by the @ref foreach_cell_bfs() macro in @ref foreach_cell_bfs.h.
     *
     * In the MPI_SINGLE_FILE case, we write this descriptor array for the tree local to each proc.
     *
     * \sa @ref hdf_get_descriptors
     * \sa @ref foreach_cell_bdf
     */
    {

#if MPI_SINGLE_FILE
        // Calculate local and global sizes for collective MPI-IO operations
        hsize_t local_size = (hsize_t)vtk_hdf_htg_data->descriptors_size;
        hsize_t global_size = local_size;

        // MPI_Allreduce: Aggregate local sizes to compute total global dataset size across all ranks
        MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        // MPI_Exscan: Calculate exclusive prefix sum to determine each rank's offset in the global dataset
        // This ensures each MPI rank writes to a non-overlapping region of the HDF5 data array
        hsize_t offset = 0;
        MPI_Exscan(&local_size, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (pid() == 0)
            offset = 0;
#else
        hsize_t local_size = (hsize_t)vtk_hdf_htg_data->descriptors_size;
        hsize_t global_size = local_size;

#endif

        // Pointer to the actual data to write
        Bit_t *descriptors = vtk_hdf_htg_data->descriptors;

        // Create a 1D array dimensions with unlimited size to allow resizing during writes
        hsize_t dims_d[1] = {global_size};
        hsize_t maxdims_d[1] = {global_size};

        // H5Screate_simple: Create a simple dataspace for the descriptors dataset
        // Allows unlimited growth, enabling scalable parallel I/O
        vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
        if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // H5Pcreate: Create a dataset creation property list
        // Used to configure parameters for the new dataset (e.g., chunking, compression)
        vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // H5Pset_chunk: Define chunked storage layout with chunk size matching total dimensions
        // Chunking is required for datasets with unlimited dimensions in HDF5
        // Improves performance for parallel I/O and enables resizing
        if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, dims_d) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // H5Tcopy: Create a copy of the 8-bit unsigned integer little-endian datatype (H5T_STD_U8LE)
        // Stores descriptor information as individual bytes with platform-independent representation
        vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_STD_U8LE);
        if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // H5Dcreate2: Create a new dataset named "Descriptors" in the VTK-HDF group
        // Associates the dataspace, datatype, and creation properties defined above
        // This dataset will hold all HyperTree descriptor data
        vtk_hdf_htg.dset_id = H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "Descriptors", vtk_hdf_htg.dset_dtype_id,
                                         vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
        if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
        // H5Pcreate: Create a dataset transfer (I/O) property list for controlling write behavior
        vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // H5Pset_dxpl_mpio: Set collective MPI-IO mode for the transfer property
        // Enables all ranks to participate in coordinated, efficient parallel writes to the same file
        // COLLECTIVE mode optimizes MPI-IO by allowing the MPI library to coordinate I/O patterns
        if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // H5Dget_space: Retrieve the file dataspace from the dataset
        // Allows selection of the specific region in the file where this rank will write its data
        vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
        if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // H5Sselect_hyperslab: Select a contiguous region (hyperslab) in the file dataspace
        // Defines where this MPI rank's data will be written: starting at 'offset', for 'local_size' elements
        // Non-overlapping hyperslabs ensure each rank writes to a unique region
        hsize_t start[1] = {offset};
        hsize_t count[1] = {local_size};
        if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        // H5Screate_simple: Create a simple memory dataspace matching the local data size
        // Describes the layout of data in the rank's local memory (a 1D array of 'local_size' elements)
        hsize_t m_dims[1] = {local_size};
        vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
        if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // H5Dwrite: Write local descriptor data to the HDF5 dataset
        // Maps data from memory (mem_space) to a specific region in the file (file_space) using collective MPI-IO
        // The xfer_plist ensures all ranks coordinate their writes for optimal performance
        // Significance: This is the actual parallel I/O operation that flushes descriptors to disk

        if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                     vtk_hdf_htg.dset_dtype_id, /* H5T_STD_U8LE: 8-bit unsigned integer type */
                     vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
                     vtk_hdf_htg.file_space,    /* file dataspace with hyperslab selected */
                     vtk_hdf_htg.xfer_plist,    /* collective MPI‐IO transfer property */
                     descriptors                /* pointer to local data */
                     ) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        // Resource cleanup: Close all HDF5 objects in reverse order of creation
        // This prevents resource leaks and ensures all data is properly flushed to disk
        H5Pclose(vtk_hdf_htg.xfer_plist);
        vtk_hdf_htg.xfer_plist = H5I_INVALID_HID;
        H5Sclose(vtk_hdf_htg.mem_space);
        vtk_hdf_htg.mem_space = H5I_INVALID_HID;
        H5Sclose(vtk_hdf_htg.file_space);
        vtk_hdf_htg.file_space = H5I_INVALID_HID;
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#else
        // Single-process write path: Use H5S_ALL to select the entire dataspace
        // Simpler than the parallel path since there's no need for hyperslab selection or MPI coordination
        if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, descriptors) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Resource cleanup: Close all HDF5 objects created for this dataset
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#endif
    }

    /*
     * Dataset: /VTKHDF/DescriptorsSize
     * Datatype: H5T_NATIVE_INT64 / int64_t
     * Dimension: {1*npe()}
     *
     * This is the number of descriptors (the individual bits, rather than the number of bytes they pack into) of the
     * tree on each process. Since we have only one tree in basilisk, the dimensions of the dataset is always
     * `{1*npe()}`.
     */
    {
#if MPI_SINGLE_FILE
        hsize_t local_size = (hsize_t)1;
        hsize_t global_size = (hsize_t)npe();
        hsize_t offset = pid() * local_size;
#else
        hsize_t local_size = (hsize_t)1;
        hsize_t global_size = (hsize_t)local_size;

#endif

        // The actual value of this locally is found from @ref vtkHDFHyperTreeGridData
        int64_t descriptors_size = vtk_hdf_htg_data->n_descriptors;

        // Set the dimension of the dataset
        hsize_t dims_d[1] = {global_size};

        // Set the maximum size of the dataset
        hsize_t maxdims_d[1] = {global_size};

        // Create the dataspace
        vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
        if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Create property list
        vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Set chunking size on the property list
        hsize_t chunk_dims[1] = {local_size};
        if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Set datatype for the dataset
        vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
        if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Create the dataset
        vtk_hdf_htg.dset_id =
            H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "DescriptorsSize", vtk_hdf_htg.dset_dtype_id,
                       vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
        if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE

        vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
        if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t start[1] = {offset};
        hsize_t count[1] = {local_size};
        if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        hsize_t m_dims[1] = {local_size};
        vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
        if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                     vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
                     vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
                     vtk_hdf_htg.file_space,    /* file dataspace with hyperslab selected */
                     vtk_hdf_htg.xfer_plist,    /* collective MPI‐IO transfer property */
                     &descriptors_size          /* pointer to local data */
                     ) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        H5Pclose(vtk_hdf_htg.xfer_plist);
        vtk_hdf_htg.xfer_plist = -1;
        H5Sclose(vtk_hdf_htg.mem_space);
        vtk_hdf_htg.mem_space = -1;
        H5Sclose(vtk_hdf_htg.file_space);
        vtk_hdf_htg.file_space = -1;
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#else
        if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &descriptors_size) <
            0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Dclose(vtk_hdf_htg.dset_id);
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
    }

    /*
     * Dataset: /VTKHDF/Mask
     * Datatype: H5T_STD_U8LE / UInt8_t
     * Dimension: {global_size}
     */
    if (vtk_hdf_htg_data->has_mask) {
#if MPI_SINGLE_FILE

        // Compute the offset by
        hsize_t local_size = (hsize_t)vtk_hdf_htg_data->mask_size;
        hsize_t global_size = local_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        hsize_t offset = 0;
        MPI_Exscan(&local_size, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (pid() == 0)
            offset = 0;

        hsize_t chunk_size = (hsize_t)(global_size / npe());
#else
        hsize_t local_size = (hsize_t)vtk_hdf_htg_data->mask_size;
        hsize_t global_size = local_size;
        hsize_t chunk_size = local_size;

#endif

        Bit_t *mask = vtk_hdf_htg_data->mask;
        hsize_t dims_d[1] = {global_size};
        hsize_t maxdims_d[1] = {global_size};

        vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
        if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t chunk_dims[1] = {chunk_size};
        if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_UINT8);
        if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_id = H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "Mask", vtk_hdf_htg.dset_dtype_id,
                                         vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
        if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#if MPI_SINGLE_FILE
        vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
        if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t start[1] = {offset};
        hsize_t count[1] = {local_size};
        if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        hsize_t m_dims[1] = {local_size};
        vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
        if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                     vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
                     vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
                     vtk_hdf_htg.file_space,
                     /* file dataspace with hyperslab selected
                      */
                     vtk_hdf_htg.xfer_plist,
                     /* collective MPI‐IO transfer property
                      */
                     mask /* pointer to local data */) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        H5Pclose(vtk_hdf_htg.xfer_plist);
        vtk_hdf_htg.xfer_plist = -1;
        H5Sclose(vtk_hdf_htg.mem_space);
        vtk_hdf_htg.mem_space = -1;
        H5Sclose(vtk_hdf_htg.file_space);
        vtk_hdf_htg.file_space = -1;
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#else
        if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, mask) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Dclose(vtk_hdf_htg.dset_id);
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
    }

    /*
     * Dataset: /VTKHDF/NumberOfCells
     * Datatype: H5T_NATIVE_INT64 / int64_t
     * Dimension: {1*npe()}
     *
     * Here we write the number of cells in the tree.
     */
    {
#if MPI_SINGLE_FILE
        hsize_t local_size = (hsize_t)1;
        hsize_t global_size = (hsize_t)npe();
        hsize_t offset = pid() * local_size;
#else
        hsize_t local_size = (hsize_t)1;
        hsize_t global_size = (hsize_t)local_size;

#endif

        int64_t *number_of_cells = &vtk_hdf_htg_data->number_of_cells;
        hsize_t dims_d[1] = {global_size};
        hsize_t maxdims_d[1] = {global_size};

        vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
        if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t chunk_dims[1] = {local_size};
        if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
        if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_id = H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "NumberOfCells", vtk_hdf_htg.dset_dtype_id,
                                         vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
        if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#if MPI_SINGLE_FILE
        vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
        if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t start[1] = {offset};
        hsize_t count[1] = {local_size};
        if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        hsize_t m_dims[1] = {local_size};
        vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
        if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                     vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
                     vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
                     vtk_hdf_htg.file_space,    /* file dataspace with hyperslab selected */
                     vtk_hdf_htg.xfer_plist,    /* collective MPI‐IO transfer property */
                     number_of_cells            /* pointer to local data */
                     ) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        H5Pclose(vtk_hdf_htg.xfer_plist);
        vtk_hdf_htg.xfer_plist = -1;
        H5Sclose(vtk_hdf_htg.mem_space);
        vtk_hdf_htg.mem_space = -1;
        H5Sclose(vtk_hdf_htg.file_space);
        vtk_hdf_htg.file_space = -1;
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#else
        if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, number_of_cells) <
            0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Dclose(vtk_hdf_htg.dset_id);
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
    }

    /*
     * Dataset: /VTKHDF/NumberOfCellsPerTreeDepth
     * Datatype: H5T_NATIVE_INT64 / int64_t
     * Dimension: {global_size}
     *
     * Here we write the number of cells at each level/depth in the tree. For example the tree
     *
     * [Level 0]                      o
     *                               / \
     *                              /   \
     * [Level 1]                   o     o
     *                            / \
     *                           /   \
     * [Level 2]                o     o
     *
     * would write the array `int64_t depth_per_tree = {1 2 2}`.
     */
    {
#if MPI_SINGLE_FILE
        hsize_t local_size = (hsize_t)vtk_hdf_htg_data->depth_per_tree;
        hsize_t global_size = local_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        hsize_t offset = 0;
        MPI_Exscan(&local_size, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (pid() == 0)
            offset = 0;
#else
        hsize_t local_size = (hsize_t)vtk_hdf_htg_data->depth_per_tree;
        hsize_t global_size = local_size;

#endif

        int64_t *number_of_cells_per_tree_depth = vtk_hdf_htg_data->number_of_cells_per_tree_depth;
        hsize_t dims_d[1] = {global_size};
        hsize_t maxdims_d[1] = {global_size};

        vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
        if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t chunk_dims[1] = {local_size};
        if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
        if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_id =
            H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "NumberOfCellsPerTreeDepth", vtk_hdf_htg.dset_dtype_id,
                       vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
        if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#if MPI_SINGLE_FILE
        vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
        if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t start[1] = {offset};
        hsize_t count[1] = {local_size};
        if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        hsize_t m_dims[1] = {local_size};
        vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
        if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Dwrite(vtk_hdf_htg.dset_id,           /* dataset handle */
                     vtk_hdf_htg.dset_dtype_id,     /* H5T_IEEE_F64LE */
                     vtk_hdf_htg.mem_space,         /* memory dataspace [local_nx] */
                     vtk_hdf_htg.file_space,        /* file dataspace with hyperslab selected */
                     vtk_hdf_htg.xfer_plist,        /* collective MPI‐IO transfer property */
                     number_of_cells_per_tree_depth /* pointer to local data */
                     ) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        H5Pclose(vtk_hdf_htg.xfer_plist);
        vtk_hdf_htg.xfer_plist = -1;
        H5Sclose(vtk_hdf_htg.mem_space);
        vtk_hdf_htg.mem_space = -1;
        H5Sclose(vtk_hdf_htg.file_space);
        vtk_hdf_htg.file_space = -1;
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#else
        if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     number_of_cells_per_tree_depth) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Dclose(vtk_hdf_htg.dset_id);
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
    }

    /*
     * Dataset: /VTKHDF/XCoordinates
     * Datatype: H5T_IEEE_F64LE / double
     * Dimension: {2 * npe()}
     *
     * In the MPI_SINGLE_FILE case we must (redudantly) write this for each proccess.
     */
    {
#if MPI_SINGLE_FILE
        hsize_t local_nx = (hsize_t)vtk_hdf_htg_data->n_x;
        hsize_t global_nx = (hsize_t)npe() * local_nx;
        hsize_t offset_nx = pid() * local_nx;
#else
        hsize_t local_nx = (hsize_t)vtk_hdf_htg_data->n_x;
        hsize_t global_nx = (hsize_t)local_nx; 
#endif

        double *xc_data = vtk_hdf_htg_data->x;
        hsize_t dims_d[1] = {global_nx};
        hsize_t maxdims_d[1] = {global_nx};

        vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
        if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        hsize_t chunk_dims[1] = {local_nx};
        if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_IEEE_F64LE);
        if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_id = H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "XCoordinates", vtk_hdf_htg.dset_dtype_id,
                                         vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
        if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
        vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
        if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t start[1] = {offset_nx};
        hsize_t count[1] = {local_nx};
        if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        hsize_t m_dims[1] = {local_nx};
        vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
        if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                     vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
                     vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
                     vtk_hdf_htg.file_space,    /* file dataspace with hyperslab selected */
                     vtk_hdf_htg.xfer_plist,    /* collective MPI‐IO transfer property */
                     xc_data                    /* pointer to local data */
                     ) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        H5Pclose(vtk_hdf_htg.xfer_plist);
        vtk_hdf_htg.xfer_plist = -1;
        H5Sclose(vtk_hdf_htg.mem_space);
        vtk_hdf_htg.mem_space = -1;
        H5Sclose(vtk_hdf_htg.file_space);
        vtk_hdf_htg.file_space = -1;
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#else
        if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, xc_data) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Dclose(vtk_hdf_htg.dset_id);
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
    }

    /*
     * Dataset: /VTKHDF/YCoordinates
     * Datatype: H5T_IEEE_F64LE / double
     * Dimension: {(2 or 1) * npe()}
     *
     * In the MPI_SINGLE_FILE case we must (redudantly) write this for each proccess.
     */
    {
#if MPI_SINGLE_FILE
        hsize_t local_ny = (hsize_t)vtk_hdf_htg_data->n_y;
        hsize_t global_ny = (hsize_t)npe() * local_ny;
        hsize_t offset_ny = pid() * local_ny;
#else
        hsize_t local_ny = (hsize_t)vtk_hdf_htg_data->n_y;
        hsize_t global_ny = (hsize_t)local_ny; 
#endif

        double *yc_data = vtk_hdf_htg_data->y;
        hsize_t dims_d[1] = {global_ny};
        hsize_t maxdims_d[1] = {global_ny};

        vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
        if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t chunk_dims[1] = {local_ny};
        if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_IEEE_F64LE);
        if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_id = H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "YCoordinates", vtk_hdf_htg.dset_dtype_id,
                                         vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
        if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
        vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
        if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t start[1] = {offset_ny};
        hsize_t count[1] = {local_ny};
        if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        hsize_t m_dims[1] = {local_ny};
        vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
        if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                     vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
                     vtk_hdf_htg.mem_space,     /* memory dataspace [local_ny] */
                     vtk_hdf_htg.file_space,    /* file dataspace with hyperslab selected */
                     vtk_hdf_htg.xfer_plist,    /* collective MPI‐IO transfer property */
                     yc_data                    /* pointer to local data */
                     ) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        H5Pclose(vtk_hdf_htg.xfer_plist);
        vtk_hdf_htg.xfer_plist = -1;
        H5Sclose(vtk_hdf_htg.mem_space);
        vtk_hdf_htg.mem_space = -1;
        H5Sclose(vtk_hdf_htg.file_space);
        vtk_hdf_htg.file_space = -1;
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#else
        if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, yc_data) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Dclose(vtk_hdf_htg.dset_id);
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
    }

    /*
     * Dataset: /VTKHDF/ZCoordinates
     * Datatype: H5T_IEEE_F64LE / double
     * Dimension: {(2 or 1) * npe()}
     *
     * In the MPI_SINGLE_FILE case we must (redudantly) write this for each proccess.
     */
    {
#if MPI_SINGLE_FILE
        hsize_t local_nz = (hsize_t)vtk_hdf_htg_data->n_z;
        hsize_t global_nz = (hsize_t)npe() * local_nz;
        hsize_t offset_nz = pid() * local_nz;
#else
        hsize_t local_nz = (hsize_t)vtk_hdf_htg_data->n_z;
        hsize_t global_nz = (hsize_t)local_nz; 
#endif

        double *zc_data = vtk_hdf_htg_data->z;
        hsize_t dims_d[1] = {global_nz};
        hsize_t maxdims_d[1] = {global_nz};

        vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
        if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t chunk_dims[1] = {local_nz};
        if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_IEEE_F64LE);
        if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_id = H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "ZCoordinates", vtk_hdf_htg.dset_dtype_id,
                                         vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
        if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
        vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
        if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t start[1] = {offset_nz};
        hsize_t count[1] = {local_nz};
        if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        hsize_t m_dims[1] = {local_nz};
        vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
        if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                     vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
                     vtk_hdf_htg.mem_space,     /* memory dataspace [local_nz] */
                     vtk_hdf_htg.file_space,    /* file dataspace with hyperslab selected */
                     vtk_hdf_htg.xfer_plist,    /* collective MPI‐IO transfer property */
                     zc_data                    /* pointer to local data */
                     ) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        H5Pclose(vtk_hdf_htg.xfer_plist);
        vtk_hdf_htg.xfer_plist = -1;
        H5Sclose(vtk_hdf_htg.mem_space);
        vtk_hdf_htg.mem_space = -1;
        H5Sclose(vtk_hdf_htg.file_space);
        vtk_hdf_htg.file_space = -1;
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#else
        if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, zc_data) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Dclose(vtk_hdf_htg.dset_id);
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
    }

    /*
     * Dataset: /VTKHDF/NumberOfDepths
     * Datatype: H5T_NATIVE_INT64 / int64_t
     * Dimension: {1 * npe()}
     *
     * The actual maximum depth per tree.
     */
    {
#if MPI_SINGLE_FILE
        hsize_t local_size = (hsize_t)1;
        hsize_t global_size = (hsize_t)npe();
        hsize_t offset = pid() * local_size;
#else
        hsize_t local_size = (hsize_t)1;
        hsize_t global_size = (hsize_t)local_size;
#endif

        int64_t *number_of_depths = &vtk_hdf_htg_data->depth_per_tree;
        hsize_t dims_d[1] = {global_size};
        hsize_t maxdims_d[1] = {global_size};

        vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
        if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t chunk_dims[1] = {local_size};
        if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
        if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_id = H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "NumberOfDepths", vtk_hdf_htg.dset_dtype_id,
                                         vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
        if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#if MPI_SINGLE_FILE
        vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
        if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t start[1] = {offset};
        hsize_t count[1] = {local_size};
        if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        hsize_t m_dims[1] = {local_size};
        vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
        if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                     vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
                     vtk_hdf_htg.mem_space,     /* memory dataspace [local_size] */
                     vtk_hdf_htg.file_space,    /* file dataspace with hyperslab selected */
                     vtk_hdf_htg.xfer_plist,    /* collective MPI‐IO transfer property */
                     number_of_depths           /* pointer to local data */
                     ) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        H5Pclose(vtk_hdf_htg.xfer_plist);
        vtk_hdf_htg.xfer_plist = -1;
        H5Sclose(vtk_hdf_htg.mem_space);
        vtk_hdf_htg.mem_space = -1;
        H5Sclose(vtk_hdf_htg.file_space);
        vtk_hdf_htg.file_space = -1;
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#else
        if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, number_of_depths) <
            0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Dclose(vtk_hdf_htg.dset_id);
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
    }

    /*
     * Dataset: /VTKHDF/NumberOfTrees
     * Datatype: H5T_NATIVE_INT64 / int64_t
     * Dimension: {1 * npe()}
     *
     * The number of trees. Since basilisk can only ever have one tree, we create a dataset with just 1.
     */
    {
#if MPI_SINGLE_FILE
        hsize_t local_size = (hsize_t)1;
        hsize_t global_size = (hsize_t)npe();
        hsize_t offset = pid() * local_size;
#else
        hsize_t local_size = (hsize_t)1;
        hsize_t global_size = (hsize_t)local_size;

#endif

        int64_t number_of_trees = vtk_hdf_htg_data->number_of_trees;
        hsize_t dims_d[1] = {global_size};
        hsize_t maxdims_d[1] = {global_size};

        vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
        if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t chunk_dims[1] = {local_size};
        if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
        if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_id = H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "NumberOfTrees", vtk_hdf_htg.dset_dtype_id,
                                         vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
        if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
        // Set up dataset transfer property list for collective I/O
        vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Create memory dataspace for a single scalar
        hsize_t m_dims[1] = {1};
        vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
        if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        // Create and select a hyperslab
        vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
        if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t start[1] = {offset};
        hsize_t count[1] = {1};
        if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        // Write
        if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                     vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
                     vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
                     vtk_hdf_htg.file_space,    /* file dataspace with hyperslab selected */
                     vtk_hdf_htg.xfer_plist,    /* collective MPI‐IO transfer property */
                     &number_of_trees           /* pointer to local data */
                     ) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        H5Pclose(vtk_hdf_htg.xfer_plist);
        vtk_hdf_htg.xfer_plist = -1;
        H5Sclose(vtk_hdf_htg.mem_space);
        vtk_hdf_htg.mem_space = -1;
        H5Sclose(vtk_hdf_htg.file_space);
        vtk_hdf_htg.file_space = -1;
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#else
        if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &number_of_trees) <
            0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Dclose(vtk_hdf_htg.dset_id);
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
    }

    /*
     * Dataset: /VTKHDF/TreeIds
     * Datatype: H5T_NATIVE_INT64 / int64_t
     * Dimension: {1 * npe()}
     *
     * The TreeIds of each tree. Again, since basilisk has only one tree, the TreeId is always zero. However, the
     * ordering in the case of a forest-of-trees has implications in informing VTK which tree belongs to which cell of
     * the regular grid described in X/Y/Z coordinates.
     */
    {
#if MPI_SINGLE_FILE
        hsize_t local_size = (hsize_t)1;
        hsize_t global_size = (hsize_t)npe();
        hsize_t offset = pid() * local_size;
#else
        hsize_t local_size = (hsize_t)1;
        hsize_t global_size = (hsize_t)local_size;
#endif

        int64_t *tree_ids = &vtk_hdf_htg_data->tree_ids;
        hsize_t dims_d[1] = {global_size};
        hsize_t maxdims_d[1] = {global_size};

        vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
        if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t chunk_dims[1] = {local_size};
        if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
        if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.dset_id = H5Dcreate2(vtk_hdf_htg.vtk_hdf.grp_vtkhdf_id, "TreeIds", vtk_hdf_htg.dset_dtype_id,
                                         vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
        if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
        vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
        if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        hsize_t start[1] = {offset};
        hsize_t count[1] = {local_size};
        if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        hsize_t m_dims[1] = {local_size};
        vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
        if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
        if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                     vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
                     vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
                     vtk_hdf_htg.file_space,    /* file dataspace with hyperslab selected */
                     vtk_hdf_htg.xfer_plist,    /* collective MPI‐IO transfer property */
                     tree_ids                   /* pointer to local data */
                     ) < 0) {
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
        }

        H5Pclose(vtk_hdf_htg.xfer_plist);
        vtk_hdf_htg.xfer_plist = -1;
        H5Sclose(vtk_hdf_htg.mem_space);
        vtk_hdf_htg.mem_space = -1;
        H5Sclose(vtk_hdf_htg.file_space);
        vtk_hdf_htg.file_space = -1;
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
        H5Dclose(vtk_hdf_htg.dset_id);
#else
        if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, tree_ids) < 0)
            vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

        H5Dclose(vtk_hdf_htg.dset_id);
        H5Tclose(vtk_hdf_htg.dset_dtype_id);
        H5Pclose(vtk_hdf_htg.dcpl_id);
        H5Sclose(vtk_hdf_htg.dset_space_id);
#endif
    }

    /*
     * Dataset: /VTKHDF/CellData/*
     * Datatype: H5T_IEEE_F32LE / float
     * Dimension: {global_size} (Scalars) or {global_size,dimension} (Vectors)
     *
     * In the remaining part of the function, we write both scalar fields and vector fields defined on the tree. These
     * must be written in breadth-first search order--the same order as the descriptors.
     *
     * In the MPI_SINGLE_FILE case, we write a PHyperTreeGrid, which neccesitates a short MPI exchange to determine
     * offsets. Since the fields are all described by the same tree, we only need to do this once for all of the vector
     * and scalar fields.
     *
     * In addition, we will use this information to choose a chunk size. This chunk size (measured in the number of
     * elements) should be sized such that we minimize the instances in which ranks/processes write into the same chunk.
     * Thus we set the chunk size to be the average of the local number of cells.
     */
    {
#if MPI_SINGLE_FILE
        hsize_t local_size = (hsize_t)vtk_hdf_htg_data->number_of_cells;
        hsize_t global_size = local_size;
        MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        hsize_t offset = 0;
        MPI_Exscan(&local_size, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        if (pid() == 0) {
            offset = 0;
        }
        hsize_t chunk_size = global_size / npe();
#else
        hsize_t local_size = (hsize_t)vtk_hdf_htg_data->number_of_cells;
        hsize_t global_size = local_size;
#endif

        /*
         * Scalars
         */
        for (scalar s in scalar_list) {
            float *s_data = malloc(sizeof(float) * local_size);
            size_t vi = 0;
            foreach_cell_BFS() {
                bool write = false;
                if (is_leaf(cell)) {
                    if (is_local(cell)) {
                        write = true;
                    } else {
                        foreach_neighbor(1) {
                          if (is_local(cell))
                            write = true;
                        }
                    }
                }
                s_data[vi++] = write ? (float)val(s) : 0.0;
            }

            hsize_t dims_d[1] = {global_size};
            hsize_t maxdims_d[1] = {global_size};

            vtk_hdf_htg.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
            if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
            if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            hsize_t chunk_dims[1] = {chunk_size};
            if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 1, chunk_dims) < 0)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if COMPRESSION
            if (H5Pset_deflate(vtk_hdf_htg.dcpl_id, COMPRESSION_LEVEL) < 0)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#endif

            vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_IEEE_F32LE);
            if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            vtk_hdf_htg.dset_id = H5Dcreate2(vtk_hdf_htg.grp_celldata_id, s.name, vtk_hdf_htg.dset_dtype_id,
                                             vtk_hdf_htg.dset_space_id, H5P_DEFAULT, vtk_hdf_htg.dcpl_id, H5P_DEFAULT);
            if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#if MPI_SINGLE_FILE
            vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
            if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            hsize_t start[1] = {offset};
            hsize_t count[1] = {local_size};
            if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start, NULL, count, NULL) < 0) {
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
            }

            hsize_t m_dims[1] = {local_size};
            vtk_hdf_htg.mem_space = H5Screate_simple(1, m_dims, NULL);
            if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
            if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
            if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                         vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F64LE */
                         vtk_hdf_htg.mem_space,     /* memory dataspace [local_nx] */
                         vtk_hdf_htg.file_space,    /* file dataspace with hyperslab selected */
                         vtk_hdf_htg.xfer_plist,    /* collective MPI‐IO transfer property */
                         s_data                     /* pointer to local data */
                         ) < 0) {
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
            }

            H5Pclose(vtk_hdf_htg.xfer_plist);
            vtk_hdf_htg.xfer_plist = -1;
            H5Sclose(vtk_hdf_htg.mem_space);
            vtk_hdf_htg.mem_space = -1;
            H5Sclose(vtk_hdf_htg.file_space);
            vtk_hdf_htg.file_space = -1;
            H5Tclose(vtk_hdf_htg.dset_dtype_id);
            H5Pclose(vtk_hdf_htg.dcpl_id);
            H5Sclose(vtk_hdf_htg.dset_space_id);
            H5Dclose(vtk_hdf_htg.dset_id);
#else
            if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, s_data) < 0)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            H5Dclose(vtk_hdf_htg.dset_id);
            H5Tclose(vtk_hdf_htg.dset_dtype_id);
            H5Pclose(vtk_hdf_htg.dcpl_id);
            H5Sclose(vtk_hdf_htg.dset_space_id);
#endif

            free(s_data);

        } /* end of “for (scalar s in scalar_list)” */

        /*
         * Vectors
         */
        for (vector v in vector_list) {

            char *vector_name;
            size_t trunc_len = (size_t)(strlen(v.x.name) - 2);
            vector_name = malloc((trunc_len + 1) * sizeof(char));
            strncpy(vector_name, v.x.name, trunc_len);
            vector_name[trunc_len] = '\0';

#if MPI_SINGLE_FILE
            hsize_t local_size = (hsize_t)vtk_hdf_htg_data->number_of_cells;
            hsize_t global_size;
            MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

            hsize_t offset = 0;
            MPI_Exscan(&local_size, &offset, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
            if (pid() == 0)
                offset = 0;
#else
            hsize_t local_size = (hsize_t)vtk_hdf_htg_data->number_of_cells;
            hsize_t global_size = local_size;

#endif

            /* 1) allocate one big flat buffer */
            float *v_data = malloc(local_size * dimension * sizeof(float));
            if (!v_data) {
                perror("malloc(v_data)");
                exit(1);
            }

            /* 2) fill it row‐major: every cell’s D components */
            size_t vi = 0;
            foreach_cell_BFS() {
                bool write = false;
                if (is_leaf(cell)) {
                    if (is_local(cell)) {
                        write = true;
                    } else {
                        foreach_neighbor(1) {
                            if (is_local(cell))
                                write = true;
                        }
                    }
                }

#if dimension == 1
                v_data[vi * dimension + 0] = write ? (float)val(v.x) : 0.0;
#elif dimension == 2
                v_data[vi * dimension + 0] = write ? (float)val(v.x) : 0.0;
                v_data[vi * dimension + 1] = write ? (float)val(v.y) : 0.0;
#else // dimension == 3
                v_data[vi * dimension + 0] = write ? (float)val(v.x) : 0.0;
                v_data[vi * dimension + 1] = write ? (float)val(v.y) : 0.0;
                v_data[vi * dimension + 2] = write ? (float)val(v.z) : 0.0;

#endif
                vi++;
            }

            /* 3) create a rank‐2 dataspace { global_size, D } */
            hsize_t dims_d[2] = {global_size, (hsize_t)dimension};
            hsize_t maxdims_d[2] = {global_size, (hsize_t)dimension};

            vtk_hdf_htg.dset_space_id = H5Screate_simple(2, dims_d, maxdims_d);
            if (vtk_hdf_htg.dset_space_id == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            /* 4) create a chunked property list, rank=2 with 2‐element chunk dims */
            vtk_hdf_htg.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
            if (vtk_hdf_htg.dcpl_id == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            hsize_t chunk_dims[2] = {chunk_size, (hsize_t)dimension};
            if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 2, chunk_dims) < 0)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            // if (H5Pset_chunk(vtk_hdf_htg.dcpl_id, 2, dims_d) < 0)
            //   vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if COMPRESSION
            if (H5Pset_deflate(vtk_hdf_htg.dcpl_id, COMPRESSION_LEVEL) < 0)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
#endif

            vtk_hdf_htg.dset_dtype_id = H5Tcopy(H5T_IEEE_F32LE);
            if (vtk_hdf_htg.dset_dtype_id == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            /* 6) create the dataset under CellData (collectively if MPI_SINGLE_FILE) */
            vtk_hdf_htg.dset_id =
                H5Dcreate2(vtk_hdf_htg.grp_celldata_id, vector_name, /* e.g. "u" or "velocity" */
                           vtk_hdf_htg.dset_dtype_id, vtk_hdf_htg.dset_space_id, H5P_DEFAULT, /* link creation */
                           vtk_hdf_htg.dcpl_id,                                               /* chunking */
                           H5P_DEFAULT /* default DCPL for filters */
                );
            if (vtk_hdf_htg.dset_id == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

#if MPI_SINGLE_FILE
            /* 7) select a 2D hyperslab [offset..offset+local_size-1] × [0..D-1] */
            vtk_hdf_htg.file_space = H5Dget_space(vtk_hdf_htg.dset_id);
            if (vtk_hdf_htg.file_space == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            hsize_t start2D[2] = {offset, 0};
            hsize_t count2D[2] = {local_size, (hsize_t)dimension};
            if (H5Sselect_hyperslab(vtk_hdf_htg.file_space, H5S_SELECT_SET, start2D, NULL, count2D, NULL) < 0) {
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
            }

            /* 8) create a 2D memory dataspace { local_size, D } */
            hsize_t m_dims[2] = {local_size, (hsize_t)dimension};
            vtk_hdf_htg.mem_space = H5Screate_simple(2, m_dims, NULL);
            if (vtk_hdf_htg.mem_space == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            /* 9) set up a collective transfer property list */
            vtk_hdf_htg.xfer_plist = H5Pcreate(H5P_DATASET_XFER);
            if (vtk_hdf_htg.xfer_plist == H5I_INVALID_HID)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
            if (H5Pset_dxpl_mpio(vtk_hdf_htg.xfer_plist, H5FD_MPIO_COLLECTIVE) < 0)
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);

            /* 10) write the slab */
            if (H5Dwrite(vtk_hdf_htg.dset_id,       /* dataset handle */
                         vtk_hdf_htg.dset_dtype_id, /* H5T_IEEE_F32LE */
                         vtk_hdf_htg.mem_space,     /* memory = {local_size,D} */
                         vtk_hdf_htg.file_space,    /* file = hyperslab of {global_size,D} */
                         vtk_hdf_htg.xfer_plist,    /* collective */
                         v_data                     /* contiguous float* */
                         ) < 0) {
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
            }

            /* 11) close everything in reverse order */
            H5Pclose(vtk_hdf_htg.xfer_plist);
            vtk_hdf_htg.xfer_plist = -1;
            H5Sclose(vtk_hdf_htg.mem_space);
            vtk_hdf_htg.mem_space = -1;
            H5Sclose(vtk_hdf_htg.file_space);
            vtk_hdf_htg.file_space = -1;
            H5Tclose(vtk_hdf_htg.dset_dtype_id);
            H5Pclose(vtk_hdf_htg.dcpl_id);
            H5Sclose(vtk_hdf_htg.dset_space_id);
            H5Dclose(vtk_hdf_htg.dset_id);
#else
            /* Serial (no MPI_SINGLE_FILE) version: write the entire dataset in one go
             */
            if (H5Dwrite(vtk_hdf_htg.dset_id, vtk_hdf_htg.dset_dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, v_data) < 0) {
                vtk_HDF_hypertreegrid_error(&vtk_hdf_htg);
            }

            H5Tclose(vtk_hdf_htg.dset_dtype_id);
            H5Pclose(vtk_hdf_htg.dcpl_id);
            H5Sclose(vtk_hdf_htg.dset_space_id);
            H5Dclose(vtk_hdf_htg.dset_id);
#endif

            free(v_data);
        } /* end of “for (vector v in vector_list)” */
    }

    vtk_hdf_hypertreegrid_data_free(vtk_hdf_htg_data);

    return vtk_hdf_htg;
}
