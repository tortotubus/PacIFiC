#include "vtkHDF.h"
#include <string.h>

typedef struct {
  /* Parent */
  vtkHDF vtk_hdf;

  hid_t grp_celldata_id;
  hid_t grp_steps_id;
  hid_t grp_celldataoffsets_id;

  hid_t attr_space_id;
  hid_t attr_dtype_id;
  hid_t dset_space_id;
  hid_t vspace_id;
  hid_t dcpl_id;
  hid_t dset_dtype_id;

  hid_t dset_id;
  hid_t attr_id;

  hid_t mem_space_id;

} vtkHDFHyperTreeGridTemporalVDS;

typedef struct {
  char  fname[1024];
  char  dname[1024];
  hid_t vspace;   /* selection in the (old) VDS */
  hid_t sspace;   /* selection in the source */
} Map;

static int append_source_to_vds(
  const char *src_file_name,
  const char *src_data_path,
  const char *dst_file_name,
  const char *dst_vdata_path
) {
  herr_t result;
  hid_t dst_file_id = -1, dst_dset_id = -1, dst_dtype_id = -1, dst_space_id = -1, dst_dcpl_id = -1;
  hid_t src_file_id = -1, src_dset_id = -1, src_dtype_id = -1, src_space_id = -1;

  Map *maps = NULL;
  size_t nsrc = 0;

  // Open the VDS and read the metadata
  dst_file_id = H5Fopen(dst_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
  dst_dset_id = H5Dopen2(dst_file_id, dst_vdata_path, H5P_DEFAULT);

  dst_dtype_id = H5Dget_type(dst_dset_id);
  dst_space_id = H5Dget_space(dst_dset_id);

  int rank = H5Sget_simple_extent_ndims(dst_space_id);
  if (rank <= 0) { fprintf(stderr, "VDS rank=%d unsupported.\n", rank); }

  hsize_t old_dims[H5S_MAX_RANK], old_max[H5S_MAX_RANK];
  result = H5Sget_simple_extent_dims(dst_space_id, old_dims, old_max);
  // if (result < 0) { assert(1==0); }

  dst_dcpl_id = H5Dget_create_plist(dst_dset_id);
  nsrc = H5Pget_virtual_count(dst_dcpl_id);

  maps = (Map*)calloc(nsrc, sizeof(Map));
  // if (!maps) { assert(1==0); }

  for (size_t i = 0; i < nsrc; i++) {
    maps[i].vspace = H5Pget_virtual_vspace(dst_dcpl_id, i);
    maps[i].sspace = H5Pget_virtual_srcspace(dst_dcpl_id, i);
  }

  // 2) Inspect the new source dataset to append
  src_file_id = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
}

void vtk_HDF_hypertreegrid_tvds_close(vtkHDFHyperTreeGridTemporalVDS *vtk_hdf_htg_tvds) {
  if (vtk_hdf_htg_tvds->grp_celldata_id >= 0)
    H5Gclose(vtk_hdf_htg_tvds->grp_celldata_id);

  if (vtk_hdf_htg_tvds->grp_steps_id >= 0)
    H5Gclose(vtk_hdf_htg_tvds->grp_steps_id);

  if (vtk_hdf_htg_tvds->grp_celldataoffsets_id >= 0)
    H5Gclose(vtk_hdf_htg_tvds->grp_celldataoffsets_id);

  vtk_HDF_close(&vtk_hdf_htg_tvds->vtk_hdf);
}

void vtk_HDF_hypertreegrid_tvds_error(vtkHDFHyperTreeGridTemporalVDS *vtk_hdf_htg_tvds) {
  vtk_HDF_hypertreegrid_tvds_close(vtk_hdf_htg_tvds);
  assert(1==2);
}

vtkHDFHyperTreeGridTemporalVDS vtk_HDF_hypertreegrid_tvds_append(
  scalar *scalar_list,
  vector *vector_list,
  double time,
  const char *src_file_name,
  const char *dst_file_name
) {
  vtkHDFHyperTreeGridTemporalVDS vtk_hdf_htg_tvds = {
    .vtk_hdf = {
      .file_id = -1,
      .fapl = -1,
      .grp_vtkhdf_id = -1
    },
    .grp_celldata_id = -1,
    .grp_steps_id = -1,
    .attr_space_id = -1,
    .attr_dtype_id = -1,
    .attr_id = -1, 
    .dset_space_id = -1,
    .dcpl_id = -1,
    .dset_dtype_id = -1,
    .dset_id = -1,
    .mem_space_id = -1,
  };

  if (pid() == 0) {
    herr_t result;

    vtk_hdf_htg_tvds.vtk_hdf.file_id = H5Fopen(dst_file_name, H5F_ACC_RDWR, H5P_DEFAULT);
    if (vtk_hdf_htg_tvds.vtk_hdf.file_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.vtk_hdf.grp_vtkhdf_id = H5Gopen2(vtk_hdf_htg_tvds.vtk_hdf.file_id, "/VTKHDF", H5P_DEFAULT);
    if (vtk_hdf_htg_tvds.vtk_hdf.grp_vtkhdf_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.grp_steps_id = H5Gopen2(vtk_hdf_htg_tvds.vtk_hdf.grp_vtkhdf_id, "Steps", H5P_DEFAULT);
    if (vtk_hdf_htg_tvds.grp_steps_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    {
      double time_val = time;  // use the argument

      // Open dataset under /VTKHDF/Steps
      vtk_hdf_htg_tvds.dset_id = H5Dopen2(vtk_hdf_htg_tvds.grp_steps_id, "Values", H5P_DEFAULT);
      if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

      // Get current size
      vtk_hdf_htg_tvds.dset_space_id = H5Dget_space(vtk_hdf_htg_tvds.dset_id);
      if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

      int rank = H5Sget_simple_extent_ndims(vtk_hdf_htg_tvds.dset_space_id);
      if (rank != 1) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

      hsize_t dims[1], maxdims[1];
      if (H5Sget_simple_extent_dims(vtk_hdf_htg_tvds.dset_space_id, dims, maxdims) < 0)
          vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

      // Extend dataset by 1 (NOTE: must be chunked & unlimited from creation time)
      hsize_t new_dims[1] = { dims[0] + 1 };
      if (H5Dset_extent(vtk_hdf_htg_tvds.dset_id, new_dims) < 0)
          vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

      // Refresh filespace after the extent change
      H5Sclose(vtk_hdf_htg_tvds.dset_space_id);
      vtk_hdf_htg_tvds.dset_space_id = H5Dget_space(vtk_hdf_htg_tvds.dset_id);
      if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

      // Select the last element
      hsize_t start[1] = { new_dims[0] - 1 };
      hsize_t count[1] = { 1 };
      if (H5Sselect_hyperslab(vtk_hdf_htg_tvds.dset_space_id, H5S_SELECT_SET, start, NULL, count, NULL) < 0)
          vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

      // Memory dataspace for one value
      vtk_hdf_htg_tvds.mem_space_id = H5Screate_simple(1, count, NULL);
      if (vtk_hdf_htg_tvds.mem_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

      // Write the new data point
      result = H5Dwrite(
          vtk_hdf_htg_tvds.dset_id,
          H5T_IEEE_F64LE,                     // memory type (matches on-disk in your schema)
          vtk_hdf_htg_tvds.mem_space_id,      // memspace (1 element)
          vtk_hdf_htg_tvds.dset_space_id,     // filespace (selected last element)
          H5P_DEFAULT,
          &time_val
      );
      if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

      // Cleanup locals
      H5Sclose(vtk_hdf_htg_tvds.mem_space_id); vtk_hdf_htg_tvds.mem_space_id = -1;
      H5Sclose(vtk_hdf_htg_tvds.dset_space_id); vtk_hdf_htg_tvds.dset_space_id = -1;
      H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    }
  }

  return vtk_hdf_htg_tvds;
}

vtkHDFHyperTreeGridTemporalVDS vtk_HDF_hypertreegrid_tvds_init(
  scalar *scalar_list,
  vector *vector_list,
  double time,
  const char *src_file_name, 
  const char *dst_file_name
) {

  vtkHDF vtk_hdf = vtk_HDF_init(dst_file_name);

  vtkHDFHyperTreeGridTemporalVDS vtk_hdf_htg_tvds = {
    .vtk_hdf = vtk_hdf,
    .grp_celldata_id = -1,
    .grp_steps_id = -1,
    .attr_space_id = -1,
    .attr_dtype_id = -1,
    .attr_id = -1, 
    .dset_space_id = -1,
    .dcpl_id = -1,
    .dset_dtype_id = -1,
    .dset_id = -1,   
    .mem_space_id = -1,
  };

  // Attribute: /VTKHDF/BranchFactor 
  {
    int64_t bf_value = 2;
    hsize_t dims_attr[1] = {1};

    vtk_hdf_htg_tvds.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
    if (vtk_hdf_htg_tvds.attr_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.attr_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.attr_id = H5Acreate2(vtk_hdf_htg_tvds.vtk_hdf.grp_vtkhdf_id, "BranchFactor", vtk_hdf_htg_tvds.attr_dtype_id, vtk_hdf_htg_tvds.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg_tvds.attr_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    if (H5Awrite(vtk_hdf_htg_tvds.attr_id, vtk_hdf_htg_tvds.attr_dtype_id, &bf_value) < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Aclose(vtk_hdf_htg_tvds.attr_id);
    H5Tclose(vtk_hdf_htg_tvds.attr_dtype_id);
    H5Sclose(vtk_hdf_htg_tvds.attr_space_id);
  }

  // Attribute: /VTKHDF/Dimensions 
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
    vtk_hdf_htg_tvds.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
    if (vtk_hdf_htg_tvds.attr_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* Use a 64-bit‐int little‐endian type */
    vtk_hdf_htg_tvds.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.attr_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.attr_id = H5Acreate2(vtk_hdf_htg_tvds.vtk_hdf.grp_vtkhdf_id, "Dimensions", vtk_hdf_htg_tvds.attr_dtype_id, vtk_hdf_htg_tvds.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg_tvds.attr_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* Now actually write dims_value[3] into the attribute */
    if (H5Awrite(vtk_hdf_htg_tvds.attr_id, vtk_hdf_htg_tvds.attr_dtype_id, dims_value) < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Aclose(vtk_hdf_htg_tvds.attr_id);
    H5Tclose(vtk_hdf_htg_tvds.attr_dtype_id);
    H5Sclose(vtk_hdf_htg_tvds.attr_space_id);
  }

  // Attribute: /VTKHDF/TransposedRootIndexing
  {
    int64_t tri_value = 0;
    hsize_t dims_attr[1] = {1};

    vtk_hdf_htg_tvds.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
    if (vtk_hdf_htg_tvds.attr_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.attr_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.attr_id = H5Acreate2(vtk_hdf_htg_tvds.vtk_hdf.grp_vtkhdf_id, "TransposedRootIndexing", vtk_hdf_htg_tvds.attr_dtype_id, vtk_hdf_htg_tvds.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg_tvds.attr_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    if (H5Awrite(vtk_hdf_htg_tvds.attr_id, vtk_hdf_htg_tvds.attr_dtype_id, &tri_value) < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Aclose(vtk_hdf_htg_tvds.attr_id);
    H5Tclose(vtk_hdf_htg_tvds.attr_dtype_id);
    H5Sclose(vtk_hdf_htg_tvds.attr_space_id);
  }

  // Attribute: /VTKHDF/Type
  {
    const char *type_str = "HyperTreeGrid";
    // hsize_t scalar_dims = 1;

    vtk_hdf_htg_tvds.attr_space_id = H5Screate(H5S_SCALAR);
    if (vtk_hdf_htg_tvds.attr_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* Create a fixed-length string datatype of length 13, null-padded, ASCII */
    vtk_hdf_htg_tvds.attr_dtype_id = H5Tcopy(H5T_C_S1);
    if (vtk_hdf_htg_tvds.attr_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    if (H5Tset_size(vtk_hdf_htg_tvds.attr_dtype_id, (size_t)13) < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    if (H5Tset_strpad(vtk_hdf_htg_tvds.attr_dtype_id, H5T_STR_NULLPAD) < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    if (H5Tset_cset(vtk_hdf_htg_tvds.attr_dtype_id, H5T_CSET_ASCII) < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.attr_id = H5Acreate2(vtk_hdf_htg_tvds.vtk_hdf.grp_vtkhdf_id, "Type", vtk_hdf_htg_tvds.attr_dtype_id, vtk_hdf_htg_tvds.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg_tvds.attr_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* Write the string (automatically null‐padded up to length 13) */
    if (H5Awrite(vtk_hdf_htg_tvds.attr_id, vtk_hdf_htg_tvds.attr_dtype_id, type_str) < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Aclose(vtk_hdf_htg_tvds.attr_id);
    H5Tclose(vtk_hdf_htg_tvds.attr_dtype_id);
    H5Sclose(vtk_hdf_htg_tvds.attr_space_id);
  }

  // Attribute: /VTKHDF/Version
  {
    int64_t vers_value[2] = {2, 4};
    hsize_t dims_attr[1] = {2};

    vtk_hdf_htg_tvds.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
    if (vtk_hdf_htg_tvds.attr_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.attr_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.attr_id = H5Acreate2(vtk_hdf_htg_tvds.vtk_hdf.grp_vtkhdf_id, "Version", vtk_hdf_htg_tvds.attr_dtype_id, vtk_hdf_htg_tvds.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg_tvds.attr_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    if (H5Awrite(vtk_hdf_htg_tvds.attr_id, vtk_hdf_htg_tvds.attr_dtype_id, vers_value) < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Aclose(vtk_hdf_htg_tvds.attr_id);
    H5Tclose(vtk_hdf_htg_tvds.attr_dtype_id);
    H5Sclose(vtk_hdf_htg_tvds.attr_space_id);
  }

  // Virtual Dataset: /VTKHDF/DepthPerTree
  {
    const char *src_dset_name = "/VTKHDF/DepthPerTree";
    const char *dst_dset_name = "/VTKHDF/DepthPerTree";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);
  }

  // Virtual Dataset: /VTKHDF/Descriptors
  {
    const char *src_dset_name = "/VTKHDF/Descriptors";
    const char *dst_dset_name = "/VTKHDF/Descriptors";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);

  }

  // Virtual Dataset: /VTKHDF/DescriptorsSize
  {
    const char *src_dset_name = "/VTKHDF/DescriptorsSize";
    const char *dst_dset_name = "/VTKHDF/DescriptorsSize";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);

  }

  // Virtual Dataset: /VTKHDF/Mask
  {
    const char *src_dset_name = "/VTKHDF/Mask";
    const char *dst_dset_name = "/VTKHDF/Mask";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);

  }

  // Virtual Dataset: /VTKHDF/NumberOfCells
  {
    const char *src_dset_name = "/VTKHDF/NumberOfCells";
    const char *dst_dset_name = "/VTKHDF/NumberOfCells";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);

  }
  
  // Virtual Dataset: /VTKHDF/NumberOfCellsPerTreeDepth
  {
    const char *src_dset_name = "/VTKHDF/NumberOfCellsPerTreeDepth";
    const char *dst_dset_name = "/VTKHDF/NumberOfCellsPerTreeDepth";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);

  }

  // Virtual Dataset: /VTKHDF/NumberOfDepths
  {
    const char *src_dset_name = "/VTKHDF/NumberOfDepths";
    const char *dst_dset_name = "/VTKHDF/NumberOfDepths";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);

  }

  // Virtual Dataset: /VTKHDF/NumberOfTrees
  {
    const char *src_dset_name = "/VTKHDF/NumberOfTrees";
    const char *dst_dset_name = "/VTKHDF/NumberOfTrees";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);

  }

  // Virtual Dataset: /VTKHDF/TreeIds
  {
    const char *src_dset_name = "/VTKHDF/TreeIds";
    const char *dst_dset_name = "/VTKHDF/TreeIds";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);

  }

  // Virtual Dataset: /VTKHDF/XCoordinates
  {
    const char *src_dset_name = "/VTKHDF/XCoordinates";
    const char *dst_dset_name = "/VTKHDF/XCoordinates";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);

  }

  // Virtual Dataset: /VTKHDF/YCoordinates
  {
    const char *src_dset_name = "/VTKHDF/YCoordinates";
    const char *dst_dset_name = "/VTKHDF/YCoordinates";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);

  }

  // Virtual Dataset: /VTKHDF/ZCoordinates
  {
    const char *src_dset_name = "/VTKHDF/ZCoordinates";
    const char *dst_dset_name = "/VTKHDF/ZCoordinates";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);

  }

  // Group: /VTKHDF/CellData 
  {
    vtk_hdf_htg_tvds.grp_celldata_id = H5Gcreate2(vtk_hdf_htg_tvds.vtk_hdf.grp_vtkhdf_id, "CellData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg_tvds.grp_celldata_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
  }

  // Virtual Dataset(s): /VTKHDF/CellData (scalar)
  for (scalar s in scalar_list)
  {
    const char base_dset_name[] = "/VTKHDF/CellData";

    char src_dset_name[100] = "";
    char dst_dset_name[100] = "";

    sprintf(src_dset_name, "%s/%s", base_dset_name, s.name);
    sprintf(dst_dset_name, "%s/%s", base_dset_name, s.name);

    // const char *src_dset_name = "/VTKHDF/CellData/";
    // const char *dst_dset_name = "/VTKHDF/CellData/";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);
  }

  // Virtual Dataset(s): /VTKHDF/CellData (vector)
  for (vector v in vector_list)
  {
    const char base_dset_name[] = "/VTKHDF/CellData";

    char* vector_name;
    size_t trunc_len = (size_t)(strlen(v.x.name) - 2);
    vector_name = malloc((trunc_len + 1) * sizeof(char));
    strncpy( vector_name, v.x.name, trunc_len );
    vector_name[trunc_len] = '\0';

    char src_dset_name[100] = "";
    char dst_dset_name[100] = "";

    sprintf(src_dset_name, "%s/%s", base_dset_name, vector_name);
    sprintf(dst_dset_name, "%s/%s", base_dset_name, vector_name);

    // const char *src_dset_name = "/VTKHDF/CellData/";
    // const char *dst_dset_name = "/VTKHDF/CellData/";

    herr_t result;

    /* open source */
    hid_t file_src = H5Fopen(src_file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t data_src = H5Dopen2(file_src, src_dset_name, H5P_DEFAULT);
    if (data_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_src = H5Dget_type(data_src);
    if (type_src < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t type_cpy = H5Tcopy(type_src); /* datatype for the VDS */
    if (type_cpy < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    hid_t src_space_id = H5Dget_space(data_src);
    if (src_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    int rank = H5Sget_simple_extent_ndims(src_space_id);
    if (rank <= 0) {
      fprintf(stderr, "Dataset %s has rank %d (unsupported)\n", src_dset_name, rank);
      vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    }

    hsize_t dims[H5S_MAX_RANK], maxdims[H5S_MAX_RANK];
    result = H5Sget_simple_extent_dims(src_space_id, dims, maxdims);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* virtual dataspace matches source extents */
    vtk_hdf_htg_tvds.vspace_id = H5Screate_simple(rank, dims, NULL);
    if (vtk_hdf_htg_tvds.vspace_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* DCPL with one virtual mapping: whole VDS -> whole source */
    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Sselect_all(vtk_hdf_htg_tvds.vspace_id);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
    
    result = H5Pset_virtual(
      vtk_hdf_htg_tvds.dcpl_id,
      vtk_hdf_htg_tvds.vspace_id,         /* VDS selection (all) */
      src_file_name,                       /* source file */
      src_dset_name,                       /* source path */
      src_space_id                         /* source selection (all) */
    );

    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* create the VDS in the destination — use the file id and absolute path */
    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.vtk_hdf.file_id,       /* loc_id (file ok for abs path) */
      dst_dset_name,                          /* name (absolute path) */
      type_cpy,                               /* datatype */
      vtk_hdf_htg_tvds.vspace_id,             /* dataspace */
      H5P_DEFAULT,                            /* lcpl */
      vtk_hdf_htg_tvds.dcpl_id,               /* dcpl (holds VDS mapping) */
      H5P_DEFAULT                             /* dapl */
    );

    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    /* --- close locals and dest temporaries for this VDS --- */
    H5Dclose(vtk_hdf_htg_tvds.dset_id); vtk_hdf_htg_tvds.dset_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id); vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.vspace_id); vtk_hdf_htg_tvds.vspace_id = -1;

    H5Sclose(src_space_id);
    H5Tclose(type_cpy);
    H5Tclose(type_src);
    H5Dclose(data_src);
    H5Fclose(file_src);

    free(vector_name);
  }

  // Group: /VTKHDF/Steps
  {
    vtk_hdf_htg_tvds.grp_steps_id = H5Gcreate2(vtk_hdf_htg_tvds.vtk_hdf.grp_vtkhdf_id, "Steps", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg_tvds.grp_steps_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
  }

  // Attribute: /VTKHDF/Steps/NSteps
  {
    herr_t result;
    int64_t nsteps_value = 1;
    hsize_t dims_attr[1] = {1};

    vtk_hdf_htg_tvds.attr_space_id = H5Screate_simple(1, dims_attr, dims_attr);
    if (vtk_hdf_htg_tvds.attr_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.attr_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.attr_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.attr_id = H5Acreate2(vtk_hdf_htg_tvds.grp_steps_id, "NSteps", vtk_hdf_htg_tvds.attr_dtype_id, vtk_hdf_htg_tvds.attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg_tvds.attr_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Awrite(vtk_hdf_htg_tvds.attr_id, vtk_hdf_htg_tvds.attr_dtype_id, &nsteps_value);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Aclose(vtk_hdf_htg_tvds.attr_id);
    H5Tclose(vtk_hdf_htg_tvds.attr_dtype_id);
    H5Sclose(vtk_hdf_htg_tvds.attr_space_id);
  }

  // Dataset: /VTKHDF/Steps/Values
  {
    herr_t result;

    double time = 0.0;

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_IEEE_F64LE);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_steps_id, "Values",
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &time
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;
  }

  // Dataset: /VTKHDF/Steps/DepthPerTreeOffsets
  {
    herr_t result;

    int64_t offset = 0;

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_steps_id, "DepthPerTreeOffsets",
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &offset
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;
  }

  // Dataset: /VTKHDF/Steps/DescriptorsOffsets
  {
    herr_t result;

    int64_t offset = 0;

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_steps_id, "DescriptorsOffsets",
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &offset
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;
  }

  // Dataset: /VTKHDF/Steps/MaskOffsets
  {
    herr_t result;

    int64_t offset = 0;

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_steps_id, "MaskOffsets",
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &offset
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;
  }

  // Dataset: /VTKHDF/Steps/NumberOfCellsPerTreeDepthOffsets
  {
    herr_t result;

    int64_t offset = 0;

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_steps_id, "NumberOfCellsPerTreeDepthOffsets",
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &offset
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;
  }

  // Dataset: /VTKHDF/Steps/PartOffsets
  {
    herr_t result;

    int64_t offset = 0;

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_steps_id, "PartOffsets",
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &offset
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;
  }

  // Dataset: /VTKHDF/Steps/NumberOfParts
  {
    herr_t result;

    int64_t number_of_parts = (int64_t) npe();

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_steps_id, "NumberOfParts",
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &number_of_parts
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;
  }

  // Dataset: /VTKHDF/Steps/TreeIdsOffsets
  {
    herr_t result;

    int64_t offset = 0;

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_steps_id, "TreeIdsOffsets",
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &offset
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;
  }

  // Dataset: /VTKHDF/Steps/XCoordinatesOffsets
  {
    herr_t result;

    int64_t offset = 0;

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_steps_id, "XCoordinatesOffsets",
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &offset
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;
  }

  // Dataset: /VTKHDF/Steps/YCoordinatesOffsets
  {
    herr_t result;

    int64_t offset = 0;

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_steps_id, "YCoordinatesOffsets",
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &offset
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;
  }

  // Dataset: /VTKHDF/Steps/ZCoordinatesOffsets
  {
    herr_t result;

    int64_t offset = 0;

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_steps_id, "ZCoordinatesOffsets",
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &offset
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;
  }

  // Group: /VTKHDF/Steps/CellDataOffsets
  {
    vtk_hdf_htg_tvds.grp_celldataoffsets_id = H5Gcreate2(vtk_hdf_htg_tvds.grp_steps_id, "CellDataOffsets", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (vtk_hdf_htg_tvds.grp_celldataoffsets_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);
  }

  // Dataset: /VTKHDF/Steps/CellDataOffsets
  for (scalar s in scalar_list)
  {
    herr_t result;

    int64_t offset = 0;

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_celldataoffsets_id, s.name,
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &offset
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;
  }

  // Dataset: /VTKHDF/Steps/CellDataOffsets
  for (vector v in vector_list)
  {
    char* vector_name;
    size_t trunc_len = (size_t)(strlen(v.x.name) - 2);
    vector_name = malloc((trunc_len + 1) * sizeof(char));
    strncpy( vector_name, v.x.name, trunc_len );
    vector_name[trunc_len] = '\0';

    herr_t result;

    int64_t offset = 0;

    hsize_t dims_d[1] = {1};
    hsize_t maxdims_d[1] = {H5S_UNLIMITED};
    hsize_t chunk[1] = {1};

    vtk_hdf_htg_tvds.dset_space_id = H5Screate_simple(1, dims_d, maxdims_d);
    if (vtk_hdf_htg_tvds.dset_space_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    if (vtk_hdf_htg_tvds.dcpl_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Pset_chunk(vtk_hdf_htg_tvds.dcpl_id, 1, chunk);
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_dtype_id = H5Tcopy(H5T_NATIVE_INT64);
    if (vtk_hdf_htg_tvds.dset_dtype_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    vtk_hdf_htg_tvds.dset_id = H5Dcreate2(
      vtk_hdf_htg_tvds.grp_celldataoffsets_id, vector_name,
      vtk_hdf_htg_tvds.dset_dtype_id, vtk_hdf_htg_tvds.dset_space_id,
      H5P_DEFAULT, vtk_hdf_htg_tvds.dcpl_id, H5P_DEFAULT
    );
    if (vtk_hdf_htg_tvds.dset_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    result = H5Dwrite(
      vtk_hdf_htg_tvds.dset_id, vtk_hdf_htg_tvds.dset_dtype_id, H5S_ALL,
      H5S_ALL, H5P_DEFAULT, &offset
    );
    if (result < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

    H5Dclose(vtk_hdf_htg_tvds.dset_id);            vtk_hdf_htg_tvds.dset_id = -1;
    H5Tclose(vtk_hdf_htg_tvds.dset_dtype_id);      vtk_hdf_htg_tvds.dset_dtype_id = -1;
    H5Pclose(vtk_hdf_htg_tvds.dcpl_id);            vtk_hdf_htg_tvds.dcpl_id = -1;
    H5Sclose(vtk_hdf_htg_tvds.dset_space_id);      vtk_hdf_htg_tvds.dset_space_id = -1;

    free(vector_name);
  }

  return vtk_hdf_htg_tvds;
}