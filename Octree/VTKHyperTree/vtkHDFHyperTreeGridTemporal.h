#include "vtkHDF.h"

typedef struct {
  /* Parent */
  vtkHDF vtk_hdf;

  hid_t grp_celldata_id;
  hid_t grp_steps_id;
  hid_t grp_celldataoffsets_id;

  hid_t attr_space_id;
  hid_t attr_dtype_id;
  hid_t dset_space_id;
  hid_t dcpl_id;
  hid_t dset_dtype_id;

  hid_t dset_id;
  hid_t attr_id;

} vtkHDFHyperTreeGridTemporalVDS;

void vtk_HDF_hypertreegrid_tvds_close(vtkHDFHyperTreeGridTemporalVDS *vtk_hdf_htg_tvds) {

}

void vtk_HDF_hypertreegrid_tvds_error(vtkHDFHyperTreeGridTemporalVDS *vtk_hdf_htg_tvds) {

}

vtkHDFHyperTreeGridTemporalVDS vtk_HDF_hypertreegrid_tvds_init(const char *src_name, const char *out_name) {
  vtkHDF vtk_hdf = vtk_HDF_init(out_name);
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

  // Group: /VTKHDF/Steps
  vtk_hdf_htg_tvds.grp_steps_id = H5Gcreate2(vtk_hdf_htg_tvds.vtk_hdf.grp_vtkhdf_id, "Steps", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (vtk_hdf_htg_tvds.grp_steps_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

  // Group: /VTKHDF/Steps/CellDataOffsets group
  vtk_hdf_htg_tvds.grp_celldataoffsets_id = H5Gcreate2(vtk_hdf_htg_tvds.grp_steps_id, "CellDataOffsets", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (vtk_hdf_htg_tvds.grp_celldataoffsets_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

  // Group: /VTKHDF/CellData 
  vtk_hdf_htg_tvds.grp_celldata_id = H5Gcreate2(vtk_hdf_htg_tvds.vtk_hdf.grp_vtkhdf_id, "CellData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (vtk_hdf_htg_tvds.grp_celldata_id < 0) vtk_HDF_hypertreegrid_tvds_error(&vtk_hdf_htg_tvds);

}