#include "io/output-common.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "io/vtk/vtkHDFPolyData.h"
#include "io/vtk/vtkPolyData.h"
#include "io/vtk/vtkTimeSeries.h"

#include "lagrangian_caps_optim/capsule-manager-ft.h"

// ============================================================================
// Function Declarations
// ============================================================================

static inline bool capsule_vtk_should_write_mesh(CapsuleMesh *mesh);
static inline double capsule_vtk_clean_double(double value);
static inline void capsule_vtk_store_coord(double *data, size_t i, coord value);
void output_hdf_capsules_series(const char *basename, int iter, double time,
                                bool overwrite);
void output_hdf_capsules(const char *basename = NULL, int iter = i,
                         double time = t, bool overwrite = true);

// ============================================================================
// Function Definitions
// ============================================================================

/*! 
 * @brief 
 */
static inline bool capsule_vtk_should_write_mesh(CapsuleMesh *mesh) {
  if (!mesh || !mesh->isactive)
    return false;
#if _MPI
  return capsule_mesh_is_local_owner(mesh);
#else
  return true;
#endif
}

/*! 
 * @brief 
 */
static inline double capsule_vtk_clean_double(double value) {
  return isfinite(value) ? value : 0.;
}

/*! 
 * @brief 
 */
static inline void capsule_vtk_store_coord(double *data, size_t i,
                                           coord value) {
  data[3 * i + 0] = capsule_vtk_clean_double(value.x);
  data[3 * i + 1] = capsule_vtk_clean_double(value.y);
  data[3 * i + 2] = capsule_vtk_clean_double(value.z);
}

/*! 
 * @brief 
 */
trace void output_hdf_capsules_series(const char *basename, int iter,
                                      double time, bool overwrite) {
  char pname[1024];
  snprintf(pname, sizeof(pname), "%s/vtkhdf-capsules-series", basename);

  char fname[1024];
  snprintf(fname, sizeof(fname), "%s/vtkhdf-capsules-series/capsules_%d.vtkhdf",
           basename, iter);

  char series_entry[1024];
  snprintf(series_entry, sizeof(series_entry),
           "vtkhdf-capsules-series/capsules_%d.vtkhdf", iter);

  char series_filename[1024];
  snprintf(series_filename, sizeof(series_filename),
           "%s/capsules.vtkhdf.series", basename);

  if (pid() == 0) {
    assert(!create_path(basename));
    assert(!create_path(pname));
    output_vtk_series(series_entry, series_filename, iter, time);
  }
#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  size_t n_points = 0;
  size_t n_polygons = 0;
  for (int c = 0; c < CAPS_COUNT(); c++) {
    CapsuleMesh *mesh = &CAPS(c);
    if (!capsule_vtk_should_write_mesh(mesh))
      continue;
    n_points += (size_t)mesh->nln;
    n_polygons += (size_t)mesh->nlt;
  }

  const size_t n_pointdata = 10;
  const size_t n_celldata = 13;
  vtkPolyData vtk_capsules = vtk_polydata_init_with_celldata(
      n_points, 0, 0, 0, n_polygons, n_pointdata, n_celldata);

  for (int c = 0; c < CAPS_COUNT(); c++) {
    CapsuleMesh *mesh = &CAPS(c);
    if (!capsule_vtk_should_write_mesh(mesh))
      continue;

    for (int n = 0; n < mesh->nln; n++) {
      coord pos = mesh->nodes[n].pos;
      vtk_polydata_add_point(&vtk_capsules,
                             (float)capsule_vtk_clean_double(pos.x),
                             (float)capsule_vtk_clean_double(pos.y),
                             (float)capsule_vtk_clean_double(pos.z));
    }
  }

  size_t point_offset = 0;
  for (int c = 0; c < CAPS_COUNT(); c++) {
    CapsuleMesh *mesh = &CAPS(c);
    if (!capsule_vtk_should_write_mesh(mesh))
      continue;

    for (int tr = 0; tr < mesh->nlt; tr++) {
      TriangleTopology *topo = LAG_TRI_TOPO(mesh, tr);
      vtk_polydata_add_triangle(
          &vtk_capsules, (int64_t)(point_offset + (size_t)topo->node_ids[0]),
          (int64_t)(point_offset + (size_t)topo->node_ids[1]),
          (int64_t)(point_offset + (size_t)topo->node_ids[2]));
    }

    point_offset += (size_t)mesh->nln;
  }

  int64_t pd_cap_id =
      vtk_polydata_add_pointdata_scalar(&vtk_capsules, "cap_id");
  int64_t pd_cap_type =
      vtk_polydata_add_pointdata_scalar(&vtk_capsules, "cap_type");
  int64_t pd_owner_pid =
      vtk_polydata_add_pointdata_scalar(&vtk_capsules, "owner_pid");
  int64_t pd_lagvel =
      vtk_polydata_add_pointdata_vector(&vtk_capsules, "lagVel", 3);
  int64_t pd_lagforce =
      vtk_polydata_add_pointdata_vector(&vtk_capsules, "lagForce", 3);
  int64_t pd_normal =
      vtk_polydata_add_pointdata_vector(&vtk_capsules, "normal", 3);
  int64_t pd_curv = vtk_polydata_add_pointdata_scalar(&vtk_capsules, "curv");
  int64_t pd_gcurv = vtk_polydata_add_pointdata_scalar(&vtk_capsules, "gcurv");
  int64_t pd_ref_curv =
      vtk_polydata_add_pointdata_scalar(&vtk_capsules, "ref_curv");
  int64_t pd_fit_iters =
      vtk_polydata_add_pointdata_scalar(&vtk_capsules, "nb_fit_iterations");

  double *pd_cap_id_data =
      vtk_polydata_get_pointdata_data(&vtk_capsules, pd_cap_id);
  double *pd_cap_type_data =
      vtk_polydata_get_pointdata_data(&vtk_capsules, pd_cap_type);
  double *pd_owner_pid_data =
      vtk_polydata_get_pointdata_data(&vtk_capsules, pd_owner_pid);
  double *pd_lagvel_data =
      vtk_polydata_get_pointdata_data(&vtk_capsules, pd_lagvel);
  double *pd_lagforce_data =
      vtk_polydata_get_pointdata_data(&vtk_capsules, pd_lagforce);
  double *pd_normal_data =
      vtk_polydata_get_pointdata_data(&vtk_capsules, pd_normal);
  double *pd_curv_data =
      vtk_polydata_get_pointdata_data(&vtk_capsules, pd_curv);
  double *pd_gcurv_data =
      vtk_polydata_get_pointdata_data(&vtk_capsules, pd_gcurv);
  double *pd_ref_curv_data =
      vtk_polydata_get_pointdata_data(&vtk_capsules, pd_ref_curv);
  double *pd_fit_iters_data =
      vtk_polydata_get_pointdata_data(&vtk_capsules, pd_fit_iters);

  size_t point_index = 0;
  for (int c = 0; c < CAPS_COUNT(); c++) {
    CapsuleMesh *mesh = &CAPS(c);
    if (!capsule_vtk_should_write_mesh(mesh))
      continue;

    for (int n = 0; n < mesh->nln; n++, point_index++) {
      CapNode *node = &mesh->nodes[n];
      pd_cap_id_data[point_index] = mesh->cap_id;
      pd_cap_type_data[point_index] = mesh->cap_type;
      pd_owner_pid_data[point_index] = mesh->owner_pid;
      capsule_vtk_store_coord(pd_lagvel_data, point_index, node->lagVel);
      capsule_vtk_store_coord(pd_lagforce_data, point_index, node->lagForce);
      capsule_vtk_store_coord(pd_normal_data, point_index, node->normal);
      pd_curv_data[point_index] = capsule_vtk_clean_double(node->curv);
      pd_gcurv_data[point_index] = capsule_vtk_clean_double(node->gcurv);
      pd_ref_curv_data[point_index] = capsule_vtk_clean_double(node->ref_curv);
      pd_fit_iters_data[point_index] = node->nb_fit_iterations;
    }
  }

  int64_t cd_cap_id = vtk_polydata_add_celldata_scalar(&vtk_capsules, "cap_id");
  int64_t cd_cap_type =
      vtk_polydata_add_celldata_scalar(&vtk_capsules, "cap_type");
  int64_t cd_owner_pid =
      vtk_polydata_add_celldata_scalar(&vtk_capsules, "owner_pid");
  int64_t cd_t1 = vtk_polydata_add_celldata_scalar(&vtk_capsules, "T1");
  int64_t cd_t2 = vtk_polydata_add_celldata_scalar(&vtk_capsules, "T2");
  int64_t cd_tmax = vtk_polydata_add_celldata_scalar(&vtk_capsules, "Tmax");
  int64_t cd_tavg = vtk_polydata_add_celldata_scalar(&vtk_capsules, "Tavg");
  int64_t cd_lambda1 =
      vtk_polydata_add_celldata_scalar(&vtk_capsules, "Lambda_1");
  int64_t cd_lambda2 =
      vtk_polydata_add_celldata_scalar(&vtk_capsules, "Lambda_2");
  int64_t cd_radius = vtk_polydata_add_celldata_scalar(&vtk_capsules, "Radius");
  int64_t cd_area = vtk_polydata_add_celldata_scalar(&vtk_capsules, "area");
  int64_t cd_normal =
      vtk_polydata_add_celldata_vector(&vtk_capsules, "normal", 3);
  int64_t cd_centroid =
      vtk_polydata_add_celldata_vector(&vtk_capsules, "centroid", 3);

  double *cd_cap_id_data =
      vtk_polydata_get_celldata_data(&vtk_capsules, cd_cap_id);
  double *cd_cap_type_data =
      vtk_polydata_get_celldata_data(&vtk_capsules, cd_cap_type);
  double *cd_owner_pid_data =
      vtk_polydata_get_celldata_data(&vtk_capsules, cd_owner_pid);
  double *cd_t1_data = vtk_polydata_get_celldata_data(&vtk_capsules, cd_t1);
  double *cd_t2_data = vtk_polydata_get_celldata_data(&vtk_capsules, cd_t2);
  double *cd_tmax_data = vtk_polydata_get_celldata_data(&vtk_capsules, cd_tmax);
  double *cd_tavg_data = vtk_polydata_get_celldata_data(&vtk_capsules, cd_tavg);
  double *cd_lambda1_data =
      vtk_polydata_get_celldata_data(&vtk_capsules, cd_lambda1);
  double *cd_lambda2_data =
      vtk_polydata_get_celldata_data(&vtk_capsules, cd_lambda2);
  double *cd_radius_data =
      vtk_polydata_get_celldata_data(&vtk_capsules, cd_radius);
  double *cd_area_data = vtk_polydata_get_celldata_data(&vtk_capsules, cd_area);
  double *cd_normal_data =
      vtk_polydata_get_celldata_data(&vtk_capsules, cd_normal);
  double *cd_centroid_data =
      vtk_polydata_get_celldata_data(&vtk_capsules, cd_centroid);

  size_t cell_index = 0;
  for (int c = 0; c < CAPS_COUNT(); c++) {
    CapsuleMesh *mesh = &CAPS(c);
    if (!capsule_vtk_should_write_mesh(mesh))
      continue;

    for (int tr = 0; tr < mesh->nlt; tr++, cell_index++) {
      TriangleState *state = LAG_TRI_STATE(mesh, tr);
      double t1 = capsule_vtk_clean_double(state->tension[0]);
      double t2 = capsule_vtk_clean_double(state->tension[1]);
      cd_cap_id_data[cell_index] = mesh->cap_id;
      cd_cap_type_data[cell_index] = mesh->cap_type;
      cd_owner_pid_data[cell_index] = mesh->owner_pid;
      cd_t1_data[cell_index] = t1;
      cd_t2_data[cell_index] = t2;
      cd_tmax_data[cell_index] = t1 > t2 ? t1 : t2;
      cd_tavg_data[cell_index] = 0.5 * (t1 + t2);
      cd_lambda1_data[cell_index] = capsule_vtk_clean_double(state->stretch[0]);
      cd_lambda2_data[cell_index] = capsule_vtk_clean_double(state->stretch[1]);
      cd_radius_data[cell_index] = capsule_vtk_clean_double(mesh->cap_radius);
      cd_area_data[cell_index] = capsule_vtk_clean_double(state->area);
      capsule_vtk_store_coord(cd_normal_data, cell_index, state->normal);
      capsule_vtk_store_coord(cd_centroid_data, cell_index, state->centroid);
    }
  }

  vtkHDFPolyData vtk_hdf_capsules =
      vtk_HDF_polydata_init_static(fname, overwrite, &vtk_capsules);
  vtk_HDF_polydata_close(&vtk_hdf_capsules);

  vtk_polydata_free(&vtk_capsules);
}

/*! 
 * @brief 
 */
trace void output_hdf_capsules(const char *basename = NULL, int iter = i,
                               double time = t, bool overwrite = true) {
  char basename_buff[1024];

  if (!basename)
    snprintf(basename_buff, sizeof(basename_buff), "%s", get_unique_basename());
  else
    snprintf(basename_buff, sizeof(basename_buff), "%s", basename);

  output_hdf_capsules_series(basename_buff, iter, time, overwrite);
}
