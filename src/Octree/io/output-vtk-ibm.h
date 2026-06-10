#include "io/output-common.h"

#include <ctype.h>
#include <string.h>
#include <stdlib.h>

#include "io/vtk/vtkHDFPolyData.h"
#include "io/vtk/vtkPolyData.h"
#include "io/vtk/vtkTimeSeries.h"

#include "ibm/IBMeshManager.h"

static inline bool output_hdf_ibm_should_write_node (IBMesh* mesh,
                                                     IBNode* node) {
  if (!mesh || !node)
    return false;
#if _MPI
  /* Sparse-managed meshes may contain local-support copies whose centers are
     owned by another rank.  Those copies are required for coupling, but they
     are not authoritative particle/IB-node records and should not appear in
     the default polydata output. */
  return node->pid == pid ();
#else
  return true;
#endif
}

#if _MPI
static inline bool output_hdf_ibm_has_sparse_managed_meshes (void) {
  foreach_ibmesh () {
    if (mesh->sparse_managed)
      return true;
  }
  return false;
}
#endif

trace void output_hdf_ibm_series (IBscalar* scalar_list = NULL,
                                 IBvector* vector_list = NULL,
                                 const char* basename,
                                 int iter,
                                 double time,
                                 bool overwrite) {

  char pname[128];
  sprintf (pname, "%s/vtkhdf-pd-series", basename);

  char fname[128];
  sprintf (fname, "%s/vtkhdf-pd-series/polydata_%d.vtkhdf", basename, iter);

  char series_entry[128];
  sprintf (series_entry, "vtkhdf-pd-series/polydata_%d.vtkhdf", iter);

  char series_filename[128];
  sprintf (series_filename, "%s/pd.vtkhdf.series", basename);

  if (pid () == 0) {
    assert (!create_path (basename));
    assert (!create_path (pname));
    output_vtk_series (series_entry, series_filename, iter, time);
  }

IBscalar* slist = scalar_list;
IBvector* vlist = vector_list;

if (!slist && !vlist) {
  foreach_ibscalar (iball) {
    if (_ibattribute[s.i].v.x.i == s.i)
      vlist = ibvectors_add (vlist, _ibattribute[s.i].v);
    else if (_ibattribute[s.i].v.x.i < 0)
      slist = iblist_add (slist, s);
  }
}


  // Init polydata

#if _MPI
  if (output_hdf_ibm_has_sparse_managed_meshes ())
    ibmeshmanager_sync_model_outputs ();
  else
    ibmeshmanager_update_pid ();
  size_t n_points_local = 0;

  foreach_ibnode_per_ibmesh () {
    if (output_hdf_ibm_should_write_node (mesh, node)) {
      n_points_local++;
    }
  }
#else
  size_t n_points_local = 0;
  foreach_ibnode_per_ibmesh () {
    n_points_local++;
    node->pid = 0.;
  }
#endif

  size_t n_points = n_points_local;
  size_t n_vertices = n_points;
  size_t n_lines = 0;
  size_t n_strips = 0;
  size_t n_polygons = 0;

  int n_scalars = iblist_len(slist);
  int n_vectors = ibvectors_len(vlist);
  const size_t n_ibm_metadata = 6;
  size_t n_pointdata = (size_t) n_scalars + (size_t) n_vectors + n_ibm_metadata;

  vtkPolyData vtk_ibm = vtk_polydata_init (
    n_points, n_vertices, n_lines, n_strips, n_polygons, n_pointdata);

  // Populate polydata
  foreach_ibnode_per_ibmesh () {
    if (output_hdf_ibm_should_write_node (mesh, node)) {
      vtk_polydata_add_point (&vtk_ibm, pos.x, pos.y, pos.z);
    }
  }

  for (size_t i = 0; i < n_points; i++) {
    vtk_polydata_add_vertex (&vtk_ibm, i);
  }

  size_t ncomp = dimension;

  double** scalar_data = NULL;
  double** vector_data = NULL;
  int64_t* scalar_ids = NULL;
  int64_t* vector_ids = NULL;

  scalar_ids = malloc (sizeof (int64_t) * n_scalars);
  scalar_data = malloc (sizeof (double*) * n_scalars);
  vector_ids = malloc (sizeof (int64_t) * n_vectors);
  vector_data = malloc (sizeof (double*) * n_vectors);

  {
    int i;

    int64_t rank_pid_id =
      vtk_polydata_add_pointdata_scalar (&vtk_ibm, "rank_pid");
    int64_t node_pid_id =
      vtk_polydata_add_pointdata_scalar (&vtk_ibm, "node_pid");
    int64_t pid_id =
      vtk_polydata_add_pointdata_scalar (&vtk_ibm, "pid");
    int64_t mesh_pid_id =
      vtk_polydata_add_pointdata_scalar (&vtk_ibm, "mesh_pid");
    int64_t mesh_gid_id =
      vtk_polydata_add_pointdata_scalar (&vtk_ibm, "mesh_gid");
    int64_t node_lid_id =
      vtk_polydata_add_pointdata_scalar (&vtk_ibm, "node_lid");

    double* rank_pid_data =
      vtk_polydata_get_pointdata_data (&vtk_ibm, rank_pid_id);
    double* node_pid_data =
      vtk_polydata_get_pointdata_data (&vtk_ibm, node_pid_id);
    double* pid_data =
      vtk_polydata_get_pointdata_data (&vtk_ibm, pid_id);
    double* mesh_pid_data =
      vtk_polydata_get_pointdata_data (&vtk_ibm, mesh_pid_id);
    double* mesh_gid_data =
      vtk_polydata_get_pointdata_data (&vtk_ibm, mesh_gid_id);
    double* node_lid_data =
      vtk_polydata_get_pointdata_data (&vtk_ibm, node_lid_id);

    i = 0;
    foreach_ibscalar (slist) {
      scalar_ids[i] = vtk_polydata_add_pointdata_scalar (&vtk_ibm, ibname (s));
      scalar_data[i] = vtk_polydata_get_pointdata_data (&vtk_ibm, scalar_ids[i]);
      i++;
    }

    i = 0;
    foreach_ibvector (vlist) {
      char* vector_name;
      size_t trunc_len = (size_t) (strlen (ibname (v.x)) - 2);
      vector_name = malloc ((trunc_len + 1) * sizeof (char));
      strncpy (vector_name, ibname (v.x), trunc_len);
      vector_name[trunc_len] = '\0';
      vector_ids[i] = vtk_polydata_add_pointdata_vector (&vtk_ibm, vector_name, dimension);
      vector_data[i] = vtk_polydata_get_pointdata_data (&vtk_ibm, vector_ids[i]);
      free (vector_name);
      i++;
    }

    i = 0;
    foreach_ibnode_per_ibmesh () {
      if (output_hdf_ibm_should_write_node (mesh, node)) {
        rank_pid_data[i] = pid ();
        node_pid_data[i] = node->pid;
        pid_data[i] = node->pid;
        mesh_pid_data[i] = mesh->pid;
        mesh_gid_data[i] = mesh->gid;
        node_lid_data[i] = node->node_lid;

        // scalar ibfields
        size_t si = 0;
        foreach_ibscalar (slist) {
          scalar_data[si][i] = ibval (s);
          si++;
        }

        // vector ibfields
        size_t vi = 0;
        foreach_ibvector (vlist) {
          vector_data[vi][i * dimension + 0] = ibval (v.x);
#if dimension >= 2
          vector_data[vi][i * dimension + 1] = ibval (v.y);
#endif
#if dimension >= 3
          vector_data[vi][i * dimension + 2] = ibval (v.z);
#endif
          vi++;
        }
 

        i++;
      }
    }
  }

  // Write vtkhdf polydata
  vtkHDFPolyData vtk_hdf_ibm =
    vtk_HDF_polydata_init_static (fname, true, &vtk_ibm);
  vtk_HDF_polydata_close (&vtk_hdf_ibm);

  // Free vtkPolyData
  vtk_polydata_free (&vtk_ibm);
}

trace void output_hdf_ibm (IBscalar* scalar_list = NULL,
                          IBvector* vector_list = NULL,
                          const char* basename = NULL,
                          int iter = i,
                          double time = t,
                          bool use_transient = false,
                          bool overwrite = true,
                          bool use_ab = false) {

  char basename_buff[64];

  if (!basename) {
    snprintf (
      basename_buff, sizeof (basename_buff), "%s", get_unique_basename ());
  } else {
    snprintf (basename_buff, sizeof (basename_buff), "%s", basename);
  }

  output_hdf_ibm_series (scalar_list, vector_list, basename_buff, iter, time, overwrite);
}
