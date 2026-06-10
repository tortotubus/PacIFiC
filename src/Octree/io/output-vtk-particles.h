#include "io/output-common.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "io/vtk/vtkHDFPolyData.h"
#include "io/vtk/vtkPolyData.h"
#include "io/vtk/vtkTimeSeries.h"

#include "particle/ParticleGrid.h"

// ============================================================================
// Function Declarations
// ============================================================================

void output_hdf_particles_series(Pscalar *scalar_list, Pvector *vector_list,
                          const char *basename, int iter, double time,
                          bool overwrite);

void output_hdf_particles(Pscalar *scalar_list, Pvector *vector_list,
                   const char *basename, int iter, double time,
                   bool use_transient, bool overwrite, bool use_ab);

// ============================================================================
// Function Definitions
// ============================================================================

/*!
 * @brief
 */
trace void output_hdf_particles_series(Pscalar *scalar_list = NULL,
                                Pvector *vector_list = NULL,
                                const char *basename, int iter, double time,
                                bool overwrite) {

  char pname[128];
  sprintf(pname, "%s/vtkhdf-particle-series", basename);

  char fname[128];
  sprintf(fname, "%s/vtkhdf-particle-series/particles_%d.vtkhdf", basename, iter);

  char series_entry[128];
  sprintf(series_entry, "vtkhdf-particle-series/particles_%d.vtkhdf", iter);

  char series_filename[128];
  sprintf(series_filename, "%s/particles.vtkhdf.series", basename);

  if (pid() == 0) {
    assert(!create_path(basename));
    assert(!create_path(pname));
    output_vtk_series(series_entry, series_filename, iter, time);
  }
#if _MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  Pscalar *slist = scalar_list;
  Pvector *vlist = vector_list;

  if (!slist && !vlist) {
    foreach_pscalar(pall) {
      if (_pattribute[s.i].v.x.i == s.i)
        vlist = pvectors_add(vlist, _pattribute[s.i].v);
      else if (_pattribute[s.i].v.x.i < 0)
        slist = plist_add(slist, s);
    }
  }

  // Init particles

  size_t n_points_local = 0;
  foreach_particle() {
#if _MPI
    if (particle->pid != pid())
      continue;
#else
    particle->pid = 0.;
#endif
    n_points_local++;
  }

  size_t n_points = n_points_local;
  size_t n_vertices = n_points;
  size_t n_lines = 0;
  size_t n_strips = 0;
  size_t n_polygons = 0;

  int n_scalars = plist_len(slist);
  int n_vectors = pvectors_len(vlist);
  size_t n_pointdata = n_scalars + n_vectors;

  vtkPolyData vtk_particles = vtk_polydata_init(n_points, n_vertices, n_lines,
                                         n_strips, n_polygons, n_pointdata);

  // Populate particles
  foreach_particle() {
    if (particle->pid == pid()) {
      vtk_polydata_add_point(&vtk_particles, pos.x, pos.y, pos.z);
    }
  }

  for (size_t i = 0; i < n_points; i++) {
    vtk_polydata_add_vertex(&vtk_particles, i);
  }

  size_t ncomp = dimension;

  double **scalar_data = NULL;
  double **vector_data = NULL;
  int64_t *scalar_ids = NULL;
  int64_t *vector_ids = NULL;

  scalar_ids = malloc(sizeof(int64_t) * n_scalars);
  scalar_data = malloc(sizeof(double *) * n_scalars);
  vector_ids = malloc(sizeof(int64_t) * n_vectors);
  vector_data = malloc(sizeof(double *) * n_vectors);

  {
    int i;

    i = 0;
    foreach_pscalar(slist) {
      scalar_ids[i] = vtk_polydata_add_pointdata_scalar(&vtk_particles, pname(s));
      scalar_data[i] = vtk_polydata_get_pointdata_data(&vtk_particles, scalar_ids[i]);
      i++;
    }

    i = 0;
    foreach_pvector(vlist) {
      char *vector_name;
      size_t trunc_len = (size_t)(strlen(pname(v.x)) - 2);
      vector_name = malloc((trunc_len + 1) * sizeof(char));
      strncpy(vector_name, pname(v.x), trunc_len);
      vector_name[trunc_len] = '\0';
      vector_ids[i] =
          vtk_polydata_add_pointdata_vector(&vtk_particles, vector_name, dimension);
      vector_data[i] = vtk_polydata_get_pointdata_data(&vtk_particles, vector_ids[i]);
      free(vector_name);
      i++;
    }

    i = 0;
    foreach_particle() {
      if (particle->pid == pid()) {

        // scalar pfields
        size_t si = 0;
        foreach_pscalar(slist) {
          scalar_data[si][i] = pval(s);
          si++;
        }

        // vector pfields
        size_t vi = 0;
        foreach_pvector(vlist) {
          vector_data[vi][i * dimension + 0] = pval(v.x);
#if dimension >= 2
          vector_data[vi][i * dimension + 1] = pval(v.y);
#endif
#if dimension >= 3
          vector_data[vi][i * dimension + 2] = pval(v.z);
#endif
          vi++;
        }

        i++;
      }
    }
  }

  // Write vtkhdf particles
  vtkHDFPolyData vtk_hdf_particles =
      vtk_HDF_polydata_init_static(fname, true, &vtk_particles);
  vtk_HDF_polydata_close(&vtk_hdf_particles);

  // Free vtkPolyData
  vtk_polydata_free(&vtk_particles);
}

/*!
 * @brief
 */
trace void output_hdf_particles(Pscalar *scalar_list = NULL,
                         Pvector *vector_list = NULL,
                         const char *basename = NULL, int iter = i,
                         double time = t, bool use_transient = false,
                         bool overwrite = true, bool use_ab = false) {

  char basename_buff[64];

  if (!basename) {
    snprintf(basename_buff, sizeof(basename_buff), "%s", get_unique_basename());
  } else {
    snprintf(basename_buff, sizeof(basename_buff), "%s", basename);
  }

  output_hdf_particles_series(scalar_list, vector_list, basename_buff, iter, time,
                       overwrite);
}
