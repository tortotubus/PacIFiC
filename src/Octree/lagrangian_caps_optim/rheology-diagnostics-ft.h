#pragma once

#include <errno.h>
#include <stdbool.h>
#include <sys/stat.h>

typedef struct CapsuleRheologyStress {
  double n1;
  double n2;
  double shear;
  double pressure;
} CapsuleRheologyStress;

typedef struct CapsuleRheologyShape {
  double count;
  double area_ratio;
  double volume_ratio;
  double td_maxmin;
  double max_radius_ratio;
  double min_radius_ratio;
  double ang_vel_z;
} CapsuleRheologyShape;

static double capsule_rheology_clean(double value) {
  return isfinite(value) ? value : 0.;
}

static void capsule_rheology_ensure_dir(const char *path) {
  if (!path || !path[0])
    return;
  if (mkdir(path, 0777) != 0 && errno != EEXIST) {
    perror("capsule_rheology_ensure_dir(): mkdir");
    assert(false);
  }
}

static bool capsule_rheology_file_needs_header(const char *path) {
  FILE *fp = fopen(path, "r");
  if (!fp)
    return true;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return true;
  }
  long size = ftell(fp);
  fclose(fp);
  return size <= 0;
}

static double capsule_rheology_fluid_wall_shear(double outside_viscosity) {
  double top_visc_stress = 0.;
  double bottom_visc_stress = 0.;

  foreach_boundary(top, reduction(+:top_visc_stress)) {
    top_visc_stress += (u.x[0, 1] - u.x[]) * Delta +
                       (u.y[1] - u.y[-1]) * 0.5 * Delta;
  }
  foreach_boundary(bottom, reduction(+:bottom_visc_stress)) {
    bottom_visc_stress += (u.x[0, 1] - u.x[]) * Delta +
                          (u.y[1] - u.y[-1]) * 0.5 * Delta;
  }

  top_visc_stress *= outside_viscosity / sq(L0);
  bottom_visc_stress *= outside_viscosity / sq(L0);
  return 0.5 * (top_visc_stress + bottom_visc_stress);
}

static CapsuleRheologyStress capsule_rheology_compute_stress(
    CapsuleMesh *mesh, double outside_viscosity, double inside_viscosity,
    bool include_viscosity_jump) {
  CapsuleRheologyStress out = {0};
  if (!mesh || !mesh->isactive || !capsule_mesh_is_local_owner(mesh))
    return out;

  double sigmaxx = 0.;
  double sigmayy = 0.;
  double sigmazz = 0.;
  double sigmaxy = 0.;
  double visc_ratio = outside_viscosity > 0.
                          ? inside_viscosity / outside_viscosity
                          : 1.;

  for (int i = 0; i < mesh->nln; i++) {
    coord r = {0};
    foreach_dimension()
      r.x = mesh->nodes[i].pos.x - mesh->centroid.x;

    sigmaxx += -mesh->nodes[i].lagForce.x * r.x;
    sigmayy += -mesh->nodes[i].lagForce.y * r.y;
    sigmazz += -mesh->nodes[i].lagForce.z * r.z;
    sigmaxy += -0.5 * (mesh->nodes[i].lagForce.x * r.y +
                       mesh->nodes[i].lagForce.y * r.x);

    if (include_viscosity_jump) {
      double nodal_area = compute_node_area(mesh, i);
      if (nodal_area <= SEPS)
        continue;

      sigmaxx += 2. * outside_viscosity * (visc_ratio - 1.) *
                 mesh->nodes[i].lagVel.x * mesh->nodes[i].normal.x *
                 nodal_area;
      sigmayy += 2. * outside_viscosity * (visc_ratio - 1.) *
                 mesh->nodes[i].lagVel.y * mesh->nodes[i].normal.y *
                 nodal_area;
      sigmazz += 2. * outside_viscosity * (visc_ratio - 1.) *
                 mesh->nodes[i].lagVel.z * mesh->nodes[i].normal.z *
                 nodal_area;
      sigmaxy += outside_viscosity * (visc_ratio - 1.) *
                 (mesh->nodes[i].lagVel.x * mesh->nodes[i].normal.y +
                  mesh->nodes[i].lagVel.y * mesh->nodes[i].normal.x) *
                 nodal_area;
    }
  }

  out.n1 = capsule_rheology_clean(sigmaxx - sigmayy);
  out.n2 = capsule_rheology_clean(sigmayy - sigmazz);
  out.shear = capsule_rheology_clean(sigmaxy);
  out.pressure = capsule_rheology_clean(-(sigmaxx + sigmayy + sigmazz) / 3.);
  return out;
}

static CapsuleRheologyShape
capsule_rheology_compute_shape(CapsuleMesh *mesh) {
  CapsuleRheologyShape out = {0};
  if (!mesh || !mesh->isactive || !capsule_mesh_is_local_owner(mesh))
    return out;

  double area = 0.;
  for (int i = 0; i < mesh->nlt; i++)
    area += LAG_TRI_STATE(mesh, i)->area;

  double rmax = 0.;
  double rmin = HUGE;
  for (int i = 0; i < mesh->nln; i++) {
    coord r = {0};
    foreach_dimension()
      r.x = mesh->nodes[i].pos.x - mesh->centroid.x;
    double radius = cnorm(r);
    if (radius > rmax)
      rmax = radius;
    if (radius < rmin)
      rmin = radius;
  }

  out.count = 1.;
  out.area_ratio =
      mesh->cap_radius > 0. ? area / (4. * pi * sq(mesh->cap_radius)) : 0.;
  out.volume_ratio =
      fabs(mesh->initial_volume) > SEPS ? mesh->volume / mesh->initial_volume
                                        : 0.;
  out.td_maxmin = (rmax + rmin) > SEPS ? (rmax - rmin) / (rmax + rmin) : 0.;
  out.max_radius_ratio = mesh->cap_radius > 0. ? rmax / mesh->cap_radius : 0.;
  out.min_radius_ratio = mesh->cap_radius > 0. ? rmin / mesh->cap_radius : 0.;
  out.ang_vel_z = mesh->ang_vel.z;
  return out;
}

static void output_capsule_rheology(const char *basename, int iter, double time,
                                    double outside_viscosity,
                                    double inside_viscosity,
                                    double shear_rate,
                                    bool include_viscosity_jump) {
  int n = CAPS_COUNT();
  double *send_stress = (double *)calloc((size_t)4 * n, sizeof(double));
  double *recv_stress = (double *)calloc((size_t)4 * n, sizeof(double));
  double *send_shape = (double *)calloc((size_t)7 * n, sizeof(double));
  double *recv_shape = (double *)calloc((size_t)7 * n, sizeof(double));
  assert(send_stress && recv_stress && send_shape && recv_shape);

  for (int c = 0; c < n; c++) {
    CapsuleMesh *mesh = &CAPS(c);
    CapsuleRheologyStress stress = capsule_rheology_compute_stress(
        mesh, outside_viscosity, inside_viscosity, include_viscosity_jump);
    CapsuleRheologyShape shape = capsule_rheology_compute_shape(mesh);

    send_stress[4 * c + 0] = stress.n1;
    send_stress[4 * c + 1] = stress.n2;
    send_stress[4 * c + 2] = stress.shear;
    send_stress[4 * c + 3] = stress.pressure;

    send_shape[7 * c + 0] = shape.count;
    send_shape[7 * c + 1] = shape.area_ratio;
    send_shape[7 * c + 2] = shape.volume_ratio;
    send_shape[7 * c + 3] = shape.td_maxmin;
    send_shape[7 * c + 4] = shape.max_radius_ratio;
    send_shape[7 * c + 5] = shape.min_radius_ratio;
    send_shape[7 * c + 6] = shape.ang_vel_z;
  }

#if _MPI
  MPI_Reduce(send_stress, recv_stress, 4 * n, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(send_shape, recv_shape, 7 * n, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
#else
  memcpy(recv_stress, send_stress, (size_t)4 * n * sizeof(double));
  memcpy(recv_shape, send_shape, (size_t)7 * n * sizeof(double));
#endif

  double fluid_wall_shear = capsule_rheology_fluid_wall_shear(outside_viscosity);

  if (pid() == 0) {
    double n1 = 0.;
    double n2 = 0.;
    double particle_shear = 0.;
    double particle_pressure = 0.;
    CapsuleRheologyShape avg = {0};

    for (int c = 0; c < n; c++) {
      n1 += recv_stress[4 * c + 0];
      n2 += recv_stress[4 * c + 1];
      particle_shear += recv_stress[4 * c + 2];
      particle_pressure += recv_stress[4 * c + 3];

      avg.count += recv_shape[7 * c + 0];
      avg.area_ratio += recv_shape[7 * c + 1];
      avg.volume_ratio += recv_shape[7 * c + 2];
      avg.td_maxmin += recv_shape[7 * c + 3];
      avg.max_radius_ratio += recv_shape[7 * c + 4];
      avg.min_radius_ratio += recv_shape[7 * c + 5];
      avg.ang_vel_z += recv_shape[7 * c + 6];
    }

    if (avg.count > 0.) {
      avg.area_ratio /= avg.count;
      avg.volume_ratio /= avg.count;
      avg.td_maxmin /= avg.count;
      avg.max_radius_ratio /= avg.count;
      avg.min_radius_ratio /= avg.count;
      avg.ang_vel_z /= avg.count;
    }

    double volume = cube(L0);
    double stress_scale = outside_viscosity * fabs(shear_rate);
    double normalized_particle_shear =
        stress_scale > 0. ? particle_shear / (stress_scale * volume) : 0.;
    double normalized_bulk_viscosity =
        stress_scale > 0.
            ? fluid_wall_shear / stress_scale + normalized_particle_shear
            : 0.;

    const char *base = (basename && basename[0]) ? basename : ".";
    capsule_rheology_ensure_dir(base);

    char path[4096];
    snprintf(path, sizeof(path), "%s/rheology.txt", base);
    bool header = capsule_rheology_file_needs_header(path);
    FILE *fp = fopen(path, "a");
    assert(fp);
    if (header) {
      fprintf(fp,
              "# iter t fluid_wall_shear particle_pressure particle_shear N1 "
              "N2 normalized_particle_shear normalized_bulk_viscosity "
              "avg_area_ratio avg_volume_ratio avg_td_maxmin "
              "avg_max_radius_ratio avg_min_radius_ratio avg_ang_vel_z "
              "owned_caps\n");
    }
    fprintf(fp,
            "%d %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g "
            "%.17g %.17g %.17g %.17g %.17g %.17g\n",
            iter, time, capsule_rheology_clean(fluid_wall_shear),
            capsule_rheology_clean(particle_pressure),
            capsule_rheology_clean(particle_shear), capsule_rheology_clean(n1),
            capsule_rheology_clean(n2),
            capsule_rheology_clean(normalized_particle_shear),
            capsule_rheology_clean(normalized_bulk_viscosity),
            capsule_rheology_clean(avg.area_ratio),
            capsule_rheology_clean(avg.volume_ratio),
            capsule_rheology_clean(avg.td_maxmin),
            capsule_rheology_clean(avg.max_radius_ratio),
            capsule_rheology_clean(avg.min_radius_ratio),
            capsule_rheology_clean(avg.ang_vel_z),
            capsule_rheology_clean(avg.count));
    fflush(fp);
    fclose(fp);
  }

  free(send_stress);
  free(recv_stress);
  free(send_shape);
  free(recv_shape);
}
