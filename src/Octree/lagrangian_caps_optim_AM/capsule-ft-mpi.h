/**
# MPI synchronization of the Lagrangian mesh

In the current implementation, each processor has its own copy of the Lagrangian
mesh. In order to update the position of all the local copies of the Lagrangian
mesh, the processors need to communicate to each other the contributions of
their flow field to the velocities of each Lagrangian nodes.
*/

bool on_face2(double p, int n, double l0) {
  if ((fabs(p/(l0/n)) - ((int)fabs(p/(l0/n)))) < 1.e-10) return true;
  else return false;
}

double correct_point_pos(double ptx, double a, double L) {
   double origin = a + L/2; // L is the length of the domain in the current dimension
    if (on_face2(ptx, N, L))
      ptx += 1.e-10;
    //FIXME: the nodes should not be sent to the other side of the domain if the boundary is not periodic...
    if (ptx > origin + L/2)
      ptx -= L;
    else if (ptx < origin - L/2)
      ptx += L;
  return ptx;
}

/*Two-sides periodicity*/
#define POS_PBC_X(X) ((u.x.boundary[left] != periodic_bc) ? (X) : correct_point_pos(X, X0, L0))
#define POS_PBC_Y(Y) ((u.x.boundary[top] != periodic_bc) ? (Y) : correct_point_pos(Y, Y0, L0*Dimensions.y/Dimensions.x))
#define POS_PBC_Z(Z) ((u.x.boundary[front] != periodic_bc) ? (Z) : correct_point_pos(Z, Z0, L0*Dimensions.z/Dimensions.x))

/*Onside periodicity, might be useful for visualization*/
#define vPOS_PBC_X(X) ((u.x.boundary[left] != periodic_bc) ? (X) : (((X - (X0 +\
  L0/2)) > L0/2.) ? (X) - L0 : (X)))
#define vPOS_PBC_Y(Y) ((u.x.boundary[top] != periodic_bc) ? (Y) : (((Y - (Y0 +\
  L0*Dimensions.y/Dimensions.x/2)) > L0*Dimensions.y/Dimensions.x/2.) ? (Y) - L0*Dimensions.y/Dimensions.x : (Y)))
#define vPOS_PBC_Z(Z) ((u.x.boundary[front] != periodic_bc) ? (Z) : (((Z - (Z0 +\
  L0*Dimensions.z/Dimensions.x/2)) > L0*Dimensions.z/Dimensions.x/2.) ? (Z) - L0*Dimensions.z/Dimensions.x : (Z)))


void reduce_lagVel(lagMesh* mesh) {
  int li = dimension; // length of one item
  double* send_data = (double*)malloc(mesh->nln*li*sizeof(double));
  double* recv_data = (double*)malloc(mpi_npe*mesh->nln*li*sizeof(double));
  for(int i=0; i<mesh->nln; i++) {
    send_data[li*i] = mesh->nodes[i].lagVel.x;
    send_data[li*i+1] = mesh->nodes[i].lagVel.y;
    #if dimension > 2
    send_data[li*i+2] = mesh->nodes[i].lagVel.z;
    #endif
  }
  MPI_Allgather(send_data, mesh->nln*li, MPI_DOUBLE, recv_data,
    mesh->nln*li, MPI_DOUBLE, MPI_COMM_WORLD);

  for(int i=0; i<mesh->nln; i++) {
    mesh->nodes[i].lagVel.x = 0.;
    mesh->nodes[i].lagVel.y = 0.;
    #if dimension > 2
    mesh->nodes[i].lagVel.z = 0.;
    #endif
  }

  for(int k=0; k<mpi_npe; k++) {
    for(int i=0; i<mesh->nln; i++) {
      mesh->nodes[i].lagVel.x += recv_data[k*li*mesh->nln+li*i];
      mesh->nodes[i].lagVel.y += recv_data[k*li*mesh->nln+li*i+1];
      #if dimension > 2
      mesh->nodes[i].lagVel.z += recv_data[k*li*mesh->nln+li*i+2];
      #endif
    }
  }
  free(send_data);
  free(recv_data);
}

void reduce_alllagVel() 
{ 
  int total_nln = 0;
  int* pos_in_pack = (int*)malloc(NCAPS*sizeof(int));
  for(int i=0; i<NCAPS; i++) 
  {
      if (CAPS(i).isactive) 
      {
        // eul2lag(&CAPS(i));
        pos_in_pack[i] = total_nln;
        total_nln += CAPS(i).nln;
      }
  }   

  int li = dimension; // length of one item
  double* send_data_pack = (double*)malloc(total_nln*li*sizeof(double));
  double* recv_data_pack = (double*)malloc(total_nln*li*sizeof(double));

  for(int i=0; i<NCAPS; i++) 
    {
      if (CAPS(i).isactive) 
      { 
        for(int j=0; j<CAPS(i).nln; j++) 
        {
            send_data_pack[pos_in_pack[i]*li + li*j] = CAPS(i).nodes[j].lagVel.x;
            send_data_pack[pos_in_pack[i]*li + li*j+1] = CAPS(i).nodes[j].lagVel.y;
            #if dimension > 2
            send_data_pack[pos_in_pack[i]*li + li*j+2] = CAPS(i).nodes[j].lagVel.z;
            #endif
        }
      }
    }   

  MPI_Allreduce(send_data_pack, recv_data_pack, total_nln*li, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(int i=0; i<NCAPS; i++) 
    {
      if (CAPS(i).isactive) 
      { 
        for(int j=0; j<CAPS(i).nln; j++) 
        {
            CAPS(i).nodes[j].lagVel.x = recv_data_pack[pos_in_pack[i]*li + li*j];
            CAPS(i).nodes[j].lagVel.y = recv_data_pack[pos_in_pack[i]*li + li*j+1];
            #if dimension > 2
            CAPS(i).nodes[j].lagVel.z = recv_data_pack[pos_in_pack[i]*li + li*j+2];
            #endif
        }
      }
    }   

  free(pos_in_pack);
  free(send_data_pack);
  free(recv_data_pack);
}


bool is_capsule_in_proc(lagMesh* mesh) {

  Cache c = {0};
  foreach()
  {
    coord checkpt={0};
    checkpt.x = x;
    checkpt.y = y;
    checkpt.z = z;
    double tetative_dist = sqrt(GENERAL_SQNORM(checkpt, mesh->centroid));
    if(tetative_dist <= mesh->circum_radius)
    {
	    cache_append( &c, point, 0 );
    }
  }

  foreach_cache(c)
  {
   if(point.level>-1)
   {
     free(c.p); //free cache
     return true;
   }
  }
  free(c.p); //free cache
  return false;
}


void compute_proc_borders(coord* proc_max, coord* proc_min)
{
  #if MULT_GRID == 1
    coord domain_length = {
      L0*L0_ratio.x,
      L0*L0_ratio.y,
      L0*L0_ratio.z
    };

    proc_min->x = X0 + mpi_coords[0]*domain_length.x/Dimensions.x;
    proc_max->x = X0 + (mpi_coords[0] + 1)*domain_length.x/Dimensions.x;
    proc_min->y = Y0 + mpi_coords[1]*domain_length.y/Dimensions.y;
    proc_max->y = Y0 + (mpi_coords[1] + 1)*domain_length.y/Dimensions.y;
    proc_min->z = Z0 + mpi_coords[2]*domain_length.z/Dimensions.z;
    proc_max->z = Z0 + (mpi_coords[2] + 1)*domain_length.z/Dimensions.z;
  #else
  double xmin = HUGE, ymin = HUGE, zmin = HUGE;
  double xmax = -HUGE, ymax = -HUGE, zmax = -HUGE;

  foreach(reduction(min:xmin) reduction(min:ymin) reduction(min:zmin)
    reduction(max:xmax) reduction(max:ymax) reduction(max:zmax))
  {
    xmin = min(xmin, x - Delta/2.);
    xmax = max(xmax, x + Delta/2.);
    ymin = min(ymin, y - Delta/2.);
    ymax = max(ymax, y + Delta/2.);
    zmin = min(zmin, z - Delta/2.);
    zmax = max(zmax, z + Delta/2.);
  }

  proc_min->x = xmin;
  proc_min->y = ymin;
  proc_min->z = zmin;
  proc_max->x = xmax;
  proc_max->y = ymax;
  proc_max->z = zmax;
  #endif
}

void gather_all_proc_borders(coord local_min, coord local_max,
  coord* all_proc_min, coord* all_proc_max)
{
  double send_data[6] = {
    local_min.x, local_min.y, local_min.z,
    local_max.x, local_max.y, local_max.z
  };
  double* recv_data = (double*)malloc(6*mpi_npe*sizeof(double));

  MPI_Allgather(send_data, 6, MPI_DOUBLE, recv_data, 6, MPI_DOUBLE,
    MPI_COMM_WORLD);

  for(int p=0; p<mpi_npe; p++) {
    all_proc_min[p].x = recv_data[6*p];
    all_proc_min[p].y = recv_data[6*p + 1];
    all_proc_min[p].z = recv_data[6*p + 2];
    all_proc_max[p].x = recv_data[6*p + 3];
    all_proc_max[p].y = recv_data[6*p + 4];
    all_proc_max[p].z = recv_data[6*p + 5];
  }

  free(recv_data);
}

coord wrap_periodic_point_in_domain(coord point)
{
  coord wrapped = point;
  coord domain_min = {X0, Y0, Z0};
  coord domain_length = {
    L0*L0_ratio.x,
    L0*L0_ratio.y,
    L0*L0_ratio.z
  };

  foreach_dimension()
  {
    if (Period.x) {
      while (wrapped.x < domain_min.x)
        wrapped.x += domain_length.x;
      while (wrapped.x >= domain_min.x + domain_length.x)
        wrapped.x -= domain_length.x;
    }
  }

  return wrapped;
}

coord clamp_nonperiodic_point_to_domain(coord point)
{
  coord clamped = point;
  coord domain_min = {X0, Y0, Z0};
  coord domain_length = {
    L0*L0_ratio.x,
    L0*L0_ratio.y,
    L0*L0_ratio.z
  };

  foreach_dimension()
  {
    if (!Period.x) {
      double lo = domain_min.x;
      double hi = domain_min.x + domain_length.x;
      double eps = 1.e-12*domain_length.x;
      if (clamped.x < lo)
        clamped.x = lo;
      if (clamped.x >= hi)
        clamped.x = hi - eps;
    }
  }

  return clamped;
}

bool point_in_box_half_open(coord point, coord box_min, coord box_max)
{
  foreach_dimension()
  {
    if (point.x < box_min.x || point.x >= box_max.x)
      return false;
  }
  return true;
}

bool point_in_box_closed(coord point, coord box_min, coord box_max)
{
  foreach_dimension()
  {
    if (point.x < box_min.x || point.x > box_max.x)
      return false;
  }
  return true;
}

int find_capsule_owner_proc(lagMesh* mesh, coord* all_proc_min, coord* all_proc_max)
{
  coord owner_point = wrap_periodic_point_in_domain(mesh->centroid);
  owner_point = clamp_nonperiodic_point_to_domain(owner_point);

  for(int p=0; p<mpi_npe; p++)
    if (point_in_box_half_open(owner_point, all_proc_min[p], all_proc_max[p]))
      return p;

  for(int p=0; p<mpi_npe; p++)
    if (point_in_box_closed(owner_point, all_proc_min[p], all_proc_max[p]))
      return p;

  return -1;
}

int estimate_owner_to_ghost_nints(lagMesh* mesh)
{
  (void) mesh;
  int nints = 4; // cap_id, cap_type, nln, nle
  #if dimension > 2
    nints++; // nlt
  #endif
  return nints;
}

int estimate_owner_to_ghost_ndoubles(lagMesh* mesh)
{
  int ndoubles = 3; // cap_es, cap_radius, circum_radius
  ndoubles += dimension; // centroid
  ndoubles += mesh->nln*dimension; // node positions
  ndoubles += mesh->nln*dimension; // node forces
  return ndoubles;
}

void pack_owner_to_ghost_capsule(lagMesh* mesh, int* int_data, int* int_pos,
  double* double_data, int* double_pos)
{
  int_data[(*int_pos)++] = mesh->cap_id;
  int_data[(*int_pos)++] = mesh->cap_type;
  int_data[(*int_pos)++] = mesh->nln;
  int_data[(*int_pos)++] = mesh->nle;
  #if dimension > 2
    int_data[(*int_pos)++] = mesh->nlt;
  #endif

  double_data[(*double_pos)++] = mesh->cap_es;
  double_data[(*double_pos)++] = mesh->cap_radius;
  double_data[(*double_pos)++] = mesh->circum_radius;
  foreach_dimension()
    double_data[(*double_pos)++] = mesh->centroid.x;
  for(int i=0; i<mesh->nln; i++)
    foreach_dimension()
      double_data[(*double_pos)++] = mesh->nodes[i].pos.x;
  for(int i=0; i<mesh->nln; i++)
    foreach_dimension()
      double_data[(*double_pos)++] = mesh->nodes[i].lagForce.x;
}

void unpack_owner_to_ghost_header(int* int_data, int* int_pos,
  double* double_data, int* double_pos, int* cap_id, int* cap_type,
  int* nln, int* nle, int* nlt, double* cap_es, double* cap_radius,
  double* circum_radius)
{
  *cap_id = int_data[(*int_pos)++];
  *cap_type = int_data[(*int_pos)++];
  *nln = int_data[(*int_pos)++];
  *nle = int_data[(*int_pos)++];
  #if dimension > 2
    *nlt = int_data[(*int_pos)++];
  #else
    *nlt = 0;
  #endif

  *cap_es = double_data[(*double_pos)++];
  *cap_radius = double_data[(*double_pos)++];
  *circum_radius = double_data[(*double_pos)++];
}

int estimate_ghost_to_owner_nints(lagMesh* mesh)
{
  (void) mesh;
  return 2; // cap_id, nln
}

int estimate_ghost_to_owner_ndoubles(lagMesh* mesh)
{
  return mesh->nln*dimension; // node velocity contribution
}

void pack_ghost_to_owner_capsule_lagVel(lagMesh* mesh, coord* lagVel,
  int* int_data, int* int_pos, double* double_data, int* double_pos)
{
  int_data[(*int_pos)++] = mesh->cap_id;
  int_data[(*int_pos)++] = mesh->nln;

  for(int i=0; i<mesh->nln; i++)
    foreach_dimension()
      double_data[(*double_pos)++] = lagVel ? lagVel[i].x :
        mesh->nodes[i].lagVel.x;
}

void pack_ghost_to_owner_capsule(lagMesh* mesh, int* int_data, int* int_pos,
  double* double_data, int* double_pos)
{
  pack_ghost_to_owner_capsule_lagVel(mesh, NULL, int_data, int_pos,
    double_data, double_pos);
}

bool sphere_intersects_box(coord center, double radius, coord box_min, coord box_max)
{
  double d2 = 0.;
  foreach_dimension()
  {
    double closest = center.x;
    if (closest < box_min.x)
      closest = box_min.x;
    else if (closest > box_max.x)
      closest = box_max.x;
    d2 += sq(center.x - closest);
  }
  return d2 <= sq(radius);
}

bool periodic_sphere_intersects_box(coord center, double radius, coord box_min, coord box_max)
{
  double xshift[3] = {0., 0., 0.};
  double yshift[3] = {0., 0., 0.};
  double zshift[3] = {0., 0., 0.};
  int nx = 1, ny = 1, nz = 1;

  coord domain_min = {X0, Y0, Z0};
  coord domain_max = {
    X0 + L0*L0_ratio.x,
    Y0 + L0*L0_ratio.y,
    Z0 + L0*L0_ratio.z
  };
  coord domain_length = {
    L0*L0_ratio.x,
    L0*L0_ratio.y,
    L0*L0_ratio.z
  };

  if (Period.x) {
    if (center.x - radius < domain_min.x) xshift[nx++] = domain_length.x;
    if (center.x + radius > domain_max.x) xshift[nx++] = -domain_length.x;
  }
  if (Period.y) {
    if (center.y - radius < domain_min.y) yshift[ny++] = domain_length.y;
    if (center.y + radius > domain_max.y) yshift[ny++] = -domain_length.y;
  }
  if (Period.z) {
    if (center.z - radius < domain_min.z) zshift[nz++] = domain_length.z;
    if (center.z + radius > domain_max.z) zshift[nz++] = -domain_length.z;
  }

  for(int ix=0; ix<nx; ix++)
    for(int iy=0; iy<ny; iy++)
      for(int iz=0; iz<nz; iz++) {
        coord image_center = {
          center.x + xshift[ix],
          center.y + yshift[iy],
          center.z + zshift[iz]
        };
        if (sphere_intersects_box(image_center, radius, box_min, box_max))
          return true;
      }

  return false;
}

bool lagmesh_bounding_sphere_intersects_box(lagMesh* mesh, coord box_min, coord box_max)
{
  return periodic_sphere_intersects_box(mesh->centroid, mesh->circum_radius,
    box_min, box_max);
}

double periodic_point_distance2(coord a, coord b)
{
  double d2 = 0.;
  coord domain_length = {
    L0*L0_ratio.x,
    L0*L0_ratio.y,
    L0*L0_ratio.z
  };

  foreach_dimension()
  {
    double d = a.x - b.x;
    if (Period.x) {
      if (d > 0.5*domain_length.x)
        d -= domain_length.x;
      else if (d < -0.5*domain_length.x)
        d += domain_length.x;
    }
    d2 += sq(d);
  }

  return d2;
}

bool bounding_spheres_intersect(coord a_center, double a_radius,
  coord b_center, double b_radius)
{
  return periodic_point_distance2(a_center, b_center) <=
    sq(a_radius + b_radius);
}

bool lagmesh_bounding_spheres_intersect(lagMesh* a, lagMesh* b)
{
  return bounding_spheres_intersect(a->centroid, a->circum_radius,
    b->centroid, b->circum_radius);
}

bool is_capsule_in_boundingbox(coord proc_max, coord proc_min, lagMesh* mesh) 
{
  return lagmesh_bounding_sphere_intersects_box(mesh, proc_min, proc_max);
}
