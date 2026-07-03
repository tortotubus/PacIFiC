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
  foreach_dimension() {
    proc_max->x = -HUGE;
    proc_min->x = HUGE;
  }

  Cache c = {0};
  foreach() cache_append( &c, point, 0 );
  foreach_cache(c)
  {
   if(point.level>-1)
   {
    foreach_dimension() {
      double cell_min = x - Delta/2.;
      double cell_max = x + Delta/2.;
      if (proc_max->x < cell_max) proc_max->x = cell_max;
      if (proc_min->x > cell_min) proc_min->x = cell_min;
    }
   }
  }
  free(c.p); //free cache
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

bool is_capsule_in_boundingbox(coord proc_max, coord proc_min, lagMesh* mesh) 
{
  return lagmesh_bounding_sphere_intersects_box(mesh, proc_min, proc_max);
}
