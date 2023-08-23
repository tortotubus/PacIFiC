/**
# MPI synchronization of the Lagrangian mesh

In the current implementation, each processor has its own copy of the Lagrangian
mesh. In order to update the position of all the local copies of the Lagrangian
mesh, the processors need to communicate to each other the contributions of
their flow field to the velocities of each Lagrangian nodes.
*/

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

void reduce_lagVel2(lagMesh* mesh) {
  int li = dimension; // length of one item
  double* send_data = (double*)malloc(mesh->nln*li*sizeof(double));
  double* recv_data = (double*)malloc(mpi_npe*mesh->nln*li*sizeof(double));
  double* recv_data2 = (double*)malloc(mesh->nln*li*sizeof(double));
      for(int i=0; i<mesh->nln; i++) {
    send_data[li*i] = mesh->nodes[i].lagVel.x;
    send_data[li*i+1] = mesh->nodes[i].lagVel.y;
    #if dimension > 2
    send_data[li*i+2] = mesh->nodes[i].lagVel.z;
    #endif
  }
  MPI_Reduce(send_data, recv_data2, mesh->nln*li, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(recv_data2, mesh->nln*li, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  for(int i=0; i<mesh->nln; i++) {
    mesh->nodes[i].lagVel.x = 0.;
    mesh->nodes[i].lagVel.y = 0.;
    #if dimension > 2
    mesh->nodes[i].lagVel.z = 0.;
    #endif
  }

  for(int i=0; i<mesh->nln; i++) {
    mesh->nodes[i].lagVel.x = recv_data2[li*i];
    mesh->nodes[i].lagVel.y = recv_data2[li*i+1];
    #if dimension > 2
    mesh->nodes[i].lagVel.z = recv_data2[li*i+2];
    #endif
  }
  free(send_data);
  free(recv_data);
  free(recv_data2);
}


bool is_capsule_in_proc(lagMesh* mesh) {

  Cache c = {0};
  double delta = (L0/(1 << grid->maxdepth));

  foreach()
  {
    coord checkpt={0};
    checkpt.x = x;
    checkpt.y = y;
    checkpt.z = z;
    double tetative_dist = sqrt(GENERAL_SQNORM(checkpt, mesh->centroid));
    if(tetative_dist <= mesh->circum_radius + delta*5)
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