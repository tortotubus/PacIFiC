/*This file contains the miscellaneous plugin functions used for the simulations */


#define rand_pos_periodic() (((double)rand() / (double)((unsigned)RAND_MAX + 1)\
  - 0.5)*L0)
#define rand_pos() (((double)rand() / (double)((unsigned)RAND_MAX + 1)\
  - 0.5)*(L0 - 2*(RADIUS + MIN_INITIAL_GAP)))


/*Judge if two capsules are too close to each other*/
#define MIN_DELTA (L0/(1 << MAXLEVEL))
#define MIN_INITIAL_GAP (3*MIN_DELTA)
#define MAX_POSITIONING_ATTEMPTS 1.e+6


bool no_intersection(lagMesh* caps, int k) {
  bool no_intersection = true;
  for(int i=0; i<k; i++) {
    foreach_dimension() {
      if (fabs(GENERAL_SQNORM(caps[i].centroid, caps[k].centroid)) < sq(2*RADIUS
        + MIN_INITIAL_GAP)) {
        no_intersection = false;
        break;
      }
    }
  }
  return no_intersection;
}


void generate_random_capsules()
{

  if (pid() == 0) {
    for(int k=0; k<NCAPS; k++) {
      bool keep_drawing_positions = true;
      int nb_attempts = 0;
      while (keep_drawing_positions && nb_attempts < MAX_POSITIONING_ATTEMPTS) {
        CAPS(k).centroid.x = rand_pos();
        CAPS(k).centroid.y = rand_pos();
        CAPS(k).centroid.z = rand_pos();
        if (no_intersection(allCaps.caps, k)) {
          for(int i=0; i<CAPS(k).nln; i++)
            foreach_dimension() CAPS(k).nodes[i].pos.x += CAPS(k).centroid.x;
          keep_drawing_positions = false;
          fprintf(stderr, "Number of attempts to insert capsule %d: %d\n", k,
            nb_attempts+1);
        }
        nb_attempts++;
      }
      if (nb_attempts == MAX_POSITIONING_ATTEMPTS) {
        fprintf(stderr, "Error: max number of attempts to insert capsule \
          %d reached.\n", k);
        assert(false); //ggd
      }
    }
  }
  #if _MPI
    /** We now inform all the other processors of the positions of the capsules */
    double centroids[3*NCAPS];
    if (pid() == 0) {
      for(int k=0; k<NCAPS; k++) {
        centroids[3*k] = CAPS(k).centroid.x;
        centroids[3*k+1] = CAPS(k).centroid.y;
        centroids[3*k+2] = CAPS(k).centroid.z;
      }
    }
    MPI_Bcast(centroids, 3*NCAPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (pid() > 0) {
      for(int k=0; k<NCAPS; k++) {
        CAPS(k).centroid.x = centroids[3*k];
        CAPS(k).centroid.y = centroids[3*k+1];
        CAPS(k).centroid.z = centroids[3*k+2];
        for(int i=0; i<CAPS(k).nln; i++)
          foreach_dimension() CAPS(k).nodes[i].pos.x += CAPS(k).centroid.x;
      }
    }
  #endif

  /** We generate stencils for the capsules **/
  for(int k=0; k<NCAPS; k++) correct_lag_pos(&CAPS(k));
    generate_lag_stencils(no_warning = true);
}


void generate_capsules_from_input()
{

   FILE *file = NULL;

 if (pid() == 0) {
    file = fopen("/media/localdata/guodong.gai/BasiliskDLM/Newtest/capsules/caps_shear/insert_position_x.dat", "r");
    for (int k = 0; k < NCAPS; k++)
    {
        if (fscanf(file, "%lf", &CAPS(k).centroid.x) != 1)
        { 
          fprintf(stderr, "Error reading number %d from the file x.\n", k+1);
          fclose(file);
          assert(false); //ggd
        }
    }
    file = NULL;

    file = fopen("/media/localdata/guodong.gai/BasiliskDLM/Newtest/capsules/caps_shear/insert_position_y.dat", "r");
    for (int k = 0; k < NCAPS; k++)
    {
        if (fscanf(file, "%lf", &CAPS(k).centroid.y) != 1)
        { 
          fprintf(stderr, "Error reading number %d from the file y.\n", k+1);
          fclose(file);
          assert(false); //ggd
        }
    }

    file = NULL;
    file = fopen("/media/localdata/guodong.gai/BasiliskDLM/Newtest/capsules/caps_shear/insert_position_z.dat", "r");
    for (int k = 0; k < NCAPS; k++)
    {
        if (fscanf(file, "%lf", &CAPS(k).centroid.z) != 1)
        { 
          fprintf(stderr, "Error reading number %d from the file z.\n", k+1);
          fclose(file);
          assert(false); //ggd
        }
    }
    fclose(file);

    for (int k = 0; k < NCAPS; k++)
    {
      for(int i=0; i<CAPS(k).nln; i++)
            foreach_dimension() CAPS(k).nodes[i].pos.x += CAPS(k).centroid.x;
    }
 }
      

  #if _MPI
    /** We now inform all the other processors of the positions of the capsules */
    double centroids[3*NCAPS];
    if (pid() == 0) {
      for(int k=0; k<NCAPS; k++) {
        centroids[3*k] = CAPS(k).centroid.x;
        centroids[3*k+1] = CAPS(k).centroid.y;
        centroids[3*k+2] = CAPS(k).centroid.z;
      }
    }
    MPI_Bcast(centroids, 3*NCAPS, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (pid() > 0) {
      for(int k=0; k<NCAPS; k++) {
        CAPS(k).centroid.x = centroids[3*k];
        CAPS(k).centroid.y = centroids[3*k+1];
        CAPS(k).centroid.z = centroids[3*k+2];
        for(int i=0; i<CAPS(k).nln; i++)
          foreach_dimension() CAPS(k).nodes[i].pos.x += CAPS(k).centroid.x;
      }
    }
  #endif

  /** We generate stencils for the capsules **/
  for(int k=0; k<NCAPS; k++) correct_lag_pos(&CAPS(k));
    generate_lag_stencils(no_warning = true);
}


void generate_capsules_ordered()
{
    int ncaps1d = cbrt(NCAPS);
    double separation = L0/ncaps1d; // Calculate separation between spheres
    double centroids[3*NCAPS];

    int sphereIndex = 0;
    for (int i = 0; i < ncaps1d; i++) {
        for (int j = 0; j < ncaps1d; j++) {
            for (int k = 0; k < ncaps1d; k++) {
                if (sphereIndex < NCAPS) {
                    centroids[3*sphereIndex] = i * separation + separation / 2. -0.5*L0;
                    centroids[3*sphereIndex+1] = j * separation + separation / 2. -0.5*L0;
                    centroids[3*sphereIndex+2] = k * separation + separation / 2. -0.5*L0;
                    sphereIndex++;
                }
            }
        }
    }

    for(int k=0; k<NCAPS; k++) 
    {
      CAPS(k).centroid.x = centroids[3*k];
      CAPS(k).centroid.y = centroids[3*k+1];
      CAPS(k).centroid.z = centroids[3*k+2];
      for(int i=0; i<CAPS(k).nln; i++)
        foreach_dimension() CAPS(k).nodes[i].pos.x += CAPS(k).centroid.x;
    }

    for (int k = 0; k < NCAPS; k++) {
        printf("ncaps1d %d, separation %lf, Sphere %d: Center(%.2f, %.2f, %.2f) \n", ncaps1d, separation, k, CAPS(k).centroid.x, CAPS(k).centroid.y, CAPS(k).centroid.z);
    }

  /** We generate stencils for the capsules **/
  for(int k=0; k<NCAPS; k++) correct_lag_pos(&CAPS(k));
    generate_lag_stencils(no_warning = true);
}



void compute_vorticity()
{
scalar omega_x[];
scalar omega_y[];
scalar omega_z[];
foreach()
{
    omega_x[] = ((u.z[0,1] - u.z[0,-1] - u.y[0,0,1] + u.y[0,0,-1])/(2.*Delta))/(SHEAR_RATE); //Dimensionless vorticity.x
    omega_y[] = ((u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta))/(SHEAR_RATE); // Dimensionless vorticity.y 
    omega_z[] = ((u.y[1] - u.y[-1] - u.x[0,1] + u.x[0,-1])/(2.*Delta))/(SHEAR_RATE); //Dimensionless vorticity.z
}
boundary({omega_x, omega_y, omega_z});

}

