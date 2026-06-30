/**
# Toolbox to perform operations on a triangulated meshes

From defining geometric computations such as normal vectors, volume and
centroid, useful macros, subdividing triangles, below is a collection of helpful
functions to deal with triangular meshes.
*/

/**
## Geometric computations

The function below computes the length of an edge. It takes as arguments
a pointer to the mesh as well as the ID of the edge of interest.
*/
double edge_length(lagMesh* mesh, int i) {
  double length = 0.;
  int v1, v2;
  v1 = LAG_EDGE_NODE_ID(mesh, i, 0);
  v2 = LAG_EDGE_NODE_ID(mesh, i, 1);
  foreach_dimension() {
    length += sq(GENERAL_1DIST(mesh->nodes[v1].pos.x, mesh->nodes[v2].pos.x, L0*L0_ratio.x));
  }
  return sqrt(length);
}

/**
The function ```compute_lengths``` below computes the lengths of all edges. It
takes as an argument a pointer to the mesh. If the optional argument
```force``` is set to ```true```, the edges' lengths are computed no matter the value of ```updated_stretches```.
*/
struct compute_lengths_type{
  lagMesh* mesh;
  bool force;
};

typedef struct compute_lengths_type _compute_lengths;

void compute_lengths(_compute_lengths p) {
  lagMesh* mesh = p.mesh;
  bool force = (p.force) ? p.force : false;
  if (force || !mesh->updated_stretches) {
    for(int i=0; i < mesh->nle; i++)
      mesh->edges[i].length = edge_length(mesh, i);
    mesh->updated_stretches = true;
  }
}

#if dimension < 3
/**
The two functions below compute the outward normal vector to all the edges of
a Lagrangian mesh, for 2D simulations.
*/
void comp_edge_normal(lagMesh* mesh, int i) {
  int node_id[2];
  for(int j=0; j<2; j++) node_id[j] = LAG_EDGE_NODE_ID(mesh, i, j);
  mesh->edges[i].normal.y = GENERAL_1DIST(mesh->nodes[node_id[0]].pos.x,
    mesh->nodes[node_id[1]].pos.x, L0*L0_ratio.x );
  mesh->edges[i].normal.x = GENERAL_1DIST(mesh->nodes[node_id[1]].pos.y,
    mesh->nodes[node_id[0]].pos.y, L0*L0_ratio.y );
  double normn = sqrt(sq(mesh->edges[i].normal.x)
    + sq(mesh->edges[i].normal.y));
  foreach_dimension() mesh->edges[i].normal.x /= normn;
}

void comp_edge_normals(lagMesh* mesh) {
  for(int i=0; i<mesh->nle; i++) comp_edge_normal(mesh, i);
}
#else // dimension > 2
/** In 3D simulations, the function below assumes that the Lagrangian mesh contains the origin and
is convex, and swaps the order of the nodes in order to compute an outward
normal vector. This only need to be performed at the creation of the mesh since
the outward property of the normal vectors won't change through the simulation.
*/
void comp_initial_area_normals(lagMesh* mesh) {
  for(int i=0; i<mesh->nlt; i++) {
    int nid[3]; // node ids
    coord centroid; /** Note: the centroid is only valid if the triangle is not
    across periodic boundaries, which is fine for this function since it is
    assumed the center of the membrane is at the origin. */
    foreach_dimension() centroid.x = 0.;
    for(int j=0; j<3; j++) {
      nid[j] = LAG_TRIANGLE_NODE_ID(mesh, i, j);
      foreach_dimension() centroid.x += mesh->nodes[nid[j]].pos.x/3;
    }
    coord normal, e[2];
    for(int j=0; j<2; j++)
      foreach_dimension()
        e[j].x = GENERAL_1DIST(mesh->nodes[nid[0]].pos.x,
          mesh->nodes[nid[j+1]].pos.x, L0*L0_ratio.x);
    foreach_dimension() normal.x = e[0].y*e[1].z - e[0].z*e[1].y;
    double norm = sqrt(sq(normal.x) + sq(normal.y) + sq(normal.z));
    double dp = 0.; // dp for "dot product"
    foreach_dimension() dp += normal.x*centroid.x;
    /** If the dot product is negative, the computed normal is inward and we
    need to swap two nodes of the triangle.*/
    if (dp < 0) {
      SET_LAG_TRIANGLE_NODE_ID(mesh, i, 1, nid[2]);
      SET_LAG_TRIANGLE_NODE_ID(mesh, i, 2, nid[1]);
      foreach_dimension() normal.x *= -1;
    }
    foreach_dimension() {
      mesh->triangles[i].centroid.x = centroid.x;
      mesh->triangles[i].normal.x = normal.x/norm;
    }
    mesh->triangles[i].area = norm/2;
  }
}

/** The two functions below compute the outward normal vector to all the
triangles of the mesh, for 3D simulations. */
void comp_triangle_area_normal(lagMesh* mesh, int i) {
  int nid[3]; // node ids
  for(int j=0; j<3; j++) nid[j] = LAG_TRIANGLE_NODE_ID(mesh, i, j);
  /** The next 15 lines compute the centroid of the triangle, making sure it is
  valid when the triangle lies across periodic boundaries. */
  foreach_dimension() mesh->triangles[i].centroid.x = 0.;

  // for(int j=0; j<3; j++) {
  //   foreach_dimension() {
  //     mesh->triangles[i].centroid.x +=
  //       ACROSS_PERIODIC(mesh->nodes[nid[j]].pos.x/3,
  //       mesh->nodes[nid[0]].pos.x/3) ? mesh->nodes[nid[j]].pos.x/3 - L0 :
  //       mesh->nodes[nid[j]].pos.x/3;
  //   }
  // }

///////ggd 
  for(int j=0; j<3; j++) {
    foreach_dimension() {
      mesh->triangles[i].centroid.x += mesh->centroid.x + GENERAL_1DIST(mesh->nodes[nid[j]].pos.x, mesh->centroid.x, L0*L0_ratio.x);
    }
  }
foreach_dimension() mesh->triangles[i].centroid.x /=3.;
////ggd 


  coord origin = {X0 + L0*L0_ratio.x/2, Y0 + L0*L0_ratio.y/2, Z0 + L0*L0_ratio.z/2}; //FIXM

  foreach_dimension() {
    if (fabs(mesh->triangles[i].centroid.x - origin.x) > L0*L0_ratio.x/2.) {
      if (mesh->triangles[i].centroid.x - origin.x > 0)
        mesh->triangles[i].centroid.x -= L0*L0_ratio.x;
      else mesh->triangles[i].centroid.x += L0*L0_ratio.x;
    }
  }
  coord normal, e[2];
  for(int j=0; j<2; j++)
    foreach_dimension()
      e[j].x = GENERAL_1DIST(mesh->nodes[nid[0]].pos.x,
        mesh->nodes[nid[j+1]].pos.x, L0*L0_ratio.x);
  foreach_dimension() normal.x = e[0].y*e[1].z - e[0].z*e[1].y;
  double norm = sqrt(sq(normal.x) + sq(normal.y) + sq(normal.z));
  foreach_dimension() mesh->triangles[i].normal.x = normal.x/norm;
  mesh->triangles[i].area = norm/2.;
}

void comp_triangle_area_normals(lagMesh* mesh) {
  for(int i=0; i<mesh->nlt; i++) comp_triangle_area_normal(mesh, i);
}
#endif

/**
If a Lagrangian node falls exactly on an edge or a vertex of the Eulerian
mesh, some issues arise when checking for periodic boundary conditions. As a
quick fix, if this is the case we shift the point position by $10^{-10}$, as
is done in the two functions below.
*/
bool on_face(double p, int n, double l0) {
  if ((fabs(p/(l0/n)) - ((int)fabs(p/(l0/n)))) < 1.e-10) return true;
  else return false;
}

void correct_node_pos(coord* node) {  
  coord origin = {X0 + L0*L0_ratio.x/2, Y0 + L0*L0_ratio.y/2, Z0 + L0*L0_ratio.z/2};
  
  foreach_dimension() {
    if (on_face(node->x, N, L0*L0_ratio.x))
      node->x += 1.e-10;
    //FIXME: the nodes should not be sent to the other side of the domain if the boundary is not periodic...
    if (node->x > origin.x + L0*L0_ratio.x/2)
      node->x -= L0*L0_ratio.x;
    else if (node->x < origin.x - L0*L0_ratio.x/2)
      node->x += L0*L0_ratio.x;
  }
}

void correct_lag_pos(lagMesh* mesh) {
  for(int i=0; i < mesh->nln; i++) {
    correct_node_pos(&mesh->nodes[i].pos);
  }
  mesh->updated_stretches = false;
  mesh->updated_normals = false;
  mesh->updated_curvatures = false;
}

/**
The function below computes the centroid of the capsule as the average of the
coordinates of all its nodes. The centroid is stored as an attribute of the
caps structure.
*/
void comp_centroid(lagMesh* mesh) {
  coord origin = {X0 + L0/2, Y0 + L0*L0_ratio.y/2, Z0 + L0*L0_ratio.z/2};
  foreach_dimension() mesh->centroid.x = 0.;
  for(int i=0; i<mesh->nln; i++)
    foreach_dimension() {
      double tentative_pos = mesh->nodes[i].pos.x - mesh->nodes[0].pos.x;
      mesh->centroid.x += (tentative_pos < origin.x - L0*L0_ratio.x/2) ?
        tentative_pos + L0*L0_ratio.x :
          ((tentative_pos > origin.x + L0*L0_ratio.x/2) ? tentative_pos - L0*L0_ratio.x :
          tentative_pos);
    }
  foreach_dimension()
    mesh->centroid.x = mesh->centroid.x/mesh->nln + mesh->nodes[0].pos.x;
  correct_node_pos(&mesh->centroid);
}

/**
The function below computes the radius of a circumscribed sphere around the
capsule, centered at the current capsule centroid. The optional padding keeps
the radius conservative for Eulerian stencil operations.
*/

void comp_circum_radius(lagMesh* mesh, double radius_padding) {
  double max_radius = 0.;
  for(int i=0; i<mesh->nln; i++) {
    double tentative_radius = sqrt(GENERAL_SQNORM(mesh->nodes[i].pos, mesh->centroid));
    max_radius = (tentative_radius > max_radius) ? tentative_radius : max_radius;
  }
  mesh->circum_radius = max_radius + radius_padding;
}

trace
void comp_volume(lagMesh* mesh) {
  coord origin = {X0 + L0*L0_ratio.x/2, Y0 + L0*L0_ratio.y/2, Z0 + L0*L0_ratio.z/2}; //FIXME
  comp_centroid(mesh);
  double volume = 0;
  for(int i=0; i<mesh->nlt; i++) {
    coord nodes[3];
    for(int j=0; j<3; j++)
      foreach_dimension() {
        double tentative_pos = mesh->nodes[LAG_TRIANGLE_NODE_ID(mesh, i, j)].pos.x
          - mesh->centroid.x;
        nodes[j].x = (tentative_pos < origin.x - L0*L0_ratio.x/2) ?
          tentative_pos + L0*L0_ratio.x :
          ((tentative_pos > origin.x + L0*L0_ratio.x/2) ? tentative_pos - L0*L0_ratio.x :
          tentative_pos);
      }
    for(int j=0; j<3; j++) {
      coord cross_product;
      foreach_dimension()
        cross_product.x = nodes[(j+1)%3].y*nodes[(j+2)%3].z -
          nodes[(j+1)%3].z*nodes[(j+2)%3].y;
      volume += cdot(nodes[j],cross_product);
    }
  }
  mesh->volume = volume/18;
}


void compute_taylor_factor(lagMesh* mesh, double* Taylor_deform, double* Inclin_angle, 
coord* rs, double* TDmaxmin, double* TDang)
{
    /* To store the components of the diagonal inertia tensor Ixx, Iyy, Izz*/
  double Ixx = 0., Iyy = 0., Izz = 0., Ixy = 0., Ixz = 0., Iyz = 0.; 
  for(int i = 0; i < mesh->nlt; i++) 
  {
    double rn = 0.;
    double rSquared = sq(GENERAL_1DIST(mesh->triangles[i].centroid.x, mesh->centroid.x, L0*L0_ratio.x)) 
    + sq(GENERAL_1DIST(mesh->triangles[i].centroid.y, mesh->centroid.y, L0*L0_ratio.y))
    + sq(GENERAL_1DIST(mesh->triangles[i].centroid.z, mesh->centroid.z, L0*L0_ratio.z)) ;     
    coord tri_vec = {0};

    foreach_dimension() tri_vec.x = GENERAL_1DIST(mesh->triangles[i].centroid.x, mesh->centroid.x, L0*L0_ratio.x);
    foreach_dimension() rn += tri_vec.x * mesh->triangles[i].normal.x;
    // rn = fabs(rn); // In case of centroid ouside of capsule, we change the sign of rn
    
    //Compute the components of inertia tensor
    Ixx += mesh->triangles[i].area / 5.0 * rn * (rSquared - tri_vec.x * tri_vec.x);
    Iyy += mesh->triangles[i].area / 5.0 * rn * (rSquared - tri_vec.y * tri_vec.y);
    Izz += mesh->triangles[i].area / 5.0 * rn * (rSquared - tri_vec.z * tri_vec.z);
    Ixy += mesh->triangles[i].area / 5.0 * rn * (0. - tri_vec.x * tri_vec.y);
    Ixz += mesh->triangles[i].area / 5.0 * rn * (0. - tri_vec.x * tri_vec.z);
    Iyz += mesh->triangles[i].area / 5.0 * rn * (0. - tri_vec.y * tri_vec.z);
 }



  double aa = 0., bb = 0., cc = 0., dd = 0.;
  aa = 1.;
  bb = -(Ixx + Iyy + Izz);
  cc = (Ixx*Iyy + Ixx*Izz + Iyy*Izz - Ixy*Ixy - Ixz*Ixz - Iyz*Iyz);
  dd = (Ixx*Iyz*Iyz + Iyy*Ixz*Ixz + Izz*Ixy*Ixy - 2*Ixy*Ixz*Iyz - Ixx*Iyy*Izz);

  // solve_cubic(aa, bb, cc, dd);
    // Normalize the coefficients
    double A = bb / aa;
    double B = cc / aa;
    double C = dd / aa;

    // Compute the discriminant
    double Q = (A*A - 3*B) / 9;
    double R = (2*A*A*A - 9*A*B + 27*C) / 54;
    double D = Q*Q*Q - R*R;

    if((D < 0) && (t > dt)) fprintf(stderr, "No roots found for inertia equivalent ellipsoid!\n");
    // Normal Case: Three real roots
    double theta = acos(R / sqrt(Q*Q*Q));
    double r = -2 * sqrt(Q);
    double lambd1 = r * cos(theta / 3) - A / 3;
    double lambd2 = r * cos((theta - 2*M_PI) / 3) - A / 3;
    double lambd3 = r * cos((theta - 4*M_PI) / 3) - A / 3;


  double lambd_tmp =  0.;
  if(lambd1 > lambd2)
  {
    lambd_tmp = lambd1;
    lambd1 = lambd2;
    lambd2 = lambd_tmp;
  }
  if(lambd2 > lambd3)
  {
    lambd_tmp = lambd2;
    lambd2 = lambd3;
    lambd3 = lambd_tmp;
  }
  if(lambd1 > lambd2)
  {
    lambd_tmp = lambd1;
    lambd1 = lambd2;
    lambd2 = lambd_tmp;
  }

  double Evx = 1.;
  double Evy = -(Ixx - lambd1) / Ixy;
  // double Evz = -(Ixy * Evx + (Iyy - lambd1) * Evy) / Iyz;

  double r1 = sqrt((5.0 / (2.0 * mesh->volume)) * (lambd2 - lambd1 + lambd3));
  double r2 = sqrt((5.0 / (2.0 * mesh->volume)) * (lambd3 - lambd2 + lambd1));
  double r3 = sqrt((5.0 / (2.0 * mesh->volume)) * (lambd1 - lambd3 + lambd2));

  *Taylor_deform = (r1 - r3)/(r1 + r3);
  *Inclin_angle = atan2(Evy, Evx);

  rs->x = r1;
  rs->y = r2;
  rs->z = r3;

//////////////////////////////////////////////////////////Taylor Maxmin
    double rmax = -HUGE;
    double rmin = HUGE;
    
    for(int i = 0; i < mesh->nln; i++) 
    {
    /** The post-processing is only carried out if we are in the shear plane */
        double projx, projy, projz;
        projx = GENERAL_1DIST(mesh->nodes[i].pos.x, mesh->centroid.x, L0*L0_ratio.x);
        projy = GENERAL_1DIST(mesh->nodes[i].pos.y, mesh->centroid.y, L0*L0_ratio.y);
        projz = GENERAL_1DIST(mesh->nodes[i].pos.z, mesh->centroid.z, L0*L0_ratio.z);
        double rad  = sqrt(sq(projx) + sq(projy) + sq(projz));
        if (rad > rmax) 
        {
          rmax = rad;
          *TDang = (fabs(projx) < 1.e-14) ? (projy>0. ? pi/2. : 3*pi/2.) : fmod(atan2(projy, projx) + pi, pi);
        }
        if (rad < rmin)
         rmin = rad;
    }

    *TDmaxmin = (rmax - rmin)/(rmax + rmin);
}


trace
void comp_capsule_geodynamics(lagMesh* mesh) {
  
  /*Compute the cell size in the grid*/
  #if MULT_GRID == 1   
    double delta = (L0/(1 << grid->maxdepth)/Dimensions.x);
  #else
    double delta = (L0/(1 << grid->maxdepth));
  #endif



  coord center_vel = {0., 0., 0.};
  for(int i=0; i<mesh->nln; i++) 
  {
    foreach_dimension()
      center_vel.x += mesh->nodes[i].lagVel.x / mesh->nln;
  }


  comp_centroid(mesh);
  comp_circum_radius(mesh, 3*delta);
  coord angvel={0.,0.,0.};


  for(int i=0; i<mesh->nln; i++) 
  {
    double tentative_radius = sqrt(GENERAL_SQNORM(mesh->nodes[i].pos, mesh->centroid));

    foreach_dimension()
      angvel.x += (GENERAL_1DIST(mesh->nodes[i].pos.y, mesh->centroid.y, L0*L0_ratio.y)*(mesh->nodes[i].lagVel.z - center_vel.z) -
        GENERAL_1DIST(mesh->nodes[i].pos.z, mesh->centroid.z, L0*L0_ratio.z)*(mesh->nodes[i].lagVel.y - center_vel.y)) /tentative_radius /tentative_radius;
 
  }
  
  /*Compute the angular velocity of the capsule*/
  foreach_dimension()
    mesh->ang_vel.x = angvel.x / mesh->nln;

  // /* To store the components of the diagonal inertia tensor Ixx, Iyy, Izz*/
  // double Ixx = 0., Iyy = 0., Izz = 0., Ixy = 0., Ixz = 0., Iyz = 0.; 
  // for(int i=0; i<mesh->nlt; i++) 
  // {
  //   double rn = 0.;
  //   double rSquared = GENERAL_SQNORM(mesh->triangles[i].centroid, mesh->centroid);     
  //   coord tri_vec = {0};

  //   foreach_dimension() tri_vec.x = GENERAL_1DIST(mesh->triangles[i].centroid.x, mesh->centroid.x);
  //   foreach_dimension() rn += tri_vec.x * mesh->triangles[i].normal.x;
  //   rn = fabs(rn); // In case of centroid ouside of capsule, we change the sign of rn
    
  //   //Compute the components of inertia tensor
  //   Ixx += mesh->triangles[i].area / 5.0 * rn * (rSquared - tri_vec.x * tri_vec.x);
  //   Iyy += mesh->triangles[i].area / 5.0 * rn * (rSquared - tri_vec.y * tri_vec.y);
  //   Izz += mesh->triangles[i].area / 5.0 * rn * (rSquared - tri_vec.z * tri_vec.z);
  //   Ixy += mesh->triangles[i].area / 5.0 * rn * (0. - tri_vec.x * tri_vec.y);
  //   Ixz += mesh->triangles[i].area / 5.0 * rn * (0. - tri_vec.x * tri_vec.z);
  //   Iyz += mesh->triangles[i].area / 5.0 * rn * (0. - tri_vec.y * tri_vec.z);
  // }

  // //printf("Ixx: %lf, Iyy: %lf, Izz: %lf Ixy: %lf, Ixz: %lf, Iyz: %lf\n", Ixx, Iyy, Izz, Ixy, Ixz, Iyz);

  // double aa = 0., bb = 0., cc = 0., dd = 0.;
  // aa = 1.;
  // bb = -(Ixx + Iyy + Izz);
  // cc = (Ixx*Iyy + Ixx*Izz + Iyy*Izz - Ixy*Ixy - Ixz*Ixz - Iyz*Iyz);
  // dd = (Ixx*Iyz*Iyz + Iyy*Ixz*Ixz + Izz*Ixy*Ixy - 2*Ixy*Ixz*Iyz - Ixx*Iyy*Izz);

  // // solve_cubic(aa, bb, cc, dd);
  //   // Normalize the coefficients
  //   double A = bb / aa;
  //   double B = cc / aa;
  //   double C = dd / aa;

  //   // Compute the discriminant
  //   double Q = (A*A - 3*B) / 9;
  //   double R = (2*A*A*A - 9*A*B + 27*C) / 54;
  //   double D = Q*Q*Q - R*R;

  //   if(D < 0) fprintf(stderr, "No roots found for inertia equivalent ellipsoid!\n");
  //   // Normal Case: Three real roots
  //   double theta = acos(R / sqrt(Q*Q*Q));
  //   double r = -2 * sqrt(Q);
  //   double lambd1 = r * cos(theta / 3) - A / 3;
  //   double lambd2 = r * cos((theta - 2*M_PI) / 3) - A / 3;
  //   double lambd3 = r * cos((theta - 4*M_PI) / 3) - A / 3;


  // double lambd_tmp =  0.;
  // if(lambd1 > lambd2)
  // {
  //   lambd_tmp = lambd1;
  //   lambd1 = lambd2;
  //   lambd2 = lambd_tmp;
  // }
  // if(lambd2 > lambd3)
  // {
  //   lambd_tmp = lambd2;
  //   lambd2 = lambd3;
  //   lambd3 = lambd_tmp;
  // }
  // if(lambd1 > lambd2)
  // {
  //   lambd_tmp = lambd1;
  //   lambd1 = lambd2;
  //   lambd2 = lambd_tmp;
  // }

  // double Evx = 1.;
  // double Evy = -(Ixx - lambd1) / Ixy;
  // double Evz = -(Ixy * Evx + (Iyy - lambd1) * Evy) / Iyz;

  // double r1 = sqrt((5.0 / (2.0 * mesh->volume)) * (lambd2 - lambd1 + lambd3));
  // double r2 = sqrt((5.0 / (2.0 * mesh->volume)) * (lambd3 - lambd2 + lambd1));
  // double r3 = sqrt((5.0 / (2.0 * mesh->volume)) * (lambd1 - lambd3 + lambd2));

  // mesh->taylor_deform = (r1 - r3)/(r1 + r3);

}




/** The function below updates the normal vectors on all the nodes as well as
the lengths and midpoints of all the edges (in 2D) or the area and centroids of
all the triangles (in 3D). */
void comp_normals(lagMesh* mesh) {
  if (!mesh->updated_normals) {
    #if dimension < 3
    compute_lengths((compute_lengths){.mesh=mesh});
    for(int i=0; i<mesh->nln; i++) {
      coord n[2];
      double l[2];
      double normn;
      for(int j=0; j<2; j++) {
        int edge_id;
        edge_id = LAG_NODE_EDGE_ID(mesh, i, j);
        l[j] = mesh->edges[edge_id].length;
        comp_edge_normal(mesh, edge_id);
        foreach_dimension() n[j].x = mesh->edges[edge_id].normal.x;
      }
      /** the normal vector at a node is the weighted average of the normal
      vectors of its edges. The average is weighted by the distance of the node
      to each of the edges' centers. */
      double epsilon = l[1]/(l[0] + l[1]);
      normn = 0.;
      foreach_dimension() {
        mesh->nodes[i].normal.x = epsilon*n[0].x + (1. - epsilon)*n[1].x;
        normn += sq(mesh->nodes[i].normal.x);
      }
      normn = sqrt(normn);
      foreach_dimension() mesh->nodes[i].normal.x /= normn;
    }
    #else // dimension == 3
    comp_triangle_area_normals(mesh);
    for(int i=0; i<mesh->nln; i++) {
      foreach_dimension() mesh->nodes[i].normal.x = 0.;
      double sw = 0.; // sum of the weights
      for(int j=0; j<LAG_NODE_NB_TRIANGLES(mesh, i); j++) {
        int tid = LAG_NODE_TRIANGLE_ID(mesh, i, j);
        double dist = 0.;
        foreach_dimension()
          dist += sq(mesh->nodes[i].pos.x - mesh->triangles[tid].centroid.x);
        dist = sqrt(dist);
        sw += dist;
        foreach_dimension()
          mesh->nodes[i].normal.x += mesh->triangles[tid].normal.x*dist;
      }
      foreach_dimension()
        mesh->nodes[i].normal.x /= sw;
      double normn = cnorm(mesh->nodes[i].normal);
      foreach_dimension() mesh->nodes[i].normal.x /= normn;
    }
    #endif
    mesh->updated_normals = true;
  }
}

#if dimension > 2

/**
## Useful macros

The macros below are useful to define an icosahedron
*/
#define GET_LD(NODE) ((fabs(fabs(NODE.pos.x) - ll) < 1.e-8) ? 0 : \
  ((fabs(fabs(NODE.pos.y) - ll) < 1.e-8 ? 1 : 2)))
#define GET_LD_SIGN(NODE) ((GET_LD(NODE) == 0) ? sign(NODE.pos.x) : \
  ((GET_LD(NODE) == 1) ? sign(NODE.pos.y) : sign(NODE.pos.z)))
#define GET_SD(NODE) ((fabs(fabs(NODE.pos.x) - sl) < 1.e-8) ? 0 : \
  ((fabs(fabs(NODE.pos.y) - sl) < 1.e-8 ? 1 : 2)))
#define GET_SD_SIGN(NODE) ((GET_SD(NODE) == 0) ? sign(NODE.pos.x) : \
  ((GET_SD(NODE) == 1) ? sign(NODE.pos.y) : sign(NODE.pos.z)))
#define GET_ZD(NODE) ((fabs(NODE.pos.x) < 1.e-8) ? 0 : \
  ((fabs(NODE.pos.y) < 1.e-8 ? 1 : 2)))

/**
## Operations on edges

The function below returns true if the two nodes $i$ and $j$ are neighbors
*/
bool is_neighbor_ft(lagMesh* mesh, int i, int j) {
  for(int k=0; k<LAG_NODE_NB_NEIGHBORS(mesh, i); k++) {
    if (LAG_NODE_NEIGHBOR_ID(mesh, i, k) == j) return true;
  }
  return false;
}

/** The function below returns true if there is an edge connecting nodes i and
j */
bool edge_exists(lagMesh* mesh, int j, int k) {
  for(int i=0; i<mesh->nle; i++) {
    if ((LAG_EDGE_NODE_ID(mesh, i, 0) == j && LAG_EDGE_NODE_ID(mesh, i, 1) == k)
      || (LAG_EDGE_NODE_ID(mesh, i, 0) == k
      && LAG_EDGE_NODE_ID(mesh, i, 1) == j)) return true;
  }
  return false;
}

/** The function below returns true if the edge $i$ is across a priodic
boundary. */
bool is_edge_across_periodic(lagMesh* mesh, int i) {
  int n[2];
  for(int k=0; k<2; k++) n[k] = LAG_EDGE_NODE_ID(mesh, i, k);
  if (ACROSS_PERIODIC(mesh->nodes[n[0]].pos.x, mesh->nodes[n[1]].pos.x, L0*L0_ratio.x)
    || ACROSS_PERIODIC(mesh->nodes[n[0]].pos.y, mesh->nodes[n[1]].pos.y, L0*L0_ratio.y)
    || ACROSS_PERIODIC(mesh->nodes[n[0]].pos.z, mesh->nodes[n[1]].pos.z, L0*L0_ratio.z))
      return true;
  return false;
}

// foreach_dimension()
// bool is_edge_across_periodic_x(lagMesh* mesh, int i) {
//   int n[2];
//   for(int k=0; k<2; k++) n[k] = mesh->edges[i].node_ids[k];
//     if (ACROSS_PERIODIC(mesh->nodes[n[0]].pos.x, mesh->nodes[n[1]].pos.x))
//       return true;
//   return false;
// }


struct write_edge_type {
  lagMesh* mesh;
  int i;
  int j;
  int k;
  bool new_mesh;
  bool overwrite;
};

typedef struct write_edge_type _write_edge;

/** The function below writes edge i, connecting nodes j and k. If the edge
exists, the function returns false (no edge creation), true otherwise (edge
creation) */
bool write_edge(_write_edge p) {
  lagMesh* mesh = p.mesh;
  int i = p.i;
  int j = p.j;
  int k = p.k;
  bool new_mesh = (p.new_mesh) ? p.new_mesh : false;
  bool overwrite = (p.overwrite) ? p.overwrite : false;
  if (!overwrite && edge_exists(mesh, j, k)) return false;
  else {
    SET_LAG_EDGE_NODE_ID(mesh, i, 0, j);
    SET_LAG_EDGE_NODE_ID(mesh, i, 1, k);
    for(int ii=0; ii<2; ii++) SET_LAG_EDGE_TRIANGLE_ID(mesh, i, ii, -1);
    if (new_mesh) {
      int nbj = LAG_NODE_NB_NEIGHBORS(mesh, j);
      int nbk = LAG_NODE_NB_NEIGHBORS(mesh, k);
      SET_LAG_NODE_NEIGHBOR_ID(mesh, j, nbj, k);
      SET_LAG_NODE_NEIGHBOR_ID(mesh, k, nbk, j);
      SET_LAG_NODE_EDGE_ID(mesh, j, nbj, i);
      SET_LAG_NODE_EDGE_ID(mesh, k, nbk, i);
      SET_LAG_NODE_NB_NEIGHBORS(mesh, j, nbj + 1);
      SET_LAG_NODE_NB_NEIGHBORS(mesh, k, nbk + 1);
    }
    return true;
  }
}


/** The function below creates a new edge between nodes i and j, and updates the
connectivity information of its nodes (but not its triangles, since they
don't exist yet). */
void new_edge(lagMesh* mesh, int i, int j) {
  int eid = mesh->nle; // id of the new edge
  int nodes[2];
  nodes[0] = i; nodes[1] = j;
  for(int k=0; k<2; k++) {
    SET_LAG_EDGE_NODE_ID(mesh, eid, k, nodes[k]);

    /** Add the edge id to the newly connected nodes */
    for(int l=0; l<LAG_NODE_NB_NEIGHBORS(mesh, nodes[k]); l++) {
      if (LAG_NODE_EDGE_ID(mesh, nodes[k], l) == -1) {
        SET_LAG_NODE_EDGE_ID(mesh, nodes[k], l, eid);
        break;
      }
    }

    /** Update the neighbors' list of the newly connected nodes */
    for(int l=0; l<LAG_NODE_NB_NEIGHBORS(mesh, nodes[k]); l++) {
      if (LAG_NODE_NEIGHBOR_ID(mesh, nodes[k], l) == -1) {
        SET_LAG_NODE_NEIGHBOR_ID(mesh, nodes[k], l, nodes[(k+1)%2]);
        break;
      }
    }

    /** The newly created edge is not yet surrounded by any triangle */
    SET_LAG_EDGE_TRIANGLE_ID(mesh, eid, k, -1);
  }
  mesh->nle++;
}

/** The function below splits an edge in two smaller edges, creating a node
at its midpoint. */
void split_edge(lagMesh* mesh, int i) {
  int nid[2];
  for(int j=0; j<2; j++) nid[j] = LAG_EDGE_NODE_ID(mesh, i, j);

  /** Create new node */
  foreach_dimension()
    mesh->nodes[mesh->nln].pos.x =
      .5*(mesh->nodes[nid[0]].pos.x + mesh->nodes[nid[1]].pos.x);
  SET_LAG_NODE_NB_NEIGHBORS(mesh, mesh->nln, 6);
  SET_LAG_NODE_NB_TRIANGLES(mesh, mesh->nln, 6);
  SET_LAG_NODE_NEIGHBOR_ID(mesh, mesh->nln, 0, nid[0]);
  SET_LAG_NODE_NEIGHBOR_ID(mesh, mesh->nln, 1, nid[1]);
  SET_LAG_NODE_EDGE_ID(mesh, mesh->nln, 0, i);
  SET_LAG_NODE_EDGE_ID(mesh, mesh->nln, 1, mesh->nle);
  for(int j=0; j<6; j++) {
    SET_LAG_NODE_TRIANGLE_ID(mesh, mesh->nln, j, -1);
    if (j>1) {
      SET_LAG_NODE_NEIGHBOR_ID(mesh, mesh->nln, j, -1);
      SET_LAG_NODE_EDGE_ID(mesh, mesh->nln, j, -1);
    }
  }

  /** Create new edge and update current one */
  write_edge((_write_edge){.mesh=mesh, .i = i, .j = nid[0], .k = mesh->nln, .overwrite = true});
  write_edge((_write_edge){.mesh=mesh, .i = mesh->nle, .j = nid[1], .k = mesh->nln});
  for (int j=0; j<2; j++) {
    SET_LAG_EDGE_TRIANGLE_ID(mesh, i, j, -1);
    SET_LAG_EDGE_TRIANGLE_ID(mesh, mesh->nle, j, -1);
  }

  /** Update node information: neighboring nodes, connecting edges */
  for(int j=0; j<LAG_NODE_NB_NEIGHBORS(mesh, nid[0]); j++)
    if (LAG_NODE_NEIGHBOR_ID(mesh, nid[0], j) == nid[1])
      SET_LAG_NODE_NEIGHBOR_ID(mesh, nid[0], j, mesh->nln);
  for(int j=0; j<LAG_NODE_NB_NEIGHBORS(mesh, nid[1]); j++)
    if (LAG_NODE_NEIGHBOR_ID(mesh, nid[1], j) == nid[0])
      SET_LAG_NODE_NEIGHBOR_ID(mesh, nid[1], j, mesh->nln);
  for(int j=0; j<LAG_NODE_NB_NEIGHBORS(mesh, nid[1]); j++)
    if (LAG_NODE_EDGE_ID(mesh, nid[1], j) == i)
      SET_LAG_NODE_EDGE_ID(mesh, nid[1], j, mesh->nle);

  mesh->nln++;
  mesh->nle++;
}

/**
## Operations on triangles

The function below returns true if the triangle connecting nodes i,j and k
already exists in the mesh. */
bool triangle_exists(lagMesh* mesh, int i, int j, int k) {
  for(int t=0; t<mesh->nlt; t++) {
    for(int a=0; a<3; a++) {
      if (LAG_TRIANGLE_NODE_ID(mesh, t, a) == i) {
        for(int b=0; b<3; b++) {
          if (b != a && LAG_TRIANGLE_NODE_ID(mesh, t, b) == j) {
            for(int c=0; c<3; c++) {
              if (c != a && c != b && LAG_TRIANGLE_NODE_ID(mesh, t, c) == k) {
                return true;
              }
            }
          }
        }
      }
    }
  }
  return false;
}

/** The function below returns true if the triangle_ids $i$ is across a periodic
boundary. */
bool is_triangle_across_periodic(lagMesh* mesh, int i) {
  int e[3];
  for(int k=0; k<3; k++) e[k] = LAG_TRIANGLE_EDGE_ID(mesh, i, k);
  if (is_edge_across_periodic(mesh, e[0])
    || is_edge_across_periodic(mesh, e[1])
    || is_edge_across_periodic(mesh, e[2]))
      return true;
  return false;
}

struct write_triangle_type {
  lagMesh* mesh;
  int tid;
  int i;
  int j;
  int k;
  bool overwrite;
};

typedef struct write_triangle_type _write_triangle;

/** The function below writes a triangle at index location tid, connecting nodes
i, j and k. It updates the connectivity information of its nodes and edges.*/
bool write_triangle(_write_triangle p) {
  lagMesh* mesh = p.mesh;
  int tid = p.tid;
  int i = p.i;
  int j = p.j;
  int k = p.k;
  bool overwrite = (p.overwrite) ? p.overwrite : false;
  if (!overwrite && triangle_exists(mesh, i, j, k)) return false;
  else {
    SET_LAG_TRIANGLE_NODE_ID(mesh, tid, 0, i);
    SET_LAG_TRIANGLE_NODE_ID(mesh, tid, 1, j);
    SET_LAG_TRIANGLE_NODE_ID(mesh, tid, 2, k);
    int c = 0;
    for(int a=0; a<3; a++) {
      int va = LAG_TRIANGLE_NODE_ID(mesh, tid, a);
      int b=(a+1)%3;
      int vb = LAG_TRIANGLE_NODE_ID(mesh, tid, b);
      for(int m=0; m<LAG_NODE_NB_NEIGHBORS(mesh, va); m++) {
        for(int n=0; n<LAG_NODE_NB_NEIGHBORS(mesh, vb); n++) {
          if (LAG_NODE_EDGE_ID(mesh, va, m) == LAG_NODE_EDGE_ID(mesh, vb, n)) {
            int eid = LAG_NODE_EDGE_ID(mesh, va, m);
            SET_LAG_TRIANGLE_EDGE_ID(mesh, tid, c, eid);
            c++;
            int p = (LAG_EDGE_TRIANGLE_ID(mesh, eid, 0)
              == -1) ? 0 : 1;
            SET_LAG_EDGE_TRIANGLE_ID(mesh, eid, p, tid);
          }
        }
      }
    }
    int nbt = LAG_NODE_NB_TRIANGLES(mesh, i);
    SET_LAG_NODE_TRIANGLE_ID(mesh, i, nbt, tid);
    SET_LAG_NODE_NB_TRIANGLES(mesh, i, nbt + 1);
    nbt = LAG_NODE_NB_TRIANGLES(mesh, j);
    SET_LAG_NODE_TRIANGLE_ID(mesh, j, nbt, tid);
    SET_LAG_NODE_NB_TRIANGLES(mesh, j, nbt + 1);
    nbt = LAG_NODE_NB_TRIANGLES(mesh, k);
    SET_LAG_NODE_TRIANGLE_ID(mesh, k, nbt, tid);
    SET_LAG_NODE_NB_TRIANGLES(mesh, k, nbt + 1);
    return true;
  }
}

struct _new_triangle {
  lagMesh* mesh;
  int i;
  int j;
  int k;
  int prev_tid;
};

/** The function below also writes a new triangle connecting nodes i, j and k,
but in the specific context of the refinement process. The variable $prev\_tid$
is used to update the connectivity information of the nodes and edges. */
void new_triangle(lagMesh* mesh, int i, int j, int k, int prev_tid) {
  assert((i < mesh->nln) && (j < mesh->nln) && (k < mesh->nln));
  int nodes[3];
  nodes[0] = i; nodes[1] = j; nodes[2] = k;

  int tid = mesh->nlt;
  for(int k=0; k<3; k++) {
    SET_LAG_TRIANGLE_NODE_ID(mesh, tid, k, nodes[k]);
    /** Specify the id of the new triangle for node k */
    bool replaced_tid = false;
    for(int l=0; l<LAG_NODE_NB_TRIANGLES(mesh, nodes[k]); l++)
      if (LAG_NODE_TRIANGLE_ID(mesh, nodes[k], l) == prev_tid) {
        SET_LAG_NODE_TRIANGLE_ID(mesh, nodes[k], l, tid);
        replaced_tid = true;
        break;
      }
    if (!replaced_tid) {
      for(int l=0; l<LAG_NODE_NB_TRIANGLES(mesh, nodes[k]); l++)
        if (LAG_NODE_TRIANGLE_ID(mesh, nodes[k], l) == -1) {
          SET_LAG_NODE_TRIANGLE_ID(mesh, nodes[k], l, tid);
          replaced_tid = true;
          break;
        }
    }
    /** Find the edges and: (i) add them to the list of edges of the new
    triangle; (ii) specify the id of the new triangle for the edges */
    /** First, we find the edge that connects node i to node i+1 */
    int ce = -1; // ce for "current edge
    for(int l=0; l<LAG_NODE_NB_NEIGHBORS(mesh, nodes[k]); l++) {
      ce = LAG_NODE_EDGE_ID(mesh, nodes[k], l);
      int cn = (LAG_EDGE_NODE_ID(mesh, ce, 0) == nodes[k]) ?
        LAG_EDGE_NODE_ID(mesh, ce, 1) : LAG_EDGE_NODE_ID(mesh, ce, 0);
      if (cn == nodes[(k+1)%3]) break;
    }
    /** Then, we add this edge to the list of edges of our new triangle */
    SET_LAG_TRIANGLE_EDGE_ID(mesh, tid, k, ce);
    /** And we have to update the triangle id of the edge */
    replaced_tid = false;
    for(int l=0; l<2; l++)
      if (LAG_EDGE_TRIANGLE_ID(mesh, ce, l) == prev_tid) {
        SET_LAG_EDGE_TRIANGLE_ID(mesh, ce, l, tid);
        replaced_tid = true;
        break;
      }
    if (!replaced_tid && LAG_EDGE_TRIANGLE_ID(mesh, ce, 0) != tid) {
      for(int l=0; l<2; l++)
        if (LAG_EDGE_TRIANGLE_ID(mesh, ce, l) == -1) {
          SET_LAG_EDGE_TRIANGLE_ID(mesh, ce, l, tid);
          replaced_tid = true;
          break;
        }
    }
  }
  mesh->nlt++;
}

void overwrite_triangle(lagMesh* mesh, int tid, int i, int j, int k) {
  int nodes[3];
  nodes[0] = i; nodes[1] = j; nodes[2] = k;
  for(int i=0; i<3; i++) {
    /** 1. Take care of the nodes of the triangle */
    SET_LAG_TRIANGLE_NODE_ID(mesh, tid, i, nodes[i]);
    /** Update the list of triangles for the node $i$ */
    bool already_there = false;
    for(int j=0; j<LAG_NODE_NB_TRIANGLES(mesh, nodes[i]); j++) {
      if (LAG_NODE_TRIANGLE_ID(mesh, nodes[i], j) == tid) already_there = true;
      else if (LAG_NODE_TRIANGLE_ID(mesh, nodes[i], j) == -1 && !already_there) {
        SET_LAG_NODE_TRIANGLE_ID(mesh, nodes[i], j, tid);
        break;
      }
    }

    /** 2. Take care of the egdes of the triangle */
    /** 2.1. Identify the edge id connecting nodes $i$, $i+1$ */
    int eid = -1; // eid for "edge id"
    for(int j=0; j<LAG_NODE_NB_NEIGHBORS(mesh, nodes[i]); j++) {
      eid = LAG_NODE_EDGE_ID(mesh, nodes[i], j);
      int cn = (LAG_EDGE_NODE_ID(mesh, eid, 0) == nodes[i]) ?
        LAG_EDGE_NODE_ID(mesh, eid, 1) : LAG_EDGE_NODE_ID(mesh, eid, 0);
      if (cn == nodes[(i+1)%3]) break;
    }
    /** 2.2 Add this edge to the triangle's list of edges */
    SET_LAG_TRIANGLE_EDGE_ID(mesh, tid, i, eid);
    /** 2.3 Add the triangle id to the edge's list of triangles */
    int index = (LAG_EDGE_TRIANGLE_ID(mesh, eid, 0) > -1 &&
      LAG_EDGE_TRIANGLE_ID(mesh, eid, 0) != tid) ? 1 : 0;
    SET_LAG_EDGE_TRIANGLE_ID(mesh, eid, index, tid);
  }
}

/** The function below returns true if node j is a vertex of triangle i */
bool is_triangle_vertex(lagMesh* mesh, int i, int j) {
  for(int k=0; k<3; k++) {
    if (LAG_TRIANGLE_NODE_ID(mesh, i, k) == j) return true;
  }
  return false;
}

/**
The function below returns the (positive) angle between the two vectors formed
by the nodes [*n1, *n2] and [*n1, *n3].
*/
double comp_angle(lagNode* n1, lagNode* n2, lagNode* n3) {
  double theta = 0.;
  foreach_dimension() {
    theta += GENERAL_1DIST(n2->pos.x, n1->pos.x, L0*L0_ratio.x)*
             GENERAL_1DIST(n3->pos.x, n1->pos.x, L0*L0_ratio.x);
  }
  double norm1 = GENERAL_SQNORM(n1->pos, n2->pos);
  double norm2 = GENERAL_SQNORM(n1->pos, n3->pos);
  theta /= sqrt(norm1)*sqrt(norm2);
  theta = acos(theta);
  return theta;
}

/**
The function below returns true if the angle of node $i$ ($0 < i < 2$) of
triangle $tid$ is greater than $\pi$ radians.
*/
bool is_obtuse_node(lagMesh* mesh, int tid, int i) {
  lagNode* n[3];
  // n[0] is the node to check if obtuse in the triangle tid
  n[0] = &(mesh->nodes[i]);
  int count_pts = 1;
  for(int j=0; j<3; j++)
  {
      if(LAG_TRIANGLE_NODE_ID(mesh, tid, j) != i)
      {
        n[count_pts] = &(mesh->nodes[LAG_TRIANGLE_NODE_ID(mesh, tid, j)]);
        count_pts ++;
      }
  }

  // Check if three points of the triangle are all found and well assigned
  assert(count_pts==3);

  if (comp_angle(n[0], n[1], n[2]) > pi/2.) return true;
  else return false;
}

/**
The function below returns true if the triangle $tid$ is obtuse at any of its
three angles.
*/
bool is_obtuse_triangle(lagMesh* mesh, int tid) {
  lagNode* n[3];
  for(int j=0; j<3; j++)
    n[j] = &(mesh->nodes[LAG_TRIANGLE_NODE_ID(mesh, tid, j)]);
  if (comp_angle(n[0], n[1], n[2]) > pi/2.) return true;
  else if (comp_angle(n[1], n[2], n[0]) > pi/2.) return true;
  else if (comp_angle(n[2], n[0], n[1]) > pi/2.) return true;
  else return false;
}

// /**
// The function below returns true if the triangle $tid$ is obtuse at any of its
// three angles.
// */
// bool is_obtuse_triangle(lagMesh* mesh, int tid) {
//   for(int i=0; i<3; i++) if (is_obtuse_node(mesh, tid, i)) return true;
//   return false;
// }

/**
## Uniform refinement of a mesh by subdividing its triangles

The function below loops through all triangles in the mesh and divide them
in four smaller ones. Keeping the correct structure of the mesh, i.e. updating
the relationship between nodes, edges, triangles and their neighbors results
in the somewhat complicated implementation below. */
void refine_mesh(lagMesh* mesh) {
  /** Perform the loop subdivision algorithm: until we reach the desired number
  of nodes, we split each triangles into four smaller ones */
  int cnt = mesh->nlt; // current number of triangles
  for(int i=0; i<cnt; i++) {
    int mid_ids[3];
    /** If not done yet, we split each edge into two */
    for(int j=0; j<3; j++) {
      int edge_id = LAG_TRIANGLE_EDGE_ID(mesh, i, j);
      if (LAG_EDGE_TRIANGLE_ID(mesh, edge_id, 1) > -1) {
        split_edge(mesh, edge_id);
        mid_ids[j] = mesh->nln-1;
      }
      else {
        mid_ids[j] = is_triangle_vertex(mesh, i,
          LAG_EDGE_NODE_ID(mesh, edge_id, 0)) ?
          LAG_EDGE_NODE_ID(mesh, edge_id, 1) :
          LAG_EDGE_NODE_ID(mesh, edge_id, 0);
      }
    }

    /** Connect the three midpoints with edges, and create corner triangles */
    for(int j=0; j<3; j++) {
      /** create edges between midpoints */
      new_edge(mesh, mid_ids[j], mid_ids[(j+1)%3]);

      /** Create the corner triangle with the new edge */
      int corner_id = -1;
      for(int k=0; k<3; k++) { // Loop over the current triangle vertices
        corner_id = LAG_TRIANGLE_NODE_ID(mesh, i, k);
        if (is_neighbor_ft(mesh, corner_id, mid_ids[j]) &&
          is_neighbor_ft(mesh, corner_id, mid_ids[(j+1)%3])) {
          break;
        }
      }
      new_triangle(mesh, mid_ids[j], mid_ids[(j+1)%3], corner_id, i);
    }
    /** Shrink the original (big) triangle into the center smaller one */
    overwrite_triangle(mesh, i, mid_ids[0], mid_ids[1], mid_ids[2]);
  }
}


/**
## Periodicity helper functions

In some situations, for instance to compute the volume of a capsule, we need to
take the dot and cross products of
neighboring nodes that can be across periodic boundaries. The next three
functions implement ``periodic-friendly" versions of the dot and cross products
that do not take into account the coordinates jump across the periodic
boundaries. The implementation relies on the assumption that a capsule is
*always* smaller in the x, y and z directions that half the domain size L0/2.
*/

/**
The function below corrects the coordinates of one node $a$ in order to
ensure it is placed on the same side of a reference coordinate (in practice,
the centroid of the capsule). The function returns the corrected node
coordinate.
*/
coord correct_periodic_node_pos(coord a, coord ref) {
    coord result;
    foreach_dimension() {
        result.x = (fabs(a.x - ref.x) < L0*L0_ratio.x/2) ? a.x :
            (a.x > ref.x) ? a.x - L0*L0_ratio.x : a.x + L0*L0_ratio.x;
    }
    return result;
}

/**
The function below corrects the coordinates of two nodes $a$ and $b$ in order to
ensure they are placed on the same side of a reference coordinate (in practice,
the centroid of the capsule). The result is stored in an array of `coord` of
length 2.
*/
void correct_periodic_nodes_pos(coord* result, coord a, coord b, coord ref) {
  foreach_dimension() {
    result[0].x = (fabs(a.x - ref.x) < L0*L0_ratio.x/2) ? a.x :
      (a.x > ref.x) ? a.x - L0*L0_ratio.x : a.x + L0*L0_ratio.x;
    result[1].x = (fabs(b.x - ref.x) < L0*L0_ratio.x/2) ? b.x :
      (b.x > ref.x) ? b.x - L0*L0_ratio.x : b.x + L0*L0_ratio.x;
  }
}

/**
The function below computes the cross product of two coordinates $a$ and $b$
that potentially lie across periodic boundaries. The coordinates are
temporarilly moved on the same side of the periodic boundary as a reference
coordinate `ref`, in practice the centroid of the capsule.
*/
foreach_dimension()
double periodic_friendly_cross_product_x(coord a, coord b, coord ref) {
  coord nodes[2];
  correct_periodic_nodes_pos(nodes, a, b, ref);
  return nodes[0].y*nodes[1].z - nodes[0].z*nodes[1].y;
}

/**
The function below computes the dot product of two coordinates $a$ and $b$
that potentially lie across periodic boundaries. The coordinates are
temporarilly moved on the same side of the periodic boundary as a reference
coordinate `ref`, in practice the centroid of the capsule.
*/
double periodic_friendly_dot_product(coord a, coord b, coord ref) {
  coord nodes[2];
  correct_periodic_nodes_pos(nodes, a, b, ref);
  return cdot(nodes[0], nodes[1]);
}





/**
The function below computes the nodal area of node $i$.
*/

#define cot(x) (cos(x)/sin(x))

trace
double compute_node_area(lagMesh* mesh, int i) {
  double area = 0.;
  for(int j=0; j<LAG_NODE_NB_TRIANGLES(mesh, i); j++) {
    int tid = LAG_NODE_TRIANGLE_ID(mesh, i, j);
    if (is_obtuse_triangle(mesh, tid)) {
      area += (is_obtuse_node(mesh, tid, i)) ? mesh->triangles[tid].area/2 :
        mesh->triangles[tid].area/4;
    }
    else {
      double voronoi_area = 0.;
      for(int k=0; k<3; k++) {
        int nid = LAG_TRIANGLE_NODE_ID(mesh, tid, k);
        if (nid != i) {
          int eid = -1; // eid for "edge id", connecting i and nid
          for (int l=0; l<3; l++) {
            int teid = LAG_TRIANGLE_EDGE_ID(mesh, tid, l); // temporary edge id
            if ((LAG_EDGE_NODE_ID(mesh, teid, 0) == i ||
              LAG_EDGE_NODE_ID(mesh, teid, 1) == i) &&
              (LAG_EDGE_NODE_ID(mesh, teid, 0) == nid ||
              LAG_EDGE_NODE_ID(mesh, teid, 1) == nid)) eid = teid;
          }
          int onid[2]; // onid for "opposite node ids"
          // find onid[0], onid[1]
          for(int l=0; l<2; l++) {
            // tneid for "triangle neighboring eid"
            int tneid = LAG_EDGE_TRIANGLE_ID(mesh, eid, l);
            for(int m=0; m<3; m++)
              if (LAG_TRIANGLE_NODE_ID(mesh, tneid, m) != i &&
                LAG_TRIANGLE_NODE_ID(mesh, tneid, m) != nid)
                onid[l] = LAG_TRIANGLE_NODE_ID(mesh, tneid, m);
          }
          // compute their angle facing the relevant triangle
          double theta[2];
          for(int l=0; l<2; l++) {
            theta[l] = comp_angle(&(mesh->nodes[onid[l]]), &(mesh->nodes[i]),
              &(mesh->nodes[nid]));
          }
          // compute the squared length of [i:nid]
          double edge_length = 0.;
          foreach_dimension()
            edge_length += sq(GENERAL_1DIST(mesh->nodes[i].pos.x,
              mesh->nodes[nid].pos.x, L0*L0_ratio.x));
          voronoi_area += (cot(theta[0]) + cot(theta[1]))*edge_length;
        }
      }
      area += voronoi_area/16.;
    }
  }
  return area;
}





#endif // dimension > 2
