/**
# Dump and restore functions for triangulated membranes
The functions below write/read the lagrangian mesh to/from a file in order
to restart simulations.

At the moment these functions are only valid for three-dimensional simulations.
*/

void dump_lagnode(FILE* fp, lagMesh* mesh, int i) {
  lagNode* node = &mesh->nodes[i];
  foreach_dimension() fwrite(&(node->pos.x), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(node->lagVel.x), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(node->normal.x), sizeof(double), 1, fp);
  fwrite(&(node->curv), sizeof(double), 1, fp);
  fwrite(&(node->gcurv), sizeof(double), 1, fp);
  fwrite(&(node->ref_curv), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(node->lagForce.x), sizeof(double), 1, fp);
  int nb_neighbors = LAG_NODE_NB_NEIGHBORS(mesh, i);
  fwrite(&nb_neighbors, sizeof(int), 1, fp);
  for(int j=0; j<6; j++) {
    int id = LAG_NODE_NEIGHBOR_ID(mesh, i, j);
    fwrite(&id, sizeof(int), 1, fp);
  }
  for(int j=0; j<6; j++) {
    int id = LAG_NODE_EDGE_ID(mesh, i, j);
    fwrite(&id, sizeof(int), 1, fp);
  }
  int nb_triangles = LAG_NODE_NB_TRIANGLES(mesh, i);
  fwrite(&nb_triangles, sizeof(int), 1, fp);
  for(int j=0; j<6; j++) {
    int id = LAG_NODE_TRIANGLE_ID(mesh, i, j);
    fwrite(&id, sizeof(int), 1, fp);
  }
  fwrite(&(node->nb_fit_iterations), sizeof(int), 1, fp);
}

void restore_lagnode(FILE* fp, lagMesh* mesh, int i) {
  lagNode* node = &mesh->nodes[i];
  foreach_dimension() fread(&(node->pos.x), sizeof(double), 1, fp);
  foreach_dimension() fread(&(node->lagVel.x), sizeof(double), 1, fp);
  foreach_dimension() fread(&(node->normal.x), sizeof(double), 1, fp);
  fread(&(node->curv), sizeof(double), 1, fp);
  fread(&(node->gcurv), sizeof(double), 1, fp);
  fread(&(node->ref_curv), sizeof(double), 1, fp);
  foreach_dimension() fread(&(node->lagForce.x), sizeof(double), 1, fp);
  int value;
  fread(&value, sizeof(int), 1, fp);
  SET_LAG_NODE_NB_NEIGHBORS(mesh, i, value);
  for(int j=0; j<6; j++) {
    fread(&value, sizeof(int), 1, fp);
    SET_LAG_NODE_NEIGHBOR_ID(mesh, i, j, value);
  }
  for(int j=0; j<6; j++) {
    fread(&value, sizeof(int), 1, fp);
    SET_LAG_NODE_EDGE_ID(mesh, i, j, value);
  }
  fread(&value, sizeof(int), 1, fp);
  SET_LAG_NODE_NB_TRIANGLES(mesh, i, value);
  for(int j=0; j<6; j++) {
    fread(&value, sizeof(int), 1, fp);
    SET_LAG_NODE_TRIANGLE_ID(mesh, i, j, value);
  }
  fread(&(node->nb_fit_iterations), sizeof(int), 1, fp);
}

void dump_edge(FILE* fp, lagMesh* mesh, int i) {
  Edge* edge = &mesh->edges[i];
  int id = LAG_EDGE_NODE_ID(mesh, i, 0);
  fwrite(&id, sizeof(int), 1, fp);
  id = LAG_EDGE_NODE_ID(mesh, i, 1);
  fwrite(&id, sizeof(int), 1, fp);
  id = LAG_EDGE_TRIANGLE_ID(mesh, i, 0);
  fwrite(&id, sizeof(int), 1, fp);
  id = LAG_EDGE_TRIANGLE_ID(mesh, i, 1);
  fwrite(&id, sizeof(int), 1, fp);
  double value = LAG_EDGE_L0(mesh, i);
  fwrite(&value, sizeof(double), 1, fp);
  fwrite(&(edge->length), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(edge->normal.x), sizeof(double), 1, fp);
}

void restore_edge(FILE* fp, lagMesh* mesh, int i) {
  Edge* edge = &mesh->edges[i];
  int value;
  fread(&value, sizeof(int), 1, fp);
  SET_LAG_EDGE_NODE_ID(mesh, i, 0, value);
  fread(&value, sizeof(int), 1, fp);
  SET_LAG_EDGE_NODE_ID(mesh, i, 1, value);
  fread(&value, sizeof(int), 1, fp);
  SET_LAG_EDGE_TRIANGLE_ID(mesh, i, 0, value);
  fread(&value, sizeof(int), 1, fp);
  SET_LAG_EDGE_TRIANGLE_ID(mesh, i, 1, value);
  double dvalue;
  fread(&dvalue, sizeof(double), 1, fp);
  SET_LAG_EDGE_L0(mesh, i, dvalue);
  fread(&(edge->length), sizeof(double), 1, fp);
  foreach_dimension() fread(&(edge->normal.x), sizeof(double), 1, fp);
}

void dump_triangle(FILE* fp, lagMesh* mesh, int i) {
  Triangle* triangle = &mesh->triangles[i];
  for(int k=0; k<3; k++) {
    int id = LAG_TRIANGLE_NODE_ID(mesh, i, k);
    fwrite(&id, sizeof(int), 1, fp);
    id = LAG_TRIANGLE_EDGE_ID(mesh, i, k);
    fwrite(&id, sizeof(int), 1, fp);
  }
  fwrite(&(triangle->area), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(triangle->normal.x), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(triangle->centroid.x), sizeof(double), 1, fp);
  foreach_dimension() {
    double value = LAG_TRIANGLE_REFSHAPE_COMPONENT(mesh, i, 0, x);
    fwrite(&value, sizeof(double), 1, fp);
    value = LAG_TRIANGLE_REFSHAPE_COMPONENT(mesh, i, 1, x);
    fwrite(&value, sizeof(double), 1, fp);
  }
  for(int k=0; k<3; k++)
    for(int l=0; l<2; l++) {
      double value = LAG_TRIANGLE_SFC(mesh, i, k, l);
      fwrite(&value, sizeof(double), 1, fp);
    }
  fwrite(&(triangle->stretch[0]), sizeof(double), 2, fp);
  fwrite(&(triangle->tension[0]), sizeof(double), 2, fp);
}

void restore_triangle(FILE* fp, lagMesh* mesh, int i) {
  Triangle* triangle = &mesh->triangles[i];
  for(int k=0; k<3; k++) {
    int value;
    fread(&value, sizeof(int), 1, fp);
    SET_LAG_TRIANGLE_NODE_ID(mesh, i, k, value);
    fread(&value, sizeof(int), 1, fp);
    SET_LAG_TRIANGLE_EDGE_ID(mesh, i, k, value);
  }
  fread(&(triangle->area), sizeof(double), 1, fp);
  foreach_dimension() fread(&(triangle->normal.x), sizeof(double), 1, fp);
  foreach_dimension() fread(&(triangle->centroid.x), sizeof(double), 1, fp);
  foreach_dimension() {
    double value;
    fread(&value, sizeof(double), 1, fp);
    SET_LAG_TRIANGLE_REFSHAPE_COMPONENT(mesh, i, 0, x, value);
    fread(&value, sizeof(double), 1, fp);
    SET_LAG_TRIANGLE_REFSHAPE_COMPONENT(mesh, i, 1, x, value);
  }
  for(int k=0; k<3; k++)
    for(int l=0; l<2; l++) {
      double value;
      fread(&value, sizeof(double), 1, fp);
      SET_LAG_TRIANGLE_SFC(mesh, i, k, l, value);
    }
  fread(&(triangle->stretch[0]), sizeof(double), 2, fp);
  fread(&(triangle->tension[0]), sizeof(double), 2, fp);
}

void dump_lagmesh(FILE* fp, lagMesh* mesh) {
  int has_storage = mesh->isactive && mesh->nodes != NULL &&
    mesh->edges != NULL;
  #if dimension > 2
    has_storage = has_storage && mesh->triangles != NULL;
  #endif
  int dump_nln = has_storage ? mesh->nln : 0;
  int dump_nle = has_storage ? mesh->nle : 0;
  #if dimension > 2
    int dump_nlt = has_storage ? mesh->nlt : 0;
  #else
    int dump_nlt = 0;
  #endif
  int dump_isactive = has_storage && mesh->isactive;

  fwrite(&(mesh->cap_id), sizeof(int), 1, fp);
  fwrite(&(mesh->cap_type), sizeof(int), 1, fp);
  fwrite(&(mesh->cap_es), sizeof(double), 1, fp);
  fwrite(&(mesh->cap_radius), sizeof(double), 1, fp);
  fwrite(&dump_nln, sizeof(int), 1, fp);
  for(int i=0; i<dump_nln; i++) dump_lagnode(fp, mesh, i);
  fwrite(&dump_nle, sizeof(int), 1, fp);
  for(int i=0; i<dump_nle; i++) dump_edge(fp, mesh, i);
  fwrite(&dump_nlt, sizeof(int), 1, fp);
  for(int i=0; i<dump_nlt; i++)
    dump_triangle(fp, mesh, i);
  foreach_dimension() fwrite(&(mesh->centroid.x), sizeof(double), 1, fp);
  foreach_dimension() fwrite(&(mesh->ang_vel.x), sizeof(double), 1, fp);
  fwrite(&(mesh->volume), sizeof(double), 1, fp);
  fwrite(&(mesh->circum_radius), sizeof(double), 1, fp);
  fwrite(&(mesh->taylor_deform), sizeof(double), 1, fp);
  double initial_volume = 0.;
  #if LAG_REF_GEOMETRY
    if (mesh->ref_geometry != NULL)
      initial_volume = LAG_INITIAL_VOLUME(mesh);
  #else
    initial_volume = LAG_INITIAL_VOLUME(mesh);
  #endif
  fwrite(&initial_volume, sizeof(double), 1, fp);
  int tmp;
  tmp = mesh->updated_stretches ? 1 : 0; fwrite(&(tmp), sizeof(int), 1, fp);
  tmp = mesh->updated_normals ? 1 : 0; fwrite(&(tmp), sizeof(int), 1, fp);
  tmp = mesh->updated_curvatures ? 1 : 0; fwrite(&(tmp), sizeof(int), 1, fp);
  tmp = dump_isactive ? 1 : 0; fwrite(&(tmp), sizeof(int), 1, fp);
}

void restore_lagmesh(FILE* fp, lagMesh* mesh) {
  mesh->topology = NULL;
  mesh->ref_geometry = NULL;
  fread(&(mesh->cap_id), sizeof(int), 1, fp);
  fread(&(mesh->cap_type), sizeof(int), 1, fp);
  fread(&(mesh->cap_es), sizeof(double), 1, fp);
  fread(&(mesh->cap_radius), sizeof(double), 1, fp);
  fread(&(mesh->nln), sizeof(int), 1, fp);
  mesh->nodes = malloc(mesh->nln*sizeof(lagNode));
  mesh->topology = allocate_lag_topology(mesh->nln, 0, 0);
  for(int i=0; i<mesh->nln; i++) restore_lagnode(fp, mesh, i);
  fread(&(mesh->nle), sizeof(int), 1, fp);
  mesh->edges = malloc(mesh->nle*sizeof(Edge));
  resize_lag_topology(mesh->topology, mesh->nln, mesh->nle, 0);
  #if LAG_REF_GEOMETRY
    mesh->ref_geometry = allocate_lag_ref_geometry(mesh->nle, 0);
  #endif
  for(int i=0; i<mesh->nle; i++) restore_edge(fp, mesh, i);
  fread(&(mesh->nlt), sizeof(int), 1, fp);
  mesh->triangles = malloc(mesh->nlt*sizeof(Triangle));
  resize_lag_topology(mesh->topology, mesh->nln, mesh->nle, mesh->nlt);
  #if LAG_REF_GEOMETRY
    resize_lag_ref_geometry(mesh->ref_geometry, mesh->nle, mesh->nlt);
  #endif
  for(int i=0; i<mesh->nlt; i++) restore_triangle(fp, mesh, i);
  foreach_dimension() fread(&(mesh->centroid.x), sizeof(double), 1, fp);
  foreach_dimension() fread(&(mesh->ang_vel.x), sizeof(double), 1, fp);
  fread(&(mesh->volume), sizeof(double), 1, fp);
  fread(&(mesh->circum_radius), sizeof(double), 1, fp);
  fread(&(mesh->taylor_deform), sizeof(double), 1, fp);
  double initial_volume;
  fread(&initial_volume, sizeof(double), 1, fp);
  SET_LAG_INITIAL_VOLUME(mesh, initial_volume);
  int tmp;
  fread(&(tmp), sizeof(int), 1, fp);
  mesh->updated_stretches = (tmp == 0) ? false : true;
  fread(&(tmp), sizeof(int), 1, fp);
  mesh->updated_normals = (tmp == 0) ? false : true;
  fread(&(tmp), sizeof(int), 1, fp);
  mesh->updated_curvatures = (tmp == 0) ? false : true;
  fread(&(tmp), sizeof(int), 1, fp);
  mesh->isactive = (tmp == 0) ? false : true;
}

/** If the simulation contains several membranes, we dump and read one mesh at
a time, but all meshes are stored in the same file. */
struct _dump_capsules {
    char* name;
    FILE* fp;
};

void dump_capsules(const char* fname, FILE* fp) {
    #if _MPI && LAG_DISTRIBUTED_CAPSULES && LAG_DISTRIBUTED_OUTPUT_GATHER_RANK0
      if (pid() != 0)
        return;
    #endif
    const char default_name[10] = "caps.dump\0";
    const char* name = fname ? fname : default_name;
    FILE* file = fp ? fp : fopen(name, "w");
    assert(file);
    for(int i=0; i<NCAPS; i++) dump_lagmesh(file, &CAPS(i));
    if(!fp) fclose(file);
}

void restore_capsules(char* filename) {
  FILE* file = fopen(filename, "r");
  assert(file);
  for(int i=0; i<NCAPS; i++) {
    restore_lagmesh(file, &CAPS(i));
    if (CAPS(i).isactive) {
      attach_shared_lag_topology(&CAPS(i));
      attach_shared_lag_ref_geometry(&CAPS(i));
      debug_lag_topology(&CAPS(i), "restore_capsules");
    }
  }
  fclose(file);
  initialize_all_capsules_stencils();
  generate_lag_stencils(true);
}
